from pathlib import Path
import csv,shutil,re,html
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path.cwd();P=R/'tests/solids/hexa8_tied';rows=list(csv.DictReader((P/'archive_comparison/comparison.csv').open()));out=R/'update/solid_benchmark/data/hexa8_tied_vs_eas9v2_archive.csv';shutil.copy2(P/'archive_comparison/comparison.csv',out)
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for ax,(case,response) in zip(axes.flat,[('5-004','IN'),('5-004','OUT'),('5-007','A'),('5-007','B')]):
 for variant,label,color in [('EAS9v2','MYSTRAN EAS9v2','#4267ac'),('TIED','MYSTRAN TIED','#d17a22'),('nastran','Nastran (existing benchmark)','#23805a')]:
  q=[x for x in rows if x['case']==case and x['response']==response and x['variant']==('EAS9v2' if variant=='nastran' else variant)];q.sort(key=lambda x:int(x['mesh'].split('x')[0]));ax.plot([int(x['mesh'].split('x')[0]) for x in q],[float(x['nastran'] if variant=='nastran' else x['mystran'])/float(x['independent']) for x in q],'-o',color=color,label=label)
 ax.axhline(1,color='black',ls='--',lw=1,label='Independent reference');ax.set(title=case+' — '+response,xlabel='Longitudinal / in-plane mesh divisions',ylabel='Displacement / independent',xticks=[6,12,24] if case=='5-004' else [4,8,16,24]);ax.grid(alpha=.2)
axes[0,0].legend(fontsize=8);fig.suptitle('CHEXA8 — TIED vs EAS9v2 on the existing beam/dome benchmark',fontsize=14)
for dest in [R/'update/solid_benchmark/figures/CHEXA8_TIED_20261018.png',R/'update/solid_benchmark/figures/CHEXA8.png',R/'update/Formulation/Solid/figures/CHEXA8_benchmark.png']:fig.savefig(dest,dpi=170)
plt.close(fig)
intro='''## CHEXA8

### TIED vs EAS9v2 on the existing benchmark

NEWSOLID CHEXA8 defaults to `PARAM,HEXATYP,TIED`; `PARAM,HEXATYP,EAS9` retains EAS9v2. The following comparison uses the **existing beam and dome benchmark decks and the previously published Nastran results**. Only MYSTRAN formulation selectors are changed. No Nastran values from the new benchmark suite are used here.

The graph compares TIED with its default ZETA policy, EAS9v2, and the original Nastran reference on identical archived meshes. The EAS9v2 rerun agrees with the published archive within output precision. Other solid-element sections are unchanged.

![CHEXA8 TIED vs EAS9v2 and existing Nastran](figures/CHEXA8_TIED_20261018.png)

[Comparison CSV and source-deck hashes](data/hexa8_tied_vs_eas9v2_archive.csv)

| Case | Mesh | Response | MYSTRAN formulation | Independent | Existing Nastran | MYSTRAN | N error % | M error % |
|---|---|---|---|---:|---:|---:|---:|---:|
'''
order=sorted([x for x in rows if x['variant'] in ['EAS9v2','TIED']],key=lambda x:(x['case'],int(x['mesh'].split('x')[0]),x['response'],x['variant']))
for x in order:intro+=f"| {x['case']} | {x['mesh']} | {x['response']} | {x['variant']} | {float(x['independent']):.9g} | {float(x['nastran']):.9g} | {float(x['mystran']):.9g} | {float(x['nastran_error_pct']):.6g} | {float(x['mystran_error_pct']):.6g} |\n"
intro+='\nOn the 24x24x1 dome, point A is 0.09328662 for TIED, 0.08710673 for EAS9v2, and 0.09335861 for the existing Nastran reference (independent: 0.094). Corresponding errors are 0.7589%, 7.3333%, and 0.6823%. On the finest beam mesh, IN error is 0.60884% for TIED versus 0.60908% for EAS9v2.\n\n'
p=R/'update/solid_benchmark/comparison_all.md';s=p.read_text(encoding='utf-8');a=s.index('## CHEXA8\n');b=s.index('## CHEXA20\n',a);p.write_bytes((s[:a]+intro+s[b:]).encode())
p=R/'README.md';s=p.read_text(encoding='utf-8');a=s.index('## CHEXA8 TIED update');b=s.index('## Solid recovery update',a)
new='''## CHEXA8 TIED update — 18 October 2026

- NEWSOLID CHEXA8 defaults to `PARAM,HEXATYP,TIED`, using mean-gradient-corrected physical shear tying with EAS9 condensation. `PARAM,HEXATYP,EAS9` retains EAS9v2.
- The [updated graph](update/solid_benchmark/comparison_all.md#chexa8) compares TIED and EAS9v2 on the existing beam/dome benchmark meshes against the previously published Nastran results.
- On the 24x24x1 dome, point A error is **0.759% for TIED**, **7.333% for EAS9v2**, and **0.682% for Nastran**. Default TIED uses ZETA; `HEXATAX` also provides explicit XI/ETA and opt-in AUTO controls.

[Windows binary and runtime](update/binaries/mystran18a-chexa8-tied-win64.zip) · [CHEXA8 formulation](update/Formulation/Solid/NEWSOLID_CONTINUUM_FORMULATIONS.md#35-chexa8-tied-default-and-retained-eas9-option) · [Benchmark comparison](update/solid_benchmark/comparison_all.md#chexa8)

'''
s=s[:a]+new+s[b:];s=s.replace('`TIED` (default, provisional)','`TIED` (default)');p.write_bytes(s.encode())
p=R/'update/solid_benchmark/README.md';s=p.read_text(encoding='utf-8');a=s.index('\n\n');s='> **18 October CHEXA8 update:** [TIED vs EAS9v2](comparison_all.md#chexa8) uses the existing beam/dome benchmark meshes and previously published Nastran results. [Comparison CSV](data/hexa8_tied_vs_eas9v2_archive.csv).'+s[a:];p.write_bytes(s.encode())
p=R/'update/Formulation/Solid/NEWSOLID_CONTINUUM_FORMULATIONS.md';s=p.read_text(encoding='utf-8').replace('`HEXATYP,TIED`, order 2; provisional','`HEXATYP,TIED`, order 2').replace(' TIED is explicitly **provisional**, pending official Nastran acceptance.','');a=s.index('Thirty-six operator comparisons');b=s.index('\n### 3.6',a);s=s[:a]+'''Thirty-six operator comparisons, six rigid modes, consistent/lumped mass checks, the full distorted patch and cyclic-axis audits pass. Explicit EAS9 matches the saved executable on 20 models. The current comparison uses TIED and EAS9v2 on the existing beam/dome benchmark decks, against the previously published Nastran results. On the 24x24x1 dome, TIED has 0.7589% point A error, EAS9v2 7.3333%, and Nastran 0.6823%. See the [comparison graph and table](../../solid_benchmark/comparison_all.md#chexa8).

'''+s[b:]
a=s.index('#### CHEXA8\n');b=s.index('#### CHEXA20\n',a);s=s[:a]+'''#### CHEXA8

![CHEXA8 TIED vs EAS9v2 and existing Nastran](figures/CHEXA8_benchmark.png)

TIED and EAS9v2 use the existing archived beam/dome meshes; Nastran values are the original benchmark results. On the 24x24x1 dome, point A is 0.09328662 (TIED), 0.08710673 (EAS9v2), and 0.09335861 (Nastran), against 0.094 independent. [Full comparison](../../solid_benchmark/comparison_all.md#chexa8).

'''+s[b:];p.write_bytes(s.encode())
# Replace only the CHEXA8 HTML section; preserve all other embedded diagrams.
p=R/'update/solid_benchmark/comparison_all.html';s=p.read_text(encoding='utf-8');a=s.index('<section data-element="CHEXA8">');b=s.index('<section data-element="CHEXA20">',a);body='<section data-element="CHEXA8"><h2>CHEXA8 — TIED vs EAS9v2</h2><p>Existing beam/dome benchmark meshes and previously published Nastran results. TIED uses default ZETA. <a href="comparison_all.md#chexa8">Full comparison table</a> · <a href="data/hexa8_tied_vs_eas9v2_archive.csv">CSV</a></p><img style="width:100%" src="figures/CHEXA8_TIED_20261018.png" alt="TIED versus EAS9v2 and existing Nastran benchmark"><p>Dome 24x24 point A: TIED 0.09328662 (0.7589% error), EAS9v2 0.08710673 (7.3333%), Nastran 0.09335861 (0.6823%), independent 0.094.</p></section>';p.write_bytes((s[:a]+body+s[b:]).encode())
p=R/'update/Formulation/Solid/NEWSOLID_CONTINUUM_FORMULATIONS.html';s=p.read_text(encoding='utf-8');s=re.sub(r'<aside style="padding:20px;border:2px solid #b78127;margin:20px">.*?</aside>','<aside style="padding:20px;border:2px solid #b78127;margin:20px"><strong>CHEXA8 TIED update:</strong> TIED is the NEWSOLID default; EAS9v2 remains selectable. <a href="NEWSOLID_CONTINUUM_FORMULATIONS.md#35-chexa8-tied-default-and-retained-eas9-option">Updated formulation</a> · <a href="../../solid_benchmark/comparison_all.md#chexa8">TIED vs EAS9v2 and existing Nastran benchmark</a>.</aside>',s,count=1);p.write_bytes(s.encode())
print('published comparison prepared: 21 MYSTRAN solves; 42 scalar responses; archived Nastran unchanged')
