from pathlib import Path
import sys,json,csv,hashlib,concurrent.futures,re
P=Path(__file__).resolve().parent;W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');sys.path.insert(0,str(W));sys.path.insert(0,str(W.parent/'codex_mod/mitc3_lee_nguyen_20261024'))
from run_suite import run
from results import parse_f06,measures
from build_suite import fields,cards
from report_tools import loaded_displacements_basic
forms={'CQUAD4_DSQK':'DSQK'}
manifest=[r for r in json.loads((W/'manifest.json').read_text()) if r['formulation'] in forms and str(r['solution'])=='101']
for form in forms:
 base=next(r for r in manifest if r['formulation']==form and r['group']=='2007')
 for n in (16,32): manifest.append(dict(base,case=f'hemisphere_n{n}',deck=str(W/f'investigations/convergence_2_007/n{n}/{form}/model.dat'),nodes=0))
jobs=[]
for rev in ['old','covariant']:
 for r in manifest:
  src=W/r['deck'];out=[]
  for line in src.read_text().splitlines():
   f=fields(line)
   if f and f[0]=='SPC1' and 'THRU' not in f and len(f)>9:
    out.extend(','.join(f[:3]+f[i:i+6]) for i in range(3,len(f),6))
   else:out.append(line)
  d=P/'runs'/rev/forms[r['formulation']]/r['case']/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text('\n'.join(out)+'\n')
  jobs.append((rev,dict(r,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest())))
def task(j):
 rev,r=j;s=run(r,P/('bin' if rev=='covariant' else 'baseline')/'mystran.exe',300,True);data=parse_f06(Path(r['deck']).with_suffix('.F06'),'MYSTRAN')
 data['stress']={sc:[v for v in vals if abs(v[1])>0] for sc,vals in data['stress'].items()}
 obs=measures(r,data)
 if r['group']=='2007' and s['status']=='completed':
  b=loaded_displacements_basic(Path(r['deck']),data)
  for m in obs:
   for card in cards(Path(r['deck']).read_text()):
    if card[0]=='FORCE':
     dof=next(i for i in range(3) if float(card[5+i])!=0)
     if m['quantity']=='loaded_U'+('X','Y','Z')[dof]: m['value']=float(b[int(card[2])][dof])
 mesh=(re.search(r'nx(\d+)',r['case']) or re.search(r'hemisphere_n(\d+)',r['case']))
 return dict(element=forms[r['formulation']],revision=rev,case=r['case'],mesh=mesh[1] if mesh else '',status=s['status'],errors=s['errors'],provenance=r['deck'],observations=obs)
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for res in pool.map(task,jobs):
  results.append(res)
  (P/'static_results.json').write_text(json.dumps(results,indent=2))
  if len(results)%20==0:print(len(results),'/',len(jobs),'completed',sum(r['status']=='completed' for r in results),flush=True)
with (P/'static_comparison.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=['element','revision','case','mesh','loadcase','observable','reference','computed','abs_error_pct','status','provenance']);writer.writeheader()
 for r in results:
  for m in r['observations']:
   ref=m['reference'];v=m['value'];err=100*abs(v/ref-1) if ref and v is not None else None
   writer.writerow(dict(element=r['element'],revision=r['revision'],case=r['case'],mesh=r['mesh'],loadcase=m['loadcase'],observable=m['quantity'],reference=ref,computed=v,abs_error_pct=err,status=r['status'] if v is not None else 'missing_observation',provenance=r['provenance']))
print('FINAL',len(results),'completed',sum(r['status']=='completed' for r in results),flush=True)
