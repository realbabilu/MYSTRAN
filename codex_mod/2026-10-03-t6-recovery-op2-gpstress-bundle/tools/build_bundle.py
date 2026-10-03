"""Create a reviewable commit bundle from the GitHub baseline and local evidence."""
from pathlib import Path
import argparse, ast, csv, difflib, hashlib, importlib.metadata, json, shutil, subprocess

ap=argparse.ArgumentParser()
ap.add_argument('--workspace',type=Path,default=Path(__file__).resolve().parents[1])
args=ap.parse_args()
workspace=args.workspace.resolve()
repo=workspace/'MYSTRAN'
name='2026-10-03-t6-recovery-op2-gpstress-bundle'
bundle=repo/'codex_mod'/name
bundle.mkdir(parents=True,exist_ok=True)

def git(*args):
    return subprocess.check_output(['git','-c','core.quotepath=false',*args],cwd=repo)
baseline=git('rev-parse','931dde4').decode().strip()
head=git('rev-parse','HEAD').decode().strip()
assert baseline==head,'Review the baseline again after committing; this is a pre-commit packager.'
branch=git('branch','--show-current').decode().strip()
changes=[line.split('\t') for line in git('diff','--name-status',baseline).decode().splitlines()]
tracked=[path for status,path in changes]
records=[]

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(src,rel,kind,provenance):
    src=Path(src)
    assert src.is_file(),src
    dest=bundle/rel
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(src,dest)
    records.append(dict(original=src.relative_to(workspace).as_posix(),bundle_path=rel,
                        kind=kind,provenance=provenance,bytes=src.stat().st_size,sha256=digest(src)))

for path in tracked:
    copy(repo/path,'files/'+path,'github-tracked-modified',f'diff from {baseline}')
(bundle/'SOURCE_FROM_931dde4.patch').write_bytes(git('diff','--binary',baseline,'--','Source'))
(bundle/'SOURCE_DIFFSTAT.txt').write_bytes(git('diff','--stat',baseline,'--','Source'))

seeds=['test_patch2001_q8_t6_v4c.py','battle_2_002_quadratic_clamped.py',
       'battle_2_003_quadratic_curved.py','battle_2_004_twistedv2.py','battle_2_006_scordelisv2.py',
       'test_quadraticv2a.py','compare_disp_v4c.py','compare_f06_vs_python.py',
       'test_result_2_001_t6_disp.py','dump_python_debug.py',
       'Simo1993_Tri6_ShellElement_v2.py','MITC6_Tri_v4.py','MacNeal_MH6T_Tri_v3.py','Rezaiee2017_Tri6_v3.py']
queue=list(seeds); seen=set(); absent=set()
backup=workspace/'test/t6_reference_alignment_backup'
verified_python=[]
while queue:
    filename=queue.pop()
    if filename in seen: continue
    src=workspace/'python'/filename
    if not src.exists(): continue
    seen.add(filename)
    old=backup/'python'/filename
    provenance='workspace reference/dependency snapshot; no GitHub baseline for this sibling file'
    if old.exists() and old.read_bytes()!=src.read_bytes():
        provenance='changed against saved pre-alignment local backup (not GitHub history)'
        verified_python.append(filename)
        patch=''.join(difflib.unified_diff(old.read_text(encoding='utf-8-sig').splitlines(True),src.read_text(encoding='utf-8-sig').splitlines(True),
                   fromfile='before/python/'+filename,tofile='after/python/'+filename))
        dest=bundle/'PYTHON_LOCAL_CHANGES'/f'{filename}.patch'
        dest.parent.mkdir(parents=True,exist_ok=True); dest.write_text(patch,encoding='utf-8')
    if filename=='dump_python_debug.py': provenance='updated during this session; no saved pre-edit baseline'
    copy(src,'validation/python/'+filename,'python-support',provenance)
    for node in ast.walk(ast.parse(src.read_text(encoding='utf-8-sig'))):
        modules=[]
        if isinstance(node,ast.Import): modules=[n.name.split('.')[0] for n in node.names]
        if isinstance(node,ast.ImportFrom) and node.level==0 and node.module:
            modules=[node.module.split('.')[0]]
        for module in modules:
            candidate=workspace/'python'/f'{module}.py'
            if candidate.exists(): queue.append(candidate.name)
            elif module in {'solver_f06_convergence','Simo1993_Q8_ShellElement_V2_MacNealPatched','Simo1993_Q8_ShellElement_v3'}:
                absent.add(module)

tests=['t6_displacement_gate.py','t6_stress_gate.py','t6_output_selection_gate.py','verify_final_t6_debug.py',
       'read_t6_op2.py','check_t6_op2_reference.py','check_duel3a_mixed_t6.py','prepare_op2_reference.py']
decks=['duel3a.dat','duel3a_simo_1.dat','duel3a_mitc6_1.dat','duel3a_mh6t_1.dat','duel3a_rezaiee.dat']
for filename in tests+decks:
    copy(workspace/'test'/filename,'validation/test/'+filename,'validation-runner' if filename in tests else 'input-deck',
         'sibling test directory outside the GitHub working tree')
copy(workspace/'test/compare_f06_vs_python.py','validation/test/compare_f06_vs_python.py','historical-helper',
     'legacy sibling copy; active comparator is validation/python/compare_f06_vs_python.py')

evidence=['T6_UPGRADE_VALIDATION_2026-10-03.md','T6_OP2_VALIDATION_2026-10-03.md','CTRIA6_SIMO_MITC_RESUME.md',
 't6_upgrade/displacement_refined/results.json','t6_upgrade/recovery_final/results.json',
 't6_upgrade/stress_recovery/results.json','t6_upgrade/output_selection/results.json',
 't6_upgrade/op2_validation.txt','t6_upgrade/debug_validation.txt',
 't6_upgrade/recovery_build.log','t6_upgrade/op2_build.log','t6_upgrade/force_header_build.log',
 'duel3a.pynastran.json','duel3a_pynastran.txt','duel3a_run_latest.compare.txt',
 'duel3a_force_header.stdout.txt','t6_upgrade/op2_reference/nastran_t6.dat',
 't6_upgrade/op2_reference/mystran_t6.dat','t6_upgrade/op2_reference/nastran_t6.pynastran.json',
 't6_upgrade/op2_reference/mystran_t6.pynastran.json','t6_upgrade/op2_reference/nastran_t6.stdout.txt']
for family in ('SIMOT6','MITC6','MH6T','REZAIEE'):
    evidence += [f't6_upgrade/stress_recovery/{family}{ext}' for ext in ('.dat','.compare.txt','.pynastran.json','.pynastran.txt')]
for path in evidence:
    copy(workspace/'test'/path,'evidence/'+path,'validation-evidence','recorded local run; not a new claim of validation')

historical=['upgrade_t6_formula.py','refine_t6_formula.py','upgrade_t6_recovery.py','finalize_t6_recovery.py',
            'remove_temporary_solver_debug.py','fix_ogs_surface_tables.py','restore_q8_scope.py']
for filename in historical:
    copy(workspace/'test'/filename,'history/one_time_migrations/'+filename,'one-time-history','migration record; do not rerun')
for filename in ('align_t6_tests.py','verify_t6_alignment.py'):
    copy(workspace/filename,'history/one_time_migrations/'+filename,'one-time-history','historical helper, superseded by final gates')
copy(Path(__file__), 'tools/build_bundle.py','packaging-tool','bundle generator')

versions={}
for package in ('numpy','scipy','matplotlib','pyNastran'):
    versions[package]=importlib.metadata.version(package)
(bundle/'validation/requirements-tested.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items())+'\n')

def describe(path):
    filename=Path(path).name
    if filename.startswith('CTRIA6_'): return 'Align formulation, signed directors, metric/drilling and seven-point recovery with final Python; preserve explicit SNORM.'
    if filename=='CQUAD8_MACQ8D.f90': return 'Existing local change: move RV/SV parameters to host scope; no new Q8 formula upgrade in this session.'
    if filename.startswith(('CQUAD4_','CQUADR_')): return 'Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta.'
    if filename in ('POLYNOM_FIT_STRE_STRN.f90','SHELL_STRESS_OUTPUTS.f90'): return 'Formatting/indentation delta; no behavior change claimed.'
    if filename in ('OFP1.f90','OFP2.f90','OFP3.f90','OFP3_ELFE_1D.f90','OFP3_STRE_PCOMP.f90','OFP3_STRN_NO_PCOMP.f90','OFP3_STRN_PCOMP.f90'):
        return 'Preserve local MAXREQ*5 output-buffer capacity/initialization change.'
    return {
      'CHK_CC_CMD_DESCRIBERS.f90':'Track CENTER/CORNER request flags independently for stress/strain/force.',
      'CC_OUTPUT_DESCRIBERS.f90':'Add CENTER/CORNER request flags.',
      'LOADB.f90':'Reserve recovery storage for seven T6 samples (MAX_STRESS_POINTS plus center slot).',
      'WRITE_ELEM_ENGR_FORCE.f90':'Native CTRIA6 type-75 OEF width 38 and corrected T6 CENTER stride/selectors.',
      'WRITE_ELEM_STRESSES.f90':'T6 F06 CENTER/six nodes; native OES width 70; signed basis GP conversion; actual OGS surface IDs and per-surface tables.',
      'ELEM_STRE_STRN_ARRAYS.f90':'Use basic-coordinate element displacement UEB for T6 recovery; remove temporary printing.',
      'OFP3_ELFE_2D.f90':'Recover/copy all seven T6 force samples; suppress empty surface force headers.',
      'OFP3_STRE_NO_PCOMP.f90':'T6 direct seven-point output; per-element metadata; native signed centroid basis.',
      'ALLOCATE_LINK9_STUF.f90':'Preserve local MAXREQ*5 allocation of OGEL and shell transform arrays.',
      'MAXREQ_OGEL.f90':'T6 row counts 7 force / 14 stress; preserve existing local output capacity policy.',
      'GPSTRESS_SURFACE_UTILS.f90':'Include TRIA6/QUAD8 in shell patches and collect all six/eight nodes.',
      'MODEL_STUF.f90':'Existing local TRIA6 NUM_SEi change 4 to 3; T6 stress/force output now explicitly reserves seven.',
      'OUTPUT2_WRITE_ELFORCE.f90':'Accept TRIA6 in OEF dispatch.',
      'OUTPUT2_WRITE_STRESS.f90':'Accept TRIA6 in OES dispatch; reopen OES after OGS closes it.'
    }.get(filename,'See baseline patch for exact change.')

rows=['# Daftar file berubah dan snapshot pendukung','',f'Baseline GitHub: `{baseline}`. Branch: `{branch}`.',
      '',f'## {len(tracked)} file tracked yang berubah','', '| File di repo | Ringkasan |','|---|---|']
rows += [f'| `{p}` | {describe(p)} |' for p in tracked]
rows += ['', '## Python dengan bukti perubahan terhadap backup lokal','',
         'Folder sibling `python` dan `test` bukan bagian clone GitHub. Diff-nya tidak disamakan dengan diff GitHub.', '']
rows += [f'- `{p}`' for p in sorted(verified_python)]
rows += ['', '`dump_python_debug.py` juga diperbarui pada sesi ini, tetapi tidak memiliki backup pra-edit.',
         '`duel3a_mht6_1.dat` diganti nama menjadi `duel3a_mh6t_1.dat`; helper perbandingan aktif disesuaikan.',
         '', '## Snapshot tambahan','',
         'Referensi T6 final, core/adaptor/dependensi Python yang diimpor, skrip validasi, deck, laporan, dan bukti hasil ikut disalin agar dapat direview. Snapshot dependensi tidak berarti file tersebut diedit pada sesi ini.',
         '', 'Lokasi dan checksum setiap salinan tercantum dalam `FILE_MANIFEST.csv` dan `manifest.json`.']
(bundle/'CHANGED_FILES.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')
with (bundle/'FILE_MANIFEST.csv').open('w',newline='',encoding='utf-8') as fh:
    writer=csv.DictWriter(fh,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
manifest=dict(date='2026-10-03',timezone='Asia/Jakarta',baseline=baseline,head=head,branch=branch,
              remote=git('remote','get-url','origin').decode().strip(),tracked_modified_count=len(tracked),
              tracked_modified=tracked,python_local_backup_verified=sorted(verified_python),
              missing_historical_imports=sorted(absent),tested_python_packages=versions,files=records)
(bundle/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(bundle=str(bundle),tracked_modified=len(tracked),python_snapshots=len(seen),
                     copied_files=len(records),missing_historical_imports=sorted(absent)),indent=2))
