from pathlib import Path
import sys,os,csv,json,subprocess,concurrent.futures,ast,hashlib
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];ARCH=ROOT.parent/'battle/solid_sap_comparison/all_solids';OUT=P/'archive_comparison';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(P));from collect_regression import read,vectors
rows=[r for r in csv.DictReader((ROOT/'update/solid_benchmark/data/comparison_all.csv').open()) if r['element']=='CHEXA8'];groups={}
for r in rows:groups.setdefault((r['case'],r['mesh']),[]).append(r)
env=os.environ.copy();env['PATH']='C:/gcc/bin;C:/gcc/openblas32/bin;'+env['PATH'];env['OPENBLAS_NUM_THREADS']='1'
def run(key,variant):
 case,mesh=key;src=ARCH/'runs'/case/mesh/'CHEXA8/mystran_NEWSOLID/model.dat';nx=src.parents[1]/'nastran/model.dat';d=OUT/case/mesh/variant;d.mkdir(parents=True,exist_ok=True)
 text=src.read_text();form='EAS9' if variant=='EAS9v2' else 'TIED';axis='AUTO' if variant=='TIED_AUTO' else 'ZETA';deck=d/'model.dat';deck.write_text(text.replace('BEGIN BULK',f'BEGIN BULK\nPARAM,HEXATYP,{form}\nPARAM,HEXATAX,{axis}'))
 # Only selector cards are added to the exact archived MYSTRAN deck.
 assert '\n'.join(l for l in deck.read_text().splitlines() if not l.startswith(('PARAM,HEXATYP','PARAM,HEXATAX')))==text.rstrip('\n')
 p=subprocess.run([str(ROOT/'Binaries/mystran.exe'),deck.name],cwd=d,env=env,capture_output=True,text=True,timeout=300);(d/'solver.stdout').write_text(p.stdout+p.stderr)
 f06=deck.with_suffix('.F06').read_text(errors='replace');assert p.returncode==0 and '*ERROR' not in f06 and 'MYSTRAN terminated normally' in p.stdout,(key,variant)
 vv=vectors(read(deck.with_suffix('.OP2')));out=[]
 for r in groups[key]:
  response=r['response'];ids=ast.literal_eval(r['point_grids'])[response];idx=0 if response in ['IN','A','B'] else 1;axis=(1 if response=='IN' else 2 if response=='OUT' else 0 if response=='A' else 1);u=vv[idx][1];v=float(np.mean([u[j][axis] for j in ids]));target=float(r['independent']);old=float(r['mystran'])
  if variant=='EAS9v2':assert abs(v-old)<max(abs(old)*2e-6,1e-10),(key,response,v,old)
  out.append(dict(case=case,mesh=mesh,nodes=r['nodes'],elements=r['elements'],response=response,variant=variant,mystran=v,archived_eas9v2=old,nastran=float(r['nastran']),independent=target,sap2000=r['sap2000'],mystran_error_pct=100*abs(v-target)/abs(target),nastran_error_pct=r['nastran_error_pct'],source_mystran_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),source_nastran_sha256=hashlib.sha256(nx.read_bytes()).hexdigest(),reference='previous all-solid benchmark; original Nastran results unchanged'))
 return out
allrows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 futures={pool.submit(run,k,v):(k,v) for k in groups for v in ['EAS9v2','TIED','TIED_AUTO']}
 for future in concurrent.futures.as_completed(futures):
  result=future.result();allrows+=result;print(futures[future],'PASS',flush=True)
  with (OUT/'comparison.csv').open('w',newline='') as f:w=csv.DictWriter(f,allrows[0]);w.writeheader();w.writerows(allrows)
