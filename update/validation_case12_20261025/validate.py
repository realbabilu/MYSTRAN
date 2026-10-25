from pathlib import Path
import json,sys,hashlib,concurrent.futures,csv
import numpy as np
P=Path(__file__).resolve().parent;W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');sys.path.insert(0,str(W))
from run_suite import run
from results import parse_f06,measures,kgg
from build_suite import fields,cards
forms=[r for r in json.loads((W/'manifest.json').read_text()) if r['group']=='2014'];jobs=[]
def selector(r,text):
 if r['formulation']=='CTRIAR_MITC3plusHB':text=text.replace('TRIARTYP,MITC3+HB','TRIARTYP,KMITC3+')
 out=[]
 for line in text.splitlines():
  f=fields(line)
  if f and f[0]=='SPC1' and 'THRU' not in f and len(f)>9:out.extend(','.join(f[:3]+f[i:i+6]) for i in range(3,len(f),6))
  else:out.append(line)
 return '\n'.join(out)+'\n'
def add(r,tag,version,text):
 f='CTRIAR_KMITC3plus' if r['formulation']=='CTRIAR_MITC3plusHB' else r['formulation'];d=P/'runs'/version/f/tag/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text(selector(r,text));jobs.append((tag,version,dict(r,formulation=f,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest())))
for r in forms:
 text=(W/r['deck']).read_text();text=text.replace('STRESS(CENTER) = ALL','STRESS(CENTER) = ALL\n    FORCE(CENTER) = ALL')
 for tag,temp,g in [('annulus',0,-100),('reverse',0,100),('zero',0,0),('double',0,-200),('mean_only',50,0),('combined',50,-100)]:
  lines=[]
  for line in text.splitlines():
   if line.startswith('TEMPP1,'):
    f=line.split(',');f[3]=str(float(temp));f[4]=str(float(g));line=','.join(f)
   lines.append(line)
  add(r,tag,'after','\n'.join(lines)+'\n')
 for label,temp,g in [('mean',100,0),('gradient',0,-100),('combined',100,-100),('reverse_gradient',0,100),('zero',0,0)]:
  tri=r['element'].startswith('CTRI');n=3 if tri else 4;nodes=['GRID,1,,0.,0.,0.','GRID,2,,2.,0.,0.','GRID,3,,0.,2.,0.'] if tri else ['GRID,1,,0.,0.,0.','GRID,2,,2.,0.,0.','GRID,3,,2.,2.,0.','GRID,4,,0.,2.,0.']
  t='SOL 101\nOUTPUT4 KGG,,,,//-1/21\nCEND\nECHO=NONE\nSPC=1\nSUBCASE 1\n TEMP(LOAD)=2\n DISPLACEMENT=ALL\n STRESS(CENTER)=ALL\n FORCE(CENTER)=ALL\nBEGIN BULK\n'
  t+=f"PARAM,{r['param']},{r['value']}\nMAT1,1,29000.,,0.3,,6.5-6,0.\nPSHELL,1,1,1.,1,,1\n"+'\n'.join(nodes)+'\n'+r['element']+',1,1,'+','.join(str(i) for i in range(1,n+1))+'\nSPC1,1,12345,1\nSPC1,1,123456,'+','.join(str(i) for i in range(2,n+1))+f'\nTEMPP1,2,1,{float(temp)},{float(g)}\nENDDATA\n'
  add(r,'restrained_'+label,'after',t)
  if label=='zero':add(r,'restrained_'+label,'before',t)
for r in json.loads((W/'manifest.json').read_text()):
 if r['group']=='2013' or r['case']=='prob_2_005_simple_mesh_06_cquadr_simo':
  for version in ['before','after']:add(r,'regression_'+r['case'],version,(W/r['deck']).read_text())
def task(j):
 tag,version,r=j;exe=P/('bin' if version=='after' else 'baseline')/'mystran.exe';s=run(r,exe,300,True);d=parse_f06(Path(r['deck']).with_suffix('.F06'),'MYSTRAN');d['stress']={sc:[v for v in vs if abs(v[1])>0] for sc,vs in d['stress'].items()}
 return dict(tag=tag,version=version,formulation=r['formulation'],status=s['status'],errors=s['errors'],deck=r['deck'],sha256=r['sha256'],displacements=d['displacements'],stress=d['stress'],forces=d['forces'],observations=measures(r,d) if tag in ['annulus','reverse','double','zero','mean_only','combined'] else [])
out=[]
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for v in pool.map(task,jobs):
   out.append(v)
   if len(out)%20==0:print(len(out),'/',len(jobs),'completed',sum(v['status']=='completed' for v in out),flush=True);(P/'results.json').write_text(json.dumps(out,indent=2))
 (P/'results.json').write_text(json.dumps(out,indent=2));print('FINAL',len(out),'completed',sum(v['status']=='completed' for v in out),flush=True)
