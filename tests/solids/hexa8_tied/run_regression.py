"""Native CSI/SAP model regression. Nastran is executed, never substituted by Python."""
from pathlib import Path
import sys,os,re,json,csv,subprocess,hashlib,time,concurrent.futures
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];WS=ROOT.parent
sys.path.insert(0,str(ROOT/'tests/solids'))
from cpylog import SimpleLogger
if not hasattr(np,'in1d'):np.in1d=lambda a,b,**kw:np.isin(np.asarray(a).ravel(),b,**kw)
if not hasattr(np,'chararray'):np.chararray=np.char.chararray
from pyNastran.op2.op2 import OP2
ENV=os.environ.copy();ENV['PATH']='C:/gcc/bin;C:/gcc/openblas32/bin;'+ENV['PATH'];ENV['OPENBLAS_NUM_THREADS']='1'
NX=Path('C:/Program Files/Siemens/Femap 2512/nastran/bin/nastran.exe');MY=ROOT/'Binaries/mystran.exe';OUT=P/'regression';OUT.mkdir(exist_ok=True)
def records(s,name):
 q=re.search(r'TABLE:\s+"'+re.escape(name)+r'"(.*?)(?=\nTABLE:|\Z)',s,re.S)
 return [dict(re.findall(r'(\w+)=("[^"]*"|\S+)',l)) for l in q.group(1).splitlines() if l.strip()] if q else []
manifest=[]
def f(v):
 t=format(float(v),'.12g')
 if len(t)>16:t=format(float(v),'.10g')
 return t if any(c in t.lower() for c in '.e') else t+'.'
def card(name,fields):
 vals=list(map(str,fields));return '\n'.join([name+','+','.join(vals[:8])]+['+,'+','.join(vals[i:i+8]) for i in range(8,len(vals),8)])
def source(path):
 s=path.read_text(errors='replace');rec=lambda n:records(s,n);num=lambda v:float(v.replace(',','.'))
 xyz={int(t['Joint']):np.array([num(t[k]) for k in ['GlobalX','GlobalY','GlobalZ']]) for t in rec('JOINT COORDINATES')}
 conn=[]
 for t in rec('CONNECTIVITY - SOLID'):
  v=[int(t[f'Joint{i}']) for i in range(1,9)];cc=[v[i] for i in [0,1,3,2,4,5,7,6]]
  x=np.array([xyz[j] for j in cc]);J=np.array([[-1,1,1,-1,-1,1,1,-1],[-1,-1,1,1,-1,-1,1,1],[-1,-1,-1,-1,1,1,1,1]])@x/8
  if np.linalg.det(J)<0:cc=[cc[i] for i in [0,3,2,1,4,7,6,5]]
  conn.append(cc)
 props={t['SolidProp']:t for t in rec('SOLID PROPERTY DEFINITIONS')};selected={t['SolidProp'] for t in rec('SOLID PROPERTY ASSIGNMENTS')};assert len(selected)==1,selected
 matname=props[next(iter(selected))]['Material'];mat=next(t for t in rec('MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES') if t['Material']==matname)
 E,nu,alpha=[num(mat.get(k,'0')) for k in ['E1','U12','A1']]
 return xyz,conn,rec,E,nu,alpha,num

def bulk(xyz,conn,E,nu,alpha=0):
 b=['PARAM,POST,-1','PARAM,AUTOSPC,N',f'MAT1,1,{f(E)},,{f(nu)},,{f(alpha)},0.','PSOLID,1,1,0,2,,1']
 b += [card('GRID',[j,'',*[f(v) for v in x]]) for j,x in xyz.items()]
 b += [card('CHEXA',[i,1,*cc]) for i,cc in enumerate(conn,1)]
 b += [card('SPC1',[1,456,*list(xyz)[i:i+6]]) for i in range(0,len(xyz),6)]
 return b

def register(name,header,b,meta):
 text='\n'.join(header+['BEGIN BULK']+b+['ENDDATA'])+'\n'
 manifest.append(dict(name=name,text=text,**meta))

def make_patch():
 sys.path.insert(0,str(P/'fixtures'));import make_5001 as p
 m=p.build('CHEXA8');xyz=m.xyz;conn=m.els;b=bulk(xyz,conn,1e6,.25)
 header=['SOL 101','CEND','SPC=1','DISPLACEMENT(PLOT)=ALL','STRESS(PRINT,PLOT)=ALL','STRAIN(PLOT)=ALL','GPSTRESS(PRINT,PLOT)=ALL']
 strains=np.eye(6)*.001;grads=[]
 for e in strains:grads.append(np.array([[e[0],e[3]/2,e[5]/2],[e[3]/2,e[1],e[4]/2],[e[5]/2,e[4]/2,e[2]]]))
 grads.append(np.array([[0,-.001,.002],[.001,0,-.003],[-.002,.003,0]]))
 for i,G in enumerate(grads,1):
  bb=list(b);hh=header+['SUBCASE 1','LOAD=1','OUTPUT(POST)','SET 101 = 1 THRU 7','VOLUME 1, SET 101, SYSTEM BASIC']
  for j in range(9,17):
   bb.append(f'SPC1,1,123,{j}')
   for d,v in enumerate(G@xyz[j],1):
    if v:bb.append(f'SPCD,1,{j},{d},{f(v)}')
  register(f'5001_patch_{i}',hh,bb,dict(case='5-001',source='nastran_test/5-001_8_Solid_Patch_Test/make_5001.py exact 7-brick distorted interior-grid patch',grads=[G.tolist()],xyz={j:x.tolist() for j,x in xyz.items()},scope='affine engineering strain '+str(i) if i<7 else 'rigid rotation'))

def make_force(name,path,case,coarsen=1):
 xyz,conn,rec,E,nu,alpha,num=source(path)
 if coarsen>1:
  n=int(round(np.sqrt(len(xyz)/2)))-1;stride=n+1;nn=n//coarsen;kept={1+k*(stride**2)+i*stride+j for k in range(2) for i in range(0,n+1,coarsen) for j in range(0,n+1,coarsen)}
  # Source coordinate grid, coarsened without moving any retained joint.
  conn=[]
  for i in range(0,n,coarsen):
   for j in range(0,n,coarsen):
    a=1+i*stride+j;c=coarsen;v=[a,a+c,a+stride**2+c,a+stride**2,a+c*stride,a+c*stride+c,a+c*stride+stride**2+c,a+c*stride+stride**2]
    x=np.array([xyz[k] for k in v]);J=np.array([[-1,1,1,-1,-1,1,1,-1],[-1,-1,1,1,-1,-1,1,1],[-1,-1,-1,-1,1,1,1,1]])@x/8
    if np.linalg.det(J)<0:v=[v[k] for k in [0,3,2,1,4,7,6,5]]
    conn.append(v)
  xyz={j:x for j,x in xyz.items() if j in kept}
 b=bulk(xyz,conn,E,nu);rest=rec('JOINT RESTRAINT ASSIGNMENTS')
 for t in rest:
  j=int(t['Joint']);d=''.join(str(i) for i in range(1,4) if t.get(f'U{i}')=='Yes')
  if j in xyz and d:b.append(f'SPC1,1,{d},{j}')
 loads=rec('JOINT LOADS - FORCE');patterns=sorted({t['LoadPat'] for t in loads});probes=[]
 for sid,pat in enumerate(patterns,1):
  ff=[]
  for t in loads:
   if t['LoadPat']!=pat:continue
   j=int(t['Joint']);assert j in xyz;v=np.array([num(t[f'F{i}']) for i in range(1,4)]);assert t['CoordSys']=='GLOBAL'
   b.append(card('FORCE',[sid,j,0,'1.',*[f(a) for a in v]]));ff.append((j,v))
  for d in range(3):
   js=[j for j,v in ff if v[d]]
   if js:probes.append(dict(sid=sid,pattern=pat,axis=d,joints=js))
 header=['SOL 101','CEND','SPC=1','DISPLACEMENT(PLOT)=ALL','STRESS(PLOT)=ALL']
 for i,pat in enumerate(patterns,1):header += [f'SUBCASE {i}',f'LOAD={i}']
 register(name,header,b,dict(case=case,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),nodes=len(xyz),elements=len(conn),coarsen=coarsen,probes=probes,scope='native retained coordinates, equivalent coarsened connectivity' if coarsen>1 else 'native geometry, material, supports and joint forces'))

def make_thermal():
 path=P/'fixtures/Example 5-011.s2k';xyz,conn,rec,E,nu,alpha,num=source(path)
 for variant in ['restrained_uniform','free_uniform','restrained_gradient']:
  b=bulk(xyz,conn,E,nu,alpha)
  if variant.startswith('restrained'):
   for t in rec('JOINT RESTRAINT ASSIGNMENTS'):
    d=''.join(str(i) for i in range(1,4) if t.get(f'U{i}')=='Yes')
    if d:b.append(f'SPC1,1,{d},{t["Joint"]}')
  else:
   # Remove exactly six rigid modes on three bottom-plane corners.
   a=next(j for j,x in xyz.items() if np.allclose(x,[0,0,0]));bnode=next(j for j,x in xyz.items() if np.allclose(x,[10,0,0]));c=next(j for j,x in xyz.items() if np.allclose(x,[0,10,0]));b += [f'SPC1,1,123,{a}',f'SPC1,1,23,{bnode}',f'SPC1,1,3,{c}']
  if variant.endswith('uniform'):b.append('TEMPD,2,100.')
  else:b += [f'TEMP,2,{j},{f(100*x[2])}' for j,x in xyz.items()]
  register('5011_'+variant,['SOL 101','CEND','SPC=1','TEMPERATURE(LOAD)=2','DISPLACEMENT(PLOT)=ALL','STRESS(PRINT,PLOT)=ALL'],b,dict(case='5-011',source=str(path),scope=variant,E=E,nu=nu,alpha=alpha))

def make_cylinder():
 for mesh in ['c','e']:
  path=P/f'fixtures/Example 5-008{mesh}.s2k';xyz,conn,rec,E,_,alpha,num=source(path)
  active=rec('ACTIVE DEGREES OF FREEDOM');assert active[0]['UZ']=='No',active
  for nu in ([.3,.49,.499,.4999] if mesh=='c' else [.49]):
   b=bulk(xyz,conn,E,nu);angles={j:0 for j in xyz}
   for t in rec('JOINT LOCAL AXES ASSIGNMENTS 1 - TYPICAL'):angles[int(t['Joint'])]=np.deg2rad(num(t['AngleA']));assert num(t['AngleB'])==num(t['AngleC'])==0
   for j in xyz:
    b.append(f'SPC1,1,3,{j}')
    if angles[j]==0:b.append(f'SPC1,1,2,{j}')
    else:b.append(f'MPC,3,{j},2,1.,{j},1,{f(-np.tan(angles[j]))}')
   pressure=rec('SOLID LOADS - SURFACE PRESSURE');assert len(pressure)==1 and pressure[0]['Face']=='4' and num(pressure[0]['Pressure'])==1
   inner=[j for j,x in xyz.items() if abs(np.linalg.norm(x[:2])-3)<.002];assert len(inner)==4
   bottom=sorted([j for j in inner if xyz[j][2]==0],key=lambda j:angles[j]);a,z=[xyz[j] for j in bottom];v=np.array([z[1]-a[1],a[0]-z[0],0.]);assert v[0]>0
   for j in inner:b.append(card('FORCE',[2,j,0,'1.',*[f(x/4) for x in v]]))
   register(f'5008_{mesh}_nu{str(nu).replace(".","p")}', ['SOL 101','CEND','SPC=1','MPC=3','LOAD=2','DISPLACEMENT(PLOT)=ALL','STRESS(PRINT,PLOT)=ALL'],b,dict(case='5-008',source=str(path),scope='CSI plane strain UZ inactive; tangential local restraints as MPC; exact consistent nodal integral of unit pressure on planar inner quad',nu=nu,E=E,inner=inner,xyz={j:x.tolist() for j,x in xyz.items()},pressure_resultant=v.tolist()))

def make_buckling():
 path=P/'fixtures/SAP5010_CHEXA8_SOL105.dat';text=path.read_text();h,b=text.split('BEGIN BULK');b=[l for l in b.splitlines() if l and l!='ENDDATA'];b=[('PSOLID,1,1,0,2,,1' if l.startswith('PSOLID,') else l) for l in b]
 register('5010_buckling',h.strip().splitlines(),b,dict(case='5-010',source=str(path),scope='matched source SOL105 one-kip preload, physical recovered stress KG and modes',buckling=True))

def large_fields(text):
 # Preserve all numeric tokens while avoiding NX 8-character free-field real truncation.
 lines=[];counter=0
 for l in text.splitlines():
  vals=l.split(',')
  if vals[0] not in ['GRID','MAT1','SPCD','FORCE','MPC','TEMP']:lines.append(l);continue
  prefix=vals[0]+'*';vals=vals[1:]
  while vals:
   counter+=1;chunk=vals[:4];vals=vals[4:];cont=f'*{counter:06d}' if vals else ''
   assert all(len(v)<=16 for v in chunk)
   lines.append(prefix.ljust(8)+''.join(v.rjust(16) for v in chunk+['']*(4-len(chunk)))+cont);prefix=cont
 return '\n'.join(lines)+'\n'

def run_one(job,solver):
 name=job['name'];folder=OUT/name/solver;folder.mkdir(parents=True,exist_ok=True);text=job['text'];exe=NX if solver=='nastran' else (P/'mystran_eas9_baseline.exe' if solver=='baseline' else MY)
 if solver!='nastran':
  form='EAS9' if solver in ['EAS9','baseline'] else 'TIED';text=text.replace('BEGIN BULK',f'BEGIN BULK\nPARAM,SOLIDTYP,NEWSOLID\nPARAM,HEXATYP,{form}\nPARAM,HEXATAX,'+('AUTO' if solver=='TIED_AUTO' else 'ZETA'))
  if solver=='baseline':text='\n'.join(l for l in text.splitlines() if not l.startswith(('PARAM,HEXATYP','PARAM,HEXATAX')))+'\n'
  if job.get('buckling'):
   bulkpart=text.split('BEGIN BULK',1)[1];text='SOL 105\nCEND\nSPC=101\nDISPLACEMENT(PLOT)=ALL\nSUBCASE 1\nLOAD=201\nSUBCASE 2\nSTATSUB=1\nMETHOD=301\nBEGIN BULK'+bulkpart
 else:text=large_fields(text)
 normalized=[]
 for line in text.splitlines():
  v=line.split(',');positions={'GRID':[3,4,5],'FORCE':[4,5,6,7]}.get(v[0],[])
  for i in positions:
   if i<len(v) and v[i]:v[i]=f(v[i])
  normalized.append(','.join(v))
 text='\n'.join(normalized)+'\n'
 deck=folder/'model.dat';deck.write_text(text);args=[str(exe),deck.name]
 if solver=='nastran':args+=['old=no','scr=yes','sdir='+str(folder)]
 started=time.time();r=subprocess.run(args,cwd=folder,env=ENV,capture_output=True,text=True,timeout=600);(folder/'solver.stdout').write_text(r.stdout+r.stderr)
 f06=folder/('model.f06' if solver=='nastran' else 'model.F06');s=f06.read_text(errors='replace') if f06.exists() else ''
 ok=r.returncode==0 and bool(s) and not re.search(r'(USER|SYSTEM)\s+FATAL|\*ERROR',s)
 if solver!='nastran':ok=ok and 'MYSTRAN terminated normally' in r.stdout
 return dict(name=name,case=job['case'],solver=solver,returncode=r.returncode,status='COMPLETED' if ok else 'SOLVER_FAILED',seconds=time.time()-started,executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())

def writecsv(path,rows):
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def main():
 make_patch();make_force('5004_native',P/'fixtures/Example 5-004.s2k','5-004')
 make_force('5007_8x8',P/'fixtures/Example 5-007a.s2k','5-007',2);make_force('5007_16x16',P/'fixtures/Example 5-007a.s2k','5-007');make_force('5007_24x24',P/'fixtures/Example 5-007b.s2k','5-007',2)
 make_thermal();make_cylinder();make_buckling()
 (OUT/'manifest.json').write_text(json.dumps([{k:v for k,v in j.items() if k!='text'} for j in manifest],indent=2))
 runs=[]
 nastran_only='--nastran-only' in sys.argv
 if nastran_only and (OUT/'runs.csv').exists():
  runs=[r for r in csv.DictReader((OUT/'runs.csv').open()) if r['solver']!='nastran']
 solvers=['nastran','EAS9','TIED','TIED_AUTO']
 if (P/'mystran_eas9_baseline.exe').exists():solvers.insert(1,'baseline')
 tasks=[(j,s) for j in manifest for s in (['nastran'] if nastran_only else solvers)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  futures={pool.submit(run_one,j,s):(j['name'],s) for j,s in tasks}
  for fut in concurrent.futures.as_completed(futures):
   try:r=fut.result()
   except Exception as e:r=dict(name=futures[fut][0],solver=futures[fut][1],status='EXCEPTION',detail=str(e))
   runs.append(r);print(r['name'],r['solver'],r['status'],flush=True);writecsv(OUT/'runs.csv',runs)
if __name__=='__main__':main()
