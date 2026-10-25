from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.sparse.linalg import spsolve
P=Path(__file__).resolve().parent;H=P/'DSQK_CovariantShear_Adaptation_20261025';W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');sys.path[:0]=[str(H),str(W)]
from core import Node
from DSQK_ConsistentSpin_ShellElement import DSQK_ConsistentSpin_ShellElement as Old
from DSQK_ConsistentSpin_CovShear_ShellElement import DSQK_ConsistentSpin_CovShear_ShellElement as New
from run_suite import run
from results import kgg,parse_f06
from assemble_reference import assemble
SHAPES={'rectangle':[[0,0,0],[2,0,0],[2,1,0],[0,1,0]],'parallelogram':[[0,0,0],[1,0,0],[1.5,1,0],[.5,1,0]],'trapezoid':[[0,0,0],[1.2,0,0],[1,.9,0],[.1,1.1,0]],'skew':[[0,0,0],[1,0,0],[1.6,1,0],[.6,1,0]],'warped':[[.1,.2,0],[1.2,.1,.3],[1.4,1.1,.2],[.3,1,.1]]}
def number(x):
 s=format(float(x),'.13g');s=s if '.' in s or 'e' in s else s+'.';return s if len(s)<=16 else format(float(x),'.8e')
def runmodel(name,text,version):
 d=P/'extended'/version/name/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text(text);s=run(dict(case=name,formulation='CQUAD4_DSQK',deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest()),P/('baseline' if version=='old' else 'bin')/'mystran.exe',240,True);return d,s

def deck(x,h=.04,conn=None,field=None):
 conn=conn or [0,1,2,3];s='SOL 101\nOUTPUT4 KGG,,,,//-1/21\nCEND\nSPC=1\nLOAD=1\nDISPLACEMENT=ALL\nSTRESS(CENTER)=ALL\nFORCE(CENTER)=ALL\nBEGIN BULK\nPARAM,QUAD4TYP,DSQK\nPARAM,AUTOSPC,N\nMAT1,1,1.E7,,.3\nPSHELL,1,1,'+number(h)+',1,,1\n';s+='\n'.join('GRID,'+str(i+1)+',,'+','.join(number(v) for v in row) for i,row in enumerate(x))+'\nCQUAD4,1,1,'+','.join(str(i+1) for i in conn)+'\n'
 if field is None:return s+'SPC1,1,123456,1,2\nFORCE,1,3,0,1.,1.,0.,1.\nENDDATA\n'
 s+='SPC1,1,12345,1\nSPC1,1,123456,2,3,4\n'
 for i,(a,b,c) in enumerate(x):
  u=np.zeros(6)
  if field=='membrane':u[0]=.001*a+.00015*b;u[1]=-.0002*b+.00015*a
  elif field=='bending':u[2]=.0005*a*a;u[4]=-.001*a
  else:u[2]=.001*a
  for j,v in enumerate(u):
   if v:s+=f'SPCD,1,{i+1},{j+1},{number(v)}\n'
 return s+'ENDDATA\n'
rows=[]
for name,pts in SHAPES.items():
 x=np.array(pts,dtype=float);Ks={}
 for version,cls in [('old',Old),('covariant',New)]:
  d,st=runmodel(name,deck(x),version);assert st['status']=='completed',st;e=cls(1,[Node(i+1,*q) for i,q in enumerate(x)],1e7,.3,.04);K=e.k_global();NK=kgg(d.with_suffix('.OU1'));err=float(np.linalg.norm(K-NK)/np.linalg.norm(K));assert err<1e-12,(name,version,err);Ks[version]=K
  residual=[]
  for j in range(6):
   u=np.zeros((4,6))
   if j<3:u[:,j]=1
   else:w=np.eye(3)[j-3];u[:,:3]=np.cross(w,x);u[:,3:]=w
   residual.append(float(np.linalg.norm(K@u.ravel())/np.linalg.norm(K)/np.linalg.norm(u)))
  if name!='warped':assert max(residual)<1e-12
  rows.append(dict(case=name,version=version,K_relative_difference=err,rigid_max=max(residual),status=st['status']))
 if name=='warped':
  vals=[]
  for j in range(3):
   w=np.eye(3)[j];u=np.c_[np.cross(w,x),np.tile(w,(4,1))].ravel();vals.append(float(np.linalg.norm((Ks['covariant']-Ks['old'])@u)/np.linalg.norm(Ks['old'])/np.linalg.norm(u)))
  rows.append(dict(case='warped_inherited_limitation',delta_K_rigid_max=max(vals),note='Both kernels violate physical rigid rotation on this warped facet; this is not certified by nullity six.'))
 if name=='warped':continue
 for field in ['membrane','bending','shear']:
  for version,cls in [('old',Old),('covariant',New)]:
   d,st=runmodel(name+'_'+field,deck(x,field=field),version);assert st['status']=='completed',st;data=parse_f06(d.with_suffix('.F06'),'MYSTRAN');e=cls(1,[Node(i+1,*q) for i,q in enumerate(x)],1e7,.3,.04);e.k_global();u=np.array([data['displacements'][1][i] for i in range(1,5)]).ravel();ul=e._build_T_5to6()@u;gp=[(a,b) for a in [-1/np.sqrt(3),1/np.sqrt(3)] for b in [-1/np.sqrt(3),1/np.sqrt(3)]];bm=sum((e._Bm(a,b) for a,b in gp))/4;bb=sum((e._Bb(a,b) for a,b in gp))/4;bs=sum((e._Bs(a,b) for a,b in gp))/4;C=1e7/(1-.3**2)*np.array([[1,.3,0],[.3,1,0],[0,0,.35]]);errs=[]
   for eid,z,sx,sy,tau in data['stress'][1]:
    if abs(z)==0:continue
    sig=C@(bm@ul-z*bb@ul);T=np.column_stack([e.t1,e.t2]);tensor=T@np.array([[sig[0],sig[2]],[sig[2],sig[1]]])@T.T;expected=np.array([tensor[0,0],tensor[1,1],tensor[0,1]]);errs.append(float(np.linalg.norm(np.array([sx,sy,tau])-expected)/max(1,np.linalg.norm(expected))))
   shear=5/6*1e7/2/1.3*.04*(bs@ul);f=data['forces'][1][0];qerr=float(np.linalg.norm(np.array(f[7:9])-shear)/max(1,np.linalg.norm(shear)));assert max(errs)<3e-6 and (version=='old' or qerr<3e-6),(name,field,version,errs,qerr,f,shear)
   if version=='covariant' and field=='bending':assert np.linalg.norm(shear)<1e-9
   rows.append(dict(case=name+'_'+field,version=version,stress_relative_error=max(errs),shear_relative_error=qerr,shear_norm=float(np.linalg.norm(shear))))
(P/'extended_matrix_recovery.json').write_text(json.dumps(rows,indent=2));print('Matrix and prescribed-field recovery checks passed; warped inherited limitation recorded',flush=True)
# Matched random distortion, exact same nodal tributary load and hard simply supported BCs as the supplied Python study.
plates=[]
for h in [.1,.001]:
 for dist in [0,.3]:
  for n in [4,8,16]:
   def nid(i,j):return i*(n+1)+j+1
   rng=np.random.default_rng(3);step=.5/n;pts={}
   for i in range(n+1):
    for j in range(n+1):
     a,b=step*i,step*j
     if dist and 0<i<n and 0<j<n:a+=dist*step*rng.uniform(-1,1);b+=dist*step*rng.uniform(-1,1)
     pts[nid(i,j)]=np.array([float(format(a,'.13g')),float(format(b,'.13g')),0.])
   t='SOL 101\nCEND\nSPC=1\nLOAD=1\nDISPLACEMENT=ALL\nBEGIN BULK\nPARAM,QUAD4TYP,DSQK\nMAT1,1,1.E7,,.3\nPSHELL,1,1,'+number(h)+',1,,1\n';t+='\n'.join(f'GRID,{k},,'+','.join(number(v) for v in v) for k,v in pts.items())+'\n';F={k:0. for k in pts};eid=0;fixed={k:set() for k in pts}
   for i in range(n):
    for j in range(n):
     ns=[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)];p=np.array([pts[k] for k in ns]);A=.5*abs(np.cross(p[2]-p[0],p[3]-p[1])[2]);eid+=1;t+=f'CQUAD4,{eid},1,'+','.join(map(str,ns))+'\n'
     for k in ns:F[k]+=A/4
   for k in range(n+1):
    for node,ds in [(nid(n,k),[3,4]),(nid(k,n),[3,5]),(nid(0,k),[1,5]),(nid(k,0),[2,4])]:fixed[node].update(ds)
   fixed[nid(0,0)].update([1,2,6])
   for k,ds in fixed.items():
    if ds:t+=f'SPC1,1,'+''.join(map(str,sorted(ds)))+f',{k}\n'
   for k,f in F.items():t+=f'FORCE,1,{k},0,{number(f)},0.,0.,1.\n'
   t+='ENDDATA\n'
   for version,cls in [('old',Old),('covariant',New)]:
    name=f'plate_h{h}_d{dist}_n{n}';d,s=runmodel(name,t,version);assert s['status']=='completed',s;K,forces,free,ids,Q,els=assemble(d,cls);u=np.zeros(K.shape[0]);u[free]=spsolve(K[free,:][:,free],forces[1][free]);data=parse_f06(d.with_suffix('.F06'),'MYSTRAN');nu=data['displacements'][1][1][2];py=u[2];err=abs(nu-py)/abs(py);assert err<3e-6,(name,version,err);D=1e7*h**3/(12*(1-.3**2));ref=.00406/D+.0737/(5/6*1e7/2/1.3*h);plates.append(dict(version=version,h=h,distortion=dist,n=n,native=nu,python=py,ratio=abs(nu)/ref,response_relative_difference=err));print(name,version,plates[-1]['ratio'],flush=True)
(P/'distorted_plates.json').write_text(json.dumps(plates,indent=2));print('24 matched native plate checks passed',flush=True)
