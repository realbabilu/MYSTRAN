from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.sparse.linalg import spsolve
P=Path(__file__).resolve().parent;W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');sys.path[:0]=[str(P/'DSQK_CovariantShear_Adaptation_20261025'),str(W)]
from DSQK_ConsistentSpin_ShellElement import DSQK_ConsistentSpin_ShellElement as Old
from DSQK_ConsistentSpin_CovShear_ShellElement import DSQK_ConsistentSpin_CovShear_ShellElement as New
from run_suite import run
from results import parse_f06,kgg
from assemble_reference import assemble
from build_suite import fields,cards

def num(v):return format(float(v),'.8e')
def native(name,text,version):
 d=P/'geometry_checks'/version/name/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text(text);s=run(dict(case=name,formulation='CQUAD4_DSQK',deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest()),P/('baseline' if version=='old' else 'bin')/'mystran.exe',240,True);return d,s
out=[]
for label,case in [('roof','prob_2_006_n4_cquadr_simo'),('hemisphere','prob_2_007_thick_cquadr_simo')]:
 src=W/f'decks/{case}/CQUAD4_DSQK/model.dat';text=src.read_text();cs=cards(text);xyz={int(c[1]):np.array(list(map(float,c[3:6]))) for c in cs if c[0]=='GRID'};rng=np.random.default_rng(20261025);moved={}
 if label=='roof':
  angles={k:np.arctan2(v[0],v[2]) for k,v in xyz.items()};mx=max(angles.values());my=max(v[1] for v in xyz.values());dx=mx/4;dy=my/4
  for k,v in xyz.items():
   a=angles[k];y=v[1]
   if 1e-6<a<mx-1e-6 and 1e-6<y<my-1e-6:
    a+=.2*dx*rng.uniform(-1,1);y+=.2*dy*rng.uniform(-1,1);r=np.hypot(v[0],v[2]);moved[k]=[r*np.sin(a),y,r*np.cos(a)]
 else:
  ang={k:(np.arctan2(v[1],v[0]),np.arcsin(v[2]/np.linalg.norm(v))) for k,v in xyz.items()};latmax=max(v[1] for v in ang.values())
  for k,v in xyz.items():
   az,lat=ang[k]
   if 1e-6<az<np.pi/2-1e-6 and 1e-6<lat<latmax-1e-6:
    az+=.2*(np.pi/2/8)*rng.uniform(-1,1);lat+=.2*(latmax/8)*rng.uniform(-1,1);r=np.linalg.norm(v);moved[k]=[r*np.cos(lat)*np.cos(az),r*np.cos(lat)*np.sin(az),r*np.sin(lat)]
 lines=[]
 for line in text.splitlines():
  f=fields(line)
  if f and f[0]=='GRID' and int(f[1]) in moved:f[3:6]=[num(v) for v in moved[int(f[1])]];line=','.join(f)
  if f and f[0]=='SPC1' and len(f)>9 and 'THRU' not in f:lines.extend(','.join(f[:3]+f[i:i+6]) for i in range(3,len(f),6))
  else:lines.append(line)
 text='\n'.join(lines)+'\n';text=text.replace('SOL 101','SOL 101\nOUTPUT4 KGG,,,,//-1/21')
 for version,cls in [('old',Old),('covariant',New)]:
  d,s=native(label+'_distorted',text,version);assert s['status']=='completed',s;K,F,free,ids,Q,els=assemble(d,cls);NK=kgg(d.with_suffix('.OU1'));ke=float(np.linalg.norm(NK-K.toarray())/np.linalg.norm(NK));u=np.zeros(K.shape[0]);u[free]=spsolve(K[free,:][:,free],F[1][free]);data=parse_f06(d.with_suffix('.F06'),'MYSTRAN');v=np.array([data['displacements'][1][k] for k in ids]).ravel();ue=float(np.linalg.norm(u-v)/np.linalg.norm(u));assert ke<1e-12 and ue<3e-6,(label,version,ke,ue);out.append(dict(case=label+'_distorted',version=version,moved_nodes=len(moved),K_relative_error=ke,u_relative_error=ue,note='Parametric surface perturbation with original boundary nodes, output frames and discrete nodal loads retained; parity diagnostic, not a new continuum load-accuracy claim.'))
for label,pts in [('concave',[[0,0,0],[2,0,0],[.3,.1,0],[0,1,0]]),('bowtie',[[0,0,0],[1,1,0],[0,1,0],[1,0,0]]),('collapsed',[[0,0,0],[2,0,0],[2,1e-15,0],[0,1e-15,0]])]:
 text='SOL 101\nCEND\nSPC=1\nLOAD=1\nBEGIN BULK\nPARAM,QUAD4TYP,DSQK\nMAT1,1,1.E7,,.3\nPSHELL,1,1,.04,1,,1\n'+'\n'.join(f'GRID,{i+1},,'+','.join(num(v) for v in p) for i,p in enumerate(pts))+'\nCQUAD4,1,1,1,2,3,4\nSPC1,1,123456,1,2\nFORCE,1,3,0,1.,0.,0.,1.\nENDDATA\n';d,s=native(label,text,'covariant');assert s['status']=='solver_error',(label,s);out.append(dict(case=label,status=s['status'],errors=s['errors'],expected_rejection=True))
# Mixed family recovery: a DSQK and a MIN3 element, prescribed shear in both.
t='SOL 101\nCEND\nSPC=1\nLOAD=1\nDISPLACEMENT=ALL\nSTRESS(CENTER)=ALL\nFORCE(CENTER)=ALL\nBEGIN BULK\nPARAM,QUAD4TYP,DSQK\nPARAM,TRIA3TYP,MIN3\nMAT1,1,1.E7,,.3\nPSHELL,1,1,.04,1,,1\n'
pts=[[0,0,0],[2,0,0],[2,1,0],[0,1,0],[10,0,0],[12,0,0],[10,1,0]]
for i,p in enumerate(pts):t+=f'GRID,{i+1},,'+','.join(num(v) for v in p)+'\n'
t+='CQUAD4,1,1,1,2,3,4\nCTRIA3,2,1,5,6,7\nSPC1,1,12345,1\nSPC1,1,123456,2,3,4,5,6,7\n'
for i,p in enumerate(pts):
 if p[0]:t+=f'SPCD,1,{i+1},3,{num(.001*p[0])}\n'
t+='ENDDATA\n';data={}
for version in ['old','covariant']:
 d,s=native('mixed_DSQK_MIN3',t,version);assert s['status']=='completed',s;data[version]=parse_f06(d.with_suffix('.F06'),'MYSTRAN')
for kind in ['stress','forces']:
 a=np.array([r for r in data['old'][kind][1] if r[0]==2 and (kind=='forces' or abs(r[1])>0)]);b=np.array([r for r in data['covariant'][kind][1] if r[0]==2 and (kind=='forces' or abs(r[1])>0)]);assert np.array_equal(a,b),(kind,a,b)
out.append(dict(case='mixed_DSQK_MIN3',MIN3_stress_and_forces_unchanged=True))
(P/'geometry_checks.json').write_text(json.dumps(out,indent=2));print('Distorted curved parity, invalid mappings and mixed family checks passed')
