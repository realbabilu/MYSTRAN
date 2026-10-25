from pathlib import Path
import sys,json,hashlib,csv
import numpy as np
from scipy.sparse.linalg import spsolve
P=Path(__file__).resolve().parent;W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');kind=sys.argv[1]
ref=P/'DSQK_CovariantShear_Adaptation_20261025';sys.path[:0]=[str(ref),str(W)]
from core import Node
from assemble_reference import assemble
from run_suite import run
from results import kgg,parse_f06
if kind=='DSQK':from DSQK_ConsistentSpin_CovShear_ShellElement import DSQK_ConsistentSpin_CovShear_ShellElement as Element
else:from DSG3_ConsistentSpin_ShellElement import DSG3_ConsistentSpin_ShellElement as Element
xyz=np.array([[0.,0.,0.],[2.,0.,0.],[2.2,1.,0.],[.1,1.2,0.]]) if kind=='DSQK' else np.array([[0.,0.,0.],[2.,0.,0.],[.3,1.,0.]])
n=len(xyz);baseK=Element(1,[Node(j+1,*x) for j,x in enumerate(xyz)],6.825e7,.3,.04).k_global();rows=[]
a=.71;b=.43;R=np.array([[np.cos(a),-np.sin(a),0],[np.sin(a),np.cos(a),0],[0,0,1]])@np.array([[1,0,0],[0,np.cos(b),-np.sin(b)],[0,np.sin(b),np.cos(b)]])
for label,coords,perm in [('skew',xyz,list(range(n))),('rotated',xyz@R.T,list(range(n))),('cyclic',xyz,list(range(1,n))+[0]),('reversed',xyz,[0]+list(range(n-1,0,-1)))]:
 coords=np.array([[float(format(v,'.13g')) for v in row] for row in coords])
 e=Element(1,[Node(j+1,*coords[j]) for j in perm],6.825e7,.3,.04);K=e.k_global();order=[6*j+d for j in perm for d in range(6)];back=np.argsort(order);K=K[np.ix_(back,back)]
 card='CQUAD4' if n==4 else 'CTRIA3';param='QUAD4TYP' if n==4 else 'TRIA3TYP'
 lines=['SOL 101','OUTPUT4 KGG,,,,//-1/21','CEND','SPC=1','LOAD=1','DISPLACEMENT=ALL','STRESS=ALL','FORCE=ALL','BEGIN BULK',f'PARAM,{param},{kind}','PARAM,AUTOSPC,N','MAT1,1,6.825E7,,.3','PSHELL,1,1,.04,1,,1']
 lines += [f'GRID,{j+1},,'+','.join((format(v,'.13g') if ('.' in format(v,'.13g') or 'e' in format(v,'.13g')) else format(v,'.13g')+'.') for v in coords[j]) for j in range(n)]
 lines += [f'{card},1,1,'+','.join(str(j+1) for j in perm),'SPC1,1,123456,1,2','FORCE,1,3,0,1.,1.,0.,0.','ENDDATA']
 d=P/'parity'/kind/label/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text('\n'.join(lines)+'\n');r=dict(case=label,formulation=kind,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest());s=run(r,P/'bin/mystran.exe',120,True);assert s['status']=='completed',s
 NK=kgg(d.with_suffix('.OU1'));err=np.linalg.norm(K-NK)/np.linalg.norm(K);assert err<1e-12,(kind,label,err)
 eig=np.linalg.eigvalsh(K);tol=max(abs(eig))*1e-9;rig=[]
 for ax in range(6):
  u=np.zeros((n,6))
  if ax<3:u[:,ax]=1
  else:w=np.eye(3)[ax-3];u[:,:3]=np.cross(w,coords);u[:,3:]=w
  rig.append(float(np.linalg.norm(K@u.ravel())/np.linalg.norm(K)/np.linalg.norm(u)))
 assert max(rig)<1e-12 and np.sum(abs(eig)<tol)==6 and min(eig)>-tol
 target=np.kron(np.eye(2*n),R)@baseK@np.kron(np.eye(2*n),R).T if label=='rotated' else baseK
 inv=np.linalg.norm(K-target)/np.linalg.norm(K);assert inv<1e-12,(label,inv)
 rows.append(dict(element=kind,case=label,stiffness_relative_difference=float(err),invariance_error=float(inv),rigid_residual_max=max(rig),nullity=int(np.sum(abs(eig)<tol))))
# Assemble exact supplied benchmark cards, retaining all load sets and output frames.
form='CQUAD4_DSQK' if n==4 else 'CTRIA3_DSG3'
for case in ['prob_2_002_nx24_cquadr_simo','prob_2_004_nx24_cquadr_simo','prob_2_007_thick_cquadr_simo']:
 src=P/'runs/covariant'/kind/case/'model.dat';d=P/'parity'/kind/case/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);text=src.read_text()
 if '2_007' in case:text=text.replace('SOL 101','SOL 101\nOUTPUT4 KGG,,,,//-1/21')
 d.write_text(text);r=dict(case=case,formulation=kind,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest());s=run(r,P/'bin/mystran.exe',120,True);assert s['status']=='completed',s
 K,forces,free,ids,Q,els=assemble(d,Element);data=parse_f06(d.with_suffix('.F06'),'MYSTRAN')
 if '2_007' in case:
  err=np.linalg.norm(K.toarray()-kgg(d.with_suffix('.OU1')))/np.linalg.norm(K.toarray());assert err<1e-12;rows.append(dict(element=kind,case=case,stiffness_relative_difference=float(err)))
 for lc,disp in data['displacements'].items():
  if lc not in forces:continue
  u=np.zeros(K.shape[0]);u[free]=spsolve(K[free,:][:,free],forces[lc][free]);v=np.array([disp[j] for j in ids]).ravel();err=np.linalg.norm(u-v)/np.linalg.norm(u);assert err<3e-6,(kind,case,lc,err)
  rows.append(dict(element=kind,case=case,loadcase=lc,relative_displacement_error=float(err)))
(P/f'parity_{kind}.json').write_text(json.dumps(rows,indent=2));print(kind,len(rows),'parity checks passed',flush=True)
