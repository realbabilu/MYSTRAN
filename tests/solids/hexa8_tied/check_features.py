from pathlib import Path
import json,re,subprocess,csv
import numpy as np
from scipy.linalg import eigh
import run_regression as r
import check_operators as o
P=r.P;OUT=P/'feature_checks';OUT.mkdir(exist_ok=True)
def run(name,text,expect_error=False):
 d=OUT/name;d.mkdir(exist_ok=True);(d/'model.dat').write_text(text)
 p=subprocess.run([str(r.MY),'model.dat'],cwd=d,env=r.ENV,capture_output=True,text=True,timeout=120);(d/'solver.stdout').write_text(p.stdout+p.stderr);f=(d/'model.F06').read_text(errors='replace')
 if not expect_error:assert p.returncode==0 and '*ERROR' not in f and 'MYSTRAN terminated normally' in p.stdout,(name,p.stdout[-500:])
 return f,d

def main():
 src=(r.OUT/'5001_patch_1/TIED/model.dat').read_text();base='\n'.join(l for l in src.splitlines() if not l.startswith(('PARAM,HEXATYP','PARAM,HEXATAX')))+'\n';rows=[]
 ff,default=run('default',base);assert 'HEXATYP=TIED' in ff
 _,explicit=run('explicit',base.replace('BEGIN BULK','BEGIN BULK\nPARAM,HEXATYP,TIED'))
 _,lower=run('lowercase',base.replace('BEGIN BULK','BEGIN BULK\nPARAM,HEXATYP,tied\nPARAM,HEXATAX,zeta'))
 assert (default/'model.OP2').read_bytes()==(explicit/'model.OP2').read_bytes()==(lower/'model.OP2').read_bytes()
 for key,value,msg in [('HEXATYP','RAW','PARAM HEXATYP MUST BE TIED OR EAS9'),('HEXATAX','MINIMUM','PARAM HEXATAX MUST BE XI, ETA, ZETA OR AUTO')]:
  ff,_=run('invalid_'+key,base.replace('BEGIN BULK',f'BEGIN BULK\nPARAM,{key},{value}'),True);assert msg in ff
 rows.append(dict(test='default_explicit_lowercase_and_invalid',status='PASS'))
 # Consistent and lumped mass through the complete solver, independent eigenvalues.
 X=o.X0;xyz={i+1:x for i,x in enumerate(X)};K=o.reference(X,np.zeros(24),np.zeros(8),np.zeros(3))[0]
 M=np.zeros((24,24))
 for x,_ in o.GP:
  for y,_ in o.GP:
   for z,_ in o.GP:
    n=o.shape(x,y,z);M+=2*np.kron(np.outer(n,n),np.eye(3))
 for coup in [0,1]:
  b=r.bulk(xyz,[list(range(1,9))],2000,.3);b=[('MAT1,1,2000.,,0.3,2.' if l.startswith('MAT1,') else l) for l in b]
  b += [f'PARAM,COUPMASS,{coup}','SPC1,1,123,1,2,3,4','EIGRL   3                       4                       10.     MAX     +E1','+E1             DGB']
  text='\n'.join(['SOL 103','CEND','SPC=1','METHOD=3','DISPLACEMENT(PLOT)=ALL','BEGIN BULK']+b+['ENDDATA'])+'\n'
  ff,_=run('modal_'+str(coup),text);tail=ff.split('R E A L   E I G E N V A L U E S')[-1];vals=np.array([float(v) for v in re.findall(r'^\s+\d+\s+\d+\s+([+-]?\d+\.\d+E[+-]\d+)',tail,re.M)])[:4]
  mass=M if coup else np.diag(M.sum(axis=1));expected=eigh(K[12:,12:],mass[12:,12:],eigvals_only=True)[:4];err=float(np.max(abs(vals/expected-1)));assert len(vals)==4 and err<1e-5,(coup,vals,expected)
  rows.append(dict(test='COUPMASS_'+str(coup),status='PASS',max_relative_eigenvalue_error=err))
 # Policy decision: near cube and severely skewed bricks fallback; clear thin bricks tie along their short direction.
 for name,X,axis in [('cube',o.X0,'ZETA'),('thin_xi',o.X0*[.02,1,1],'XI'),('thin_eta',o.X0*[1,.02,1],'ETA'),('thin_zeta',o.X0*[1,1,.02],'ZETA'),('skewed',o.X0@np.array([[1,0,0],[.9,.2,0],[0,0,.02]]),'ZETA')]:
  a=o.probe(X,np.zeros(24),np.zeros(8),np.zeros(3),axis='AUTO')[0];b=o.probe(X,np.zeros(24),np.zeros(8),np.zeros(3),axis=axis)[0];assert np.allclose(a,b,rtol=1e-10,atol=1e-8)
  rows.append(dict(test='AUTO_'+name,status='PASS',selected_axis=axis))
 r.writecsv(OUT/'feature_matrix.csv',rows);print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
