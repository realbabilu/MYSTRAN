from pathlib import Path
import sys,os,subprocess,json,csv
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];WS=ROOT.parent
REF=P/'reference'
sys.path.insert(0,str(REF))
from chexa8_tied import CHEXA8TiedEAS,_perm_for
from chexa8_eas_family import GP,Bmat,Mnat,T0_matrix,shape
from core import Node
ENV=os.environ.copy();ENV['PATH']='C:/gcc/bin;C:/gcc/openblas32/bin;'+ENV['PATH'];ENV['OPENBLAS_NUM_THREADS']='1'
X0=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)
def probe(X,u,temp,point,form='TIED',axis='ZETA',nu=.3,baseline=False):
 inp=f'{form} {axis} 2000 {nu}\n'+'\n'.join(' '.join(map(str,row)) for row in X)+'\n'+' '.join(map(str,np.r_[u,temp,point]))+'\n'
 exe=P/('baseline_probe.exe' if baseline else 'probe.exe');out=P/'state.bin'
 r=subprocess.run([str(exe),str(out)],input=inp,cwd=ROOT,env=ENV,capture_output=True,text=True)
 assert r.returncode==0,(r.returncode,r.stderr)
 raw=np.fromfile(out);res=[];off=0
 for sh in [(24,24),(24,),(6,24),(6,),(6,),(24,24)]:
  n=np.prod(sh);res.append(raw[off:off+n].reshape(sh,order='F'));off+=n
 return res

def reference(X,u,temp,point,nu=.3):
 el=CHEXA8TiedEAS(1,[Node(i+1,*x) for i,x in enumerate(X)],2000,nu)
 J0,d0,_=el._jac(0,0,0);T=T0_matrix(J0);D=el._D
 t=J0[0].copy();n=J0[2]/np.linalg.norm(J0[2]);t-=n*(t@n);t/=np.linalg.norm(t);R=np.column_stack((t,np.cross(n,t),n))
 Px=el._proj(R,(0,2));Py=el._proj(R,(1,2))
 def at(p):
  _,det,g=el._jac(*p);B=Bmat(g);xi,eta,zeta=p
  Bt=B+Px@(el._Bat((0,eta,zeta))-B)+Py@(el._Bat((xi,0,zeta))-B)
  return det,g,B,Bt,(d0/det)*(T@Mnat(el.modes,*p)),shape(*p)
 pts=[at((x,y,z)) for x,_ in GP for y,_ in GP for z,_ in GP]
 mean=sum(det*(bt-b) for det,g,b,bt,m,n in pts)/sum(p[0] for p in pts)
 k=np.zeros((24,24));ku=np.zeros((24,9));ka=np.zeros((9,9));f=np.zeros(24);fa=np.zeros(9);alp=np.array([1,1,1,0,0,0])*1e-5
 for det,g,b,bt,m,n in pts:
  b=bt-mean;eth=alp*(n@temp);k+=det*b.T@D@b;ku+=det*b.T@D@m;ka+=det*m.T@D@m;f+=det*b.T@D@eth;fa+=det*m.T@D@eth
 sol=np.linalg.solve(ka,ku.T);eta=np.linalg.solve(ka,fa-ku.T@u);k-=ku@sol;f-=ku@np.linalg.solve(ka,fa)
 det,g,b,bt,m,n=at(point);b=bt-mean;beff=b-m@sol;e=b@u+m@eta;s=D@(e-alp*(n@temp));kg=np.zeros((24,24))
 for det,g,b,bt,m,n in pts:
  v=D@((bt-mean)@u+m@eta-alp*(n@temp));S=np.array([[v[0],v[3],v[5]],[v[3],v[1],v[4]],[v[5],v[4],v[2]]]);kg+=det*np.kron(g.T@S@g,np.eye(3))
 return [k,f,beff,e,s,kg]
def main():
 rng=np.random.default_rng(618);Q,_=np.linalg.qr(rng.normal(size=(3,3)));rows=[];u=rng.normal(size=24)*1e-4;temp=rng.uniform(0,50,8);point=np.array([.2,-.3,.4])
 shapes={'regular':X0*[1,.8,.5],'distorted':X0*[1,.8,.5]+rng.normal(size=(8,3))*.08,'thin':X0*[1,.8,.02]}
 shapes.update({k+'_rotated':x@Q.T for k,x in list(shapes.items())})
 for nu in [.3,.4999]:
  for name,X in shapes.items():
   old=probe(X,u,temp,point,'EAS9',nu=nu,baseline=True);new=probe(X,u,temp,point,'EAS9',nu=nu)
   assert all(np.array_equal(a,b) for a,b in zip(old,new)),('legacy changed',name)
   for axis,ax in [('XI',0),('ETA',1),('ZETA',2)]:
    perm=_perm_for(ax);ix=np.array([[3*i+d for d in range(3)] for i in perm]).ravel()
    expected=reference(X[perm],u[ix],temp[perm],point[[(ax+1)%3,(ax+2)%3,ax]],nu)
    # Natural modes are cyclically invariant; evaluate directly with cyclic physical projectors.
    got=probe(X,u,temp,point,axis=axis,nu=nu)
    mapped=[got[0][np.ix_(ix,ix)],got[1][ix],got[2][:,ix],got[3],got[4],got[5][np.ix_(ix,ix)]]
    errors=[float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30)) for a,b in zip(mapped,expected)]
    eig=np.linalg.eigvalsh(got[0]);zero=int(np.sum(abs(eig)<max(abs(eig))*1e-11));neg=int(np.sum(eig < -max(abs(eig))*1e-11))
    assert max(errors)<1e-8,(name,axis,nu,errors)
    assert zero==6 and neg==0,(name,axis,nu,eig)
    # Uniform stress patch and stress-free thermal affine field.
    ut=(X*1e-4).ravel();thermal=probe(X,ut,np.ones(8)*10,point,axis=axis,nu=nu)
    terr=float(np.linalg.norm(thermal[4])/max(np.linalg.norm(2000*1e-4),1e-30))
    assert terr<1e-7,(name,axis,terr)
    rows.append(dict(geometry=name,axis=axis,nu=nu,max_operator_relative_error=max(errors),zero_modes=zero,negative_modes=neg,thermal_stress_relative=terr,legacy_bitwise_equal=True))
 with (P/'operator_matrix.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)
 print('PASS',len(rows),'operator/axis/Poisson cases; legacy bitwise equal')
if __name__=='__main__':main()
