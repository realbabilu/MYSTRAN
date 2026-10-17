"""Identical FIXED_QUADRATIC DAT patch, Python first, then Fortran kernels."""
import csv
import os
import subprocess
import shutil
import sys
from pathlib import Path
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]
WORK=SOURCE.parent
REF=WORK/'python-solids/handoff/python/reference_solids'
sys.path.insert(0,str(REF))
from core import Node
from solid3d_cpenta15 import CPENTA15,shape_penta15,dshape_penta15,TRI_GAUSS7,LINE_GAUSS3
from solid3d_chexa20 import CHEXA20,shape_chexa20,dshape_chexa20,GAUSS3
from solid3d_cpenta15 import b_matrix_penta15
from solid3d_chexa20 import b_matrix_chexa20
from solid3d_chexa8_eas import isotropic_D

DECK=WORK/'nastran_test/5-001_GPSTRESS_BASIC/quadratic_audit/fixed_quadratic_normalized.dat'
D=isotropic_D(1e6,.25)
G0=np.array([[2,1,1],[1,2,1],[1,1,2]],float)/2000


def states():
    out=[('manual',G0)]
    for name,i,j in [('exx',0,0),('eyy',1,1),('ezz',2,2),('gxy',0,1),('gyz',1,2),('gxz',0,2)]:
        g=np.zeros((3,3));g[i,j]=.001 if i==j else .0005;g[j,i]=g[i,j]
        out.append((name,g))
    out.append(('rotation',np.array([[0,-.001,0],[.001,0,0],[0,0,0]])))
    return out


def write_csv(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def load_case(pid):
    model=BDF(debug=False);model.read_bdf(str(DECK),xref=False)
    els=[e for e in model.elements.values() if e.pid==pid]
    ids=sorted({n for e in els for n in e.node_ids});index={n:i for i,n in enumerate(ids)}
    xyz=np.array([model.nodes[n].xyz for n in ids]);origin=xyz.min(axis=0)
    constraints={}
    for card in model.spcs[1]:
        for nid in card.nodes:
            if nid not in index:continue
            for c in str(card.components):
                if c in '123':constraints[3*index[nid]+int(c)-1]=0.
    for card in model.loads[1]:
        if card.type!='SPCD':continue
        for nid,c,value in zip(card.nodes,card.components,card.enforced):
            if nid in index:constraints[3*index[nid]+int(c)-1]=float(value)
    fixed=np.array(sorted(constraints));free=np.setdiff1d(np.arange(3*len(ids)),fixed)
    boundary={i for i,x in enumerate(xyz-origin) if np.any(np.isclose(x,0,atol=1e-12)) or np.any(np.isclose(x,1,atol=1e-12))}
    assert set(fixed)=={3*i+j for i in boundary for j in range(3)},'BC differs from full exterior exact displacement'
    manual=(xyz-origin)@G0.T
    np.testing.assert_allclose(np.array([constraints[i] for i in fixed]),manual.ravel()[fixed],atol=1e-15,rtol=0)
    cls=CPENTA15 if pid==4 else CHEXA20
    objects=[]
    for e in els:
        nodes=[Node(n,*model.nodes[n].xyz) for n in e.node_ids]
        obj=cls(e.eid,nodes,1e6,.25)
        dofs=np.array([3*index[n]+j for n in e.node_ids for j in range(3)])
        objects.append((e,obj,dofs))
    return model,ids,xyz,origin,fixed,free,objects


def quadrature(pid):
    if pid==4:return [(np.array([r,s,t]),wt*wz) for r,s,wt in TRI_GAUSS7 for t,wz in LINE_GAUSS3]
    return [(np.array([r,s,t]),wr*ws*wt) for r,wr in GAUSS3 for s,ws in GAUSS3 for t,wt in GAUSS3]


def solve_and_check(K,objects,xyz,origin,ids,fixed,free,pid,phase):
    rows=[];nodal=[];stress_rows=[]
    kff=K[np.ix_(free,free)]
    eig=np.linalg.eigvalsh(kff)
    assert eig[0]>0,(pid,phase,eig[0])
    for name,g in states():
        exact=((xyz-origin)@g.T).ravel();u=exact.copy()
        u[free]=np.linalg.solve(kff,-K[np.ix_(free,fixed)]@exact[fixed])
        scale=max(np.max(abs(exact)),1e-15)
        uerr=np.max(abs(u-exact))/scale
        eps=np.array([g[0,0],g[1,1],g[2,2],g[0,1]+g[1,0],g[1,2]+g[2,1],g[0,2]+g[2,0]])
        expected=D@eps;sscale=max(np.max(abs(expected)),1e6*np.max(abs(g)))
        serr=0.
        for e,obj,dofs in objects:
            for q,w in quadrature(pid):
                sig=obj.strains_stresses_at(u[dofs],*q)[1]
                serr=max(serr,float(np.max(abs(sig-expected))/sscale))
                stress_rows.append(dict(phase=phase,type=obj.__class__.__name__,state=name,eid=e.eid,r=q[0],s=q[1],t=q[2],
                                        **dict(zip(['sxx','syy','szz','sxy','syz','szx'],sig))))
        j8=ids.index(pid*100000+8)*3
        rows.append(dict(phase=phase,type=objects[0][1].__class__.__name__,state=name,nodes=len(ids),elements=len(objects),
                         free_dofs=len(free),u_rel=uerr,stress_rel=serr,symmetry=np.linalg.norm(K-K.T)/np.linalg.norm(K),
                         min_constrained_eigenvalue=eig[0],joint8_ux=u[j8],joint8_uy=u[j8+1],joint8_uz=u[j8+2]))
        for nid,actual,wanted in zip(ids,u.reshape(-1,3),exact.reshape(-1,3)):
            nodal.append(dict(phase=phase,type=objects[0][1].__class__.__name__,state=name,grid=nid,
                              ux=actual[0],uy=actual[1],uz=actual[2],max_abs_error=max(abs(actual-wanted))))
    return rows,nodal,stress_rows


def python_audit():
    rows=[];nodal=[];stress=[];shape_rows=[]
    for pid in (4,6):
        model,ids,xyz,origin,fixed,free,objects=load_case(pid)
        K=np.zeros((3*len(ids),3*len(ids)))
        for e,obj,dofs in objects:K[np.ix_(dofs,dofs)]+=obj.k_local()
        np.save(ROOT/f'python_K_pid{pid}.npy',K)
        r,n,s=solve_and_check(K,objects,xyz,origin,ids,fixed,free,pid,'Python');rows+=r;nodal+=n;stress+=s
        shape,derivative=(shape_penta15,dshape_penta15) if pid==4 else (shape_chexa20,dshape_chexa20)
        for e,obj,dofs in objects:
            if pid==4:
                corners=np.array([(0,0,-1),(1,0,-1),(0,1,-1),(0,0,1),(1,0,1),(0,1,1)],float)
                edges=[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]
            else:
                corners=np.array([(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)],float)
                edges=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
            natural=np.vstack([corners,[(corners[a]+corners[b])/2 for a,b in edges]])
            kronecker=np.max(abs(np.array([shape(*q) for q in natural])-np.eye(len(natural))))
            midside=max(np.linalg.norm(obj._xyz[len(corners)+k]-(obj._xyz[a]+obj._xyz[b])/2) for k,(a,b) in enumerate(edges))
            assert kronecker<1e-12 and midside<1e-12
            partition=gradient=deriv_error=0.;min_det=float('inf')
            for q,w in quadrature(pid):
                N=shape(*q);dn=derivative(*q);h=2e-5
                numerical=np.array([(shape(*(q-2*h*np.eye(3)[j]))-8*shape(*(q-h*np.eye(3)[j]))+
                                     8*shape(*(q+h*np.eye(3)[j]))-shape(*(q+2*h*np.eye(3)[j])))/(12*h) for j in range(3)])
                partition=max(partition,abs(sum(N)-1));gradient=max(gradient,np.max(abs(dn.sum(axis=1))))
                deriv_error=max(deriv_error,np.max(abs(dn-numerical)))
                min_det=min(min_det,np.linalg.det(dn@obj._xyz))
            shape_rows.append(dict(type=obj.__class__.__name__,eid=e.eid,partition_error=partition,
                                   derivative_sum_error=gradient,numerical_derivative_error=deriv_error,min_detj=min_det,
                                   quadrature_weight_sum=sum(w for q,w in quadrature(pid)),
                                   kronecker_error=kronecker,midside_coordinate_error=midside))
        print(pid,'Python worst u/stress',max(x['u_rel'] for x in r),max(x['stress_rel'] for x in r),flush=True)
    write_csv('python_patch.csv',rows);write_csv('python_nodal.csv',nodal);write_csv('python_gauss_stress.csv',stress)
    write_csv('python_shape_checks.csv',shape_rows)
    assert max(max(r['u_rel'],r['stress_rel']) for r in rows)<1e-7


def fortran_audit(tag='before'):
    env=os.environ.copy();env['PATH']='C:/gcc/bin;'+str(SOURCE/'Binaries')+os.pathsep+env['PATH']
    summary=[];patch=[]
    for pid in (4,6):
        model,ids,xyz,origin,fixed,free,objects=load_case(pid)
        for reduced in ('N','Y'):
            K=np.zeros((3*len(ids),3*len(ids)))
            for e,obj,dofs in objects:
                n=len(e.node_ids);nd=3*n
                perm=list(range(9))+[12,13,14,9,10,11] if pid==4 else list(range(12))+list(range(16,20))+list(range(12,16))
                inverse=np.argsort((np.array(perm)[:,None]*3+np.arange(3)).ravel())
                data=f'{"PENTA15" if pid==4 else "HEXA20"} NEWSOLID {n} {reduced}\n'
                data+='\n'.join(' '.join(format(v,'.17g') for v in row) for row in [*obj._xyz[perm],*D])+'\n'
                result=subprocess.run([str(ROOT/'quadratic_probe.exe'),'stiffness'],input=data,text=True,
                                      capture_output=True,cwd=ROOT,env=env,timeout=30)
                assert result.returncode==0,(e.eid,result.stderr,result.stdout[-1000:])
                raw=np.fromstring(result.stdout,sep=' ')[:nd*nd].reshape(nd,nd)
                ke=raw[np.ix_(inverse,inverse)]
                full=obj.k_local()
                py=full if reduced=='N' else selective_reference(obj,pid)
                np.savez(ROOT/f'{tag}_pid{pid}_eid{e.eid}_reduced{reduced}.npz',python=py,fortran=ke)
                summary.append(dict(stage=tag,type=obj.__class__.__name__,eid=e.eid,reduced=reduced,
                                    matrix_relative=np.linalg.norm(ke-py)/np.linalg.norm(py),
                                    full_integration_difference=np.linalg.norm(ke-full)/np.linalg.norm(full),
                                    symmetry=np.linalg.norm(ke-ke.T)/np.linalg.norm(ke),
                                    null_modes=int(sum(abs(np.linalg.eigvalsh(ke))<np.linalg.norm(ke,2)*1e-9))))
                K[np.ix_(dofs,dofs)]+=ke
            r,_,_=solve_and_check(K,objects,xyz,origin,ids,fixed,free,pid,f'Fortran_{tag}_reduced{reduced}')
            patch+=r
            print(pid,reduced,'matrix worst',max(x['matrix_relative'] for x in summary if x['type']==obj.__class__.__name__ and x['reduced']==reduced),
                  'patch worst',max(x['u_rel'] for x in r),max(x['stress_rel'] for x in r),flush=True)
    write_csv(f'fortran_matrices_{tag}.csv',summary);write_csv(f'fortran_patch_{tag}.csv',patch)
    if tag!='before':
        assert max(x['matrix_relative'] for x in summary)<1e-7
        assert all(x['null_modes']==6 for x in summary)
        assert max(max(x['u_rel'],x['stress_rel']) for x in patch)<1e-7


def selective_reference(obj,pid):
    normal=D.copy();normal[3:,:]=0;normal[:,3:]=0;shear=D-normal
    fn=b_matrix_penta15 if pid==4 else b_matrix_chexa20
    k=np.zeros((3*len(obj.nodes),)*2)
    for q,w in quadrature(pid):
        det,_,b=fn(obj._xyz,*q);k+=b.T@normal@b*det*w
    g=1/np.sqrt(3)
    low=[((r,s,t),1/6) for r,s in [(1/6,1/6),(2/3,1/6),(1/6,2/3)] for t in [-g,g]] if pid==4 else [((r,s,t),1) for r in [-g,g] for s in [-g,g] for t in [-g,g]]
    for q,w in low:
        det,_,b=fn(obj._xyz,*q);k+=b.T@shear@b*det*w
    return k


def baseline_matrices_from_cache():
    rows=[]
    for pid in (4,6):
        objects=load_case(pid)[-1]
        for reduced in ('N','Y'):
            for e,obj,dofs in objects:
                ke=np.load(ROOT/f'before_pid{pid}_eid{e.eid}_reduced{reduced}.npz')['fortran']
                full=obj.k_local();ref=full if reduced=='N' else selective_reference(obj,pid)
                rows.append(dict(stage='before',type=obj.__class__.__name__,eid=e.eid,reduced=reduced,
                                 matrix_relative=np.linalg.norm(ke-ref)/np.linalg.norm(ref),
                                 full_integration_difference=np.linalg.norm(ke-full)/np.linalg.norm(full),
                                 symmetry=np.linalg.norm(ke-ke.T)/np.linalg.norm(ke)))
    write_csv('fortran_matrices_before.csv',rows)


def probe_input(obj,pid,reduced='Y'):
    perm=list(range(9))+[12,13,14,9,10,11] if pid==4 else list(range(12))+list(range(16,20))+list(range(12,16))
    data=f'{"PENTA15" if pid==4 else "HEXA20"} NEWSOLID {len(perm)} {reduced}\n'
    data+='\n'.join(' '.join(format(v,'.17g') for v in row) for row in [*obj._xyz[perm],*D])+'\n'
    return data,perm


def geometry_audit():
    env=os.environ.copy();env['PATH']='C:/gcc/bin;'+str(SOURCE/'Binaries')+os.pathsep+env['PATH']
    rows=[];cache={};ordering=[];boundary=[]
    for pid in (4,6):
        model,ids,xyz,origin,fixed,free,objects=load_case(pid)
        shape,derivative=(shape_penta15,dshape_penta15) if pid==4 else (shape_chexa20,dshape_chexa20)
        for e,obj,dofs in objects:
            data,perm=probe_input(obj,pid)
            inv=np.argsort(perm);idof=np.argsort((np.array(perm)[:,None]*3+np.arange(3)).ravel());n=len(perm)
            for fi,py in enumerate(perm):
                for axis in range(3):
                    ordering.append(dict(type=obj.__class__.__name__,eid=e.eid,fortran_node=fi+1,python_node=py+1,
                                         grid=e.node_ids[py],axis=axis+1,fortran_translation_dof=3*fi+axis+1,
                                         fortran_six_slot_dof=6*fi+axis+1,python_dof=3*py+axis+1))
            Bs=[]
            for q,w in quadrature(pid):
                result=subprocess.run([str(ROOT/'quadratic_probe.exe'),'shape'],input=data+' '.join(map(str,q))+'\n',
                                      text=True,capture_output=True,env=env,cwd=ROOT,check=True)
                a=np.fromstring(result.stdout,sep=' ')
                N=a[:n][inv];dn=a[n:4*n].reshape(3,n)[:,inv];det=a[4*n]
                b=a[4*n+1:].reshape(6,3*n)[:,idof]
                expected=shape(*q);dexpected=derivative(*q)
                rows.append(dict(type=obj.__class__.__name__,eid=e.eid,r=q[0],s=q[1],t=q[2],
                                 partition_error=abs(sum(N)-1),gradient_sum=np.max(abs(dn.sum(axis=1))),
                                 shape_python_difference=np.max(abs(N-expected)),derivative_python_difference=np.max(abs(dn-dexpected)),
                                 determinant=det,determinant_python_difference=abs(det-np.linalg.det(dexpected@obj._xyz))))
                Bs.append(b)
            cache[e.eid]=np.array(Bs)
        nc=6 if pid==4 else 8
        corner_ids={n for e,obj,dofs in objects for n in e.node_ids[:nc]}
        for i,nid in enumerate(ids):
            boundary.append(dict(type=objects[0][1].__class__.__name__,grid=nid,x=xyz[i,0],y=xyz[i,1],z=xyz[i,2],
                                 midside=nid not in corner_ids,constrained_translations=sum(3*i+j in fixed for j in range(3))))
    write_csv('fortran_shape_checks.csv',rows)
    write_csv('node_dof_ordering.csv',ordering);write_csv('boundary_constraints.csv',boundary)
    np.savez(ROOT/'fortran_gauss_B.npz',**{str(k):v for k,v in cache.items()})
    assert max(r['shape_python_difference'] for r in rows)<1e-12
    assert max(r['derivative_python_difference'] for r in rows)<1e-12
    assert min(r['determinant'] for r in rows)>0
    result=subprocess.run([str(ROOT/'quadratic_probe.exe'),'weights'],input=data,text=True,
                          capture_output=True,env=env,cwd=ROOT,check=True)
    a=np.fromstring(result.stdout,sep=' ');line=a[:6].reshape(3,2);tria=a[6:].reshape(7,3)
    write_csv('quadrature_weights.csv',
              [dict(rule='line3',point=i+1,r=x,s=0,weight=w) for i,(x,w) in enumerate(line)]+
              [dict(rule='triangle7',point=i+1,r=x,s=y,weight=w) for i,(x,y,w) in enumerate(tria)])
    print('Fortran shapes/Jacobians checked:',len(rows),'line weight sum',sum(line[:,1]),'triangle weight sum',sum(tria[:,2]),flush=True)


def solver_audit():
    sys.path.insert(0,str(SOURCE/'tests/gpstress'))
    from read_ogs import records,read_volumes
    env=os.environ.copy();env['PATH']='C:/gcc/bin;'+str(SOURCE/'Binaries')+os.pathsep+env['PATH']
    cache=np.load(ROOT/'fortran_gauss_B.npz')
    summaries=[];nodes=[];stress=[]
    for pid in (4,6):
        model,ids,xyz,origin,fixed,free,objects=load_case(pid)
        perm=probe_input(objects[0][1],pid)[1]
        for name,g in states():
            bdf=BDF(debug=False);bdf.sol=101;bdf.case_control_deck=CaseControlDeck([])
            bdf.add_param('POST',-1);bdf.add_param('SOLIDTYP','NEWSOLID')
            bdf.add_mat1(1,1e6,None,.25);bdf.add_psolid(pid,1,integ=3)
            for nid,x in zip(ids,xyz):bdf.add_grid(nid,x)
            for e,obj,dofs in objects:
                ns=[e.node_ids[p] for p in perm]
                (bdf.add_cpenta if pid==4 else bdf.add_chexa)(e.eid,pid,ns)
            boundary=sorted({ids[i//3] for i in fixed})
            bdf.add_spc1(1,'123',boundary)
            exact=((xyz-origin)@g.T).ravel()
            for dof in fixed:
                if exact[dof]!=0:bdf.add_spcd(1,[ids[dof//3]],[str(dof%3+1)],[float(exact[dof])])
            filename=ROOT/f'pid{pid}_{name}.dat'
            bdf.write_bdf(str(filename),size=16,enddata=True)
            text=filename.read_text();bulk=text[text.index('BEGIN BULK'):]
            cc=['SOL 101','CEND','SPC = 1','LOAD = 1','DISPLACEMENT(PLOT) = ALL','STRESS(PLOT) = ALL',
                'GPSTRESS(PRINT,PLOT) = ALL','OUTPUT(POST)',f'SET 101 = {objects[0][0].eid} THRU {objects[-1][0].eid}',
                f'VOLUME {pid}, SET 101, SYSTEM BASIC']
            filename.write_text('\n'.join(cc)+'\n'+bulk)
            result=subprocess.run([str(ROOT/'mystran_audit.exe'),filename.name],cwd=ROOT,env=env,capture_output=True,text=True,timeout=60)
            filename.with_suffix('.stdout').write_text(result.stdout+result.stderr)
            f06=filename.with_suffix('.F06').read_text()
            assert 'MYSTRAN terminated normally' in result.stdout,(pid,name,result.stdout[-500:])
            assert '*ERROR' not in f06 and '*FATAL' not in f06,(pid,name,'solver diagnostic')
            rec=records(filename.with_suffix('.L5A'))
            raw=rec[-6*len(ids):]
            assert all(len(x)==8 for x in raw)
            import struct
            u=np.array([struct.unpack('<d',x)[0] for x in raw]).reshape(-1,6)[:,:3].ravel()
            err=np.max(abs(u-exact))/max(np.max(abs(exact)),1e-15)
            eps=np.array([g[0,0],g[1,1],g[2,2],g[0,1]+g[1,0],g[1,2]+g[2,1],g[0,2]+g[2,0]])
            expected=D@eps;scale=max(np.max(abs(expected)),1e6*np.max(abs(g)))
            serr=0.
            for e,obj,dofs in objects:
                for q,B in zip(quadrature(pid),cache[str(e.eid)]):
                    sigma=D@B@u[dofs];python=obj.strains_stresses_at(u[dofs],*q[0])[1]
                    serr=max(serr,float(np.max(abs(sigma-expected))/scale))
                    stress.append(dict(type=obj.__class__.__name__,state=name,eid=e.eid,r=q[0][0],s=q[0][1],t=q[0][2],
                                       python_difference=np.max(abs(sigma-python)),
                                       **dict(zip(['sxx','syy','szz','sxy','syz','szx'],sigma))))
            table=read_volumes(filename.with_suffix('.OP2'))[0]
            gp_error=np.max(abs(np.array(table['values'])[:,:6]-expected))/scale
            j8=ids.index(pid*100000+8)*3
            summaries.append(dict(type=objects[0][1].__class__.__name__,state=name,u_rel=err,stress_rel=serr,gpstress_rel=gp_error,
                                  joint8_ux=u[j8],joint8_uy=u[j8+1],joint8_uz=u[j8+2],fatal_errors=0))
            for nid,actual,wanted in zip(ids,u.reshape(-1,3),exact.reshape(-1,3)):
                nodes.append(dict(type=objects[0][1].__class__.__name__,state=name,grid=nid,ux=actual[0],uy=actual[1],uz=actual[2],
                                  displacement_error=max(abs(actual-wanted))))
            print(pid,name,'solver',err,serr,'GPSTRESS',gp_error,flush=True)
    write_csv('mystran_solver_patch.csv',summaries);write_csv('mystran_solver_nodal.csv',nodes);write_csv('mystran_solver_gauss.csv',stress)
    assert max(max(r['u_rel'],r['stress_rel'],r['gpstress_rel']) for r in summaries)<1e-7


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='baseline-cache':baseline_matrices_from_cache()
    elif len(sys.argv)>1 and sys.argv[1]=='geometry':geometry_audit()
    elif len(sys.argv)>1 and sys.argv[1]=='solver':solver_audit()
    elif len(sys.argv)>1:fortran_audit(sys.argv[1])
    else:python_audit()
