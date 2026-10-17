"""Native OES/OSTR recovery regression on the unchanged distorted 5-001 mesh.

Run: python -m unittest discover -s tests/gpstress -p test_solid_element_recovery.py -v
Requires a built executable and the project's NumPy/pyNastran environment.
"""
from pathlib import Path
import json,os,re,subprocess,tempfile,unittest
import numpy as np
from read_ogs import read_volumes

REPO=Path(__file__).resolve().parents[2]
SOURCE=REPO/'update/solid_benchmark/5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH_FIXED_QUADRATIC_MYSTRAN_EAS9.dat'
EXE=REPO/'Binaries/mystran.exe'
ENV=os.environ.copy();ENV['PATH']=str(EXE.parent)+';C:/gcc/bin;'+ENV['PATH'];ENV['OPENBLAS_NUM_THREADS']='1'
COUNTS={1:168,2:168,3:14,4:14,5:7,6:7,7:42,8:42}

def read_op2(p):
    if not hasattr(np,'in1d'):np.in1d=lambda a,b,**kw:np.isin(np.asarray(a).ravel(),b,**kw)
    if not hasattr(np,'chararray'):np.chararray=np.char.chararray
    from pyNastran.op2.op2 import OP2
    from pyNastran.op2.op2_interface.msc_tables import MSC_ELEMENTS
    from cpylog import SimpleLogger
    MSC_ELEMENTS.setdefault(255,'CPYRAM')
    m=OP2(debug=False,log=SimpleLogger(level='error'),mode='msc')
    m.set_additional_result_tables_to_read({b'PVT0':False})
    m.set_results(['stress','strain','displacements'])
    m.read_op2(str(p),build_dataframe=False)
    return m

def tables(model,kind):
    result=[]
    for prefix in ['ctetra','cpenta','chexa','cpyram']:
        result.extend(getattr(getattr(model.op2_results,kind),prefix+'_'+kind).values())
    return result

class SolidRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='solid_recovery_',dir=REPO/'build')
        cls.work=Path(cls.tmp.name);cls.text=SOURCE.read_text();cls.cards=[]
        for raw in cls.text.splitlines():
            line=raw.split('$')[0].strip()
            if ',' not in line:continue
            f=[v.strip() for v in line.split(',')]
            if f[0]=='+':cls.cards[-1].extend(f[1:])
            else:cls.cards.append(f)
        cls.xyz={int(f[1]):np.array(list(map(float,f[3:6]))) for f in cls.cards if f[0]=='GRID'}
        cls.fixed={(int(n),int(d)-1) for f in cls.cards if f[0]=='SPC1' for n in f[3:] if n for d in f[2] if d in '123'}
        mu=4e5;lam=4e5;cls.D=np.diag([2*mu]*3+[mu]*3);cls.D[:3,:3]+=lam

    @classmethod
    def tearDownClass(cls):
        assert cls.work.resolve().is_relative_to((REPO/'build').resolve())
        cls.tmp.cleanup()

    def run_deck(self,name,text,allow_errors=False):
        d=self.work/name;d.mkdir(exist_ok=True);deck=d/'model.dat';deck.write_text(text)
        p=subprocess.run([str(EXE),deck.name],cwd=d,env=ENV,capture_output=True,text=True,timeout=120)
        (d/'solver.stdout').write_text(p.stdout+p.stderr)
        self.assertEqual(p.returncode,0,p.stdout[-1000:]+p.stderr[-1000:])
        f06=deck.with_suffix('.F06').read_text(errors='replace')
        if not allow_errors:self.assertFalse([l for l in f06.splitlines() if '*ERROR' in l][:12])
        return read_op2(deck.with_suffix('.OP2')),deck,f06

    def affine(self,g):
        def u(n):
            x=self.xyz[n].copy();x[0]-=2*(n//100000-1)
            return g@x
        lines=[l for l in self.text.splitlines() if not l.startswith('SPCD,')]
        i=lines.index('ENDDATA')
        def real(value):
            value=format(float(value),'.12g')
            return value if any(c in value.upper() for c in '.ED') else value+'.'
        lines[i:i]=[f'SPCD,1,{n},{d+1},{real(u(n)[d])}' for n,d in sorted(self.fixed)]
        return '\n'.join(lines)+'\n',u

    def check_field(self,model,eps):
        expected={'strain':eps,'stress':self.D@eps}
        for kind,target in expected.items():
            ids=set()
            for a in tables(model,kind):
                ids.update(map(int,a.element_node[:,0]))
                error=np.max(abs(a.data[0,:,:6]-target))/max(float(np.max(abs(target))),1.)
                self.assertLess(error,1e-7,(kind,error))
            self.assertEqual(len(ids),462,kind)
            for group,count in COUNTS.items():self.assertEqual(sum(n//100000==group for n in ids),count)

    def test_supplied_patch_full_stress_and_strain(self):
        m,p,_=self.run_deck('supplied',self.text)
        self.check_field(m,np.full(6,.001))
        for code in [27,28]:
            volumes=read_volumes(p.with_suffix('.OP2'),code)
            self.assertEqual(len(volumes),8);self.assertEqual(sum(len(v['grids']) for v in volumes),204)
        for a in tables(m,'stress'):
            # Every standard solid record has one center and the physical corners.
            count={'CTETRA':5,'CPENTA':7,'CHEXA':9,'CPYRAM':6}[a.element_name]
            for eid in np.unique(a.element_node[:,0]):
                rows=a.element_node[a.element_node[:,0]==eid,1]
                self.assertEqual(len(rows),count);self.assertEqual(int(rows[0]),0)
                self.assertEqual(len(set(map(int,rows))),count)

    def test_six_constant_strain_states_and_rotation(self):
        for i in range(7):
            with self.subTest(state=i):
                eps=np.zeros(6);g=np.zeros((3,3))
                if i<3:g[i,i]=eps[i]=.001
                elif i<6:
                    a,b=[(0,1),(1,2),(2,0)][i-3];g[a,b]=g[b,a]=.0005;eps[i]=.001
                else:g[0,1]=-.001;g[1,0]=.001
                text,u=self.affine(g);m,_,_=self.run_deck('state_'+str(i),text)
                self.check_field(m,eps)
                a=next(iter(m.displacements.values()))
                target=np.array([u(int(n)) for n in a.node_gridtype[:,0]])
                self.assertLess(np.max(abs(a.data[0,:,:3]-target)),1e-8)

    def test_valid_solid_batch_survives_another_recovery_failure(self):
        # LEGACY retains direct corner evaluation: 500002 fails, but the six
        # valid bricks must still be written and the diagnostic must name 500002.
        text=self.text.replace('PARAM,SOLIDTYP,NEWSOLID','PARAM,SOLIDTYP,LEGACY')
        m,_,f06=self.run_deck('legacy_bad_corner',text,allow_errors=True)
        ids={int(e) for a in tables(m,'stress') for e in a.element_node[:,0] if 500001<=int(e)<=500007}
        self.assertEqual(ids,{500001,500003,500004,500005,500006,500007})
        error_lines=[l for l in f06.splitlines() if '*ERROR  9201' in l]
        self.assertTrue(error_lines)
        self.assertTrue(any('500002' in l and 'HEXA8' in l for l in error_lines))
        self.assertFalse(any('800042' in l for l in error_lines))

    def test_uniform_thermal_expansion_has_zero_stress(self):
        text,_=self.affine(np.eye(3)*.001)
        text=text.replace('SUBCASE 1','TEMPERATURE(LOAD) = 2\nSUBCASE 1')
        text=text.replace('MAT1,1,1.0E6,,0.25,1.0','MAT1,1,1.0E6,,0.25,1.0,0.00001,0.')
        text=text.replace('ENDDATA','TEMPD,2,100.\nENDDATA')
        m,p,_=self.run_deck('thermal_expansion',text)
        for a in tables(m,'strain'):
            np.testing.assert_allclose(a.data[0,:,:6],np.broadcast_to([.001,.001,.001,0,0,0],a.data[0,:,:6].shape),atol=1e-10)
        for a in tables(m,'stress'):self.assertLess(np.max(abs(a.data[0,:,:6])),1e-5)
        for v in read_volumes(p.with_suffix('.OP2')):self.assertLess(np.max(abs(np.array(v['values'])[:,:6])),1e-5)

    def test_pyramid_nonconstant_strain_and_temperature(self):
        probe=self.work/'probe.exe'
        cmd=['C:/gcc/bin/gfortran.exe','-static','-fcheck=all','-ffree-line-length-none',
             '-I',str(REPO/'Binaries/mod'),str(Path(__file__).with_name('solid_recovery_probe.f90')),
             str(REPO/'build/CMakeFiles/mystran.dir/objects.a'),
             '@'+str(REPO/'build/CMakeFiles/mystran.dir/linkLibs.rsp'),'-o',str(probe)]
        subprocess.run(cmd,cwd=REPO,env=ENV,check=True,capture_output=True,text=True)
        xyz=np.array([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0],[0,0,1],
                      [0,-1,0],[1,0,0],[0,1,0],[-1,0,0],
                      [-.5,-.5,.5],[.5,-.5,.5],[.5,.5,.5],[-.5,.5,.5]],float)
        targets=np.vstack([[0,0,.25],xyz[:5]])
        for n in [5,13]:
            for field in ['mechanical','thermal']:
                with self.subTest(nodes=n,field=field):
                    x=xyz[:n];u=np.zeros((n,3));t=np.zeros(n);alpha=np.zeros(6)
                    if field=='mechanical':
                        # Quadratic pyramid reproduces u_x=x*y. The linear
                        # pyramid's nodal field is x*y/(1-z), whose collapsed
                        # gradients are b, a, a*b (apex output is the mean).
                        u[:,0]=x[:,0]*x[:,1]
                        eps=np.zeros((6,6));eps[:,0]=targets[:,1];eps[:,3]=targets[:,0]
                        if n==5:eps[:,5]=targets[:,0]*targets[:,1]
                        stress=eps@self.D.T
                    else:
                        t=20+5*x[:,0]+7*x[:,1]+9*x[:,2];alpha[:3]=1e-5
                        eps=np.zeros((6,6));stress=-(20+5*targets[:,0]+7*targets[:,1]+9*targets[:,2])[:,None]*(self.D@alpha)
                    lines=[('PYRAM5' if n==5 else 'PYRAM13')+' '+str(n)]
                    lines+=[' '.join(map(str,row)) for row in [*x,*self.D,u.ravel(),t,alpha]]
                    p=subprocess.run([str(probe)],input='\n'.join(lines)+'\n',cwd=self.work,env=ENV,capture_output=True,text=True)
                    self.assertEqual(p.returncode,0,p.stderr)
                    actual=np.loadtxt(p.stdout.splitlines())
                    np.testing.assert_allclose(actual[:,:6],eps,atol=1e-12)
                    np.testing.assert_allclose(actual[:,6:],stress,rtol=1e-10,atol=1e-8)

if __name__=='__main__':unittest.main()
