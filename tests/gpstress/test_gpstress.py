"""End-to-end solid volume regression tests; run with Python unittest."""
import csv
import os
import re
import subprocess
import unittest
from pathlib import Path
import numpy as np
from generate_patch import ROOT, generate
from read_ogs import read_volumes, pynastran_volumes, records, export

EXE=ROOT.parents[1]/'Binaries'/'mystran.exe'
EXPECTED=np.array([2000.,2000.,2000.,400.,400.,400.])
COUNTS=[4,4,6,6,8,8,5,5]


def run(name,text):
    deck=ROOT/f'{name}.dat'
    deck.write_text(text)
    env=os.environ.copy()
    env['PATH']=str(EXE.parent)+os.pathsep+'C:/gcc/bin'+os.pathsep+env['PATH']
    result=subprocess.run([str(EXE),deck.name],cwd=ROOT,env=env,
                          capture_output=True,text=True,timeout=60)
    (ROOT/f'{name}.stdout').write_text(result.stdout+result.stderr)
    f06=deck.with_suffix('.F06').read_text(errors='replace')
    return deck,f06,result.stdout


def f06_rows(text):
    result={}
    subcase=volume=0
    active=False
    for line in text.splitlines():
        if 'G R I D   P O I N T   S T R E S S E S' in line:
            active=True
        if not active:
            continue
        if line.strip().startswith('SUBCASE ='):
            subcase=int(line.split('=')[1])
        if line.strip().startswith('VOLUME  ='):
            volume=int(line.split('=')[1])
        parts=line.split()
        if len(parts)==7 and parts[0].isdigit():
            try:
                values=[float(x) for x in parts[1:]]
            except ValueError:
                continue
            result[(subcase,volume,int(parts[0]))]=values
    return result


class VolumeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.templates={m:generate(m).read_text() for m in ('LEGACY','NEWSOLID')}

    def check_patch(self,mode):
        deck,f06,stdout=run(f'regular_{mode.lower()}',self.templates[mode])
        self.assertIn('MYSTRAN terminated normally',stdout)
        tables=export(deck.with_suffix('.OP2'),deck.with_suffix('.csv'))
        self.assertEqual([t['volume'] for t in tables],list(range(1,9)))
        frows=f06_rows(f06)
        for table,count in zip(tables,COUNTS):
            self.assertEqual(len(table['grids']),count)
            self.assertEqual(table['grids'],list(range(table['volume']*100000+1,table['volume']*100000+count+1)))
            values=np.array(table['values'])
            self.assertTrue(np.isfinite(values).all())
            np.testing.assert_allclose(values[:,:6],np.broadcast_to(EXPECTED,(count,6)),rtol=1e-5)
            np.testing.assert_allclose(values[:,6],-2000.,rtol=1e-5)
            np.testing.assert_allclose(values[:,7],1200.,rtol=1e-5)
            for grid,value in zip(table['grids'],values):
                np.testing.assert_allclose(value[:6],frows[(1,table['volume'],grid)],rtol=1e-6)
        model=pynastran_volumes(deck.with_suffix('.OP2'))
        actual=model.grid_point_stresses_volume_direct
        self.assertEqual(len(actual),8)
        for key,result in actual.items():
            table=next(t for t in tables if t['volume']==result.ogs_id)
            np.testing.assert_array_equal(result.node,table['grids'])
            np.testing.assert_allclose(result.data[0],table['values'])

    def test_principal_volume_payload(self):
        deck,_,stdout=run('principal_volume',self.templates['NEWSOLID'])
        self.assertIn('MYSTRAN terminated normally',stdout)
        direct=read_volumes(deck.with_suffix('.OP2'))
        principal=read_volumes(deck.with_suffix('.OP2'),28)
        self.assertEqual(len(principal),8)
        for d,p in zip(direct,principal):
            self.assertEqual(d['grids'],p['grids'])
            self.assertEqual(p['header'][:12],(13,28,p['volume'],1,1,0,0,0,1,15,1,3))
            for tensor,row in zip(d['values'],p['values']):
                xx,yy,zz,xy,yz,zx=tensor[:6]
                a=np.array([[xx,xy,zx],[xy,yy,yz],[zx,yz,zz]])
                eig=np.linalg.eigvalsh(a)[[2,0,1]]
                np.testing.assert_allclose(row[:3],eig,rtol=1e-5)
                axes=np.array(row[3:12]).reshape(3,3)
                for j in range(3):
                    if np.linalg.norm(axes[:,j])>0:
                        np.testing.assert_allclose(a@axes[:,j],row[j]*axes[:,j],rtol=1e-5,atol=.01)
                np.testing.assert_allclose(row[12:],tensor[6:],rtol=1e-6)
        text=self.templates['NEWSOLID'].replace('SYSTEM BASIC','SYSTEM BASIC DIRECT')
        deck,_,_=run('direct_only',text)
        self.assertEqual(len(read_volumes(deck.with_suffix('.OP2'),28)),0)
        text=self.templates['NEWSOLID'].replace('SYSTEM BASIC','SYSTEM BASIC PRINCIPAL')
        deck,_,_=run('principal_only',text)
        self.assertEqual(len(read_volumes(deck.with_suffix('.OP2'))),0)
        self.assertEqual(len(read_volumes(deck.with_suffix('.OP2'),28)),8)

    def test_patch_legacy(self): self.check_patch('LEGACY')
    def test_patch_newsolid(self): self.check_patch('NEWSOLID')

    def test_single_volume(self):
        text=self.templates['LEGACY']
        for v in range(2,9):
            text=re.sub(rf'^VOLUME {v},.*\n','',text,flags=re.M)
        deck,f06,stdout=run('single_volume',text)
        self.assertIn('MYSTRAN terminated normally',stdout)
        self.assertEqual(len(read_volumes(deck.with_suffix('.OP2'))),1)

    def test_rectangular_transform(self):
        text=self.templates['LEGACY'].replace('SYSTEM BASIC','SYSTEM CORD 42')
        # Local x=basic y, local y=-basic x, local z=basic z.
        text=text.replace('BEGIN BULK','BEGIN BULK\nCORD2R,42,,0.,0.,0.,0.,0.,1.\n,0.,1.,0.')
        deck,f06,stdout=run('rotated_volume',text)
        self.assertIn('MYSTRAN terminated normally',stdout)
        for table in read_volumes(deck.with_suffix('.OP2')):
            self.assertEqual(table['cid'],42)
            expected=[2000,2000,2000,-400,-400,400]
            np.testing.assert_allclose(np.array(table['values'])[:,:6],np.broadcast_to(expected,(len(table['grids']),6)),rtol=1e-5)

    def test_subcase_scope(self):
        text=self.templates['LEGACY'].replace('SPC = 1','SPC = 1',1)
        text=text.replace('GPSTRESS(PRINT,PLOT) = ALL\nOUTPUT(POST)',
            'SUBCASE 10\n GPSTRESS(PRINT,PLOT) = ALL\nSUBCASE 20\n GPSTRESS = NONE\nOUTPUT(POST)')
        deck,f06,stdout=run('subcase_scope',text)
        self.assertIn('MYSTRAN terminated normally',stdout)
        self.assertEqual({t['subcase'] for t in read_volumes(deck.with_suffix('.OP2'))},{10})

    def test_selected_domains(self):
        text=self.templates['LEGACY'].replace('GPSTRESS(PRINT,PLOT) = ALL','GPSTRESS(PRINT,PLOT) = 99')
        text=text.replace('OUTPUT(POST)','OUTPUT(POST)\nSET 99 = 1, 3 THRU 4')
        deck,f06,stdout=run('selected_domains',text)
        self.assertIn('MYSTRAN terminated normally',stdout)
        self.assertEqual([t['volume'] for t in read_volumes(deck.with_suffix('.OP2'))],[1,3,4])

    def test_many_volumes(self):
        extra='\n'.join(f'SET {1000+v} = 100001\nVOLUME {v}, SET {1000+v}, SYSTEM BASIC' for v in range(9,270))
        text=self.templates['LEGACY'].replace('BEGIN BULK',extra+'\nBEGIN BULK')
        deck,f06,stdout=run('many_volumes',text)
        self.assertIn('MYSTRAN terminated normally',stdout)
        self.assertEqual(len(read_volumes(deck.with_suffix('.OP2'))),269)

    def test_errors(self):
        cases=[('missing_element','SET 101 = 100001','SET 101 = 123456','MISSING ELEMENT ID'),
               ('missing_range','SET 101 = 100001','SET 101 = 100001 THRU 100003','MISSING IDS'),
               ('undefined_set','SET 101,','SET 999,','UNDEFINED OR EMPTY ELEMENT SET'),
               ('invalid_cid','SYSTEM BASIC','SYSTEM CORD 876','UNDEFINED OUTPUT COORDINATE'),
               ('duplicate_volume','VOLUME 2,','VOLUME 1,','DUPLICATE GPSTRESS VOLUME'),
               ('empty_volume','SET 101 = 100001','SET 101 = 100001 THRU 100001 EXCEPT 100001','EMPTY VOLUME'),
               ('unsupported_family','SET 101 = 100001','SET 101 = 999999','UNSUPPORTED ELEMENT FAMILY')]
        for name,old,new,message in cases:
            with self.subTest(name=name):
                _,f06,stdout=run(name,self.templates['LEGACY'].replace(old,new))
                self.assertIn(message,f06)
                self.assertNotIn('MYSTRAN terminated normally',stdout)

    def test_existing_output_unchanged(self):
        gp=self.templates['LEGACY']
        no=gp.replace('GPSTRESS(PRINT,PLOT) = ALL\n','')
        no=re.sub(r'OUTPUT\(POST\).*?BEGIN BULK','BEGIN BULK',no,flags=re.S)
        a,_,sa=run('control_no_gp',no)
        b,_,sb=run('control_with_gp',gp)
        self.assertIn('MYSTRAN terminated normally',sa)
        self.assertIn('MYSTRAN terminated normally',sb)
        def payloads(path):
            rec=records(path)
            out=[]
            in_ogs=False
            for item in rec:
                if len(item)==8 and re.fullmatch(b'[A-Z0-9 ]{8}',item):
                    in_ogs=item==b'OGS1    '
                if not in_ogs and len(item)>32:
                    out.append(item)
            return out
        self.assertEqual(payloads(a.with_suffix('.OP2')),payloads(b.with_suffix('.OP2')))

    def test_attached_basic_structure_and_existing_payloads(self):
        import hashlib
        fixture=ROOT/'5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH.dat'
        self.assertEqual(hashlib.sha256(fixture.read_bytes()).hexdigest(),
                         'b4efe70cd735d82bd037c6f93307e4117137ecf084a746d52203cdcddec04074')
        gp=fixture.read_text()
        no=gp.replace('  GPSTRESS(PRINT,PLOT) = ALL\n','')
        no=re.sub(r'OUTPUT\(POST\).*?BEGIN BULK','BEGIN BULK',no,flags=re.S)
        a,_,_=run('attached_no_gp',no)
        b,f06,_=run('attached_with_gp',gp)
        tables=read_volumes(b.with_suffix('.OP2'))
        self.assertEqual([len(t['grids']) for t in tables],[47,47,16,16,16,16,23,23])
        self.assertEqual([t['volume'] for t in tables],list(range(1,9)))
        self.assertTrue(all(t['cid']==0 for t in tables))
        def existing(path):
            out=[]; in_ogs=False
            for item in records(path):
                if len(item)==8 and re.fullmatch(b'[A-Z0-9 ]{8}',item):
                    in_ogs=item==b'OGS1    '
                if not in_ogs and len(item)>32: out.append(item)
            return out
        self.assertEqual(existing(a.with_suffix('.OP2')),existing(b.with_suffix('.OP2')))
        # Do not certify this deck as a patch pass: preserve its existing errors.
        self.assertIn('HAS JACOBIAN LESS THAN OR EQUAL TO ZERO',f06)


if __name__=='__main__': unittest.main(verbosity=2)
