"""Generate eight independent affine solid patches with valid regular geometry.

All solid translations are prescribed. A separate spring supplies one free DOF
so the normal SOL 101 pipeline is exercised without a zero-size solve.
"""
from pathlib import Path
import numpy as np
from pyNastran.bdf.bdf import BDF

ROOT = Path(__file__).resolve().parent
HEX = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
TET = [HEX[i] for i in (0,1,3,4)]
PENT = [HEX[i] for i in (0,1,3,4,5,7)]
PYR = HEX[:4]+[(.5,.5,1)]
TE = [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
PE = [(0,1),(1,2),(2,0),(0,3),(1,4),(2,5),(3,4),(4,5),(5,3)]
HE = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,5),(2,6),(3,7),(4,5),(5,6),(6,7),(7,4)]
YE = [(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)]


def generate(mode='LEGACY', name=None):
    name = name or f'regular_{mode.lower()}'
    b = BDF(debug=False)
    b.sol = 101
    cc = ['SPC = 1', 'LOAD = 1', 'DISPLACEMENT(PLOT) = ALL',
          'STRESS(PLOT) = ALL', 'STRAIN(PLOT) = ALL',
          'GPSTRESS(PRINT,PLOT) = ALL', 'OUTPUT(POST)']
    for v in range(1,9):
        cc += [f'SET {100+v} = {v*100000+1}', f'VOLUME {v}, SET {100+v}, SYSTEM BASIC']
    from pyNastran.bdf.case_control_deck import CaseControlDeck
    # pyNastran BDF does not parse VOLUME; insert Case Control after writing bulk.
    b.case_control_deck = CaseControlDeck([])
    b.add_param('POST', -1)
    b.add_param('SOLIDTYP', mode)
    b.add_mat1(1, 1e6, None, .25)
    all_grids = []
    shapes = [(TET,[],b.add_ctetra),(TET,TE,b.add_ctetra),
              (PENT,[],b.add_cpenta),(PENT,PE,b.add_cpenta),
              (HEX,[],b.add_chexa),(HEX,HE,b.add_chexa),
              (PYR,[],b.add_cpyram),(PYR,YE,b.add_cpyram)]
    strain = np.full((3,3), .0005)
    np.fill_diagonal(strain, .001)
    for v,(corners,edges,add) in enumerate(shapes,1):
        xyz = [np.array(x,dtype=float) for x in corners]
        xyz += [(xyz[i]+xyz[j])/2 for i,j in edges]
        nids = [v*100000+i+1 for i in range(len(xyz))]
        b.add_psolid(v,1,integ=2 if v % 2 else 3)
        add(v*100000+1,v,nids)
        for nid,x in zip(nids,xyz):
            b.add_grid(nid,x)
            for c,u in enumerate(strain@x,1):
                if u != 0:
                    b.add_spcd(1,[nid],[str(c)],[float(u)])
        all_grids += nids
    b.add_spc1(1,'123',all_grids)
    b.add_grid(999999,[10.,0.,0.])
    b.add_grid(999998,[11.,0.,0.])
    b.add_spc1(1,'123456',[999998])
    b.add_celas2(999999,1.,[999999,999998],c1=1,c2=1)
    b.add_force(1,999999,1.,[1.,0.,0.])
    bulk = ROOT/f'{name}.dat'
    b.write_bdf(str(bulk),size=16)
    text = bulk.read_text()
    start = text.index('BEGIN BULK')
    bulk.write_text('SOL 101\nCEND\n'+'\n'.join(cc)+'\n'+text[start:]+'\nENDDATA\n')
    return bulk


if __name__ == '__main__':
    for mode in ('LEGACY','NEWSOLID'):
        print(generate(mode))
