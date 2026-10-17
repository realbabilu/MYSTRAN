"""Compare actual Simcenter 2512 and MYSTRAN OGS1 by subcase/volume/GRID."""
import csv
from pathlib import Path
import numpy as np
from read_ogs import read_volumes, pynastran_volumes

ROOT=Path(__file__).resolve().parent
WORKSPACE=ROOT.parents[2]
NAME='5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH'
REFERENCE=WORKSPACE/'nastran_test/5-001_GPSTRESS_BASIC/nastran'/f'{NAME}.op2'
MYSTRAN=ROOT/f'{NAME}.OP2'
EXPECTED=np.array([2000.,2000.,2000.,400.,400.,400.])


def compare():
    reference=read_volumes(REFERENCE)
    actual=read_volumes(MYSTRAN)
    independent=pynastran_volumes(REFERENCE).grid_point_stresses_volume_direct
    assert len(reference)==len(actual)==len(independent)==8
    summary=[]
    with (ROOT/'simcenter_2512_comparison.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['subcase','volume','grid','component','simcenter_2512','mystran','analytical','mystran_minus_simcenter'])
        for ref in reference:
            my=next(t for t in actual if (t['subcase'],t['volume'])==(ref['subcase'],ref['volume']))
            indexed=dict(zip(my['grids'],my['values']))
            common=np.array([indexed[g] for g in ref['grids']])[:,:6]
            target=np.array(ref['values'])[:,:6]
            decoded=next(v for v in independent.values() if v.ogs_id==ref['volume'])
            np.testing.assert_array_equal(decoded.node,ref['grids'])
            np.testing.assert_allclose(decoded.data[0],ref['values'])
            for grid,sc,ms in zip(ref['grids'],target,common):
                for c,x,y,z in zip(('sxx','syy','szz','sxy','syz','szx'),sc,ms,EXPECTED):
                    writer.writerow([ref['subcase'],ref['volume'],grid,c,x,y,z,y-x])
            sc_error=float(np.max(abs(target-EXPECTED)/EXPECTED))
            my_error=float(np.max(abs(common-EXPECTED)/EXPECTED))
            delta=float(np.max(abs(common-target)/EXPECTED))
            summary.append([ref['volume'],len(ref['grids']),len(my['grids']),sc_error,my_error,delta,
                            'PASS' if sc_error<=1e-5 else 'FAIL'])
    with (ROOT/'simcenter_2512_summary.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['volume','simcenter_grids','mystran_grids','simcenter_analytical_maxrel',
                         'mystran_analytical_maxrel_at_common_grids','mystran_simcenter_maxrel','simcenter_patch'])
        writer.writerows(summary)
    for row in summary: print(*row,sep=', ')


if __name__=='__main__': compare()
