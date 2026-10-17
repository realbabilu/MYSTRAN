"""Audit the supplied remeshed 5-001 without changing its geometry or loads.

Exports tensor comparisons and affine-displacement residuals. Independent
Python reference elements additionally check CPENTA6/15 and CHEXA20 recovery
using the displacements actually solved by MYSTRAN (not imposed ideal values).
"""
import csv
import struct
import sys
from pathlib import Path
import numpy as np
from pyNastran.bdf.bdf import BDF
from read_ogs import records, export, pynastran_volumes
from test_gpstress import EXPECTED, f06_rows
from generate_patch import PE, HE, PENT, HEX

ROOT=Path(__file__).resolve().parent
WORKSPACE=ROOT.parents[2]


def displacements(path):
    rec=records(path)
    for i,r in enumerate(rec):
        if r!=b'OUGV1   ': continue
        j=i+1
        while len(rec[j])!=584: j+=1
        j+=1
        while len(rec[j])==4: j+=1
        return {n//10:np.array(a[:3]) for n,t,*a in struct.iter_unpack('<ii6f',rec[j])}
    raise ValueError('no OUGV1 displacements')


def basis(points,family):
    r,s,t=np.asarray(points).T
    if family=='wedge':
        return np.array([np.ones_like(r),r,s,t,r*t,s*t]).T
    return np.array([np.ones_like(r),r,s,t,r*s,r*t,s*t,r*s*t]).T


def independent_reference(model,u):
    source=WORKSPACE/'python-solids/cpenta6_handoff_v2/CPENTA6_MYSTRAN_CODEX_HANDOFF/python/source'
    reference=WORKSPACE/'python-solids/handoff/python/reference_solids'
    sys.path.insert(0,str(reference))
    sys.path.insert(0,str(source))
    from core import Node
    from solid3d_cpenta6_ans_eas import CPENTA6ANSEAS
    from solid3d_chexa8_eas import isotropic_D
    sys.path.insert(0,str(reference))
    from solid3d_cpenta15 import b_matrix_penta15
    from solid3d_chexa20 import b_matrix_chexa20
    D=isotropic_D(1e6,.25)
    g=1/np.sqrt(3)
    result={}
    for volume in (3,4,6):
        accum={}
        for element in sorted(model.elements.values(),key=lambda e:e.eid):
            if element.pid!=volume: continue
            ids=element.node_ids
            xyz=np.array([model.nodes[n].xyz for n in ids])
            ue=np.array([u[n] for n in ids])
            if volume in (3,4):
                points=np.array([(r,s,t) for t in (-g,g) for r,s in ((1/6,1/6),(2/3,1/6),(1/6,2/3))])
                target=np.array([(0,0,-1),(1,0,-1),(0,1,-1),(0,0,1),(1,0,1),(0,1,1)],dtype=float)
                if volume==4:
                    target=np.vstack([target,[(target[a]+target[b])/2 for a,b in PE]])
                family='wedge'
                if volume==3:
                    obj=CPENTA6ANSEAS(element.eid,[Node(n,*x) for n,x in zip(ids,xyz)],1e6,.25,families='ABE')
                    samples=np.array([obj.strains_stresses_at(ue.ravel(),*p)[1] for p in points])
                else:
                    # Python reference orders top-face midsides before verticals.
                    order=list(range(9))+[12,13,14,9,10,11]
                    samples=np.array([D@b_matrix_penta15(xyz[order],*p)[2]@ue[order].ravel() for p in points])
            else:
                target=np.array([(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),
                                 (-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)],dtype=float)
                points=target*g
                target=np.vstack([target,[(target[a]+target[b])/2 for a,b in HE]])
                family='hex'
                order=list(range(12))+list(range(16,20))+list(range(12,16))
                samples=np.array([D@b_matrix_chexa20(xyz[order],*p)[2]@ue[order].ravel() for p in points])
            coefficients=np.linalg.solve(basis(points,family),samples)
            nodal=basis(target,family)@coefficients
            for nid,stress in zip(ids,nodal):
                accum.setdefault(nid,[]).append(stress)
        result[volume]={nid:np.mean(v,axis=0) for nid,v in accum.items()}
    return result


def audit(deck):
    deck=Path(deck)
    model=BDF(debug=False)
    model.read_bdf(str(deck),xref=False)
    tables=export(deck.with_suffix('.OP2'),ROOT/'attached_gpstress.csv')
    rows=f06_rows(deck.with_suffix('.F06').read_text(errors='replace'))
    u=displacements(deck.with_suffix('.OP2'))
    eps=np.full((3,3),.0005)
    np.fill_diagonal(eps,.001)
    reference=independent_reference(model,u)
    summary=[]
    with (ROOT/'attached_comparison.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['volume','grid','component','mystran','analytical','absolute_error','relative_error','python_reference'])
        for table in tables:
            volume=table['volume']
            a=np.array(table['values'])[:,:6]
            delta=np.array([u[n]-eps@model.nodes[n].xyz for n in table['grids']])
            # Independent patches have different origins; remove rigid translation.
            affine_error=np.max(np.ptp(delta,axis=0))
            max_python=0.
            for nid,stress in zip(table['grids'],a):
                np.testing.assert_allclose(stress,rows[(1,volume,nid)],rtol=2e-7,atol=1e-5)
                py=reference.get(volume,{}).get(nid)
                if py is not None:
                    max_python=max(max_python,float(np.max(abs(stress-py))))
                    # Reference uses single-precision OP2 displacement values;
                    # derivatives amplify that input rounding on small elements.
                    np.testing.assert_allclose(stress,py,rtol=1e-5,atol=2e-3)
                for c,value,expected in zip(('sxx','syy','szz','sxy','syz','szx'),stress,EXPECTED):
                    j=('sxx','syy','szz','sxy','syz','szx').index(c)
                    writer.writerow([volume,nid,c,value,expected,abs(value-expected),abs(value-expected)/abs(expected),'' if py is None else py[j]])
            max_rel=float(np.max(abs(a-EXPECTED)/abs(EXPECTED)))
            summary.append([volume,len(a),max_rel,affine_error,'' if volume not in reference else max_python,'PASS' if max_rel<=1e-5 else 'FAIL'])
    with (ROOT/'attached_summary.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['volume','grids','max_relative_stress_error','affine_displacement_range','max_python_stress_difference','patch_status'])
        writer.writerows(summary)
    model=pynastran_volumes(deck.with_suffix('.OP2'))
    assert len(model.grid_point_stresses_volume_direct)==8
    print('volume, grids, maximum relative stress error, affine residual, Python difference, patch status')
    for row in summary: print(*row,sep=', ')
    return summary


if __name__=='__main__':
    audit(sys.argv[1] if len(sys.argv)>1 else ROOT/'5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH.dat')
