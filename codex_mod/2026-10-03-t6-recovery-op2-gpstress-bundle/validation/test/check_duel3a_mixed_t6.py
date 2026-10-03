"""Use the T6 comparator on the T6 subset of the mixed duel3a deck."""
from pathlib import Path
import sys, contextlib, json
import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python'))
import compare_f06_vs_python as c
import print_pernode_stress_2001_v3 as v3
from dump_python_debug import resolve_reference
D = c.parse_mystran(str(root/'test/duel3a.F06'))
_,entry = resolve_reference(v3,'MITC6')
cls, xy, conn, bnd, x0, y0, tag = entry
v3.FLIP_SHEAR_IF_NEG_E3 = False
metrics = {}
for sc,mode in ((1,'membrane'),(2,'bending')):
    fem,elems = v3.build_and_solve(cls,xy,conn,bnd,x0,y0,mode)
    disp_errors=[]; stress_errors=[]; moment_errors=[]; missing=[]
    for nid in xy:
        actual = D['disp'].get(sc,{}).get(nid)
        if actual is None:
            missing.append(['displacement',nid]); continue
        ref = fem.u[fem.nodes[nid].dofs]
        disp_errors.append(float(np.max(np.abs(actual-ref))/max(np.max(np.abs(ref)),1e-8)))
    for eid,e in elems.items():
        p = v3.collect_element(e,fem,conn[eid],mode,tag)[0]
        if sc==1:
            actual = D['elstr'].get(sc,{}).get(eid)
            if actual is None: missing.append(['CENTER stress',eid]); continue
            ref = np.array(v3.fibers(p['sig_loc'],p['mom_loc']))
            stress_errors.append(float(np.max(np.abs(actual-ref))/max(np.max(np.abs(ref)),1e-8)))
        else:
            actual = D['elfrc'].get(sc,{}).get(eid,{}).get('CENTER')
            if actual is None: missing.append(['CENTER moment',eid]); continue
            ref = p['mom_loc']
            moment_errors.append(float(np.max(np.abs(actual[3:6]-ref))/max(np.max(np.abs(ref)),1e-8)))
    metrics[mode] = dict(nodes=len(disp_errors),centers=len(stress_errors),moments=len(moment_errors),
                        max_displacement_rel=max(disp_errors,default=0),max_stress_rel=max(stress_errors,default=0),
                        max_moment_rel=max(moment_errors,default=0),missing=missing)
with (root/'test/duel3a_run_latest.compare.txt').open('w',encoding='utf-8') as out, contextlib.redirect_stdout(out):
    print('Mixed deck: compare T6 subset only (element IDs 51–60, grid IDs 301–325).')
    print(json.dumps(metrics,indent=2))
    c.report_python('MITC6',D,'duel3a')
print(json.dumps(metrics,indent=2))
