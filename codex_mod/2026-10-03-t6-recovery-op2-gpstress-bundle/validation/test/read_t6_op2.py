"""Read native CTRIA6 OP2 tables with pyNastran and check against F06.

Usage: python test/read_t6_op2.py test/duel3a.OP2
"""
from pathlib import Path
import sys, json
import numpy as np
from pyNastran.op2.op2 import OP2
from cpylog import SimpleLogger
root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python'))
from compare_f06_vs_python import parse_mystran
path = Path(sys.argv[1]).resolve()
model = OP2(debug=False, log=SimpleLogger(level='warning'))
model.read_op2(str(path),build_dataframe=False)
F = parse_mystran(str(path.with_suffix('.F06')))
summary = {}
for table, is_force in ((model.op2_results.stress.ctria6_stress,False),
                        (model.op2_results.force.ctria6_force,True)):
    for key,obj in table.items():
        sc = obj.isubcase
        assert obj.element_type == 75
        errors=[]
        for i,(eid,gid) in enumerate(obj.element_node):
            if is_force:
                actual = F['elfrc'].get(sc,{}).get(int(eid),{}).get('CENTER' if gid==0 else int(gid))
                if actual is not None:
                    np.testing.assert_allclose(obj.data[0,i],actual,rtol=1e-4,atol=1e-10)
                    errors.append(float(np.max(np.abs(obj.data[0,i]-actual))))
            else:
                actual = F['elstr'].get(sc,{}).get(int(eid)) if gid==0 else F['elgrd'].get(sc,{}).get(int(eid),{}).get(int(gid))
                if actual is not None:
                    ref = actual[i%2]
                    np.testing.assert_allclose(obj.data[0,i,1:4],ref,rtol=1e-4,atol=1e-8)
                    errors.append(float(np.max(np.abs(obj.data[0,i,1:4]-ref))))
        record=dict(element_type=obj.element_type,num_wide=obj.num_wide,shape=list(obj.data.shape),
                    elements=np.unique(obj.element_node[:,0]).tolist(),f06_samples=len(errors),
                    max_f06_abs_error=max(errors,default=0))
        summary[f'{"force" if is_force else "stress"}_SC{sc}']=record
        if '--reference' not in sys.argv:
            assert errors, f'No F06 samples compared for subcase {sc}'
assert summary, 'Native CTRIA6 results missing'
if '--reference' not in sys.argv:
    for key,obj in model.grid_point_surface_stresses.items():
        if obj.ogs_id != 600:
            continue
        errors=[]
        for i,(gid,eid) in enumerate(obj.node_element):
            fiber=obj.location[i]
            if isinstance(fiber,bytes): fiber=fiber.decode()
            fiber=str(fiber).strip()
            ref=F['gp'].get(obj.isubcase,{}).get(int(gid),{}).get(fiber)
            assert ref is not None, f'Missing F06 GPSTRESS {gid}/{fiber}'
            np.testing.assert_allclose(obj.data[0,i,:3],ref,rtol=1e-4,atol=1e-8)
            errors.append(float(np.max(np.abs(obj.data[0,i,:3]-ref))))
        summary[f'GPSTRESS_600_SC{obj.isubcase}']=dict(shape=list(obj.data.shape),f06_samples=len(errors),
                                                   max_f06_abs_error=max(errors,default=0))
print('READ OP2 + CTRIA6 CHECKS PASS:',path)
print(model.get_op2_stats(short=True))
print(json.dumps(summary,indent=2))
path.with_suffix('.pynastran.json').write_text(json.dumps(summary,indent=2))
