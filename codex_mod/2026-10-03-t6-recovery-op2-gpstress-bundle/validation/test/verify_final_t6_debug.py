"""Verify final class identity, element labels and nodal shape interpolation."""
from pathlib import Path
import sys, subprocess
import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python'))
import dump_python_debug as dump
import print_pernode_stress_2001_v3 as v3
out = root/'test/t6_upgrade'
expected = {
    'SIMOT6':'Simo1993_Tri6_ShellElement_v2',
    'MITC6':'MITC6_Tri_v4',
    'MH6T':'MacNeal_MH6T_Tri_v3',
    'REZAIEE':'Rezaiee2017_Tri6_v3',
}
for name, module in expected.items():
    canonical, entry = dump.resolve_reference(v3,name)
    assert canonical == name and entry[0].__module__ == module
    for i,(r,s) in enumerate(v3.T6_NODE_RS):
        np.testing.assert_allclose(v3.t6_N(r,s),np.eye(6)[i],atol=1e-14)
    subprocess.run([sys.executable,str(root/'python/dump_python_debug.py'),name,'--elem','51','--attrs'],cwd=out,check=True)
    for mode in ('membrane','bending'):
        text = (out/f'debug_{name}_{mode}.txt').read_text()
        assert f'# reference: {module}.' in text and 'E?' not in text and '/BS gagal:' not in text
        vals = dump.load_dump(out/f'debug_{name}_{mode}.txt')
        for pt in ('CTR','C1','C2','C3','M12','M23','M31'):
            for matrix in ('BM','BB','BS'):
                assert f'E51/{pt}/{matrix}/r01' in vals
    print(name,'final reference and both dumps PASS')
for alias, target in {'SimoT6':'SIMOT6','MITC6v4':'MITC6','MHT6v3':'MH6T','Rezaiee_v3':'REZAIEE'}.items():
    assert dump.resolve_reference(v3,alias)[0] == target
print('All T6 dump checks PASS')
