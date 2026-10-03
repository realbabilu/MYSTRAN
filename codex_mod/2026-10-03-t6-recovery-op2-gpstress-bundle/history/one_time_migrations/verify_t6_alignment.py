from pathlib import Path
import ast
import importlib.util
import importlib
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'python'))
names = ['test_patch2001_q8_t6_v4c.py', 'battle_2_002_quadratic_clamped.py',
         'battle_2_003_quadratic_curved.py', 'battle_2_004_twistedv2.py',
         'battle_2_006_scordelisv2.py', 'test_quadraticv2a.py']
for name in names:
    path = ROOT / 'python' / name
    tree = ast.parse(path.read_text(encoding='utf-8-sig'))
    missing = sorted({node.module for node in ast.walk(tree)
                      if isinstance(node, ast.ImportFrom) and node.module
                      and importlib.util.find_spec(node.module) is None})
    print(name, 'syntax OK; missing imports:', missing)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in ('T6_ELEMENTS', 'T6_ELEMS') for t in node.targets):
            print('  T6:', ast.unparse(node.value))

import test_patch2001_q8_t6_v4c as patch
for mode in ('membrane', 'bending'):
    for name in ('SIMOT6', 'MITC6', 'MH6T', 'REZAIEE'):
        cls, stress, tag, nodes, conn, boundary, x0, y0 = patch.REGISTRY[name]
        ok, disp_error, stress_error, *_ = patch.build_and_solve(
            cls, nodes, conn, boundary, x0, y0, mode, stress, tag=tag)
        print(f'PATCH {mode} {name}: pass={ok} displacement_error={disp_error:.8e} stress_or_moment_error={stress_error:.8e}')

clamped = importlib.import_module('battle_2_002_quadratic_clamped')
curved = importlib.import_module('battle_2_003_quadratic_curved')
twisted = importlib.import_module('battle_2_004_twistedv2')
roof = importlib.import_module('battle_2_006_scordelisv2')
general = importlib.import_module('test_quadraticv2a')
for module in (clamped, curved, twisted, roof, general):
    registry = getattr(module, 'T6_ELEMENTS', getattr(module, 'T6_ELEMS', {}))
    assert set(registry) == {'SIMOT6', 'MITC6', 'MH6T', 'REZAIEE'}, module.__name__
    print('IMPORT/REGISTRY OK:', module.__name__)
roof.N_ARC = roof.N_LONG = 2
for name, cls in clamped.T6_ELEMENTS.items():
    for lc in (1, 2):
        print('SMOKE 2-002', name, lc, clamped.solve_one(2, 'T6', cls, lc))
        print('SMOKE 2-003', name, lc, curved.solve_one(name, 2, lc))
        print('SMOKE 2-004', name, lc, twisted.run_one(name, 2, lc))
    model = roof.build_t6_model(cls)
    model.build()
    model.solve_static()
    print('SMOKE 2-006', name, model.get_displacement(roof.nid(2, 2))[2])
