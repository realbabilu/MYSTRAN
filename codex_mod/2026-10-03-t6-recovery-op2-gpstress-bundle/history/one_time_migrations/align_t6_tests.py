from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parent
PY = ROOT / 'python'
IMPORTS = '''from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from MITC6_Tri_v4 import MITC6_Tri_v4
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
'''
REGISTRY = '''{
    "SIMOT6": Simo1993_Tri6_ShellElement_v2,
    "MITC6": MITC6_Tri_v4,
    "MH6T": MacNeal_MH6T_Tri_v3,
    "REZAIEE": Rezaiee2017_Tri6_ShellElement_v3,
}'''
FILES = [
    'battle_2_002_quadratic_clamped.py',
    'battle_2_003_quadratic_curved.py',
    'battle_2_004_twistedv2.py',
    'battle_2_006_scordelisv2.py',
    'test_quadraticv2a.py',
]
BACKUP = ROOT / 'test' / 't6_reference_alignment_backup'
BACKUP.mkdir(exist_ok=True)

def save(path, source):
    target = BACKUP / path.relative_to(ROOT)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        shutil.copy2(path, target)
    path.write_text(source, encoding='utf-8')

for name in FILES:
    path = PY / name
    source = path.read_text(encoding='utf-8-sig')
    # Remove superseded direct T6 imports; Q8 definitions stay independent.
    source = re.sub(r'^from (?:Simo1993_Tri6\w*|MacNeal_MH6T\w*|MITC6_Tri\w*|Rezaiee2017_Tri6\w*)\s+import[^\n]*\n', '', source, flags=re.M)
    source = source.replace('from core import Model\n', 'from core import Model\n' + IMPORTS, 1)
    key = 'T6_ELEMS' if name == 'test_quadraticv2a.py' else 'T6_ELEMENTS'
    source, count = re.subn(r'\b' + key + r'\s*=\s*\{[^}]*\}', key + ' = ' + REGISTRY, source, count=1)
    assert count == 1, name
    save(path, source)

# Use one final implementation per T6 family in the patch-test registry.
path = PY / 'test_patch2001_q8_t6_v4c.py'
source = path.read_text(encoding='utf-8-sig')
start = source.index("    'SIMOT6': (") if "    'SIMOT6': (" in source else source.index("    'MITC6': (")
end = source.index('# ═', start)
entries = []
for label, cls in [('SIMOT6', 'Simo1993_Tri6_ShellElement_v2'),
                   ('MITC6', 'MITC6_Tri_v4'),
                   ('MH6T', 'MacNeal_MH6T_Tri_v3'),
                   ('REZAIEE', 'Rezaiee2017_Tri6_ShellElement_v3')]:
    entries.append(f"    '{label}': (\n        {cls}, stress_t6, 'T6',\n        NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6),\n")
source = source[:start] + '\n'.join(entries) + '}\n\n' + source[end:]
save(path, source)

# Keep old solver output as historical data; only rename the input deck.
old = ROOT / 'test' / 'duel3a_mht6_1.dat'
new = ROOT / 'test' / 'duel3a_mh6t_1.dat'
if old.exists():
    assert not new.exists(), new
    shutil.copy2(old, BACKUP / old.name)
    old.rename(new)
print('Aligned six test registries; renamed MH6T input deck. Backups:', BACKUP)

for rel in ['python/compare_disp_v4c.py', 'python/compare_f06_vs_python.py',
            'python/test_result_2_001_t6_disp.py', 'test/compare_f06_vs_python.py']:
    path = ROOT / rel
    source = path.read_text(encoding='utf-8-sig')
    save(path, source.replace('duel3a_mht6_1', 'duel3a_mh6t_1'))

# These two historical Q8 modules are absent locally. Keep them available
# when installed, and report the omission without blocking all T6 runs.
path = PY / 'battle_2_006_scordelisv2.py'
source = path.read_text(encoding='utf-8-sig')
source = source.replace('from Simo1993_Q8_ShellElement_V2_MacNealPatched import (\n    Simo1993_Q8_ShellElement_v2_MacNealPatched,\n)\n', '')
source = source.replace('from Simo1993_Q8_ShellElement_v3 import Simo1993_Q8_ShellElement_v3\n', '')
source = source.replace('    "SimoQ8v2": Simo1993_Q8_ShellElement_v2_MacNealPatched,\n', '')
source = source.replace('    "SimoQ8v3": Simo1993_Q8_ShellElement_v3,\n', '')
marker = 'T6_ELEMENTS = '
optional = '''# Optional historical Q8 variants; missing files must not block T6 tests.
from importlib import import_module
for _label, _module, _class in (
    ("SimoQ8v2", "Simo1993_Q8_ShellElement_V2_MacNealPatched", "Simo1993_Q8_ShellElement_v2_MacNealPatched"),
    ("SimoQ8v3", "Simo1993_Q8_ShellElement_v3", "Simo1993_Q8_ShellElement_v3"),
):
    try:
        Q8_ELEMENTS[_label] = getattr(import_module(_module), _class)
    except ModuleNotFoundError as exc:
        if exc.name != _module:
            raise
        print(f"SKIP {_label}: {_module}.py is unavailable")

'''
if '# Optional historical Q8 variants' not in source:
    source = source.replace(marker, optional + marker, 1)
save(path, source)
