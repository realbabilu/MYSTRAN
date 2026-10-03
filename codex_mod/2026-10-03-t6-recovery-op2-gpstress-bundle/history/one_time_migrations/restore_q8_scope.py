from pathlib import Path
import re
root = Path(__file__).resolve().parents[1]
path = root/'python/battle_2_006_scordelisv2.py'
original = (root/'test/t6_reference_alignment_backup/python/battle_2_006_scordelisv2.py').read_text(encoding='utf-8-sig')
current = path.read_text(encoding='utf-8-sig')
registry = re.search(r'T6_ELEMENTS\s*=\s*\{[^}]*\}', current).group()
imports = '\n'.join(line for line in current.splitlines() if line.startswith(('from Simo1993_Tri6_ShellElement_v2 ', 'from MacNeal_MH6T_Tri_v3 ', 'from MITC6_Tri_v4 ', 'from Rezaiee2017_Tri6_v3 ')))
original = original.replace('from core import Model\n', 'from core import Model\n'+imports+'\n')
original = re.sub(r'T6_ELEMENTS\s*=\s*\{[^}]*\}', registry, original, count=1)
path.write_text(original,encoding='utf-8')
print('Q8 imports and registry restored; final T6 registry retained.')
