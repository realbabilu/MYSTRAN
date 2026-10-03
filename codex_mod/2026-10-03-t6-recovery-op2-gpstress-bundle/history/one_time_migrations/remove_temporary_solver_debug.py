"""Remove temporary matrix/recovery printing, preserving solver diagnostics."""
from pathlib import Path
import re, shutil
root = Path(__file__).resolve().parents[1]
backup = root/'test/t6_upgrade/source_before_debug_cleanup'
backup.mkdir(parents=True,exist_ok=True)

def edit(rel, fn):
    path = root/rel
    saved = backup/path.name
    if not saved.exists():
        shutil.copy2(path,saved)
    source = path.read_text()
    cleaned = fn(source)
    if source != cleaned:
        path.write_text(cleaned)
        print(path.name)

def remove_block(s, marker):
    while marker in s:
        pos = s.index(marker)
        start = s.rfind('\n',0,pos)+1
        end = s.index('ENDIF',pos)+len('ENDIF')
        end = s.index('\n',end)+1
        s = s[:start]+s[end:]
    return s

def kernel(s):
    s = remove_block(s,'IF ((DEBUG(246) > 0) .AND. (EID == 1)) THEN')
    return re.sub(r'^\s*USE DEBUG_PARAMETERS[^\n]*\n','',s,flags=re.M)
for family in ('SIMO1993','MITC6'):
    edit(f'MYSTRAN/Source/EMG/EMG4/CTRIA6_{family}.f90',kernel)

def recovery(s):
    s = re.sub(r"^.*WRITE\(ERR.*' DEBUG: ELEM_STRE_STRN 2D TYPE='[^\n]*\n",'',s,flags=re.M)
    s = remove_block(s,"IF ((TYPE(1:5) == 'TRIA6') .AND. (DEBUG(249) > 0) .AND. (EID <= 60)) THEN")
    return s.replace('EID, SHELL_T, ZS','EID, SHELL_T')
edit('MYSTRAN/Source/LK9/L92/ELEM_STRE_STRN_ARRAYS.f90',recovery)

def writer(s):
    s = remove_block(s,'IF (DEBUG(249) > 0) THEN')
    return re.sub(r'^\s*USE DEBUG_PARAMETERS[^\n]*\n','',s,flags=re.M)
edit('MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRESSES.f90',writer)

def force(s):
    pos = s.index("'DEBUG_Q8 EID=21")
    start = s.rfind('IF (EID == 21) THEN',0,pos)
    return remove_block(s,s[start:pos].split('\n')[0])
edit('MYSTRAN/Source/LK9/L92/OFP3_ELFE_2D.f90',force)
