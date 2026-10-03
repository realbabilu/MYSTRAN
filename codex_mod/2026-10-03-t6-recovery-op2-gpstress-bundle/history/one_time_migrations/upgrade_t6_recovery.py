"""Align T6 recovery ordering with final Python: center, then six nodes."""
from pathlib import Path
import shutil
root = Path(__file__).resolve().parents[1]
backup = root / 'test/t6_upgrade/source_before_recovery'

def edit(rel, fn):
    path = root / rel
    saved = backup / path.name
    saved.parent.mkdir(parents=True, exist_ok=True)
    if saved.exists():
        raise RuntimeError(f'Already backed up: {path}')
    shutil.copy2(path, saved)
    source = path.read_text()
    path.write_text(fn(source))

for family in ('SIMO1993', 'MITC6', 'MH6T', 'REZAIEE'):
    def kernel(s):
        s = s.replace('RN6(6), SN6(6)', 'RN7(7), SN7(7)')
        if 'RN7(7)' not in s:
            s = s.replace('W6(6)', 'W6(6), RN7(7), SN7(7)', 1)
        s = s.replace('      RN6 = (/ONE/3.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0, 0.5D0/)\n', '')
        s = s.replace('      SN6 = (/ONE/3.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0/)\n', '')
        pos = s.index("      IF (OPT(1) == 'Y')")
        s = s[:pos] + ('      RN7 = (/ONE/3.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0, 0.5D0, 0.0D0/)\n'
                       '      SN7 = (/ONE/3.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.5D0, 0.5D0/)\n\n') + s[pos:]
        start = s.index("      IF (OPT(3) == 'Y')")
        end = s.index("      IF (OPT(4) == 'Y')", start)
        block = s[start:end]
        block = block.replace('DO I=1,3', 'DO I=1,7').replace('DO I=1,6', 'DO I=1,7', 1)
        block = block.replace('RN6(I)', 'RN7(I)').replace('SN6(I)', 'SN7(I)')
        block = block.replace('R3(I)', 'RN7(I)').replace('S3(I)', 'SN7(I)')
        block = block.replace('I <= MAX_STRESS_POINTS', 'I <= MAX_STRESS_POINTS+1')
        return s[:start]+block+s[end:]
    edit(f'MYSTRAN/Source/EMG/EMG4/CTRIA6_{family}.f90', kernel)

def allocation(s):
    end = s.index('      END SUBROUTINE CALC_MAX_STRESS_POINTS')
    start = s.rfind('      ENDDO', 0, end)+len('      ENDDO')
    return s[:start]+'\n      IF (NCTRIA6 > 0) MAX_STRESS_POINTS = MAX(MAX_STRESS_POINTS,6_LONG)'+s[start:]
edit('MYSTRAN/Source/LK1/L1A/LOADB.f90', allocation)

def stress(s):
    s = s.replace('NUM_PTS_ELEM = NUM_SEi(I) + 1  ! TRIA6 CORNER: CENTER + NUM_SEi corners',
                  "NUM_PTS_ELEM = 1\n                           IF (STRE_CORNER_REQ .OR. GPSTRESS_REQ .OR. &\n                               (STRE_LOC == 'CORNER  ') .OR. (STRE_LOC == 'GAUSS   ')) NUM_PTS_ELEM = 7")
    start = s.index('                   ! For CORNER request on TRIA6')
    end = s.index("                   IF (TYPE == 'BEAM    ')", start)
    s = s[:start]+s[end:]
    s = s.replace("                  IF (TYPE == 'BEAM    ') THEN\n                     STRESS_OUT(:,:) = STRESS_RAW(:,:)",
                  "                  IF ((TYPE == 'BEAM    ') .OR. (TYPE(1:5) == 'TRIA6')) THEN\n                     STRESS_OUT(:,:) = STRESS_RAW(:,:)")
    return s
edit('MYSTRAN/Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90', stress)

def writer(s):
    s = s.replace('NUM_GRID_PTS = 3  ! 3 corner nodes for CORNER request',
                  'NUM_GRID_PTS = MIN(6_LONG,NUM_PTS-1)  ! All six T6 nodes')
    start = s.index('      SUBROUTINE BUILD_SURFACE_PATCH_OUTPUT')
    end = s.index('      END FUNCTION GET_FIRST_PATCH_ELEM_FOR_GRID', start)
    block = s[start:end].replace('DO J=1,3', "DO J=1,MERGE(6,3,TYPE(1:5) == 'TRIA6')")
    return s[:start]+block+s[end:]
edit('MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRESSES.f90', writer)
print('T6 recovery updated; originals saved in', backup)
