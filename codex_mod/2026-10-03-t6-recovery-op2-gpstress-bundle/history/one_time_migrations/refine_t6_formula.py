from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
emg=root/'MYSTRAN/Source/EMG/EMG4'
for family in ('SIMO1993','MITC6','MH6T','REZAIEE'):
    path=emg/f'CTRIA6_{family}.f90'
    src=path.read_text(encoding='utf-8')
    if family in ('SIMO1993','MITC6'):
        src=src.replace('-NVAL(II)*NORMS_LOC(II,1)','-NVAL(II)*E3(1)').replace('-NVAL(II)*NORMS_LOC(II,2)','-NVAL(II)*E3(2)').replace('-NVAL(II)*NORMS_LOC(II,3)','-NVAL(II)*E3(3)')
    if family=='MH6T':
        src=re.sub(r'      IF \(\(DABS\(E3F\(1\)\).*?      ENDIF\n      END SUBROUTINE FIXED_FRAME_T6',
                   '      END SUBROUTINE FIXED_FRAME_T6',src,flags=re.S)
    path.write_text(src,encoding='utf-8')
print('Pointwise drilling normal and MH6T centroid frame aligned.')
