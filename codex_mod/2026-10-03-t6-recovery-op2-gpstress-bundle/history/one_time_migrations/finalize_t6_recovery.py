"""Preserve explicit SNORM cards and remove temporary MH6T diagnostics."""
from pathlib import Path
root = Path(__file__).resolve().parents[1]
for family in ('SIMO1993','MITC6','MH6T','REZAIEE'):
    path = root/f'MYSTRAN/Source/EMG/EMG4/CTRIA6_{family}.f90'
    s = path.read_text()
    s = s.replace('MAX_STRESS_POINTS, SOL_NAME', 'MAX_STRESS_POINTS, SOL_NAME, NSNORM',1)
    s = s.replace('ELGP, GRID_SNORM,', 'ELGP, GRID_SNORM, GRID_ID, SNORM,',1)
    s = s.replace('INTEGER(LONG) :: II, BIDX', 'INTEGER(LONG) :: II, BIDX, ISN',1)
    start = s.index('      SUBROUTINE CALC_NODAL_NORMALS_T6')
    end = s.index('      END SUBROUTINE CALC_NODAL_NORMALS_T6',start)
    block = s[start:end]
    insert = '''! Explicit SNORM overrides geometric directors; automatically averaged
! GRID_SNORM values do not replace the final Python geometric defaults.
         IF (NSNORM > 0 .AND. ALLOCATED(GRID_SNORM) .AND. ALLOCATED(SNORM)) THEN
            BIDX = BGRID(II)
            IF (BIDX > 0) THEN
               DO ISN=1,NSNORM
                  IF (SNORM(ISN,1) /= GRID_ID(BIDX)) CYCLE
                  SN = GRID_SNORM(BIDX,:)
                  NM = VNORM(SN)
                  IF (NM > 1.0D-15) THEN
                     SN = SN/NM
                     IF (DOT_PRODUCT(SN,NORMS(II,:)) < ZERO) SN = -SN
                     NORMS(II,:) = SN
                  ENDIF
                  EXIT
               ENDDO
            ENDIF
         ENDIF
'''
    pos = block.rfind('      ENDDO')
    block = block[:pos]+insert+block[pos:]
    s = s[:start]+block+s[end:]
    if family == 'MH6T':
        s = '\n'.join(line for line in s.split('\n') if "'T6K250 '" not in line)
        start = s.index('      IF ((DEBUG(250) > 0) .AND. (EID == 1)) THEN')
        end = s.index('      ENDIF',start)+len('      ENDIF\n')
        s = s[:start]+s[end:]
        s = s.replace('      USE DEBUG_PARAMETERS, ONLY      :  DEBUG\n','')
    if family == 'MITC6':
        s = s.replace('!   D:\\18a\\bending_only\\Shell\\gemini2\\shit\\validation\\q8\\MITC6_Tri_v1.py',
                      '!   python/MITC6_Tri_v4.py')
        s = s.replace('assumed covariant strains reconstructed from the Python v1 port.',
                      'assumed covariant strains and bending from the final Python v4 reference.')
    path.write_text(s)
print('Explicit SNORM support preserved; temporary MH6T diagnostics removed.')
