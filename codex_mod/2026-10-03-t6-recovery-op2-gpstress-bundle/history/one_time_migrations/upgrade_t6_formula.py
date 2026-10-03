"""One-time T6-only port corrections to the user's final Python formulas."""
from pathlib import Path
import re, shutil
root=Path(__file__).resolve().parents[1]
emg=root/'MYSTRAN/Source/EMG/EMG4'
backup=root/'test/t6_upgrade/source_before_formula'
backup.mkdir(parents=True,exist_ok=True)
def routine(src,name):
    return re.search(r'      SUBROUTINE '+name+r'\b.*?      END SUBROUTINE '+name,src,re.S).group()
def replace(src,name,new):
    return src.replace(routine(src,name),new)
simo=(emg/'CTRIA6_SIMO1993.f90').read_text(encoding='utf-8-sig')
mitc=(emg/'CTRIA6_MITC6.f90').read_text(encoding='utf-8-sig')
mh=(emg/'CTRIA6_MH6T.f90').read_text(encoding='utf-8-sig')
norm=routine(simo,'CALC_NODAL_NORMALS_T6')
# Python references use element pointwise normals unless supplied explicitly.
norm=re.sub(r'         IF \(ALLOCATED\(GRID_SNORM\)\) THEN.*?         ENDIF\n      ENDDO',
            '      ENDDO',norm,flags=re.S)
drill=routine(mitc,'BDRILL_T6_AT')
for family in ('SIMO1993','MITC6','MH6T','REZAIEE'):
    path=emg/f'CTRIA6_{family}.f90'
    src=path.read_text(encoding='utf-8-sig')
    target=backup/path.name
    if not target.exists(): shutil.copy2(path,target)
    # Pointwise coherent physical basis, as _local_basis_at_point in Python.
    basis=routine(src,'SURFACE_BASIS_T6')
    basis=basis.replace('G3(3), TMP(3), NM','G3(3), TMP(3), NC(3), NM, SGN')
    anchor='      NM = VNORM(G1)'
    orient='''      CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
      NC = MATMUL(DN(1,:), XYZN)
      TMP = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(NC, TMP, G3)
      SGN = ONE
      IF (G3(3) < -1.0D-6*VNORM(G3)) SGN = -ONE
      E3 = SGN*E3
'''
    assert anchor in basis
    basis=basis.replace(anchor,orient+anchor,1)
    src=replace(src,'SURFACE_BASIS_T6',basis)
    cov=routine(src,'COV_MAP_T6')
    cov=cov.replace('DOT_PRODUCT(E1F, GC','DOT_PRODUCT(E1, GC').replace('DOT_PRODUCT(E2F, GC','DOT_PRODUCT(E2, GC')
    src=replace(src,'COV_MAP_T6',cov)
    src=replace(src,'CALC_NODAL_NORMALS_T6',norm)
    if family in ('MITC6','REZAIEE'):
        src=replace(src,'BB_T6_AT',routine(mh,'BB_T6_AT'))
    if family in ('MH6T','REZAIEE'):
        d=drill.replace('-NVAL(II)*NORMS_LOC(II,1)','-NVAL(II)*E3(1)').replace('-NVAL(II)*NORMS_LOC(II,2)','-NVAL(II)*E3(2)').replace('-NVAL(II)*NORMS_LOC(II,3)','-NVAL(II)*E3(3)')
        src=replace(src,'BDRILL_T6_AT',d)
    if family=='REZAIEE':
        # v3 uses edge membrane samples and RC_SHEAR for the off-edge sample.
        old=routine(src,'BM_REZAIEE_AT')
        decl=old[:old.index('      CALL FIXED_FRAME_T6')]
        body='''      CALL FIXED_FRAME_T6(XYZN, E1F, E2F, E3F)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, ZERO, RR11, SS11, RS11)
      CALL NAT_MEM_ROWS_T6(XYZN, R2MITC, ZERO, RR21, SS21, RS21)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, 1.0D0/DSQRT(3.0D0), RRC, SSC, RSC)
      VALS3(1,:) = RR11
      VALS3(2,:) = RR21
      VALS3(3,:) = RRC
      PTS3(1,:) = (/R1MITC, ZERO/)
      PTS3(2,:) = (/R2MITC, ZERO/)
      PTS3(3,:) = (/R1MITC, 1.0D0/DSQRT(3.0D0)/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 1, A1, B1, C1V)
      CALL NAT_MEM_ROWS_T6(XYZN, ZERO, R1MITC, RR11, SS11, RS11)
      CALL NAT_MEM_ROWS_T6(XYZN, ZERO, R2MITC, RR12, SS12, RS12)
      CALL NAT_MEM_ROWS_T6(XYZN, 1.0D0/DSQRT(3.0D0), R1MITC, RRC, SSC, RSC)
      VALS3(1,:) = SS11
      VALS3(2,:) = SS12
      VALS3(3,:) = SSC
      PTS3(1,:) = (/ZERO, R1MITC/)
      PTS3(2,:) = (/ZERO, R2MITC/)
      PTS3(3,:) = (/1.0D0/DSQRT(3.0D0), R1MITC/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 1, A2, B2, C2V)
      CALL NAT_MEM_ROWS_T6(XYZN, R2MITC, R1MITC, RR21, SS21, RS21)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, R2MITC, RR12, SS12, RS12)
      CALL NAT_MEM_ROWS_T6(XYZN, R1MITC, R1MITC, RR11, SS11, RS11)
      VALS3(1,:) = RR21 + SS21 - TWO*RS21
      VALS3(2,:) = RR12 + SS12 - TWO*RS12
      VALS3(3,:) = RR11 + SS11 - TWO*RS11
      PTS3(1,:) = (/R2MITC, R1MITC/)
      PTS3(2,:) = (/R1MITC, R2MITC/)
      PTS3(3,:) = (/R1MITC, R1MITC/)
      CALL FIT_AFFINE36_T6(VALS3, PTS3, 2, A3, B3, C3V)
      BERR = A1 + B1*R + C1V*S
      BESS = A2 + B2*R + C2V*S
      BEQQ = A3 + B3*R + C3V*(ONE - R - S)
      BERS = 0.5D0*(BERR + BESS - BEQQ)
'''
        tail=old[old.index('      CALL COV_MAP_T6'):]
        src=replace(src,'BM_REZAIEE_AT',decl+body+tail)
        shear=routine(mitc,'BS_MITC6_AT').replace('BS_MITC6_AT','BS_REZAIEE_AT').replace('RCMITC','RCSHEAR')
        shear=shear.replace('XYZN, NORMS, R1MITC, R1MITC, BRT','XYZN, NORMS, ZERO, R1MITC, BRT').replace('XYZN, NORMS, R1MITC, R2MITC, BRT, BST, JAC)\n      RHS(4','XYZN, NORMS, ZERO, R2MITC, BRT, BST, JAC)\n      RHS(4')
        shear=shear.replace('ONE, R1MITC, R1MITC, -(R1MITC*R1MITC), -(R1MITC*R1MITC)','ONE, ZERO, R1MITC, ZERO, ZERO').replace('ONE, R1MITC, R2MITC, -(R1MITC*R1MITC), -(R1MITC*R2MITC)','ONE, ZERO, R2MITC, ZERO, ZERO')
        src=replace(src,'BS_REZAIEE_AT',shear)
        # Python v3 common shear point is RC=1/3, not RC_SHEAR.
        src=src.replace('RCSHEAR = ONE/SQRT3','RCSHEAR = ONE/3.0D0')
        nat=routine(src,'NAT_SHEAR_ROWS_T6')
        nat=re.sub(r'      NM = VNORM\(T0\)\n      IF \(NM > 1.0D-15\) T0 = T0/NM\n','',nat)
        src=replace(src,'NAT_SHEAR_ROWS_T6',nat)
    path.write_text(src,encoding='utf-8')
print('Applied T6 formula port corrections. Backups:',backup)
