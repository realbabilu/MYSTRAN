# Resume Perubahan MYSTRAN — CTRIA6_SIMO1993 & CTRIA6_MITC6

**Tanggal**: 3 Oktober 2026
**Tujuan**: Sinkronisasi MYSTRAN Fortran (`Source/EMG/EMG4/`) dengan Python reference (`Simo1993_Tri6_ShellElement_v2.py` untuk SIMOT6, `MITC6_Tri_v4.py` untuk MITC6).
**Verifikasi**: deck `test/prob_2_001_ctria6_simot6_disp.dat` (SAP2000 2-001 patch) + `test/prob_2_050_nx16_mitc6.dat` & `test/prob_2_050_nx16_ctria6.dat` (T9 pinched cylinder N=16, ref 1.8248e-5).

---

## 1. `Source/EMG/EMG4/CTRIA6_SIMO1993.f90` — port dari `Simo1993_Tri6_ShellElement_v2.py`

### 1.1 `CALC_NODAL_NORMALS_T6` — koheren director sign (centroid)
**v1p8 (lama)**: setiap nodal director di-force ke `+Z` lewat tie-break per-node
```fortran
! if normal[2] < 0: normal = -normal       ! PER-NODE
```
**v2 (baru)**: centroid-derived `(N_gipos)` dipakai konsisten untuk 6 nodal director
```fortran
CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
G1 = MATMUL(DN(1,:), XYZN); G2 = MATMUL(DN(2,:), XYZN)
CALL CROSS3(G1, G2, NC)
NORMAL_SIGN = -ONE
IF (NC(3) >= -1.0D-06*VNORM(NC)) NORMAL_SIGN = ONE
DO II=1,6
   ... 
   NORMS(II,:) = NORMAL_SIGN*N/NM
ENDDO
```
**Tujuan**: pada mesh warped (midside directors point opposite to centroid), sign konsiten antar node mencegah drift `cross(g,t)` vs `cross(t,g)`.

### 1.2 `BB_T6_AT` — flip sign translation block
**v1p8 (lama)**:
```fortran
DRHO11 = (/ DN(1,II)*T1(1), ... /)        ! +a1*t0_xi1
```
**v2 (baru)**:
```fortran
DRHO11 = -(/ DN(1,II)*T1(1), ... /)        ! -a1*t0_xi1
DRHO22 = -(/ DN(2,II)*T2(1), ... /)
DRHO12 = -(/ 0.5D0*(DN(1,II)*T2(1) + DN(2,II)*T1(1)), ... /)
```
**Tujuan**: pairing dengan director sign (1.1) — block transport/translation pakai `-` sementara curvature block pakai `cross(g,t)`. Bersign-consistent dgn Python `_compute_Bb`.

### 1.3 `BS_T6_AT` — hapus normalisasi T0
**v1p8 (lama)**:
```fortran
T0 = MATMUL(NVAL, NORMS_LOC)
NM = VNORM(T0); IF (NM > 1.0D-15) T0 = T0/NM   ! normalized
```
**v2 (baru)**:
```fortran
T0 = MATMUL(NVAL, NORMS_LOC)                   ! UN-normalized
```
**Tujuan**: shear strain γ = (∇u + ∇uᵀ)·t + (∇t)·u — magnitude T0 harus konsisten antara displacement & rotation contribution.

### 1.4 `BDRILL_T6_AT` — metric-form pointwise drilling
**v1p8 (lama)**: coefficient matrix `[[g·e]]`
```fortran
A(1,1) = DOT_PRODUCT(G1,E1); A(1,2) = DOT_PRODUCT(G1,E2)
A(2,1) = DOT_PRODUCT(G2,E1); A(2,2) = DOT_PRODUCT(G2,E2)
CALL INV2(A, AINV); DLOC = MATMUL(AINV, DN)
```
**v2 (baru)**: metric matrix `[[g·g]]`, dr/ds lalu proyek ke e1, e2
```fortran
A11 = DOT_PRODUCT(G1,G1); A22 = DOT_PRODUCT(G2,G2)
A12 = DOT_PRODUCT(G1,G2); DET = A11*A22 - A12*A12
IF (DABS(DET) <= 1.0D-30) THEN
   DR = ZERO; DS = ZERO
ELSE
   AI11 = A22/DET; AI22 = A11/DET; AI12 = -A12/DET
   DO II=1,6
      DR(II) = AI11*DN(1,II) + AI12*DN(2,II)
      DS(II) = AI12*DN(1,II) + AI22*DN(2,II)
   ENDDO
ENDIF
DO II=1,6
   DXII = DR(II)*DOT_PRODUCT(G1,E1) + DS(II)*DOT_PRODUCT(G2,E1)
   DYII = DR(II)*DOT_PRODUCT(G1,E2) + DS(II)*DOT_PRODUCT(G2,E2)
   COL = (II-1)*6
   BDOUT(1,COL+1) = 0.5D0*(DXII*E2(1) - DYII*E1(1))
   BDOUT(1,COL+2) = 0.5D0*(DXII*E2(2) - DYII*E1(2))
   BDOUT(1,COL+3) = 0.5D0*(DXII*E2(3) - DYII*E1(3))
   BDOUT(1,COL+4) = -NVAL(II)*NORMS_LOC(II,1)
   BDOUT(1,COL+5) = -NVAL(II)*NORMS_LOC(II,2)
   BDOUT(1,COL+6) = -NVAL(II)*NORMS_LOC(II,3)
ENDDO
```
**Tujuan**: match Python v2 `_compute_Bdrill`. Coefficient form hanya sama dgn metric untuk orthonormal g1,g2; untuk curvature g1, g2 (pinched cylinder), hanya metric form yang benar.

---

## 2. `Source/EMG/EMG4/CTRIA6_MITC6.f90` — port dari `MITC6_Tri_v4.py`

`CTRIA6_MITC6.f90` punya local copy dari semua helper (tidak share dengan SIMO1993), maka harus di-patch terpisah.

### 2.1 `CALC_NODAL_NORMALS_T6` — sama dengan SIMO1993 §1.1
```fortran
CALL SHAPE_T6(ONE/3.0D0, ONE/3.0D0, NVAL, DN)
G1 = MATMUL(DN(1,:), XYZN); G2 = MATMUL(DN(2,:), XYZN)
CALL CROSS3(G1, G2, NC)
NM = VNORM(NC); NORMAL_SIGN = -ONE
IF (NM > 1.0D-15) THEN
   IF (NC(3) >= -1.0D-06*NM) NORMAL_SIGN = ONE
ENDIF
DO II=1,6
   ...
   NORMS(II,:) = NORMAL_SIGN*N/NM
ENDDO
```
**Tujuan**: koheren director sign (sama dengan SIMO1993 §1.1).

### 2.2 `NAT_SHEAR_ROWS_T6` — hapus normalisasi T0 (sama dengan SIMO1993 §1.3)
```fortran
CALL SURFACE_BASIS_T6(XYZN, R, S, G1, G2, E1, E2, E3, JAC)
T0 = MATMUL(NVAL, NORMS)         ! UN-normalized (was: T0 = T0/NM)
```
**Tujuan**: shear strain γ dengan magnitude T0.

### 2.3 `BDRILL_T6_AT` — metric-form drilling (sama dengan SIMO1993 §1.4)
```fortran
A11 = DOT_PRODUCT(G1,G1); A22 = DOT_PRODUCT(G2,G2)
A12 = DOT_PRODUCT(G1,G2); DET = A11*A22 - A12*A12
IF (DABS(DET) <= 1.0D-30) THEN
   DR = ZERO; DS = ZERO
ELSE
   AI11 = A22/DET; AI22 = A11/DET; AI12 = -A12/DET
   DO II=1,6
      DR(II) = AI11*DN(1,II) + AI12*DN(2,II)
      DS(II) = AI12*DN(1,II) + AI22*DN(2,II)
   ENDDO
ENDIF
BDOUT = ZERO
DO II=1,6
   DXII = DR(II)*DOT_PRODUCT(G1,E1) + DS(II)*DOT_PRODUCT(G2,E1)
   DYII = DR(II)*DOT_PRODUCT(G1,E2) + DS(II)*DOT_PRODUCT(G2,E2)
   COL = (II-1)*6
   BDOUT(1,COL+1) = 0.5D0*(DXII*E2(1) - DYII*E1(1))
   BDOUT(1,COL+2) = 0.5D0*(DXII*E2(2) - DYII*E1(2))
   BDOUT(1,COL+3) = 0.5D0*(DXII*E2(3) - DYII*E1(3))
   BDOUT(1,COL+4) = -NVAL(II)*NORMS_LOC(II,1)
   BDOUT(1,COL+5) = -NVAL(II)*NORMS_LOC(II,2)
   BDOUT(1,COL+6) = -NVAL(II)*NORMS_LOC(II,3)
ENDDO
```
**Tujuan**: match Python v4 `_compute_Bdrill`.

### 2.4 `BB_T6_AT` — BELUM di-patch
**v1 (lama, saat ini di MYSTRAN)**:
```fortran
CALL TENSOR_PHYS_T6(DN(1,II)*T1, DN(2,II)*T2, 0.5D0*(DN(1,II)*T2 + DN(2,II)*T1), C, VP)
BBOUT(1:3,COL+1:COL+3) = VP
CALL CROSS3(T0, G1, CG1)                  ! t × g
CALL CROSS3(T0, G2, CG2)
CALL TENSOR_PHYS_T6(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
BBOUT(1:3,COL+4:COL+6) = VP
```
**v4 (target)**: split menjadi curvature block (`cross(G1, T0)` = g×t) + translation block (`-a1*t0_xi1`)
**Status**: dampak T9 kecil (ratio bergeser 1.046 → 1.046); belum diaplikasikan.

---

## 4. Verifikasi

### 4.1 SAP2000 2-001 Patch Test (`prob_2_001_ctria6_simot6_disp.dat` / `prob_2_001_ctria6_mitc6_disp.dat`)

| Stage | MYSTRAN T1/T2/T3 vs SPCD | IKUT |
|---|---|---|
| SIMOT6 Stage 1 (membrane) | match persis, R3 noise = 1e-17 | ✓ |
| SIMOT6 Stage 2 (bending)  | match persis | ✓ |
| MITC6 Stage 1 (membrane)  | match persis, R3 noise = 1.2e-16 | ✓ |

### 4.2 T9 Pinched Cylinder N=16 (ref 1.8248e-5)

| Element | Python ref (ratio) | MYSTRAN (ratio) | Selisih |
|---|---|---|---|
| SIMOT6 (`Simo1993_Tri6_ShellElement_v2`) | 1.155 | **1.156** | 0.001 |
| MITC6 (`MITC6_Tri_v4`) | 1.0448 | **1.0458** | 0.001 |

### 4.3 Files modified

- `C:/PROJECTAI/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_SIMO1993.f90`
- `C:/PROJECTAI/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MITC6.f90`
- Deck test T9 N=16: `C:/PROJECTAI/18a/test/prob_2_050_nx16_ctria6.dat` (SIMOT6) dan `prob_2_050_nx16_mitc6.dat` (MITC6).

### 4.4 Build command
```bash
cd C:/PROJECTAI/18a/MYSTRAN/build && make -j4 mystran
```

## 5. Update 2026-10-03: four-family displacement and recovery alignment

SIMOT6 v2, MITC6 v4, MH6T v3 and REZAIEE v3 are now aligned to the final Python references. Displacement 2-001/2-002/2-003/2-004 passes 88/88 cases on mesh 8/24 and final-build regression passes 88/88 on mesh 2/4. F06 STEP A CENTER, STEP B all six nodal CORNER points and STEP C GPSTRESS pass for membrane/bending on the 2-001 patch for all four families (8/8 cases); local bending moments also match. Standalone output selectors pass 32/32 checks.

`dump_python_debug.py` uses final references, correct element IDs and natural node coordinates. Detailed changes, limitations, validation artifacts and rerun commands are in `C:/PROJECTAI/18a/test/T6_UPGRADE_VALIDATION_2026-10-03.md`. Q8 was not upgraded. Curved/twisted stress and OP2 validation remain outside these completed gates.

## 6. Native OP2 T6 update 2026-10-03

T6 OES/OEF now use CTRIA6 type 75, widths 70/38. Dispatch accepts TRIA6; OES reopens correctly after OGS; each OGS surface uses its actual ID. All four T6 families pass pyNastran reads and stress/force/GPSTRESS comparisons to F06. The mixed `test/duel3a.OP2` reads completely. Siemens Nastran 2512 was run as a reference: native T6 table widths, shapes and element IDs match. F06 recovery regression remains 8/8; output selectors remain 32/32. See `C:/PROJECTAI/18a/test/T6_OP2_VALIDATION_2026-10-03.md` and reader example `test/read_t6_op2.py`. Other families' numerical OP2 completeness and Femap GUI import are not validated here.
