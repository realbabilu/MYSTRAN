# Roadmap Integrasi DKMQ24 / DKMT18 / MITC3+ ke MYSTRAN

Dokumen ini merangkum rencana teknis untuk memasukkan elemen shell baru ke MYSTRAN/Nastran-style workflow:

- `CQUADR` / `CQUAD4-like` → DKMQ24 / DKMQ-based quadrilateral shell
- `CTRIAR` / triangular higher shell → DKMT18
- `CTRIA3` option → MITC3+ / MITC3+ SNORM
- Composite branch nanti: `PCOMP` / `PCOMPG` → ABD/As laminate

Catatan penting:

> `MITC3+` untuk tahap awal hanya dipakai untuk static/modal.  
> `MITC3+` **tidak dimasukkan ke buckling dulu**, karena geometric stiffness / prestress recovery-nya belum divalidasi seperti DKMT18 dan DKMQ24.

---

## 1. Strategi Branch

Disarankan pakai dua branch besar agar aman untuk MYSTRAN existing code.

```text
feature/dkmq-dkmt-isotropic
    Fokus:
      - CQUADR / DKMQ24 isotropic
      - CTRIAR / DKMT18 isotropic
      - CTRIA3 / MITC3+ static/modal only
      - MAT1 + PSHELL
      - static benchmarks
      - modal consistent mass
      - buckling operator untuk DKMQ24/DKMT18 only

feature/dkmq-dkmt-composite
    Fokus:
      - PCOMP / PCOMPG parser
      - CLT ABD/As laminate
      - rho per ply
      - ply stress recovery
      - composite static/modal/buckling
      - DKMQ24/DKMT18 composite
```

Jangan langsung mengubah default existing element MYSTRAN. Aktifkan lewat parameter atau option.

Contoh parameter internal/testing:

```nastran
PARAM,SHELLFORM,NEW
PARAM,CQUADRFORM,DKMQ24
PARAM,CTRIARFORM,DKMT18
PARAM,CTRIA3FORM,MITC3PLUS
PARAM,BUCKOP,ARPACK
```

---

## 2. Mapping Elemen

### 2.1 `CQUADR` / Quadrilateral Shell

Target:

```text
CQUADR / CQUAD4-like
    -> DKMQ24
```

Fungsi:

```text
Static:
    yes

Modal:
    yes, pakai consistent mass

Buckling:
    yes, pakai operator K^-1 KG

Composite:
    nanti, branch composite
```

Material path tahap isotropic:

```text
MAT1 + PSHELL
    -> A, B, D, As equivalent
    -> B = 0 untuk isotropic
```

Formula isotropic equivalent:

```text
Q = E/(1-nu^2) * [[1, nu, 0],
                  [nu, 1, 0],
                  [0, 0, (1-nu)/2]]

A  = Q * t
B  = 0
D  = Q * t^3 / 12
As = kappa * G * t * I
```

### 2.2 `CTRIAR`

Target:

```text
CTRIAR
    -> DKMT18
```

Fungsi:

```text
Static:
    yes

Modal:
    yes, pakai consistent mass

Buckling:
    yes, pakai operator K^-1 KG

Composite:
    nanti, branch composite
```

`CTRIAR` lebih cocok untuk DKMT18 dibanding `CTRIA3` karena DKMT18 adalah triangular shell yang lebih kaya field interpolation.

### 2.3 `CTRIA3`

Target:

```text
CTRIA3
    -> existing CTRIA3 default atau MITC3+ optional
```

Untuk MITC3+:

```text
Static:
    yes

Modal:
    experimental yes setelah consistent mass stabil

Buckling:
    no, jangan dulu
```

Alasan MITC3+ buckling ditunda:

```text
1. KG / geometric stiffness belum divalidasi.
2. Prestress resultant recovery belum dibuktikan sekuat DKMQ24 dan DKMT18.
3. MITC3+ pada patch test shear punya karakter khusus.
4. Untuk buckling, coupling rotasi-transverse harus sangat hati-hati.
```

Jika user meminta buckling dengan MITC3+:

```text
WARNING:
MITC3+ buckling is not implemented/validated.
Use DKMT18/CTRIAR or DKMQ24/CQUADR for buckling.
```

---

## 3. Material dan Density

### 3.1 Non-composite `PSHELL + MAT1`

Untuk isotropic/non-composite:

```text
rhoh = rho * t + NSM
rhoI = rho * t^3 / 12
```

`rho` diambil dari `MAT1`.

Jika `PSHELL` punya `NSM`, tambahkan ke mass per area:

```text
rhoh_total = rho * t + NSM
```

Untuk rotary inertia tahap awal:

```text
rhoI = rho * t^3 / 12
```

NSM rotary inertia bisa diabaikan dulu atau diperlakukan sesuai existing MYSTRAN convention.

### 3.2 Composite `PCOMP`

Untuk composite:

```text
Setiap ply bisa punya MID berbeda.
Setiap MID bisa punya E, nu, G, rho berbeda.
```

Per ply:

```text
E1_i, E2_i, nu12_i, G12_i, G13_i, G23_i
rho_i
theta_i
t_i
z_bot_i, z_top_i
```

Mass:

```text
rhoh = NSM + Σ rho_i * (z_top_i - z_bot_i)

rhoI = Σ rho_i * (z_top_i^3 - z_bot_i^3) / 3
```

Stiffness:

```text
A = Σ Qbar_i * (z_top_i - z_bot_i)

B = 1/2 Σ Qbar_i * (z_top_i^2 - z_bot_i^2)

D = 1/3 Σ Qbar_i * (z_top_i^3 - z_bot_i^3)
```

Shear:

```text
As = Σ kappa_i * Qsbar_i * t_i
```

---

## 4. PSHELL `MID1/MID2/MID3/MID4`

`MID1/MID2/MID3/MID4` adalah PSHELL Nastran-way ala MSC/NX, bukan PCOMP utama.

Konsep PSHELL:

```text
MID1:
    membrane material

MID2:
    bending material

MID3:
    transverse shear material

MID4:
    membrane-bending coupling material
```

Untuk tahap awal MYSTRAN:

```text
MID1:
    wajib untuk PSHELL isotropic path

MID2:
    boleh support untuk bending stiffness D

MID3:
    boleh support untuk shear As

MID4:
    jangan default dulu
    experimental only
```

Default aman:

```text
MID2 blank:
    inherit MID1

MID3 blank:
    inherit MID1

MID3 = 0:
    treat as high shear stiffness / zero shear flexibility jika mengikuti Nastran convention yang dipilih

MID4 blank:
    B = 0
```

Untuk composite `PCOMP`, jangan pakai `MID4` untuk coupling. Coupling `B` harus datang dari stacking sequence.

---

## 5. Static Stiffness

Unified stiffness formula untuk DKMQ24/DKMT18 composite-capable:

```text
K = ∫ (
      Bm^T A  Bm
    + Bm^T B  Bb
    + Bb^T B  Bm
    + Bb^T D  Bb
    + Bs^T As Bs
) dA
+ K_drill
```

Untuk isotropic/non-composite:

```text
B = 0
A, D, As dari PSHELL/MAT1
```

Untuk composite:

```text
A, B, D, As dari PCOMP/CLT
```

Secara branch implementasi boleh dipisah:

```text
branch isotropic:
    B = 0 only

branch composite:
    full A/B/D/As
```

---

## 6. Modal

### 6.1 Consistent Mass

Untuk shell:

```text
M_trans = ∫ rhoh * N^T N dA
          untuk ux, uy, uz

M_rot   = ∫ rhoI * N^T N dA
          untuk rx, ry, rz
```

Untuk isotropic:

```text
rhoh = rho * t + NSM
rhoI = rho * t^3 / 12
```

Untuk composite:

```text
rhoh = NSM + Σ rho_i * t_i
rhoI = Σ rho_i * (z_top_i^3 - z_bot_i^3) / 3
```

### 6.2 Modal Equation

```text
K phi = omega^2 M phi
```

Output:

```text
freq_Hz = sqrt(omega^2) / (2*pi)
```

### 6.3 Element Support

```text
DKMQ24:
    modal yes

DKMT18:
    modal yes

MITC3+:
    modal experimental

MITC3+ buckling:
    no
```

---

## 7. Buckling

### 7.1 Jangan Pakai Schur Dense untuk Production

Schur benchmark:

```text
Keff = Kaa - Kai Kii^-1 Kia
Keff phi_w = lambda KGww phi_w
```

Ini benar untuk validasi kecil, tapi berat untuk model besar.

### 7.2 Production Route: Operator `K^-1 KG`

Gunakan buckling equation:

```text
K phi = lambda KG phi
```

Ubah menjadi reciprocal eigenproblem:

```text
OP(x) = K^-1 KG x

OP phi = mu phi

lambda = 1 / mu
```

Workflow solver:

```text
1. Assemble K full shell
2. Assemble KG dari recovered membrane resultants
3. Apply boundary conditions
4. Factor K sekali
5. ARPACK reverse communication:
       y = KG * x
       z = solve(K, y)
       return z
6. ARPACK cari mu positif terbesar
7. lambda_crit = 1 / mu_max
```

Ini cocok untuk:

```text
DKMQ24
DKMT18
```

Tidak untuk MITC3+ tahap awal.

### 7.3 Recover Membrane Resultants untuk KG

Dari prebuckling displacement:

```text
eps0, kappa -> element strain/curvature
```

Resultants:

```text
N = A eps0 + B kappa
M = B eps0 + D kappa
```

Untuk isotropic:

```text
B = 0
N = A eps0
```

Untuk composite unsymmetric:

```text
B != 0
N = A eps0 + B kappa
M = B eps0 + D kappa
```

Geometric stiffness plate/shell awal:

```text
KG = ∫ G_w^T N G_w dA
```

Tanda kompresi harus konsisten:

```text
compression Nxx < 0
KG untuk eigen K = lambda KG memakai -N
```

---

## 8. Benchmark Wajib Sebelum Merge

### 8.1 Static Isotropic

Untuk DKMQ24/DKMT18:

```text
patch_tests.py
Scordelis-Lo
Raasch hook
LE5 Z cantilever
MacNeal shell tests jika tersedia
```

Untuk MITC3+:

```text
patch_tests.py
Scordelis-Lo
Raasch hook
LE5, jika stress recovery valid
```

MITC3+ buckling tidak dites dulu.

### 8.2 Static Composite

Untuk branch composite:

```text
composite_patch_tests.py
composite_field_patch_tests.py
composite_plate_bending.py
```

Expected:

```text
isotropic reproduction PASS
symmetric laminate B ≈ 0
unsymmetric laminate B != 0
Q4/T3 convergence OK
```

### 8.3 Modal

Non-composite:

```text
square plate modal isotropic
DKMQ24 vs DKMT18 convergence
```

Composite:

```text
sym_crossply
unsym_crossply
sym_angle
```

Mass:

```text
consistent mass default
lumped mass optional/debug
```

### 8.4 Buckling

Non-composite and composite:

```text
DKMQ24
DKMT18
```

Method:

```text
operator K^-1 KG
```

Do not use:

```text
active w-only eigen without condensation
```

MITC3+:

```text
buckling disabled
```

---

## 9. Warning / Error Behavior

### 9.1 MITC3+ Buckling

If user requests buckling with MITC3+:

```text
WARNING:
MITC3+ buckling is not implemented or validated.
Use DKMT18/CTRIAR or DKMQ24/CQUADR for buckling.
```

### 9.2 Composite with Missing Density

If ply material has no density:

```text
WARNING:
Ply MID xxxx has no density. Modal mass may be incomplete.
```

Option:

```text
- reject modal run
- or assume rho=0 and warn
```

For static, missing rho is not fatal.

### 9.3 PSHELL MID4

If `MID4` is present:

```text
WARNING:
PSHELL MID4 membrane-bending coupling is experimental in this shell formulation.
```

Or disable until validated.

---

## 10. Implementation Order

### Phase 1: Isotropic Branch

```text
1. Add DKMQ24/CQUADR path.
2. Add DKMT18/CTRIAR path.
3. Add MITC3+/CTRIA3 optional static path.
4. Add PSHELL/MAT1 -> A,D,As equivalent.
5. Add static benchmarks.
6. Add consistent mass.
7. Add modal benchmarks.
8. Add buckling operator K^-1 KG for DKMQ24/DKMT18 only.
```

### Phase 2: Composite Branch

```text
1. Add PCOMP/PCOMPG parser.
2. Add ply material orthotropic Qbar.
3. Add A/B/D/As CLT.
4. Add rhoh/rhoI per ply.
5. Add composite static benchmarks.
6. Add composite modal.
7. Add composite buckling operator.
8. Add ply stress/strain recovery.
9. Add failure index later.
```

---

## 11. Current Prototype Conclusion

### DKMQ24 / DKMT18 Static Composite

```text
PASS
```

Prototype evidence:

```text
isotropic reproduction PASS
ABD/B coupling PASS
unsymmetric static PASS
plate bending convergence good
```

### Modal Consistent Mass

```text
PASS
```

Observed:

```text
DKMQ24 vs DKMT18 mode-1 convergence good
consistent mass stable
```

### Buckling Operator

```text
PASS for DKMQ24/DKMT18 benchmark
```

Important lesson:

```text
Wrong:
    solve only active KG / w DOF

Correct:
    OP(x) = K^-1 KG x
    lambda = 1/mu
```

### MITC3+

```text
Static:
    usable / promising

Modal:
    experimental

Buckling:
    not implemented yet
```

---

## 12. Final Recommendation

For MYSTRAN integration:

```text
Use separate development branches:
    1. isotropic DKMQ/DKMT/MITC3+
    2. composite PCOMP DKMQ/DKMT

Do not enable MITC3+ buckling yet.

Use operator buckling for DKMQ24/DKMT18:
    OP = K^-1 KG

Use density per ply for PCOMP:
    rhoh = Σ rho_i t_i + NSM
    rhoI = Σ rho_i ∫z^2 dz
```

Best short-term default:

```text
CQUADR -> DKMQ24
CTRIAR -> DKMT18
CTRIA3 -> existing default, with optional MITC3+ for testing

Buckling:
    DKMQ24/DKMT18 only

Composite:
    branch 2 after isotropic branch stable
```
