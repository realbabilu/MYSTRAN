# CTRIA6 Formulation Smoke Compare

Date: August 15, 2026

## Scope

This note captures a quick numerical comparison of the four active `CTRIA6` shell formulations now exposed through `PARAM,TRIA6TYP`:

- `SIMOT6`
- `MITC6`
- `MH6T`
- `REZAIEE`

The purpose here is not full verification against NASTRAN or the Python reference yet. This is a first pass to see whether the four branches:

- compile
- assemble and solve
- produce comparable global displacement levels on common smoke decks

## Decks used

### Problem 2-002

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_simot6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_mitc6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_mh6t.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_rezaiee.dat`

### Problem 2-003

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_simot6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mitc6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mh6t.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_rezaiee.dat`

## Comparison basis

For this smoke pass, the comparison uses the `ABS*` row from the MYSTRAN displacement table in each subcase:

- `T1 T2 T3 R1 R2 R3`

This is intentionally coarse, but it is fast and good enough to flag major formulation divergence.

## Problem 2-002 summary

### ABS displacement rows

`SIMOT6`

- subcase 1: `2.582878E-04  3.970085E-03  6.672183E-03  6.345939E-03  1.296330E-03  1.124343E-03`
- subcase 2: `1.334486E-02  2.380737E-01  4.006895E-01  3.816465E-01  6.683837E-02  6.746176E-02`
- subcase 3: `1.064105E-02  2.384441E-01  2.378962E-01  3.819856E-01  3.980149E-02  6.744032E-02`
- subcase 4: `1.114615E-06  1.198609E-06  2.740125E-05  2.536528E-06  1.822182E-05  3.054211E-07`
- subcase 5: `3.124588E-03  6.671059E-02  6.671012E-02  6.351866E-02  1.113357E-02  2.024593E-02`
- subcase 6: `2.864824E-03  5.894556E-02  5.888389E-02  2.445876E-02  9.841259E-03  2.207884E-02`

`MITC6`

- subcase 1: `2.583020E-04  3.969966E-03  6.672865E-03  6.343973E-03  1.295443E-03  1.124372E-03`
- subcase 2: `1.334525E-02  2.380697E-01  4.007117E-01  3.814870E-01  6.681970E-02  6.746430E-02`
- subcase 3: `1.064076E-02  2.384341E-01  2.378983E-01  3.818120E-01  3.980575E-02  6.743912E-02`
- subcase 4: `1.116046E-06  1.201455E-06  2.745023E-05  2.610195E-06  1.816747E-05  3.055053E-07`
- subcase 5: `3.124564E-03  6.670888E-02  6.671500E-02  6.349493E-02  1.113754E-02  2.024565E-02`
- subcase 6: `2.864779E-03  5.894410E-02  5.888940E-02  2.444717E-02  9.849804E-03  2.207829E-02`

`MH6T`

- subcase 1: `2.583294E-04  3.971179E-03  6.673775E-03  6.351127E-03  1.296120E-03  1.124380E-03`
- subcase 2: `1.334701E-02  2.381546E-01  4.007573E-01  3.819854E-01  6.679977E-02  6.746603E-02`
- subcase 3: `1.064308E-02  2.385314E-01  2.379760E-01  3.823397E-01  3.980253E-02  6.744476E-02`
- subcase 4: `1.117191E-06  1.233917E-06  2.737549E-05  2.693117E-06  1.822016E-05  3.162305E-07`
- subcase 5: `3.124753E-03  6.671997E-02  6.671877E-02  6.359382E-02  1.113930E-02  2.024591E-02`
- subcase 6: `2.864762E-03  5.894860E-02  5.888781E-02  2.456081E-02  9.853193E-03  2.207839E-02`

`REZAIEE`

- subcase 1: `2.583068E-04  3.970239E-03  6.672060E-03  6.340441E-03  1.296679E-03  1.124439E-03`
- subcase 2: `1.334476E-02  2.380679E-01  4.006988E-01  3.818299E-01  6.681311E-02  6.746175E-02`
- subcase 3: `1.064127E-02  2.384457E-01  2.378832E-01  3.824953E-01  3.978898E-02  6.744094E-02`
- subcase 4: `1.115357E-06  1.185842E-06  2.732784E-05  1.634156E-06  1.823563E-05  2.797326E-07`
- subcase 5: `3.124686E-03  6.671119E-02  6.670833E-02  6.357314E-02  1.113601E-02  2.024600E-02`
- subcase 6: `2.865098E-03  5.894862E-02  5.887726E-02  2.447570E-02  9.849509E-03  2.207976E-02`

### Reading

- all four branches are very close on `prob_2_002`
- no branch shows a gross stiffness or kinematic anomaly on this deck
- `MH6T` and `REZAIEE` are numerically distinct from `SIMOT6` and `MITC6`, but still remain in the same displacement band

## Problem 2-003 summary

### ABS displacement rows

`SIMOT6`

- subcase 1: `4.296845E-04  2.271656E-03  0.0  0.0  0.0  1.283410E-03`
- subcase 2: `0.0  0.0  7.121750E-02  2.151988E-02  2.315182E-02  0.0`

`MITC6`

- subcase 1: `5.117961E-05  8.985788E-04  0.0  0.0  0.0  2.004154E-04`
- subcase 2: `0.0  0.0  1.064637E-01  6.336533E-03  3.438242E-02  0.0`

`MH6T`

- subcase 1: `8.723350E-04  6.362418E-03  0.0  0.0  0.0  3.053693E-03`
- subcase 2: `0.0  0.0  7.275874E-02  2.081868E-03  3.210449E-02  0.0`

`REZAIEE`

- subcase 1: `5.117961E-05  8.985788E-04  0.0  0.0  0.0  2.004154E-04`
- subcase 2: `0.0  0.0  2.440827E-01  5.594281E-02  5.434726E-02  0.0`

### Loaded-node spot check

For the loaded edge nodes `GRID 5` and `GRID 6`:

`SIMOT6`

- subcase 1, grid 5: `4.296845E-04  2.192847E-03  0.0  0.0  0.0  1.283410E-03`
- subcase 1, grid 6: `1.876295E-04  2.271656E-03  0.0  0.0  0.0  1.230420E-03`
- subcase 2, grid 5: `0.0  0.0  7.121750E-02 -2.079928E-02 -2.315182E-02  0.0`
- subcase 2, grid 6: `0.0  0.0  7.004582E-02 -2.151988E-02 -2.257854E-02  0.0`

`MITC6`

- subcase 1, grid 5: `5.117961E-05  8.856538E-04  0.0  0.0  0.0  9.880795E-05`
- subcase 1, grid 6: `1.881765E-05  8.985788E-04  0.0  0.0  0.0  1.977852E-04`
- subcase 2, grid 5: `0.0  0.0  1.047864E-01  7.835649E-04 -3.438242E-02  0.0`
- subcase 2, grid 6: `0.0  0.0  1.064637E-01  7.001298E-04 -3.227312E-02  0.0`

`MH6T`

- subcase 1, grid 5: `8.723350E-04  6.133758E-03  0.0  0.0  0.0  2.919274E-03`
- subcase 1, grid 6: `3.151255E-04  6.362418E-03  0.0  0.0  0.0  3.053693E-03`
- subcase 2, grid 5: `0.0  0.0  7.067975E-02 -1.861569E-03 -3.210449E-02  0.0`
- subcase 2, grid 6: `0.0  0.0  7.275874E-02 -1.953174E-03 -3.180556E-02  0.0`

`REZAIEE`

- subcase 1, grid 5: `5.117961E-05  8.856538E-04  0.0  0.0  0.0  9.880795E-05`
- subcase 1, grid 6: `1.881765E-05  8.985788E-04  0.0  0.0  0.0  1.977852E-04`
- subcase 2, grid 5: `0.0  0.0  2.384165E-01 -5.177489E-02 -2.747098E-02  0.0`
- subcase 2, grid 6: `0.0  0.0  2.386270E-01 -5.594281E-02 -4.482430E-02  0.0`

### Reading

- `prob_2_003` is not a small-perturbation compare like `prob_2_002`; the four formulations split noticeably
- `MITC6` and `REZAIEE` match exactly in subcase 1 for the extracted displacement rows and loaded nodes
- this exact in-plane match is consistent with the Python source design, where `REZAIEE` intentionally keeps the same assumed membrane interpolation as `MITC6`
- `MH6T` is much more flexible than `SIMOT6` and `MITC6` in subcase 1
- `REZAIEE` is dramatically more flexible in subcase 2 than the other three branches
- because these differences are large, this deck should be treated as an audit trigger, not yet as a passed equivalence check

## Immediate conclusion

- formulation plumbing is working for all four `CTRIA6` branches
- `prob_2_002` supports the view that the new ports are numerically healthy at least on a simple smoke deck
- `prob_2_003` shows that `MH6T` still needs deeper formulation-level audit in MYSTRAN
- `prob_2_003` also shows that `REZAIEE` needs a targeted out-of-plane and shear audit, but its in-plane match with `MITC6` is currently consistent with the intended Python formulation

## Direct Python vs MYSTRAN check for `prob_2_003`, `nx=2`

To avoid mixing different harness assumptions, a direct Python tip-displacement check was run on August 15, 2026 using:

- `battle_2_003_curved_plot.py`
- `build_t6(..., alternate_diag=False)`
- the same four triangle formulations:
  - `SIMOT6`
  - `MITC6`
  - `MH6T`
  - `REZAIEE`

The MYSTRAN side used:

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_simot6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mitc6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mh6t.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_rezaiee.dat`

The MYSTRAN decks were checked and they do **not** blanket-lock `RZ` on all nodes. They only clamp the two root nodes:

- `SPC1,1,123456,1`
- `SPC1,1,123456,2`

### Tip comparison

The values below compare Python average tip displacement against the average of `GRID 5` and `GRID 6` from the MYSTRAN `F06` displacement table.

| Formulation | Python `Uy` | MYSTRAN `Uy` | MY/PY | Python `Uz` | MYSTRAN `Uz` | MY/PY |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `SIMOT6` | `1.208532E-02` | `2.232252E-03` | `0.184708` | `4.076211E-01` | `7.063166E-02` | `0.173278` |
| `MITC6` | `1.090429E-03` | `8.921163E-04` | `0.818134` | `2.915596E-01` | `1.056251E-01` | `0.362276` |
| `MH6T` | `5.451729E-02` | `6.248088E-03` | `0.114607` | `4.678699E-01` | `7.171925E-02` | `0.153289` |
| `REZAIEE` | `1.090429E-03` | `8.921163E-04` | `0.818134` | `4.261961E-01` | `2.385218E-01` | `0.559653` |

### Reading

- `MITC6` and `REZAIEE` are both reasonably close to Python in the in-plane `Uy` case
- `SIMOT6` and `MH6T` are much stiffer than Python in the in-plane `Uy` case
- all four MYSTRAN branches remain stiffer than Python in the out-of-plane `Uz` case
- `REZAIEE` is the closest branch to Python in the out-of-plane comparison, but is still substantially stiffer
- this means the large differences seen in `prob_2_003` are not explained by the old blanket-`RZ` bug in the legacy Python harness

## Recommended next audit

1. compare the same `prob_2_003` deck directly against the Python reference outputs for `SIMOT6`, `MITC6`, `MH6T`, and `REZAIEE`
2. isolate whether the divergence is dominated by membrane, bending, or transverse shear
3. verify local frame construction and shear interpolation for `MH6T`
4. compare the out-of-plane `REZAIEE` response directly against the Python reference before treating the large shear-driven split as a bug

## Follow-up on August 15, 2026

Two additional geometry/basis hypotheses were tested in the Fortran `CTRIA6` family:

1. align `SURFACE_BASIS_T6` with Python so that `e1 = g1 / |g1|`
2. align `COV_MAP_T6` with Python so that the covariant-to-physical map uses the pointwise local basis `E1,E2` instead of the centroid fixed frame

These edits were applied to:

- [CTRIA6_SIMO1993.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_SIMO1993.f90)
- [CTRIA6_MITC6.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MITC6.f90)
- [CTRIA6_MH6T.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MH6T.f90)
- [CTRIA6_REZAIEE.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_REZAIEE.f90)

After rebuild and rerun of the four `prob_2_003_nx02_ctria6_*` decks, the reported `GRID 5/6` displacements were unchanged to the printed precision.

### Reading

- the large `prob_2_003` mismatch is **not** explained by the current `SURFACE_BASIS_T6` definition
- it is also **not** explained by whether `COV_MAP_T6` uses a centroid fixed frame or the pointwise local basis
- the next high-probability audit target should therefore move deeper into the membrane and assumed-strain construction itself rather than the local coordinate mapping layer
