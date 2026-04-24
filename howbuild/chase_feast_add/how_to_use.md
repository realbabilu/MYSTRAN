# CHASE/FEAST - How To Use

This note is for runtime usage after build integration.

## 1) Pick one example deck

Pass-focused quick set (Eigen8):
- `examples/midas_eigen08_cbar_axis_debug.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbar_axis_debug_feast.dat`
- `examples/midas_eigen08_cbar_axis_debug_chase.dat`
- `examples/midas_eigen08_cbeam_axis_debug.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbeam_axis_debug_feast.dat`
- `examples/midas_eigen08_cbeam_axis_debug_chase.dat`
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli_feast.dat`
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli_chase.dat`

Validation-focused shell set (Eigen13):
- `examples/eigen13_arpack.dat`
- `examples/eigen13_chase.dat`
- `examples/eigen13_feast_true.dat`

## 2) How to switch method in deck
Use `PARAM,LANCMETH,...`:
- `PARAM,LANCMETH,ARPACK` (default/stable path)
- `PARAM,LANCMETH,FEAST`
- `PARAM,LANCMETH,CHASE`

Important:
- `LANCMETH` is used in the **LANCZOS path** (`EIGRL`), not in `EIGR/MGIV`.
- If your deck still uses `EIGR,...,MGIV`, summary will remain MGIV and FEAST/CHASE dispatch is not exercised.
- New guard message confirms this explicitly:
  - `*WARNING 4910: PARAM LANCMETH=... IS IGNORED FOR EIG METHOD "MGIV ..."`

Recommended solver pairing for modal runs:
- `PARAM,SOLLIB,SPARSE`
- `PARAM,SPARSEFLAVOR,SUPERLU`

## 3) Run command
Example:

```powershell
E:\mystran17\mystran\Binaries\mystran.exe E:\mystran17\mystran\howbuild\chase_feast_add\examples\midas_eigen08_cbeam_axis_debug_feast.dat
```

## 4) What to check
- `*.ERR`: warnings/fallback info
- `*.F06`: eigenvalue table / frequencies
- Compare FEAST/CHASE results against ARPACK baseline deck with same model.

## 5) Current behavior note
`LANCMETH=FEAST/CHASE` is integrated as experimental path with guarded behavior in LINK4 flow.  
This package includes Eigen8 (quick) and Eigen13 (stress validation) examples.

Current fallback behavior:
- `LANCMETH=CHASE`:
  - if external backend is **not linked**, fallback is **ARPACK LANCZOS** (`*WARNING 4912`).
  - if external backend is linked but native call fails, fallback is **ARPACK LANCZOS** (`*WARNING 4917`).
- `LANCMETH=FEAST`:
  - if external backend is **not linked** or native solve fails, fallback is **ARPACK LANCZOS** (`*WARNING 4911/4916`).

## 6) FEAST caution and Eigen13 reference

FEAST is sensitive to interval and matrix conditioning. Use this checklist:
1. Ensure interval deck (`EIGRL,V1,V2`) is present for FEAST native path.
2. Confirm FEAST native actually runs (`WARNING 4918` in `.ERR`).
3. Confirm it did not fallback (`WARNING 4916` absent).

Reference validation result:
- `validation_eigen13_proof.md`
- `lancmeth_compare_eigen13.csv`
