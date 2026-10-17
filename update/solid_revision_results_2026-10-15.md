> Subsequent CHEXA8 EAS9 v2 feature update: see [CHEXA8 results](hexa8_eas9_v2_results_2026-10-15.md). Condensed thermal RHS and enhanced center/corner recovery now pass; physical GP-stress Kg passes independent and full SOL 105 checks. Earlier CHEXA8 thermal/corner failures below describe the preceding increment. CPENTA6 failures, solid pressure input and output-format/Nastran limitations remain open.

> Subsequent CPYRAM13 v3 update: see [v3 results](cpyram13_v3_results_2026-10-15.md). Current consistent mass is Jacobi27, and physical GP initial-stress Kg8/Kg27 plus thermal 5-011 pass Python/analytic checks. Earlier mass64 and untested-Kg descriptions below describe the preceding increment. New SAP/Nastran thermal/modal/buckling acceptance remains pending; other listed failures remain active.

# Solid revision results — 2026-10-15

Status: tested stiffness candidates; **not full production certification**.
Base: branch `v18.00.a`, commit `9e904aaf`. No commit or push was made.
Audit/change plan: `solid_revision_plan_2026-10-15.md`.

## Implemented behavior

* `CPYRAM` accepts exactly 5 or 13 GRID IDs in free, small and large fields.
  Five-node `CPYRA` remains a compatibility spelling. Nonstandard names and
  12/14-node connectivity are rejected. CPYRA14/Liu and pyramid EAS54 code
  were removed. Reversed winding permutes physical midside edge pairs too.
* PYRAM13 uses Bedrosian analytic interpolation in collapsed coordinates and
  exactly 2x2 Legendre x 2 Jacobi(alpha=2,beta=0) stiffness quadrature.
  Collapsed Jacobian weighting is accounted for once. No stabilization.
  It is available under both LEGACY and NEWSOLID.
* PYRAM5 stays standard in both modes. LEGACY remains the overall default.
* NEWSOLID PENTA6/HEXA8 use **EAS9 v2**, not the older handoff modes. The
  enhanced natural covariant strains are mapped using the inverse center
  Jacobian, with engineering-shear factors. Nine internal parameters are
  condensed locally: `Kuu - Kua solve(Kaa,Kua^T)`. No global enrichment.
  Failed local solves now give an explicit fatal diagnostic, not a silent
  fallback. PENTA retains the exact tiny regularizer in the v2 reference.
* Local matrices remain 18x18, 24x24 and 39x39. Translational indices map
  into the existing six-component GRID storage. Sparse assembly is unchanged.
* PYRAM13 mass uses independent 4x4x4 quadrature: reusing the eight stiffness
  points produced a deficient mass rank 24/39. Consistent mass now has rank
  39/39 in the tested regular, warped and rotated cases. Modal row-sum mass
  is rejected with a request for `PARAM,COUPMASS,1`.
* Pyramid PSOLID data are copied explicitly to prevent stale material/coordinate
  settings from a preceding element. General material-coordinate parity is
  still NOT VALIDATED.
* A small numeric-locale guard before BDF real conversion addresses the
  observed MinGW decimal-comma input failure. The user's regional setting is
  now en-US with decimal point. This is independent of the element formulas.

## Reference correction and provenance

The original handoff includes `solid3d_cpenta6_eas.py` and
`solid3d_chexa8_eas.py`, whose Cartesian-mode formulations do **not** reproduce
the v2 twisted-beam evidence. Unchanged v2 files were recovered from the
supplied local archive
`nastran_test/6/SAP5004_CPENTA6_CPYRA5_Nastran_Python_Comparison_v1.zip`:

* `solid_5004/solid3d_eas9_frame_v2.py`
* `solid_5004/solid3d_cpyra5_eas_v2.py` (the `natural_to_cartesian` helper)

Python references are retained under `tests/solids/references`; hashes and
sources are in `tests/solids/provenance.json`. No experimental pyramid EAS
or 14-node formulation was ported. No empirical tuning was performed.

## Required v2 gates

These are actual Fortran element matrices, not a separate Fortran reimplementation
inside the tests. 5-001 uses the supplied reconstructed seven-hexa topology,
converted to 14 wedges or 42 quadratic pyramids. The rotation check rotates
both geometry and an independent non-affine nodal force vector.

| Formulation | 5-001 displacement relative error | Rotated-load solution relative error | Scaled free eigenvalue minimum |
|---|---:|---:|---:|
| CHEXA8 EAS9 v2 | 2.43e-16 | 4.92e-16 | 0.5895 |
| CPENTA6 EAS9 v2 | 2.65e-16 | 1.10e-15 | 0.4947 |
| CPYRAM13 | 9.78e-16 | 1.66e-15 | 0.08740 |

Full matrices also agree with the unchanged v2 Python reference on regular,
warped and rotated elements to numerical roundoff. Symmetry, nonnegative
spectrum within roundoff, mass, affine recovery and uniform thermal expansion
were checked separately. `nu=0.499` passes the same kernel screen; this is
**not** the complete 5-008 cylinder benchmark.

### 5-004 reproduction

Ratios below divide displacements by **IN=0.005429, OUT=0.001749**.
They are not Python-minus-Nastran differences.

| Element / mesh | MYSTRAN IN displacement | MYSTRAN OUT displacement | IN ratio | OUT ratio |
|---|---:|---:|---:|---:|
| CPENTA6 v2, 12x2x2 | 0.00516571756 | 0.00166267063 | 0.951504 | 0.950641 |
| CPENTA6 v2, 24x4x2 | 0.00521154655 | 0.00168459013 | **0.959946** | **0.963173** |
| CHEXA8 v2, 12x2x2 | 0.00535674673 | 0.00173048885 | 0.986691 | 0.989416 |
| CHEXA8 v2, 24x4x2 | 0.00539593305 | 0.00174618303 | **0.993909** | **0.998389** |

Full displacement fields match Python within approximately 5.4e-8 relative
error in float32 OP2. The wedge test sets PSOLID IN=2, ISOP=1 explicitly;
LEGACY selective shear integration is unchanged and is a separate option.

Supplied Nastran standard-wedge tip means on 24x4x2 are 0.00519550126 and
0.00167075777. EAS v2 is approximately +0.309% / +0.828% more flexible,
consistent with the supplied archived evidence; these are different formulations.
No matching HEXA8 5-004 Nastran OP2 is supplied, so its Nastran gate is
**NOT VALIDATED**. Its HEXA8 mesh is reconstructed from the same structured
GRID/BC/load model, not asserted to be an exact official SAP2000 mesh.

### CPYRAM13 Nastran displacement gates

5-004 full-field MYSTRAN/Python differences are below 4.4e-8; supplied
Nastran OP2 differences are below 3.1e-6 relative. Maximum tip-mean difference
is approximately 0.000293%, reproducing the handoff.

| 5-007 L1 mesh | Full-field MYSTRAN/Nastran relative error | Gate 1e-4 |
|---|---:|---|
| 4x4 | 1.4337e-4 | **FAIL** |
| 8x8 | 1.3607e-5 | PASS |
| 16x16 | 2.0103e-5 | PASS |
| 24x24 | 8.6851e-6 | PASS |

MYSTRAN/Python full-field agreement is below 5e-7 on these meshes. Python
global residuals reach approximately 3.4e-6 in the thin dome. The coarse
Nastran gate is retained as failed; there was no stabilization, quadrature
retuning or tolerance change to hide it. L2 meshes remain NOT VALIDATED.

## Independent feature acceptance ledger

| Feature | CPYRAM13 | CPENTA6 / CHEXA8 EAS9 v2 |
|---|---|---|
| Stiffness | Python matrix/patch/rotation pass; Nastran 5-004 and larger L1 dome pass; coarse dome gate fails | Python v2 matrix/patch/rotation/5-004 pass; wedge NX comparison reported; HEXA8 NX NOT VALIDATED |
| Consistent mass | Independent quadrature, Python matrix and positive-rank checks pass; Nastran modal/dynamic NOT VALIDATED | Python mass matrices pass; Nastran dynamics NOT VALIDATED |
| Lumped mass | Row-sum has 15 negative translational entries in regular element; modal use explicitly rejected | Independent dynamic/lumping parity NOT VALIDATED |
| Pressure | Base-face PLOAD4 G1/G3 test rejected as unimplemented; no pyramid pressure kernel | Independent face/sign/traction parity NOT VALIDATED |
| Thermal | Uniform expansion/affine test passes; gradients/Nastran NOT VALIDATED | Uniform expansion passes; enhanced gradient RHS **KNOWN_FAIL** |
| Geometric stiffness / buckling | **NOT VALIDATED** | **NOT VALIDATED** |
| Stress / strain | Affine center recovery passes; full corner/coordinate/format parity NOT VALIDATED | Affine center recovery passes; enhanced corner recovery **KNOWN_FAIL** |
| OP2 stress/strain | Center-only result cannot fill standard six-location layout; export suppressed with explicit warning, F06 center retained | Enhanced recovery not accepted for general results |

The negative EAS feature tests independently compare the actual thermal RHS
and corner strain-displacement operators with their Schur-consistent forms:
gradient thermal relative errors are 0.003823 (wedge), 0.030612 (hexa);
corner operator relative errors are 0.45079 and 0.46101. These are documented
failures, **not accepted outputs**, and were not corrected by tuning stiffness.

## KNOWN_FAIL / KNOWN_RISK register

* **KNOWN_FAIL: THIN_DOME_LOCKING** — standard CPYRAM5 remains the baseline;
  the historical 24x24 5-007 normalized displacement is ~0.00077869. EAS18
  ~0.00080589 and nodal smoothing alpha=.8 ~0.00296729 do not cure it.
* **KNOWN_RISK: LOCAL_HOURGLASS** — regular CPYRAM13 eight-point stiffness
  rank 32/39 means six rigid motions plus one integration mode. The warped
  test has rank 33, so rank 32 is not asserted for every geometry. Positive
  tested patch assembly does not certify arbitrary meshes.
* **KNOWN_FAIL: COARSE_DOME_NASTRAN_GATE** — 4x4 L1 exceeds the retained 1e-4
  full-field comparison tolerance. Rank-one stabilization and 27-point
  stiffness remain excluded because of documented bending degradation.
* **KNOWN_FAIL: UNCONDENSED_ENHANCED_THERMAL_RHS** and
  **MISSING_ENHANCED_CORNER_STRAIN** — general EAS thermal/recovery not accepted.
* **KNOWN_FAIL: PYRAM13_ROW_SUM_MODAL_MASS** — negative masses; explicit guard.
  Full consistent mass avoids this deficiency but modal Nastran parity is pending.
* **KNOWN_FAIL: PYRAMID_PLOAD4_FACE_OPTION** — pressure input path is unimplemented.
* **KNOWN_RISK: CENTER_ONLY_PYRAM13_OUTPUT** — F06 center only; incompatible
  OP2 recovery records are suppressed rather than populated with adjacent rows.
* **KNOWN_RISK: DISTORTED_CORNER_JACOBIAN** — the supplied irregular HEXA8 patch
  contains a negative recovery-corner Jacobian despite valid volume points;
  stiffness-only patch success cannot validate corner output.

## Next acceptance work

Preserve the reproduced v2 stiffness. Implement and independently regress
enhanced thermal RHS and recovery, pyramid face pressure, geometric stiffness,
modal/mass parity, output formats/coordinates and mixed-topology transitions.
Supply or generate matched Nastran cases for these features, HEXA8 5-004 and
the complete 5-008 screen. Until then, cells marked NOT VALIDATED/KNOWN_FAIL
must not be promoted to production acceptance.

Build and exact rerun commands are in `tests/solids/README.md`. Executed JSON
tables are adjacent to that README; generated solver files are ignored.
