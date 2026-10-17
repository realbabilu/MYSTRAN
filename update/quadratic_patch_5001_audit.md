# NEWSOLID CPENTA15 / CHEXA20 — 5-001 quadratic audit

## Result

**Both quadratic elements pass** the supplied FIXED_QUADRATIC geometry and
prescribed-displacement patch in Python and the complete MYSTRAN SOL 101
pipeline, after explicitly reconciling node ordering. Manual loading, six
independent constant-strain states and one rigid-body rotation were tested.
Every displacement, Gauss-stress and GPSTRESS error is below **1e-7**.

| Complete MYSTRAN solve: worst over eight states | CPENTA15 | CHEXA20 |
|---|---:|---:|
| Relative nodal displacement error | 1.13722e-10 | 3.46945e-15 |
| Relative Gauss stress error | 8.37835e-10 | 4.36557e-14 |
| Relative GPSTRESS error | 4.17345e-10 | 1.16043e-14 |
| Fatal/error diagnostics in isolated quadratic decks | 0 | 0 |
| Minimum integration-point detJ | 0.0124390601 | 0.00248365932 |

This supersedes the earlier quadratic patch failures in
`gpstress_volume_5001.md`. In the complete eight-solid deck, all quadratic
groups now pass; CPENTA6 remains outside this repair and still fails.

## Input and reference provenance

Primary input, retained unchanged:
`nastran_test/mystran/5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH_FIXED_QUADRATIC.dat`.
SHA256: `2AC0CBDA32170D10AA4951FC601E3B837DDD2FDFFB44159CE2F7B1FE4CE565E3`.

Python references used:

- `python-solids/solid3d_patch_test_macneal.py` and its report: prescribed
  affine-field method, six unit strain states and rigid rotation. The original
  Python builder uses X2=0.626; the user's latest instruction instead selects
  the supplied DAT with **X2=0.826**, so this audit keeps the DAT geometry.
- `python-solids/handoff/python/reference_solids/solid3d_cpenta15.py`, SHA256
  `D6C5602C899FBC2BC4AEC30E73F26E6E0B176194190467A2D4CDFCC2CAF931A7`.
- `python-solids/handoff/python/reference_solids/solid3d_chexa20.py`, SHA256
  `E9BB1E07E11B98F251798A816EDFAD696D3A02EAB354092D792D740E71B43592`.

The FIXED DAT's CHEXA parent lines contain too many free fields for normal BDF
parsing. A separate formatting-only normalized copy retains every intended
node. The Python audit then reads all 20 CHEXA nodes correctly. User input
files were not edited.

Material is E=1e6 psi, nu=0.25. CPENTA15 uses 14 elements / 56 grids;
CHEXA20 uses seven elements / 48 grids. Physical connectivity, GRID IDs,
coordinates, material and exterior displacement constraints come from the
FIXED DAT. No midside coordinate was moved and X2 remains 0.826.

## Node and DOF reconciliation

The supplied Python kernels/FIXED DAT use bottom, top, vertical midsides.
MYSTRAN's existing shape functions use bottom, vertical, top midsides. The
explicit permutations below map FIXED/Python positions to MYSTRAN positions:

- CPENTA15: `[1..9, 13..15, 10..12]`.
- CHEXA20: `[1..12, 17..20, 13..16]`.

Matrices are permuted in both row and column directions before comparison.
Three translation DOFs per grid are compared; existing six-slot element
storage uses slots 1,2,3 of each grid. No rotation or internal global DOF was
introduced. `node_dof_ordering.csv` records every element, physical GRID ID,
node position and DOF index used in the mapping.

`5-001_FIXED_QUADRATIC_MYSTRAN_ORDER.dat` is an explicitly named solver-order
equivalent, not a silent overwrite of the FIXED DAT. Its physical edge-node
assignments are the same as those tested by Python. Running the literal
FIXED DAT without normalizing its format and reconciling its convention is
not claimed to pass MYSTRAN.

## Python first

Python assembled the full reference element matrices and solved interior
translations with the supplied SPC1/SPCD boundary values. Every exterior
translation, including midside translations, was checked against the exact
field evaluated at that grid's own coordinate. Interior translations were
solved, not imposed. There are 102 free translations for CPENTA15 and 72 for
CHEXA20; the constrained matrices are positive definite.

| Python worst relative error | CPENTA15 | CHEXA20 |
|---|---:|---:|
| Displacement | 2.71051e-15 | 2.71051e-15 |
| Gauss stress | 4.15884e-14 | 2.97860e-14 |

The eight states are the manual tensor, exx, eyy, ezz, gxy, gyz, gxz and
rotation about z. Unit strain amplitudes are 0.001; engineering shear is
split equally between the off-diagonal displacement-gradient entries.
Rigid rotation uses an antisymmetric gradient and expects zero stress.

Errors use maximum absolute component error divided by the maximum exact
field component. For zero-stress rotation, stress is scaled by E times the
rotation amplitude, rather than dividing by zero. Component/nodal CSVs
retain absolute errors as well.

## Root cause and isolation

With reduced shear disabled, Fortran element matrices already agree with
the reference full-integration matrices and pass the patch. Shape functions,
node mapping and the basic isoparametric Jacobian were not the cause of the
quadratic solution failure.

The failing high-order selective-integration path overwrote the first entries
of the high-order `B` and `DETJ` arrays with low-order shear samples, then
integrated those arrays using high-order indices and weights. Thus some
normal-strain contributions also received the wrong Jacobians, while shear
samples were mixed between incompatible rules.

| Matching selective-integration matrix comparison | CPENTA15 | CHEXA20 |
|---|---:|---:|
| Worst relative element matrix difference before repair | 0.412532 | 0.369437 |
| Worst relative element matrix difference after repair | 1.08914e-9 | 5.20924e-15 |

Before repair, the kernel-assembled patch displacement errors reached
0.112124 for CPENTA15 and 0.0793202 for CHEXA20; stress errors reached
0.705828 and 1.14289 respectively.

The comparison distinguishes full integration from selective integration:
their non-affine element matrices are intentionally different. The selective
Fortran matrix is compared against an independent Python selective-rule
matrix, not required to equal the full-rule matrix. Both rules preserve the
constant-strain patch after the repair.

## Implementation correction

Only NEWSOLID high-order quadratic paths were changed:

- `Source/EMG/EMG5/PENTA.f90`: retain the complete 7-triangle x 3-line B/DETJ
  data. Replace only the pure shear energy integral using a separate
  3-triangle x 2-line rule with its own Jacobians and weights.
- `Source/EMG/EMG5/HEXA.f90`: retain the complete 3x3x3 data and integrate the
  shear block separately on 2x2x2 points. This keeps normal/shear sampling
  and weights consistent. Direct/shear constitutive cross terms stay on the
  full rule.
- NEWSOLID CHEXA20 corner stress/strain output is extrapolated from its valid
  3x3x3 integration field using tensor-product Lagrange interpolation. It no
  longer attempts to invert a negative geometric corner Jacobian merely to
  write recovered corner output. Center recovery remains a direct evaluation.

This preserves standard displacement-based quadratic elements and the
existing selective-integration option. It does not replace them with an
enhanced or stabilized formulation. No coordinate, load, constraint, assembly,
global-DOF, shell, beam, or linear-solid formulation code was changed.

The straight-edge quadratic geometry has the same physical corner-Jacobian
sign issue as the underlying distorted linear hexahedron at one corner.
All integration points used here have positive detJ, and no invalid corner
Jacobian is inverted in the quadratic recovery path. This is not a claim
that detJ is positive everywhere in the entire natural-coordinate domain.

## Shape, derivative and quadrature verification

483 element/Gauss locations were compared against the actual Fortran shape
functions, derivatives, B matrices and Jacobians after node permutation.

- Maximum partition-of-unity error: 2.23e-16 for CPENTA15 and 4.45e-16 for
  CHEXA20.
- Fortran versus Python derivative differences: <=1.12e-16.
- Python analytic derivatives were independently checked with a five-point
  numerical derivative.
- Line-rule weights sum to 2.000000000000001.
- The current Fortran seven-point triangle coefficients sum to
  0.4999999995 instead of 0.5. Their decimal truncation explains the remaining
  approximately 1e-9 CPENTA15 matrix discrepancy. This is below the requested
  tolerance; the shared global quadrature routine was not changed.
- Both full and selective matrices have exactly six element rigid-body
  modes. No unintentional element zero modes were found.
- Minimum constrained global eigenvalues for the corrected selective rule
  are approximately 218016.67 (CPENTA15) and 237421.10 (CHEXA20), both positive.
  Matrix symmetry is retained to roundoff.

## Complete solver, assembly, constraints and recovery

The complete SOL 101 runs use generated single-family decks with the exact
physical meshes/material/boundary data from FIXED, the declared MYSTRAN node
permutation, and the normal global assembly/constraint routines. They pass
with no ERROR/FATAL messages for all eight states per family.

The audit reads actual double-precision solved displacements from LINK5A,
retained using an audit-local MYSTRAN.INI beside an audit executable copy.
It does not substitute ideal nodal values or differentiate rounded OP2
displacements. Actual solved displacements are applied to the measured
Fortran Gauss-point B matrices and constitutive tensor, and compared against
the Python kernels. F06/OGS1 recovered tensors are also checked independently.

Exact manual-field joint 8 displacement is
`[0.0008885, 0.0011785, 0.0011570]` inches. The existing Python report prints
these rounded as `[0.00089, 0.00118, 0.00116]`. Both solvers reproduce the
unrounded affine field, not an adjusted rounded target.

## Eight-solid deck and remaining scope

The full solver-order equivalent of FIXED now gives analytical GPSTRESS
passes for volumes 1,2,4,5,6,7,8. In particular all four quadratic groups
CTETRA10, CPENTA15, CHEXA20 and CPYRAM13 pass. CPENTA6, volume 3, retains its
previous approximately 0.512528 maximum relative component error.

The combined eight-solid run also retains the existing CHEXA8 corner-output
diagnostic. Those linear-element paths were intentionally not modified by
this quadratic repair. Use the isolated CPENTA15/CHEXA20 decks for the clean
zero-error acceptance runs; do not interpret the combined deck as an all-eight
acceptance pass.

## Artifacts and reproduction

Source audit: `source/tests/quadratic_patch/`.
User run package: `nastran_test/mystran/quadratic_patch_audit/`.

```powershell
H:/Projects/MYSTRAN/.venv/Scripts/python.exe source/tests/quadratic_patch/run_audit.py
```

The audit exports:

- `python_patch.csv`, `python_shape_checks.csv`;
- `fortran_matrices_before.csv`, `fortran_matrices_after.csv` and saved raw
  element-matrix pairs (`.npz`) in the source audit directory;
- `fortran_patch_before.csv`, `fortran_patch_after.csv`;
- `mystran_solver_patch.csv`, `mystran_solver_nodal.csv`,
  `mystran_solver_gauss.csv`;
- `node_dof_ordering.csv`, `boundary_constraints.csv`, `quadrature_weights.csv`.

`run_mystran.bat` in the user run package executes the two clean manual
quadratic decks. Six-strain and rotation decks are also reproducible through
the audit runner. The main rebuilt executable is `source/Binaries/mystran.exe`.

Regression results: all ten GPSTRESS tests pass; the existing solid-kernel
suite passes for regular, distorted and rotated low-order solids in LEGACY
and NEWSOLID, including its mass/thermal checks. No commits were pushed.
