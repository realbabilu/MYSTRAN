# NEWSOLID formulation selection, 2026-10-16

User-selected formulations:

| BDF element | NEWSOLID formulation | Current acceptance |
| --- | --- | --- |
| CHEXA8 | EAS9 v2 covariant center-frame transformation and local Schur condensation | Previous independent operators and SOL 101/103/105 manufactured checks pass; matched Nastran features, pressure and stress OP2 limits remain open |
| CPENTA6 | EAS9 v2 from solid3d_eas9_frame_v2.py, including the same local regularization | Stiffness reproduction, patch/frame and 5-004 regression pass; buckling NOT FINAL PASS; general enhanced thermal/recovery remain documented failures |
| CTETRA4 | Nodal-smoothed blend alpha=0.9 | Full solver static agreement with supplied Python passes; buckling NOT FINAL PASS |
| CPYRAM5 | Standard five-node pyramid | Default preserved; buckling NOT FINAL PASS |

Use `PARAM,SOLIDTYP,NEWSOLID`. The default SOLIDTYP remains LEGACY. Connectivity and global GRID DOF architecture are unchanged. Nodal smoothing adds patch stiffness couplings between existing GRID translational DOFs; it introduces no new DOFs.

## Archive audit

README_V4.md, CODEX_NEXT.md, README_V5.md and v4 comparison.csv were read completely. The supplied v4/v5 runners execute standard formulations, not CPENTA6 EAS9 v2, CHEXA8 EAS9 v2 or nodal-smoothed CTETRA4. Consequently their standard-element results must not be labeled as v2/smoothed results.

V4 standard CTETRA4 prescribed-prestress screen: weak 68.479219 kip, strong 271.933752 kip. Standard CPYRAM5: weak 35.314902 kip, strong 273.988511 kip, far above the stated Euler weak/strong references 0.888/127.91 kip. No executed CPYRAM5 EAS/smoothing comparison in these packages establishes a best five-node variant. Preserve standard CPYRAM5 and its locking warning; do not replace it or tune K to fit this screen. CTETRA10/CPENTA15 are much closer to Euler in the supplied standard-family screen, but they do not replace the user's selected linear-element formulations.

V4 includes `solid3d_ctetra4_smooth.py` with the alpha=0.9 blend reference. Its source describes a diagnostic linear/eigen trial, not final buckling certification. V5 includes the correct EAS9 v2 reference but explicitly executes no EAS9 computations. Both packages prescribe uniform prestress for 5-010; the actual SAP applied load is absent, so they do not establish two-stage SAP parity.

## Implementation and validation

The existing CTETRA4 smoother already used ALPHA=0.9 and the NEWSOLID selector. CPENTA6 already used the correct natural-to-Cartesian EAS9 v2 transform and local Kaa solve, including Python's tiny trace-scaled diagonal regularization. No stiffness tuning or alternate EAS source was introduced.

Audit identified an unconditional smoothing assembly call in ESP during buckling step 2. The patch correction alpha*(K_smooth-K_standard) is an elastic K term and must not be added to prestress Kg. `CTETRA4S_SMOOTH_ASSEMBLY` now returns before mesh/patch assembly when SOL_NAME=BUCKLING and LOAD_ISTEP>1, for COUNT and ADD alike. This prevents an invalid mixing of elastic stiffness into the geometric matrix; it does not certify the remaining smoothed stress/Kg formulation.

Actual full MYSTRAN SOL 101 runs use the supplied v4 5-010 tetrahedral mesh (78 nodes, 144 tetrahedra, 224 free translational DOFs) and an independent random nodal load. Comparison to unchanged supplied Python:

| Mode | Maximum displacement relative error |
| --- | ---: |
| LEGACY standard CTETRA4 | 4.840e-8 |
| NEWSOLID nodal-smoothed alpha=0.9 | 3.210e-8 |

These errors include float32 OP2 output. A production-linked regression probe verifies that both smoothing COUNT/ADD return without touching mesh arrays during buckling step 2. Existing solid kernel checks were rerun to verify the CPENTA6 v2 selection independently of the smoothing path. The existing 5-004 CPENTA regression was rerun separately.

## Retained limits

Buckling remains NOT FINAL PASS for CTETRA4, CPENTA6 EAS9 v2 and standard CPYRAM5. CTETRA4 smoothed stress recovery, thermal consistency, geometric-stiffness derivation and matched Nastran/SAP two-stage validation remain pending. The current patch smoother leaves multi-material patches unchanged; heterogeneous interfaces and nonbasic GRID output coordinates are not certified by this homogeneous test. CPENTA6 general enhanced thermal RHS and corner recovery remain KNOWN_FAIL as previously measured. CPYRAM5 locking remains KNOWN_FAIL. No passing stiffness/static result promotes those features to acceptance.

Reproduce: `tests/solids/check_tetra_smooth_selection.py`, `check_kernels.py`, `compare_penta_beam.py`. Results: `tests/solids/tetra_smooth_selection_results.json`. V4 references and SHA256 provenance are preserved under `tests/solids/references/multitype_v4/`; the v5 package is extracted read-only for inspection under the workspace python-solids directory.
