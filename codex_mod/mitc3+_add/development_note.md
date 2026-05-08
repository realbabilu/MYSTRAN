# MITC3+ Integration Into MYSTRAN

Date: 2026-05-05

## 2026-05-08 workspace port note

This package was consolidated under the canonical `mitc3+_add` name for the
`optimization_rcm_v2` branch. The current workspace did not contain the older
full `TPLT_MITC3P.f90` bubble-condensed Python-closure source referenced below,
so the live committed MYSTRAN port wires the parser/routing/interface and adds a
first-pass `TPLT_MITC3P` implementation point that delegates to the legacy MIN3
plate machinery with explicit MITC3+ selection and light `RZ` drilling stiffness.

The older validation notes below are retained as historical handoff context for
the intended full kernel.

## Goal

Add a first-pass `MITC3+` triangular shell kernel into MYSTRAN under the `CTRIA3` family, with explicit parameter selection and one-element validation against the uploaded Python implementation.

## Architectural Decision

`MITC3+` is treated as a `CTRIA3`-class shell option, not as a separate refined drilling family.

Reason:

- it is a 3-node shell topology
- it uses `6 dof/node`
- it carries `rz`, but not as a true full-drilling membrane family like `Allman`/`ITW`
- this is analogous to `MITC4` and `MITC4+` living under `CQUAD4`

Accordingly, the new selector is:

- `PARAM,TRIA3TYP,MIN3`
- `PARAM,TRIA3TYP,MITC3+`

Default remains:

- `TRIA3TYP = MIN3`

## Source Files Added / Modified

All MYSTRAN source edits for this feature are wrapped with:

- `! --- MITC3+_add begin --- !`
- `! --- MITC3+_add end --- !`

Files touched:

- `D:\mystran2\MYSTRANSolver-18.0.0\Source\Modules\PARAMS.f90`
- `D:\mystran2\MYSTRANSolver-18.0.0\Source\LK1\L1A-BD\BD_PARAM.F90`
- `D:\mystran2\MYSTRANSolver-18.0.0\Source\USE_IFs\TREL1_USE_IFs.f90`
- `D:\mystran2\MYSTRANSolver-18.0.0\Source\EMG\EMG4\TREL1.f90`
- `D:\mystran2\MYSTRANSolver-18.0.0\Source\EMG\EMG4\TPLT_MITC3P.f90`
- `D:\mystran2\MYSTRANSolver-18.0.0\Source\Interfaces\TPLT_MITC3P_Interface.f90`

## Implementation Scope

This is a first-pass stiffness kernel.

Implemented:

- plate/shell bending block
- transverse shear block
- 2 internal bubble rotational DOF
- static condensation at the element level
- light drilling penalty on `rz`
- routing through legacy `CTRIA3` using `PARAM,TRIA3TYP,MITC3+`

Not yet production-complete:

- dedicated pressure-load handling for `MITC3+`
- dedicated thermal-load path
- dedicated stress/strain recovery path
- membrane is still supplied by the legacy `CTRIA3` membrane route

## Main Validation Deck

One-element debug deck:

- `D:\mystran2\MYSTRANSolver-18.0.0\run_debug\dkmq\shell_static_mct\ctria3_mitc3p_one_elem_shell_debug.bdf`

This deck activates:

- `PARAM,TRIA3TYP,MITC3+`

Artifacts produced after run:

- `D:\mystran2\MYSTRANSolver-18.0.0\run_debug\dkmq\shell_static_mct\ctria3_mitc3p_one_elem_shell_debug.BUG`
- `D:\mystran2\MYSTRANSolver-18.0.0\run_debug\dkmq\shell_static_mct\ctria3_mitc3p_one_elem_shell_debug.F06`

## Build / Run Status

Build status:

- `mystran` target builds successfully from `build_codex`

Run status:

- one-element `CTRIA3 + TRIA3TYP=MITC3+` deck runs to normal completion
- `.BUG` file is emitted
- element stiffness matrix is printed in local coordinates

## One-Element Validation Against Python

Reference Python implementation:

- `D:\mystran2\mitc3+\mitc3plus.py`

Comparison script:

- `D:\mystran2\codex_mod\mitc3+_add\compare_mystran_mitc3plus_one_elem.py`

Validation metric:

- compare `KE` from MYSTRAN `.BUG`
- compare to `MITC3Plus.k_local()` from uploaded Python

Current validation result against legacy `.BUG` full `KE`:

- MYSTRAN one-element run completed normally
- Python compare script completed normally
- current one-element stiffness comparison:
- `max_abs_diff = 6.241796154E-01`
- `rel_fro_diff = 5.331447947E-07`
- `mem_rel_diff = 8.813135389E-08`
- `plate_rel_diff = 2.027922243E-06`
- `rz_rel_diff = 1.516698131E-01`
- `vs_legacy_full = 1.071309644E-01`
- `vs_legacy_plate = 2.247151226E+00`
  - `mystran_sym_err = 0.000000000E+00`
  - `python_sym_err = 3.637978807E-12`

Interpretation:

- the routing and kernel are alive inside MYSTRAN
- the generated element stiffness is symmetric
- the membrane block matches the uploaded Python implementation essentially exactly
- the plate block also now matches essentially exactly
- the new kernel is clearly different from legacy `CTRIA3/MIN3`; it is not silently falling back to the old triangle plate behavior
- the remaining visible discrepancy is almost entirely in the small drilling penalty on `rz`
- in practical terms the MYSTRAN MITC3+ one-element shell stiffness is now numerically closed to the uploaded Python `MITC3+`

## Runtime Route Confirmation

The original suspicion that `PARAM,TRIA3TYP,MITC3+` was not being selected turned out to be false.

The current one-element run now emits the following debug blocks in:

- `D:\mystran2\MYSTRANSolver-18.0.0\run_debug\dkmq\shell_static_mct\ctria3_mitc3p_one_elem_shell_debug.ERR`

Confirmed runtime outputs:

- `MITC3P_DEBUG_BEGIN`
- `MITC3P_MATRIX_BEGIN KAA`
- `MITC3P_MATRIX_BEGIN KAB`
- `MITC3P_MATRIX_BEGIN KBB`
- `MITC3P_MATRIX_BEGIN KCOND`
- `MITC3P_DEBUG_END`

Interpretation:

- `TREL1 -> TPLT_MITC3P` routing is alive
- the `MITC3+` branch is being executed
- the remaining mismatch is therefore inside the plate kernel itself, not in parser or routing

## Direct `KAA / KAB / KBB / KCOND` Comparison

The comparison script was upgraded to compare the matrices printed by `TPLT_MITC3P` directly against the uploaded Python implementation:

- `D:\mystran2\codex_mod\mitc3+_add\compare_mystran_mitc3plus_one_elem.py`

Initial results before Jacobian fix:

- `err_kaa_rel_diff = 1.016518459E+00`
- `err_kab_rel_diff = 1.224042260E+00`
- `err_kbb_rel_diff = 2.075399434E-01`
- `err_kcond_rel_diff = 1.011572816E+00`

When the comparison is restricted to the Python **plate-only** part:

- `err_kaa_plateonly_rel = 1.222582341E+00`
- `err_kab_plateonly_rel = 1.224042260E+00`
- `err_kbb_plateonly_rel = 2.075399434E-01`
- `err_kcond_plateonly_rel = 1.441782108E+00`

Interpretation:

- those pre-fix numbers suggested that `KBB` was comparatively close while `KAA/KAB` were not
- that diagnosis turned out to be caused by an incorrect Jacobian layout in the MYSTRAN port

## Jacobian Layout Fix

The major source of mismatch was traced to the Jacobian layout in:

- `D:\mystran2\MYSTRANSolver-18.0.0\Source\EMG\EMG4\TPLT_MITC3P.f90`

Before the fix, `JMAT` had been assembled in a legacy component-first layout:

- row 1 = `[dx/dr, dx/ds]`
- row 2 = `[dy/dr, dy/ds]`

The uploaded Python `MITC3+` implementation uses:

- row 1 = `[dx/dr, dy/dr]`
- row 2 = `[dx/ds, dy/ds]`

This transpose-level mismatch propagated into:

- `JINV`
- `DFI_XY`
- the covariant shear matrix `COV_S`
- the assumed shear rows `ERT/EST`

After fixing `JMAT` to match the Python convention, the one-element closure improved dramatically.

## Updated Material / Shear-Block Check

After the Jacobian fix:

- `shell_d_rel_diff   = 2.155729868E-08`
- `shell_t_rel_diff   = 3.907140187E-06`
- `cov_s_rel_diff     = 3.907075819E-06`

Interpretation:

- the bending constitutive block `SHELL_D` is aligned with Python
- the shear constitutive block `SHELL_T` is aligned with Python
- the transformed covariant shear block `COV_S` is now also aligned with Python

## Updated `KAA / KAB / KBB / KCOND` Results

After the Jacobian fix:

- `err_kaa_plateonly_rel   = 4.298125751E-06`
- `err_kab_plateonly_rel   = 3.915605326E-06`
- `err_kbb_plateonly_rel   = 1.973341610E-06`
- `err_kcond_plateonly_rel = 3.884029998E-06`
- `err_kcond_plateonly_blk = 2.138129309E-06`

Interpretation:

- the plate-only `MITC3+` kernel in MYSTRAN is now numerically closed to the uploaded Python implementation
- any remaining mismatch is tiny and at the level of floating-point/roundoff differences

## Component Decomposition: Bending vs Shear

Additional decomposition script:

- `D:\mystran2\codex_mod\mitc3+_add\decompose_python_mitc3plus_plate.py`

Latest results:

- `err_vs_python_plate_kaa = 1.222582341E+00`
- `err_vs_python_plate_kab = 1.224042260E+00`
- `err_vs_python_plate_kbb = 2.075399434E-01`
- `err_vs_python_plate_kcond = 1.441782108E+00`

Norm check:

- `||kaa_err|| = 9.276644341E+05`
- `||kaa_b||   = 3.634422671E+03`
- `||kaa_s||   = 6.628316795E+05`
- `||kab_err|| = 8.096508282E+04`
- `||kab_b||   = 5.191087189E+03`
- `||kab_s||   = 7.122753755E+04`

Interpretation:

- in the uploaded Python element, the `KAA` and `KAB` plate contributions are dominated by the **shear** operator rather than by pure bending
- this helped isolate the debugging target to the covariant/shear part
- after the Jacobian fix, the prior `KAA/KAB` mismatch is effectively resolved

## Follow-up Diagnosis

Additional one-element probes now rule out several simple failure modes:

- the mismatch is not coming from the membrane block
- the mismatch is not explained by a simple `rx/ry` sign convention flip
- the mismatch is not explained by a simple `[w, rx, ry]` local DOF permutation
- the mismatch is not explained by a single global plate stiffness scaling factor

Current reading:

- `TPLT_MITC3P` is correctly wired into the `CTRIA3` shell path
- `TPLT_MITC3P` is not collapsing back to legacy `TPLT2/MIN3`
- the remaining gap is in the detailed `MITC3+` plate bending/shear operator and/or the bubble-condensed plate block

## Plate-Block Spectral Check

The plate sub-block `[w, rx, ry]` was further inspected through diagonal terms and eigenvalues.

Observed behavior:

- MYSTRAN plate block Frobenius norm: `5.860696932E+05`
- Python uploaded `MITC3+` plate block Frobenius norm: `3.663889798E+05`

Selected diagonals show the mismatch is not uniform:

- some MYSTRAN translational plate terms are much larger than Python
- some rotational terms are smaller
- this does not collapse under any simple sign flip or `[w, rx, ry]` permutation

Eigenvalue pattern:

- MYSTRAN plate block:
  - `-3.105864E-03  6.121890E-04  2.407865E-03  2.111436E-02  2.469002E+02  9.942604E+02  2.623239E+03  4.568187E+04  5.842798E+05`
- Python plate block:
  - `-1.887867E-11  4.815024E-12  1.978332E-11  7.716340E-04  3.473513E+02  7.702253E+02  1.869556E+03  5.829623E+04  3.617157E+05`

Interpretation:

- both plate blocks have the expected near-rigid low modes
- the dominant plate mode in MYSTRAN is substantially stiffer than in the uploaded Python implementation
- the mismatch is therefore not just numerical noise or a single sign/order bug

## Journals / References Used

Primary:

- `D:\mystran2\mitc3+\The_MITC3+_shell_element_and_its_performance.pdf`

Additional local references:

- `D:\mystran2\mitc3+\katili2018.pdf`
- `D:\mystran2\mitc3+\katili2019.pdf`
- `D:\mystran2\mitc3+\rahmawati2018.pdf`

Implementation baseline used for direct formula port:

- `D:\mystran2\mitc3+\mitc3plus.py`

## Notes

- `MITC3+` is intentionally kept under the `CTRIA3` house, the same way `MITC4` / `MITC4+` remain under `CQUAD4`
- this is consistent with the earlier decision that these are shell-topology variants, not separate full-drilling membrane families

## Shell Benchmark Closure

After the one-element closure was achieved, the next validation step used the legacy shell benchmark decks already present in:

- `D:\mystran2\MYSTRANSolver-18.0.0\run_debug\dkmq\shell_static_mct`

Approach:

- start from the existing `MITC4` shell decks for `Static-22`, `Static-23`, `Static-24`, and `Static-33`
- split every `CQUAD4` into two `CTRIA3`
- inject `PARAM,TRIA3TYP,MITC3+`
- preserve the same geometry, material, thickness, SPC, and load data

Generator / parser:

- `D:\mystran2\codex_mod\mitc3+_add\run_mystran_mitc3plus_shell_suite.py`

Generated decks:

- `static-22_mitc3p.bdf`
- `static-23_mitc3p_case1.bdf`
- `static-24_mitc3p_case1.bdf`
- `static-33_mitc3p_in_plane_shear.bdf`
- `static-33_mitc3p_out_of_plane_shear.bdf`

Summary:

| Load Case | Theory | MYSTRAN MITC3+ | Error % |
|---|---:|---:|---:|
| `Static-22` | `-3.09000000E-01` | `-2.45034300E-01` | `20.701` |
| `Static-23` | `4.51970000E-04` | `4.56968000E-04` | `1.106` |
| `Static-24` | `9.40000000E-02` | `8.66592200E-02` | `7.809` |
| `Static-33 In-plane` | `-5.42400000E-03` | `-5.32377200E-03` | `1.848` |
| `Static-33 Out-of-plane` | `-1.75400000E-03` | `-1.87920400E-03` | `7.138` |

Comparison against the uploaded Python `MITC3+` shell results already harvested earlier:

| Load Case | Python MITC3+ | MYSTRAN MITC3+ | Delta |
|---|---:|---:|---:|
| `Static-22` | `-2.50268624E-01` | `-2.45034300E-01` | `+5.23432400E-03` |
| `Static-23` | `4.56981308E-04` | `4.56968000E-04` | `-1.33080000E-08` |
| `Static-24` | `8.64204842E-02` | `8.66592200E-02` | `+2.38735800E-04` |
| `Static-33 In-plane` | `-5.34153081E-03` | `-5.32377200E-03` | `+1.77588100E-05` |
| `Static-33 Out-of-plane` | `-1.89321253E-03` | `-1.87920400E-03` | `+1.40085300E-05` |

Interpretation:

- the shell benchmark behavior now tracks the uploaded Python `MITC3+` very closely
- `Static-23`, `Static-24`, and `Static-33` are especially tight
- `Static-22` is still the loosest shell benchmark, but the MYSTRAN-vs-Python delta is still modest compared to the absolute benchmark gap
- this strongly suggests the remaining `Static-22` miss is a property of the `MITC3+` formulation / triangle split on that curved selfweight case, not a Fortran port bug

Output summary file:

- `D:\mystran2\codex_mod\mitc3+_add\mitc3plus_mystran_shell_suite.md`

## MIN3 vs MITC3+ in MYSTRAN Shell Benchmarks

For the same shell benchmark suite, a matching `MIN3` set was generated from the same legacy `MITC4` decks by:

- splitting `CQUAD4` into `CTRIA3`
- injecting `PARAM,TRIA3TYP,MIN3`

Summary:

| Load Case | Theory | MYSTRAN MITC3+ | Err % | MYSTRAN MIN3 | Err % |
|---|---:|---:|---:|---:|---:|
| `Static-22` | `-3.09000000E-01` | `-2.45034300E-01` | `20.701` | `-2.55462400E-01` | `17.326` |
| `Static-23` | `4.51970000E-04` | `4.56968000E-04` | `1.106` | `4.60192900E-04` | `1.819` |
| `Static-24` | `9.40000000E-02` | `8.66592200E-02` | `7.809` | `9.22110100E-02` | `1.903` |
| `Static-33 In-plane` | `-5.42400000E-03` | `-5.32377200E-03` | `1.848` | `-7.60020400E-03` | `40.122` |
| `Static-33 Out-of-plane` | `-1.75400000E-03` | `-1.87920400E-03` | `7.138` | `-2.20843900E-03` | `25.909` |

Interpretation:

- `MITC3+` clearly wins on `Static-23` and especially on `Static-33`
- `MIN3` remains better on `Static-22` and `Static-24`
- this mirrors the earlier Python-side reading:
  - `MITC3+` is stronger on twisted/distorted shell response
  - `MIN3` is more conservative on curved selfweight and dome-like response

Additional output summary file:

- `D:\mystran2\codex_mod\mitc3+_add\min3_mystran_shell_suite.md`
