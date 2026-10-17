# CPENTA6 ABE default / EAS9_V2 alternate — implementation and validation

Status: **INCOMPLETE NATIVE VALIDATION**. Both parameter choices compile and run independently. This is not a claim that all requested benchmark conversions or formulation acceptance gates are complete.

## Implementation

- NEWSOLID CPENTA6 defaults to `PARAM,PENTATYP,ABE`. `PARAM,PENTATYP,EAS9` selects the frozen wedge EAS9_V2. Selection applies to static, modal and buckling phases without automatic switching.
- `PARAM,PENTABUC,...` now produces an explicit removal error. ANSEAS/EASANS are not accepted PENTATYP values. LEGACY ignores PENTATYP with an F06 notice. Invalid new selectors are rejected. The trailing restart field stores PENTATYP; old ANSEAS restart values are rejected.
- ABE implements edge-midpoint tied shear and the exact nine A/B/E modes, including the frozen center covariant transformation and determinant scaling. No projection, empirical stiffness factor or fitted constant was added.
- Both variants use the same condensed state interface for stiffness, enhanced thermal RHS, total strain/stress recovery and physical initial-stress Kg. The local 9x9 block uses LAPACK DPOTRF/DPOTRS with DPOCON monitoring. Nonpositive blocks or reciprocal condition below 1e-12 fail explicitly. Python's trace regularization is not copied into Fortran.
- Kg integrates recovered Gauss-point stress with standard displacement gradients. It is an initial-stress interface, not a derived nonlinear ANS/EAS tangent. Existing CPENTA6 mass and surface-load paths are retained. Other solid dispatch is unchanged.

Relevant source: `Source/Modules/EAS9_FRAME.f90`, `Source/EMG/EMG5/PENTA.f90`, `Source/LK9/L92/ELEM_STRE_STRN_ARRAYS.f90`, `Source/Modules/PARAMS.f90`, `Source/LK1/L1A-BD/BD_PARAM.F90`, and restart/interface imports.

## Frozen references and reproducibility

The latest input is `python-solids/MYSTRAN_CPENTA6_ABE_DEFAULT_CODEX_HANDOFF_20261016_v2_WITH_5001.zip`, extracted separately under `python-solids/cpenta6_handoff_v2/`. README, CODEX_IMPLEMENTATION and VALIDATION_POLICY were read. The ABE candidate and `reference/ABE_CLAUDE_FINAL_ORIGINAL.py` have the same SHA256:

`bca87d8424529e91efa6b0b0031013564b84170305f42b5d3cdfe34a43047dfd`.

Retained deliverables are in `tests/solids/cpenta6_abe_validation/`. `provenance.json` hashes canonical inputs, source, executable and DATs. `results/mystran/manifest.json` identifies the actual F06/OP2 files under `tests/solids/generated/cpenta6_abe/`, with hashes; native files are retained there rather than duplicated. Python matrix/solution arrays are identified in `results/python/manifest.json`. Native DAT pairs are byte-identical after replacing the one PENTATYP line (37 pairs).

Run from `H:\Projects\MYSTRAN` with `.venv\Scripts\python.exe`:

```powershell
C:\gcc\bin\cmake.exe --build source\build -j 12
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_abe.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_body.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_s2k.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_plate.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_buckling.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_pressure.py
.\.venv\Scripts\python.exe source/tests/solids/check_penta6_winkler.py 100
.\.venv\Scripts\python.exe source/tests/solids/report_penta6_abe.py
```

The runners require the supplied v2 handoff and existing reference dependency files. They do not modify original S2Ks or frozen Python sources.

## Results and limits

Accuracy uses independent references; SAP only sets `max(1.2*SAP_error, SAP_error+0.5)` percent tolerance. NEAR adds one percentage point. Parity, validity and accuracy are reported separately. Counts in `summary.json` count metrics, not independent models. An accuracy PASS does not imply acceptance of untested native features.

| Test | Result |
|---|---|
| Production K, thermal RHS, Kg, center/corner stress/strain, regular/distorted/rotated | 30 operator parity checks PASS, relative error <1e-10 on these test geometries |
| Default vs explicit ABE SOL101/SOL105; explicit EAS9; LEGACY isolation; invalid/removed selector | PASS on synthetic column |
| 5-001 v2 MacNeal, 14 wedges, all eight states | Frozen ABE worst u=0.0793928278, stress=0.4621245461: FAIL. EAS9 frozen patch: PASS. Native joint 8/displacement reproduces Python. Native corner recovery: INVALID for nonpositive corner Jacobians in elements 4/9; separate displacement decks retained. Exact frozen interior stress samples are retained in Python CSV and are not misrepresented as native stress parity. |
| 5-002 A/B/C, six load cases each | 36 ABE/EAS9 native solves; exact BODY1 RBE2 matches Python rigid transformation; full displacement parity PASS. Accuracy and REVIEW flags are in the joined CSV. |
| 5-003 A/B/C, IN/OUT | Both variants run on all three original meshes. Full displacement parity PASS. Accuracy failures remain visible. |
| 5-004 IN/OUT | ABE mean 0.00525163278 / 0.00144334067 in: NEAR / FAIL. EAS9 mean 0.00259941044 / 0.000793123089 in: FAIL / FAIL. Both full displacement parity PASS. |
| 5-005 B/E/G/J, UNIF/CONC | Both variants run with original 0.01 in thickness and consistent triangular face-equivalent pressure loads. Displacement parity PASS. Independent acceptance is PROVISIONAL because manual 0.0001 in thickness conflicts with S2K. No rescaling. The triangular consistent RHS differs from historical equal quad-corner loading; results are not expected to be bit-identical to that older load discretization. Native PLOAD4 parity is not asserted. |
| 5-007 A/B dome | Both variants run; full displacement parity PASS; severe locking/accuracy FAIL remains. |
| 5-010, two separately tagged prestress methods | Prescribed uniform stress is a Python screening only. Native SOL105 uses the SOL101 recovered field from a unit resultant top traction. ABE weak 0.8912347 kip: PASS; strong 269.7684 vs 127.91 kip: FAIL. EAS9 weak/strong accuracy FAIL. Native factor parity <2.3e-7, MAC minimum >0.99999994 over 14 matched modes. Axis and signs are checked by physical translation directions. Mode arrays retained. Canonical wedge split differs from the historical converter, so historical EAS9 values are not substituted. |
| 5-011 uniform 100 F | Both variants yield xx=yy=-733.3333 psi, uniform accuracy PASS; displacement parity PASS. Distorted thermal gradient is covered by operator tests, not by a separate full native distorted thermal benchmark. |
| 5-008 C/E | Gauge-tagged native audit runs for both variants; local angular restraints and flat-face pressure resultant are preserved, with displacement parity PASS. Accuracy remains INVALID because original S2K has free axial rigid translation and reference axial stress/plane-strain/output conventions are unresolved. A single axial gauge is documented; the historical all-node axial SPC surrogate is not used. |
| 5-012 B/D | Full native 100x100x1 model (20,402 solid nodes, 20,000 wedges) preserves all five S2K anchors. Exact bilinear face stiffness is represented by an algebraically equivalent CELAS2 network plus a fixed ground node. Pair coefficients are negative off-diagonal representations, not empirical nodal springs; the assembled foundation remains the positive-definite integrated face matrix. Small-mesh matrix reconstruction and native parity pass. Full native displacement parity PASS. ABE B=-0.914253 in: PASS; D=-0.1739468 in: NEAR. EAS9 B=-0.5461586 / D=-0.1062044 in: FAIL. |
| 5-013 | INVALID: native pore pressure/effective stress parity not implemented/verified; buoyancy substitute excluded |

5-004 F06 resultant reaction comparison exceeds 1e-6 in some cases (about 4.6e-6 on a cancellation-sensitive printed sum). It is recorded as `FAIL_F06_PRECISION_UNRESOLVED`, with no tolerance relaxation. Full reaction parity remains open despite displacement parity passing.

## Regression

Build passed with gfortran/OpenBLAS. `check_kernels.py`, `check_kernels.py 0.499`, `check_patch_frame.py`, `check_eas_feature_limits.py`, `check_hexa8_v2.py`, `check_cpyram13_v3.py`, `check_tetra_smooth_selection.py`, and `compare_penta_beam.py` passed. Historical wedge kernel/frame/beam comparisons explicitly select EAS9; ABE is tested separately. The former negative CPENTA6 enhanced-thermal/recovery gate now expects those features to pass. Logs are `tests/solids/abe_*.log`.

The earlier `check_penta6_anseas.py` and `penta6_anseas_results.json` describe the retired single-mode buckling experiment and are historical; they are superseded by `check_penta6_abe.py`. Their success is not claimed for the new parameter interface.

## Open formulation and integration issues

1. ABE fails the distorted/MacNeal affine patch. Baseline is preserved without tuning; the default change is a user selection, not certification.
2. Strong-axis buckling and dome locking remain unacceptable. The canonical topology must stay explicit in comparisons.
3. The 5-008 audit verifies the actual flat-face resultant and local restraint transformation, but the source axial gauge and benchmark stress conventions remain unresolved. Native PLOAD4 is still unverified; equivalent nodal pressure loads do not certify PLOAD4.
4. The full Winkler model now uses exact integrated face stiffness and original anchors. General nonlinear foundation/contact behavior is outside this linear benchmark; no historical anchor modification is used.
5. Native pore pressure and effective-stress location/averaging need implementation for 5-013.
6. Native stress OP2 readback fails on an OES1X1 record-marker mismatch; F06 output also needs element-ID/point mapping audit. Uniform thermal values and production recovery operators do not certify general native stress export.
7. MacNeal corner detJ<=0 is reported INVALID; no extrapolation or hidden geometry repair is applied. Native parity at the four frozen interior stress samples remains unfinished.
8. Full reaction precision/parity remains open. General pressure, modal/mass solver parity and full distorted thermal benchmark are not implied by the checks above.

The requested full port validation is therefore **not declared complete**.
