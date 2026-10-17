# CHEXA8 NEWSOLID EAS9 v2 audit and results, 2026-10-15

## Audit and incremental implementation

The current CHEXA8 NEWSOLID stiffness already matched the supplied `solid3d_eas9_frame_v2.py`: nine covariant natural strain modes transformed by the center Jacobian, detJ0/detJ scaling, and element-local Schur condensation Kuu-Kua inv(Kaa) Kau. The audit found that thermal RHS, corner recovery and differential stiffness were still conventional element paths. These were not accepted as complete v2 behavior.

The implementation retains that validated stiffness and applies the same local equations to the remaining linear operators:

1. Integrate Fu = integral(B^T D eth) and Fa = integral(M^T D eth), then use Fth = Fu-Kua inv(Kaa) Fa.
2. Recover local enhanced parameters eta = inv(Kaa)(Fa-Kau u), total strain = B u+M eta, and stress = D(total strain-eth). Nodal temperature is interpolated at each integration/recovery point. Existing output-axis transformations run after this recovery.
3. Update displacement-only BE/SE center/corner operators to B-M inv(Kaa) Kau. Full recovery also includes the thermal enhanced contribution, which cannot be represented by a displacement-only operator.
4. Replace only CHEXA8 NEWSOLID Kg by the physical initial-stress integral grad(Ni)^T sigma grad(Nj) times I3, using enhanced stresses at each of eight integration points from the actual static displacement and thermal subcase. Compression is negative.
5. Verify independent operators and production ELMDIS/OPT(6)/ELEM_STRE recovery, then run patch/frame 5-001, twisted beam 5-004, and full SOL 101/103/105 cases.

No global enrichment DOFs, connectivity changes, stabilization or tuning were introduced. CHEXA8 LEGACY, CHEXA20, CPENTA6 and CPYRAM formulations are preserved by explicit guards. NEWSOLID CHEXA8 requires PSOLID IN=2 (default for linear CHEXA8); IN=3 now fails explicitly rather than silently selecting a conventional non-v2 element.

`EAS9_FRAME.f90` contains shared Cartesian mapping and the v2 state operator. Existing production HEXA stiffness is retained; tests also compare it with the shared operator, which uses an independent multi-RHS local solve. Both are checked against unchanged Python v2.

## Independent results

Regular, warped and rotated geometries were tested for mechanical, uniform thermal and gradient thermal cases at nu=.22 and .499. K, consistent mass, condensed thermal RHS, initial-stress Kg and all nine center/corner stress/strain locations were compared independently. Maximum relative error across these operator checks: **1.494e-13**, below the 1e-9 gate. Tests require stiffness rank 18 (six rigid modes) and mass rank 24.

* Gradient thermal RHS error: 3.314e-16; enhanced corner strain operator error: 5.471e-17. Earlier CHEXA8 thermal/corner KNOWN_FAIL entries are resolved by these tests. CPENTA6's corresponding failures remain open.
* Geometry/displacement/temperature frame covariance, compression/tension sign and translation-null Kg checks pass.
* Batch tests exercise actual HEXA OPT(6), ELMDIS and ELEM_STRE_STRN_ARRAYS; they do not validate only a detached helper.
* Consistent mass uses ordinary mechanical interpolation with no enhanced inertia. Mass is compared separately from K and Kg; no local condensation is applied to mechanical mass.
* 5-001 patch and geometry/load rotation tests remain unchanged and pass.
* 5-004 full MYSTRAN CHEXA8 runs reproduce Python v2 on 12x2x2 and 24x4x2 meshes. Fine-mesh displacement/reference ratios remain 0.993909 and 0.998389. These are ratios to the benchmark reference, not Python/Nastran differences. No matching CHEXA8 Nastran OP2 was supplied.

## Full solver checks

| Check | Result |
| --- | --- |
| SOL 101 uniform free thermal expansion | OP2 displacement error 2.526e-12; F06 center stress magnitude 4.337e-14 |
| SOL 103 consistent modal mass | First four eigenvalues match independent K/M eigenproblem within 1.005e-7 relative printed precision |
| SOL 105 actual loaded static solution + buckling | First factors 2060.925, 2060.925, 18151.3, 18151.3 |
| Independent 10-element static/Kg/eigen assembly | First critical load 2060.924549; free-DOF static residual 6.045e-15 |
| Independent 20-element static/Kg/eigen assembly | First critical load 2050.359843; free-DOF static residual 1.219e-14 |

The square cantilever has length 10, area 1, E=1e6 and nu=.25, with a unit compressive end load. The Euler comparison is 2056.167584. The 10/20-element ratios are 1.002314/0.997175. These are manufactured verification cases, not substitutes for matched SAP/Nastran acceptance.

## KNOWN_FAIL / KNOWN_RISK

* **KNOWN_FAIL: SOLID_PLOAD4_G1_G3_INPUT**. A separate full BDF pressure check still fails in the existing global BD_PLOAD4 parser with `CODE NOT WRITTEN YET FOR PLOAD4 OPTION USING G1 AND G3/G4`. Pressure was not modified or marked accepted based on other passing operators.
* **KNOWN_FAIL: CHEXA_STRESS_REQUEST_OP2_READBACK**. The full thermal stress-output OP2 cannot be read by installed pyNastran (FatalError after OES1X1). The same test also fails with LEGACY, confirming this is not introduced by EAS recovery. Displacement-only OP2 and F06 center/corner stress checks are separate. Full stress/strain OP2 and Femap/Patran output compatibility remain pending.
* **KNOWN_RISK: NASTRAN_FEATURE_PARITY_PENDING**. Matching Nastran CHEXA8 EAS9 v2 thermal, modal, buckling and 5-004 outputs are not supplied. The independent and full MYSTRAN checks above pass; blanket Nastran acceptance is not claimed.
* **KNOWN_RISK: LINEAR_INITIAL_STRESS_SCOPE**. Kg is the conventional physical initial-stress operator based on enhanced small-strain prestress. This does not establish a complete finite-strain nonlinear EAS material/geometric tangent or follower-pressure tangent. No such tangent is claimed by the supplied v2 linear stiffness reference.
* **KNOWN_RISK: CORNER_JACOBIAN**. Distorted elements may have valid Gauss Jacobians but invalid recovery-corner Jacobians; these fail explicitly. Passing patch stiffness does not certify every recovery corner on arbitrary meshes.
* Retained CPENTA6 thermal/corner failures and CPYRAM pressure/hourglass/coarse-dome/modal-output limits are unchanged; see their result reports.

Reproduce with the project Python: `tests/solids/check_kernels.py`, `check_eas_feature_limits.py`, `check_hexa8_v2.py`, `check_patch_frame.py` and `compare_hexa_beam.py` after a normal CMake build. Machine-readable details: `tests/solids/hexa8_v2_results.json`, `eas_feature_limits.json`, `hexa_beam_results.json`. Generated BDF/probe/output files stay ignored.
