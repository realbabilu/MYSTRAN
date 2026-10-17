# CPYRAM13 thermal/buckling v3 results, 2026-10-15

The supplied v3 ZIP is implemented incrementally in the shared PYRAM13 kernel used by both LEGACY and NEWSOLID. CPYRAM5, EAS9 and other element formulations were not changed by this increment. Standard 13-node connectivity and the 39 translational element DOFs remain unchanged; no global enrichment DOFs were added.

## Operators

* K remains the previously validated 8-point Bedrosian stiffness, with no stabilization or parameter tuning.
* Consistent mass now uses independent **27-point Jacobi quadrature**, replacing the earlier 64-point Legendre choice. Regular, warped and rotated matrices agree with v3 Python within 4e-16 relative error and have rank 39. Consistent mass is exercised in MODES with `PARAM,COUPMASS,1`; the existing MYSTRAN mass expander uses row sums for static solutions. Indefinite row-sum modal mass remains explicitly rejected.
* Thermal load remains the direct 8-point integral of B-transpose D alpha (N T - Tref). Uniform and nonuniform nodal temperature tests agree with Python. This verifies the thermal RHS separately from recovery.
* Integration-point recovery evaluates sigma = D (B u - alpha (N T - Tref)) with engineering shear strains. Its 8/27-point results agree with Python independently of the thermal RHS.
* CPYRAM13 Kg now integrates grad(Ni)-transpose sigma grad(Nj) times I3. Compression has negative physical stress. The former center-stress rotation-only formulation remains untouched for CPYRAM5, but is no longer used for CPYRAM13.
* Default CPYRAM13 Kg uses 8 points. `PSOLID` IN=3 (for example `PSOLID,1,1,,3`) selects 27-point Kg. IN=2 or default selects 8. **K remains 8 points for either setting**; this fixed-K8 behavior is specific to this formulation and must not be interpreted as general Nastran integration-order equivalence. Mass remains 27 points independently.
* Production `PYRAM OPT(6)` uses `ELMDIS` and the actual thermal subcase to recover stress at each Kg point. Tests exercise this path and compare KED against independently integrated GP stress for every batch element. Compression/tension sign, translation null modes and rotation covariance are checked separately (Kg rotation error 5.1e-16).

## Thermal 5-011

Supplied S2K-derived mesh: 600 CPYRAM13, 1703 GRID points, 4629 free translational DOFs. E=1e6 psi, nu=.25, alpha=5.5e-6, temperature rise=100. Reference basic-axis stress is (-733.333333333, -733.333333333, 0, 0, 0, 0) psi.

Assembly of actual Fortran operators and SciPy global solution gives:

| Independent check | Result |
| --- | ---: |
| Maximum displacement error vs analytic, in | 2.494e-17 |
| Maximum in-plane normal stress error, psi | 5.230e-12 |
| Maximum out-of-plane normal stress, psi | 1.056e-11 |
| Maximum shear stress, psi | 4.626e-12 |
| Relative free-DOF equilibrium residual | 1.247e-15 |

A separate full MYSTRAN SOL 101 BDF run passes parsing, thermal assembly, global solve and output. Maximum OP2 displacement error is 2.437e-11 in (float32 output); F06 center principal stress error is 3.334e-5 psi (printed precision). Principal stresses are used because local element axes differ across pyramids. This output check does not replace the integration-point checks above.

## Buckling 5-010: screening only

Supplied mesh: 144 CPYRAM13, 457 GRID points, 1352 free translational DOFs. E=3600 ksi, nu=.2. Prescribed uniform nominal sigma_zz=-1/12 kip/in2 represents a 1 kip compression reference.

| Kg points | Weak load, kip | Strong load, kip |
| --- | ---: | ---: |
| 8 | 0.8897927712 | 126.5347054 |
| 27 | 0.8897927763 | 126.5347041 |

Largest relative discrepancy against supplied Python eigenvalues is 2.46e-8. K8 is unchanged between these calculations. The generalized eigenproblem is solved using G=-Kg, consistently with compression-negative stress. These are assembled Fortran operators with a SciPy eigensolve, **not a full MYSTRAN SOL 105 run**.

## Regressions and acceptance limits

Build succeeds. Existing kernel checks and 5-001 patch/frame tests pass. CPYRAM13 twisted beam 5-004 passes on 12x2x2 and 24x4x2 meshes against Python and supplied Nastran OP2; largest full-field relative differences are 4.39e-8 and 3.01e-6 respectively.

* **KNOWN_RISK: SCREEN_ONLY_BUCKLING**: the supplied S2K has no applied compression load. A real loaded static solution followed by stress recovery and full solver buckling is still required for SAP/Nastran certification. Uniform assumed prestress does not establish that parity.
* **KNOWN_RISK: V3_NASTRAN_PARITY_PENDING**: new mass/modal, thermal 5-011 and buckling 5-010 do not have matched Nastran output in this ZIP. Python and analytic checks pass, but these features are not accepted as Nastran-validated.
* **KNOWN_RISK: GENERAL_CENTER_THERMAL_RECOVERY**: new integration-point recovery uses N*T correctly, including gradients. Existing center F06 thermal recovery still uses its historical nodal average; nonuniform-temperature center output is not certified by the uniform 5-011 check.
* **KNOWN_RISK: CENTER_ONLY_PYRAM13_OUTPUT**: GP recovery is available internally for Kg/tests; F06 remains center-only and incompatible OP2 stress/strain output remains suppressed. No corner/format certification is claimed.
* **KNOWN_FAIL: PYRAMID_PLOAD4_FACE_OPTION** remains unimplemented and independently documented in the earlier negative card test. Pressure was not altered or promoted by thermal/Kg success.
* Existing **LOCAL_HOURGLASS**, **COARSE_DOME_NASTRAN_GATE**, **PYRAM13_ROW_SUM_MODAL_MASS**, **THIN_DOME_LOCKING**, EAS enhanced thermal/corner failures and distorted-corner risks remain in `solid_revision_results_2026-10-15.md`. No tuning was used to mask them.

Run `tests/solids/check_cpyram13_v3.py` with the project Python after building. It compiles a probe linked to actual production objects, verifies the operators independently, runs the full thermal BDF and reproduces both buckling screens. Detailed machine-readable results: `tests/solids/cpyram13_v3_results.json`. Unchanged v3 references, original meshes/S2K and SHA256 provenance are under `tests/solids/references/cpyram13_v3/`. Generated BDF/output/probe files remain ignored.
