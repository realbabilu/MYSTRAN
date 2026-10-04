# Q8/T6 checkpoint — 4 October 2026

Branch: `v18.00.a`. Changes after the T6 recovery commit `4c3fba46` include the previously unpushed Q8 formulation commit `4bce7f97` and the present native shell upgrades. This is a validation checkpoint with the limitations below, not a certification of every shell configuration.

## Behavior changes

- Native Q8 formulations and point recovery use their supplied Python references (SIMOQ8, ANS8BDG6, MITC8, HBQ8; MacNeal uses MacNealQ8_1992_native_v3_drill.py).
- Shared quadratic surface pressure retains PLOAD4 corner intensities, connectivity normals, CID=0 directions and native displacement fields. PLOAD2 and continuation/THRU/LOAD protocols are supported; unsupported options fail explicitly.
- Shared translational surface mass supports RHO*t plus PSHELL NSM. COUPMASS=1 is consistent; default COUPMASS=0 uses a positive scaled diagonal rather than equal-node allocation. This changes default modal frequencies. No rotary/drilling inertia is added.
- Native membrane geometric stiffness and Q8 buckling are enabled. ARPACK buckling uses the reciprocal problem with positive K metric, positive reciprocal eigenvalue selection and explicit eigenpair residual checks.
- Native uniform-temperature and TEMPP1-gradient loads use the active Bm/Bb operators in basic coordinates. Recovery reads the gradient after all nodal mean temperatures and uses physical fiber strain eps0-z*kappa.
- Native T6 transverse shear recovery bypasses legacy PHI_SQ scaling. Nastran-aligned engineering moments use +D*kappa while retaining the physical fiber recovery convention. Legacy Q4/T3 output is preserved.
- MacNeal recovery samples CENTER plus all eight nodes directly, fills each point basis and fixes negative-fiber CENTER Sxy. F06 surface coverage is complete for the flat signed patch.
- Native Q8/T6 FORCE F06 selectors filter CENTER, CORNER/BILIN, combined and default CENTER rows, including MAX/MIN/ABS. Default FORCE is independent of STRESS(CORNER). Complete native OEF layout is retained, as confirmed against Siemens Nastran.

## Validation evidence

These earlier suites were run on their dated checkpoints. The last recovery/output changes do not edit stiffness, mass or applied-load integration. Full local solver files and reference harnesses remain in C:/PROJECTAI/18a/test and codex_mod, outside this Git repository; compact latest recovery/output results are committed alongside this note. No binaries or large solver archives are committed.

| Area | Recorded validation |
|---|---|
| Pressure | 108 element-vector cases, 27 benchmark cases, 76 protocol cases, 11 affine checks and 16 recovery regressions |
| Mass/modal | 36 native eigen cases, 18 equivalent CONM2 cases, 72 protocol cases; independent Nastran comparisons recorded |
| Buckling | 189 native cases; 30 Siemens jobs and 81 cross-solver comparisons, with formulation differences retained |
| Thermal | 45 Python/solver and 45 reaction cases, nine annulus recovery audits; ten Siemens jobs |
| Latest FORCE selectors | 60/60 cases across nine families and mixed MacNeal field; four additional Siemens selector jobs |
| Latest MacNeal recovery | 6/6 signed, shear, FORCE-only and mixed-field cases |
| Latest T6 shear | 24/24 cases |
| Engineering moments | 24 other-family cases plus six dedicated MacNeal cases |
| Latest legacy/native regression | 45 displayed force rows; native default now CENTER, legacy Q4/T3 unchanged |

The latest selector gates compare selected F06 values exactly to the saved pre-filter recovery and assert complete OP2 values/point IDs unchanged. Surface and stress values remain unchanged when requested. Native OP2 Q8 carries CENTER+4 corners; T6 CENTER+3 corners for CENTER and CORNER requests. Four Siemens 2512 engines independently confirmed these layouts. The rebuilt source has passed compilation; git diff --check is required before this commit.

Validated local executable SHA256: `b4b9765366e36f833488f0781301e10e793bc7b4c0923913715e67cdc4dbc3fd`.

## Remaining work

- General reversed-normal mechanical recovery, curved/warped recovery, non-basic coordinate systems, offsets, variable thickness and composites need further validation. Signed pressure and thermal-gradient checks do not certify all mechanical normal/fiber behavior.
- Q8 shear magnitude differs from Siemens CQUAD8 despite the documented sign/analytic patch checks; no universal numerical equality is claimed.
- Modified curved Q8 pressure fields have a physical first-moment residual around 0.028% of force times length scale.
- Buckling agreement with the port does not imply all families meet the independent plate-theory accuracy screen; spectral/shift ranges and thermal prestress remain open.
- Mass scope is translational; rotary/drilling inertia and general curved/full-six-DOF mode comparisons remain open.
- FORCE-only without GPSTRESS has no surface table. Different location selectors per subcase still use existing first-subcase/global case-control handling.

See the local dated validation bundles for full inputs, operators, results and tolerances. Historical thermal/mass notes describe their earlier recovery state; the sign, MacNeal and FORCE behavior documented here supersedes those historical output statements.
