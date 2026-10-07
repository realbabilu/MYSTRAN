# MITC6 MYSTRAN upgrade to enhanced Python V5e

Selector: `PARAM,TRIA6TYP,MITC6`.
Production source: `MYSTRAN/Source/EMG/EMG4/CTRIA6_MITC6.f90`.
Reference: `python/MITC6_Tri_v5e.py`, identical to the enhanced-elements copy;
pressure helper is `python/enhanced_Q8_T6_pressure_roof/enhanced_elements/shell_pressure_enhanced.py`.
Original Python files were preserved. The native baseline implements V4.

## Preserved evidence

Workspace directory: `test/mitc6_upgrade/v5e_checkpoint`.
Seven native runs were completed before source edits. `baseline/` contains
full DAT/F06/OP2/ERR/NEU/stdout, executable, source, pressure helper/interface
and Python references. `v5e/` contains upgraded snapshots and matched reruns.
`comparison.json` retains all first-three-mode factors. All seven before/after
deck hashes match. Original corrected decks and prior benchmark outputs remain.

Baseline executable SHA256:
`80cccf87736c26a192beb675ca2f4e2cd8450b6c40c65e22af77e53bc7f0e714`.
Upgraded executable SHA256:
`13ce3108c61f9ac75c1c92ca576caebf58b3524bbabbd6d658d8faa287d28cae`.

2-005: corrected `prob_2_005_simple_uniform_nx16_ctria6.dat`, existing pressure
harness serialization and explicit AUTOSPC=N. 2-016: canonical 16x16 panels,
shear/bending at t=.01 and t=1. 2-017: column meshes 2x12 and 4x24 at t=1.
Roof checks reuse the captured baseline executable, independently retaining
old and new solver output under `roof/`; they were run after the main port.

## Ported formulation

- Membrane: edge-tangential natural-strain samples on s=0, r=0 and r+s=1,
  plus centroid affine closure, matching enhanced V5e default membrane_tying=edge.
- Shear: e_st samples moved from interior (R1,R1)/(R1,R2) to (0,R1)/(0,R2).
  The interpolation matrix rows change with the sampling points.
- Existing bending, director sign and drilling beta=.02 retained.
- Consistent/positive HRZ mass uses Python's six-point triangle quadrature,
  with rotary inertia for all three rotations. PSHELL NSM is translational only.
- Initial-stress stiffness uses six-point triangle quadrature and the upgraded
  membrane recovery, matching Python kg_global, with the existing native
  nondegenerate-flat-geometry guard retained.
- Normal/directional PLOAD pressure uses standard T6 shape functions and
  eight-point-per-direction Duffy quadrature, matching the enhanced helper.
  Positive normal pressure follows connectivity, fixed direction is normalized
  in basic axes, and intensity is per actual surface area. No nodal couples.
- Shared QUADRATIC_SURFACE_PRESSURE gains an optional HIGH_ORDER argument;
  only MITC6 requests it. All other callers retain their existing five-point rule.
- Thermal uses the upgraded Bm and the physical connectivity-normal TEMPP1
  convention. Existing physical fiber stress and generalized-moment signs retained.

Python's optional interior membrane configuration and variable-thickness arrays
are not exposed as new MYSTRAN PARAM options. The native selector uses V5e defaults.

## Measured before / after

| Case | Existing V4 | Enhanced V5e | Assessment |
|---|---:|---:|---|
| 2-005 center Uz | -13.06639671 | -12.97978592 | error vs -12.97: +.743228% to +.075450% |
| 2-016 shear, t=.01, first factor | .0004285837 | .0004234285 | stress error +.840770% to -.372189% |
| 2-016 bending, t=.01, first factor | .001173875 | .001158739 | stress error +.769538% to -.529789% |
| 2-016 shear, t=1 | 398.8755 | 397.7278 | changed; agrees with actual Python V5e |
| 2-016 bending, t=1 | 1080.449 | 1077.180 | changed; agrees with actual Python V5e |
| 2-017 column 2x12 | 126.5009 | 126.5009 | unchanged at F06 precision |
| 2-017 column 4x24 | 126.4852 | 126.4852 | unchanged at F06 precision |
| 2-006 linear-geometry roof Uz | -.30377501 | -.30441040 | error vs -.3086: -1.563509% to -1.357616% |
| 2-006 cylinder roof Uz | -.28382701 | -.30164146 | error vs -.3086: -8.027542% to -2.254872% |

Thin-plate critical stress=factor/t for q=1. References: shear .042501034,
bending .116491057. This reproduces UPDATED_ELEMENT_COMPARISON.md's MITC6
V4-to-V5e errors +.841% to -.372%, and +.770% to -.530%. Thick factors are not
compared against exact thin-plate stress values. No accuracy improvement is
claimed solely because the thick factors decrease. All three column factors
are unchanged at F06 precision.

Enhanced_shell_pressure_report.md's roof before value -7.787% uses Python
pre-enhancement V5e/interior sampling; the current MYSTRAN baseline is V4,
hence its -8.028% before value differs. The enhanced after value -2.255% agrees.

## Validation completed

- Production Bm/Bb/Bs/Bdrill and both mass modes against actual Python V5e:
  15 geometry/thickness configurations x 11 points = 165 points, including
  curved, irregular curved and reversed connectivity. Maximum B difference
  6.141e-15; mass difference 1.666e-16. `operator_validation/results.json`.
- Pressure plate plus both roof geometries: loads, displacement and equivalent
  FORCE solves PASS. Plate Python Uz=-12.9797581879, native -12.9797859192;
  full-displacement per-component relative difference <=1.297e-5, reflecting
  thin-model numerical sensitivity. Roof nodal differences <=1.502e-7.
  `python_pressure/results.json`, `roof/results.json`.
- Six native SOL105 comparisons PASS; first three factors agree with Python
  within 4.222e-7, translation MAC >.999, internal double residual <=3.420e-8.
  Evidence: `test/buckling_quadratic/solver_results/mitc6_v5e_python/mystran/summary.json`
  and `mitc6_v5e_python_column/mystran/summary.json` in the same parent directory.
  The API adapter only maps k_geometric_global to Python's actual kg_global;
  it does not substitute the geometric-stiffness algorithm. Raw float32 OP2
  mode reconstructed residuals are large for thin models; internal residuals
  and mode MAC provide the gate.
- Six SOL103 density/NSM/mixed cases x consistent/HRZ PASS. Five frequencies
  agree within 3.670e-7. `modal/results.json`. The supplied Python HRZ function
  divides 0/0 for rho=0; the validation driver explicitly uses the exact zero
  material mass in the NSM-only case. Native zero-density mass stays finite.
- Five thermal cases (free/restrained expansion, both gradient windings,
  annulus) PASS against actual Python V5e and independent physical expectations.
  Adapter maps the alpha/dT keyword convention to V5e's native constructor-alpha
  and dT/dT_grad API without changing operators or force signs.
  `thermal_python.json` and
  `test/Q8_T6_thermal_validation/solver_results/mitc6_v5e_actual_python/mystran/summary.json`.
- 36 single-element pressure checks PASS: normal, fixed direction, varying
  pressure on flat/curved/reversed/rotated geometry, MITC6 plus ANS8BDG6 and
  SIMOT6 default-order regression. `pressure_single/results.json`.
- CMake build and git diff --check succeeded. Existing compiler warnings remain.

## Reproduction from workspace root

```text
python test/mitc6_upgrade/run_v5e_checkpoint.py v5e
python test/verify_mitc6_v5e_operators.py
python test/mitc6_upgrade/validate_v5e_solver.py pressure
python test/mitc6_upgrade/validate_v5e_solver.py buckling
python test/mitc6_upgrade/validate_v5e_solver.py buckling --only MITC6 --meshes 1 2 --cases 017 --thickness 1 --tag mitc6_v5e_python_column
python test/mitc6_upgrade/validate_v5e_mass.py
python test/mitc6_upgrade/validate_v5e_thermal.py
python test/mitc6_upgrade/validate_v5e_roof.py
```

Baseline capture refuses to overwrite a completed baseline manifest.
Scope: isotropic constant-thickness linear surface extension, not proof of a
complete published MITC6 formulation. Curved elastic pressure response is
validated on the stated specimens; curved-shell or rotational geometric
stiffness, follower pressure, composites and variable thickness are not ported.
