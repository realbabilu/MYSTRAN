# HBQ8 MYSTRAN: preserved V2 baseline and Python V3 upgrade

Selector: `PARAM,QUAD8TYP,HBQ8`.
Source: `MYSTRAN/Source/EMG/EMG4/CQUAD8_HBQ8.f90`.
Reference: supplied `python/HBQ8_Kikuchi_v3.py` and its real HBQ8 base,
`HBQ8_Kikuchi_v2.py`, `KikuchiMacNeal_HBQ8_v1.py`, `HBQ8_ShellElement.py`,
`macneal_kikuchi_shapefun.py`, and `shell_q8_common.py`.

## Preserved results

Workspace evidence directory: `test/hbq8_upgrade/v3_checkpoint`.
Before source edits, seven existing HBQ8 models were run using the existing
MYSTRAN executable. `baseline/` retains full DAT/F06/OP2/ERR/NEU/stdout,
source, executable, Python dependencies, manifest and results.json.
`v3/` retains matching upgraded runs, source and executable. All seven deck
SHA256 hashes match before and after. The original corrected source deck and
prior benchmark output directories were preserved.

Baseline executable SHA256:
`6670f4a0140dcbc17ac2c3658f8c60525bf3593b236568cebf00ca29bf9dadac`.
Upgraded executable SHA256:
`80cccf87736c26a192beb675ca2f4e2cd8450b6c40c65e22af77e53bc7f0e714`.

2-005 uses `test/pressure_validation/corrected/prob_2_005_simple_uniform_nx16_cquad8.dat`:
16x16 elements, 833 nodes. Its solver copy adds the HBQ8 selector, explicit
free-field continuation serialization, OLOAD/SPCFORCES output and AUTOSPC=N.
2-016 uses canonical 16x16 shear/bending panels at t=.01 and t=1.
2-017 uses canonical 2x12 and 4x24 column meshes at t=1.

## Formulation changes

- Kikuchi interpolation retained, with standard Q8 geometry area and modified
  field tangent metric for operators.
- Shear blend: Bs=(1-w)*Bs_tensor6+w*Bs_bdg4;
  w=max(0,1-15/s), s=sqrt(element area)/h, area evaluated with 3x3 quadrature.
- Membrane tying automatically enabled when nodal-normal angle exceeds .01 rad
  or relative midside offset exceeds .001. Flat straight-sided elements retain
  direct membrane. HBQ8's own membrane tying uses +/-1/sqrt(3) for the linear
  direction (6+6+4 samples), unlike the ANS8 normal-membrane tying locations.
- Drilling unchanged: beta_drill=1e-3, penalty coefficient beta_drill*G*h.
- Consistent/positive HRZ mass now uses 3x3 active-field interpolation and rotary
  inertia in all three directions. PSHELL NSM contributes to translations only;
  rotary inertia comes from material density.
- Recovery, thermal and initial-stress reference state use upgraded B operators.
  Existing native 5x5 normal-pressure integration and flat 4x4 geometric
  stiffness retained. Physical TEMPP1 gradient follows connectivity normal;
  its Bb thermal force has the existing negative sign.

This port exposes V3 defaults through HBQ8. It adds no HBQ8 PARAM switches for
Python's optional shear pattern, membrane tying, MacNeal field or rotary mode.
The global ANS-specific PARAM controls are not used by this HBQ8 routine.

## Before / after

| Case | Existing | V3 | Assessment |
|---|---:|---:|---|
| 2-005 center Uz | -13.51560497 | -12.97469330 | error vs -12.97: +4.2067% to +0.0362% |
| 2-016 shear, t=.01, first factor | .0004593225 | .0004242379 | critical-stress error +8.0733% to -0.1817% |
| 2-016 bending, t=.01, first factor | .001249732 | .001161134 | critical-stress error +7.2814% to -0.3242% |
| 2-016 shear, t=1 | 400.8203 | 400.8203 | unchanged at F06 precision |
| 2-016 bending, t=1 | 1086.480 | 1086.480 | unchanged at F06 precision |
| 2-017 column 2x12 | 126.5028 | 126.5028 | unchanged at F06 precision |
| 2-017 column 4x24 | 126.4854 | 126.4854 | unchanged at F06 precision |

Thin 2-016 critical stress=factor/t for these q=1 decks. References:
shear .042501034; bending .116491057. The measured improvement reproduces the
user's Python V2-to-V3 thin-plate table. Thick-plate factors are not evaluated
against thin-plate exact references. All first-three-factor results and deck
hashes are saved in `comparison.json`; every column factor is unchanged at F06
precision. These results use the actual available HBQ8 base, not a surrogate.

## Validation

- Production Bm/Bb/Bs/Bdrill and consistent/HRZ mass against Python V3:
  15 geometry/thickness configurations, 18 points each (270 points), including
  curved, irregular-curved and reversed connectivity. Maximum B error
  3.775e-15; mass error 6.662e-16. `operator_validation/results.json`.
- 2-005 native pressure load, displacement and equivalent FORCE comparison:
  PASS. Python Uz=-12.9746301928; native Uz=-12.9746932983; maximum displacement
  per-component relative difference 4.900e-6 (thin-model sensitivity).
  `python_pressure/results.json`.
- Six SOL105 comparisons against actual Python V3: PASS. First three factors
  agree within 4.114e-7, translation MAC >.999, internal double-precision
  residual <=2.663e-8. Evidence in
  `test/buckling_quadratic/solver_results/hbq8_v3_python/mystran/summary.json`
  and `hbq8_v3_python_column/mystran/summary.json` in the same parent directory.
  Raw float32 OP2-mode reconstructed residuals are much larger on thin models;
  the gate uses solver internal residuals and mode MAC.
- Six SOL103 cases: consistent/HRZ x material density/NSM/mixed mass, PASS.
  Five frequencies agree with actual Python V3 K/M within 3.413e-7.
  `modal/results.json`.
- Curved prescribed-displacement stress/force recovery: center/corners/global
  GP stresses and N/M/Q agree with actual Python V3, maximum relative stress
  error 3.483e-6; PASS. `curved_recovery/results.json`.
  Physical fiber strain is Bm*u-z*Bb*u, so Z1=-h/2 uses plus generalized moment;
  the legacy plotting helper's opposite fiber ordering is corrected in this
  validation script without changing production stress conventions.
- Five thermal cases: expansion, restraint, both gradient windings and annulus:
  PASS. Actual Python V3 reference CSVs were generated under HBQ8v3 labels,
  preserving V2 references; physical gradient patch expectations are independent.
  `thermal_python.json` and
  `test/Q8_T6_thermal_validation/solver_results/hbq8_v3_actual_python/mystran/summary.json`.
- CMake build succeeded; git diff --check passed. Existing compiler dependency,
  allocatable-assignment and unused-helper warnings remain.

## Reproduction from workspace root

```text
python test/hbq8_upgrade/run_v3_checkpoint.py v3
python test/verify_hbq8_v3_operators.py
python test/hbq8_upgrade/validate_v3_solver.py pressure
python test/hbq8_upgrade/validate_v3_solver.py buckling
python test/hbq8_upgrade/validate_v3_solver.py buckling --only HBQ8 --meshes 1 2 --cases 017 --thickness 1 --tag hbq8_v3_python_column
python test/hbq8_upgrade/validate_v3_mass.py
python test/hbq8_v3_curved_stress_gate.py
python test/hbq8_upgrade/validate_v3_thermal.py
```

The baseline checkpoint refuses to overwrite a completed baseline manifest.
The V2 operator test uses frozen pre-upgrade source to remain reproducible.

Scope: experimental linear surface extension, not a full published HBQ8
reproduction. Flat membrane geometric stiffness only; no rotational or curved
stress tangent, composites or variable-thickness extension. Normal-pressure
first-moment limitations on general curved geometry remain; Enhanced report
3D pressure corrections are absent from the supplied HBQ8 source and are not
part of this port. These benchmarks do not establish universal accuracy.
