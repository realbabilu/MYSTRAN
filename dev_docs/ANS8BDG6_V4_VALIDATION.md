# MYSTRAN ANS8BDG6: saved baseline and Python V4 port

The active selector is `PARAM,QUAD8TYP,ANS8BDG6`. Source:
`MYSTRAN/Source/EMG/EMG4/CQUAD8_ANS8BDG6.f90`. The rebuilt executable is
`MYSTRAN/Binaries/mystran.exe`.

## Preserved baseline

Before any ANS8BDG6 source edit, seven decks were run with the original executable.
`baseline/` contains the source, executable, Python reference dependencies,
manifest, results, and full DAT/F06/OP2/ERR/NEU/stdout outputs for every case.
`v4/` contains the corresponding final source/executable snapshot and reruns.
The deck hashes match before/after for all seven cases. Existing benchmark
outputs and the user-supplied corrected deck were not replaced.

Baseline executable SHA256:
`f31ebbe2f9fdbd5f6ffc0de417883c1394201ce84b44926605c9674b8076d91f`.
Final executable SHA256:
`6670f4a0140dcbc17ac2c3658f8c60525bf3593b236568cebf00ca29bf9dadac`.

2-005 uses the supplied corrected 16x16 CQUAD8 model, 833 nodes and 256 elements.
Its copy adds the ANS8BDG6 selector, OLOAD/SPCFORCES output, explicit free-field
continuations and `PARAM,AUTOSPC,N`, following the existing pressure harness.
2-016 uses the saved canonical shear and bending decks, 16x16 panels, t=0.01
and t=1. 2-017 uses saved canonical column decks with 2x12 and 4x24 panels, t=1.

## Changes ported from supplied V4

- Default Kikuchi field retained. The port follows `ANS8_BDG6_v4.py`
  and `shell_q8_common.py`, not the separate Enhanced pressure report.
- Shear `auto`: Bs=(1-w)*Bs_tensor6+w*Bs_bdg4,
  w=max(0,1-15/s), s=sqrt(standard-geometry area)/thickness.
- Membrane tying `auto`: 6+6+4 covariant samples for curved/warped geometry;
  direct membrane operator for flat straight-sided elements. Normal-angle
  threshold 0.01 rad, relative midpoint-deviation threshold 0.001.
- ANS membrane auto, Kikuchi field, and ANS drilling penalty retained unchanged: kt=0.1 and beta_drill=kt*h.
- Consistent and positive HRZ mass use the active field and 3x3 area integration,
  with V4's default rotary_inertia='all'. PSHELL NSM remains translational
  midsurface mass; rotary mass comes from material density only.
- Stiffness, output recovery, thermal loading and prestress recovery use the
  same upgraded B operators. Native 5x5 pressure and 4x4 initial-stress
  integration are retained.

MYSTRAN TEMPP1 uses the physical connectivity-normal gradient. Its existing
negative Bb thermal term is deliberately retained; the supplied Python
f_thermal has the opposite gradient argument convention. This is a card/API
mapping, not a change to the elastic V4 operators. The default is now PARAM,ANSSHEAR,AUTO. Explicit TENSOR6 and BDG4 remain supported;
DIRECT is also accepted, matching Python. PARAM,ANSFIELD,KIKUCHI/STANDARD and
PARAM,ANSMEM,AUTO/ON/OFF retain their existing meanings. Python's optional
MacNeal-Harder field is not an accepted native ANSFIELD value. Explicit TENSOR6
reproduced the pre-upgrade 2-005 displacement exactly at OP2 precision.


## Matched before/after results

| Case | Baseline | V4 | Comparison |
|---|---:|---:|---|
| 2-005 center Uz (inch) | -13.51560497 | -12.97469330 | signed error vs -12.97: +4.2067% to +0.0362% |
| 2-016 shear, t=.01, first factor | .0004593225 | .0004242379 | critical stress error: +8.0733% to -0.1817% |
| 2-016 bending, t=.01, first factor | .001249732 | .001161134 | critical stress error: +7.2814% to -0.3242% |
| 2-016 shear, t=1, first factor | 400.8203 | 400.8203 | unchanged at F06 precision |
| 2-016 bending, t=1, first factor | 1086.480 | 1086.480 | unchanged at F06 precision |
| 2-017 column 2x12, first factor | 126.5028 | 126.5028 | unchanged at F06 precision |
| 2-017 column 4x24, first factor | 126.4854 | 126.4854 | unchanged at F06 precision |

For the thin 2-016 decks, reference edge line load is q=1, so critical stress
is factor/t. The existing driver references are shear 0.042501034 and bending
0.116491057. The factors and critical stresses must not be interchanged.
The thick plate is not assessed against an exact thin-plate reference here.
Full three-factor comparisons and deck hashes are in `comparison.json`.

## Validation completed on final executable

- Production Bm/Bb/Bs/Bdrill routines (AUTO, TENSOR6, BDG4, DIRECT; KIKUCHI/STANDARD; membrane AUTO/OFF/ON): 60 geometry/thickness/parameter configurations,
  18 evaluation points each (1080 points), against supplied Python V4.
  Maximum absolute error 3.109e-15. Curved, irregular-curved and reversed
  connectivity included. Consistent/HRZ mass in all configurations agrees
  within 8.882e-16. Evidence: `operator_validation/results.json`.
- 2-005 pressure loads, displacement and equivalent nodal FORCE solve passed.
  Python V4 center Uz=-12.9746301928; native Uz=-12.9746932983. Full displacement
  maximum per-component relative difference below 4.900e-6, consistent with
  thin-model conditioning. Evidence: `python_pressure/results.json`.
- Six SOL105 comparisons to actual Python V4 passed: four 16x16 plate cases
  and two column meshes. Maximum first-three-factor relative difference
  4.114e-7 (F06 rounding); translation mode MAC exceeds .999. Solver internal
  double-precision residuals are below 2.663e-8. Raw float32 OP2 modes have
  much larger reconstructed residuals on thin models; the gate evaluates
  internal residuals and mode MAC. Evidence:
  `../../buckling_quadratic/solver_results/ans8_v4_python/mystran/summary.json`
  and `../../buckling_quadratic/solver_results/ans8_v4_python_column/mystran/summary.json`.
- Six native SOL103 checks against Python V4 K/M passed: consistent/HRZ,
  material density, NSM only and mixed density/NSM. Maximum five-frequency
  relative difference 3.413e-7. Evidence: `modal/results.json`.
- Native thermal free/restrained expansion, signed gradient patches in both
  windings and annulus passed. Patches use independent physical expectations;
  annulus also agrees with the retained corrected thermal reference.
  Evidence: `../../Q8_T6_thermal_validation/solver_results/ans8_v4/mystran/summary.json`.
- CMake build succeeded. Existing LAPACK duplicate dependency warnings,
  unused helper warnings and an allocatable-assignment warning remain.

## Reproduction

```
python test/verify_ans8_v4_operators.py
python test/ans8_upgrade/run_v4_checkpoint.py v4
python test/ans8_upgrade/validate_v4_solver.py pressure
python test/ans8_upgrade/validate_v4_solver.py buckling
python test/ans8_upgrade/validate_v4_solver.py buckling --only ANS8BDG6 --meshes 1 2 --cases 017 --thickness 1 --tag ans8_v4_python_column
python test/ans8_upgrade/validate_v4_mass.py
python test/Q8_T6_thermal_validation/solver_thermal_gate.py --only ANS8BDG6 --gradient --annulus --tag ans8_v4
```

The baseline checkpoint refuses to overwrite its completed manifest.
`verify_ans8_v3_operators.py` uses the frozen V3 source after production
advances to V4, so the old operator reference remains reproducible.

Limitations: flat membrane initial-stress stiffness only, no rotational or
curved-shell stress tangent; no composites or variable-thickness extension.
The Enhanced report's 3D geometry-reproduction correction is not part of the
supplied V4 source and has not been ported. Modified-field pressure retains
the previously identified physical first-moment limitation on general curved
geometry. These comparisons do not establish universal mesh/geometry accuracy.

## Agreement with Python improvement reports

UPDATED_ELEMENT_COMPARISON.md reports ANS8 Kikuchi v3 to v4 thin-plate errors
+8.073% to -0.182% (shear) and +7.281% to -0.324% (bending). These native runs
reproduce both improvements. Column and thick-plate results remain unchanged.
Enhanced_shell_pressure_report.md reports +0.0358% for the flat 2-005 plate;
native +0.0362% agrees within thin-model numerical sensitivity. This flat
comparison does not validate the report's curved-shell pressure enhancements.
The full baseline executable and solver outputs are retained independently.
