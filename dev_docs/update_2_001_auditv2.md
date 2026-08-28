# update_2-001_auditv2

## Scope

This audit bundle records the 2-001 patch-test element updates made after the
general stress/force recovery audit.

Changed MYSTRAN paths:

- `Source/EMG/EMG4/CQUADR_MITC4PHB.f90`
- `Source/EMG/EMG4/TPLT_MITC3P.f90`
- `Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90`

## MITC4PD

`MITC4PD` recovery now delegates its stress-recovery pass to the validated
`CQUADR_DKMQ24R` recovery path. Its GPSTRESS output also keeps the fixed
element-frame convention to avoid a second point-local basis transform.

In the refreshed 2-001 audit, `MITC4PD` matches the Python-style patch target
well for LC1 and LC2.

## MITC3+ CTRIA3

Important naming note:

- `MITC3+` is the `CTRIA3` version, through `TPLT_MITC3P.f90`.
- `MITC3+HB` is the separate `CTRIAR` version, through `CTRIAR_MITC3PHB.f90`.

For `MITC3+` / `CTRIA3`, the transverse-shear stiffness path was aligned more
closely with `MITC3p_ShellElement_6dof.py`:

- use the same 3-point triangular integration rule
- use the direct Python-style MITC3+ shear interpolation operator
- use `SHELL_T` directly for transverse shear
- remove the older rotational sign-flip and transform wrapping from the
  condensed CTRIA3 MITC3+ stiffness path

This fixes LC3 constant shear for `MITC3+` / `CTRIA3`:

```text
PYTHON MITC3+6d LC3 w_internal = 0.0053485522866, avg_err = 2.753594788%
MYSTRAN MITC3+ LC3 w_internal = 0.005348552,     avg_err = 2.753600000%
```

## Known Remaining Failure

`MITC3+` / `CTRIA3` LC2 bending in problem 2-001 is still a known fail versus
Python `MITC3+6d`.

Current refreshed CSV comparison:

```text
PYTHON MITC3+6d LC2:
  mxx = +1.111111111e-07
  myy = +1.111111111e-07
  mxy = +3.333333333e-08
  avg_err = 0.040034%

MYSTRAN MITC3+ CTRIA3 LC2:
  mxx = -1.109710500e-07
  myy = -1.115107788e-07
  mxy = +3.335987625e-08
  avg_err = 4.506543%
```

The sign of `mxx` and `myy` follows the MYSTRAN/NASTRAN bending convention used
by the 2-001 comparison. The remaining issue is the LC2 moment spread/uniformity
in the `CTRIA3` `MITC3+` path, not a parser problem and not the `CTRIAR`
`MITC3+HB` element.

Future work should upgrade `TPLT_MITC3P.f90` so the `CTRIA3` `MITC3+` LC2
bending result behaves like Python `MITC3+6d` while preserving the now-fixed LC3
constant-shear behavior.

## Verification

Commands run:

```text
cmake --build D:\18a\MYSTRAN\build -j4
D:\18a\MYSTRAN\Binaries\mystran.exe D:\18a\python\linear\working_mystran_2_001\prob_2_001_mh_ctria3_mitc3plus.dat
D:\18a\MYSTRAN\Binaries\mystran.exe D:\18a\python\linear\working_mystran_2_001\prob_2_001_patchshear_ctria3_mitc3plus.dat
python D:\18a\python\linear\battle_2_001_patch.py
```

The build completed. The benchmark refreshed
`D:\18a\python\linear\value_summary_2_001_patch.csv`.
