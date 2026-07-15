# MIN4 Basic Element Stress Patch Test

## Purpose

This note records the intended validation path for SAP2000 shell patch test `2-001` when using ordinary MYSTRAN element stress output.

The target comparison is element stress in the basic/global coordinate system, not grid-point stress recovery.

## Recommended deck controls

Use ordinary shell element stress output with `PARAM,STR_CID,0`:

```nastran
SUBCASE 1
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL
SUBCASE 2
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL
BEGIN BULK
PARAM,POST,-1
PARAM,POSTEXT,YES
PARAM,STR_CID,0
```

For OP2-only validation, the F06 detailed table does not need to be authoritative. The OP2 `OES1X1` shell stress table is the target.

## Why `PARAM,STR_CID,0`

SAP2000 problem `2-001` publishes global engineering component values:

```text
Sxx = 1333.3333
Syy = 1333.3333
Sxy = 400.0
```

Default Nastran/MYSTRAN shell stress output is element-local. In a distorted patch, elements whose local axes are not aligned with global X/Y can show different `XX/YY/XY` values while preserving the same physical stress tensor.

`PARAM,STR_CID,0` requests MYSTRAN to transform shell stress and strain tensors to the basic coordinate system before output. This is the correct path for this SAP-style patch test.

## Element scope

This note applies to the ordinary element stress path for:

- `CQUAD4` using MIN4/MITC4+ branch behavior.
- `CTRIA3` using MIN3/MITC3+ branch behavior.
- `CQUADR` and `CTRIAR` where the same output-family assumptions are valid.

`CQUAD8` can use the same coordinate-output rule, but it is not the primary patch-test target here.

## Validator path

Use ordinary shell stress paths:

```text
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/XX
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/YY
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/XY
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/MAJOR
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/MINOR
SC/#/SHELLSTRESSES/EID/#/CORNER/#/ZMID/VONMISES
```

These paths intentionally map to OP2 `OES1X1`, not OP2 `OGS1`.

`op2_query.py` now exposes a `ZMID` midpoint alias from OP2 plate `Z1/Z2` rows so the same case-file paths can be used for F06 and OP2 validation.

## Current local decks

The local validation decks for this path are:

```text
D:\18a\MYSTRAN_Validation-main\decks\shell\prob_2_001_thin.dat
D:\18a\MYSTRAN_Validation-main\decks\shell\prob_2_001_thick.dat
```

The MSC comparison decks used during investigation were kept separate from the MYSTRAN decks, for example:

```text
D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_3_msc*.dat
D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_msc*.dat
```

## Pass criterion

For the membrane subcase of SAP2000 `2-001`, the center/mid-plane element stress values should match the published global component table when `PARAM,STR_CID,0` is present.

Do not switch this case to `GPSTRESS` unless the validator is explicitly changed to compare grid-point surface stress results.

