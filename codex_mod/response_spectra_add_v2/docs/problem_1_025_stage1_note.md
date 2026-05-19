# Problem 1-025 Stage 1 Note

## What 1-025 really is

Problem 1-025 is **not** a true two-direction response spectrum benchmark.
It is a **single excitation direction** benchmark with **torsional coupling**.

From `Problem 1-025.s2k`:
- RS cases: `RS-CQC`, `RS-SRSS`, `RS-ABS`
- all use `LoadName=U1`
- all use `CoordSys=GLOBAL`
- all use `Angle=0`
- all use `Function=ELCEN`
- all use `TransAccSF=386.4`

So the comparison is still:
- modal superposition method: `CQC`, `SRSS`, `ABS`
- under a single global X excitation

## Why it is harder than 1-024

Unlike Problem 1-024, this case includes rotational inertia at the diaphragm masters.

Joint added masses in SAP:
- Joint 49: `Mass1=1.24224`, `Mass2=1.24224`, `MMI3=174907.4`
- Joint 50: `Mass1=1.24224`, `Mass2=1.24224`, `MMI3=174907.4`
- Joint 51: `Mass1=1.24224`, `Mass2=1.24224`, `MMI3=174907.4`

This means the stage mass source is:
- translational mass in X and Y
- rotational mass inertia about diaphragm Z

## Correct MYSTRAN mass representation for Stage 1

For Problem 1-025 Stage 1, the preferred MYSTRAN concentrated-mass representation is:
- `CONM2` at joints `49`, `50`, `51`

Reason:
- `CMASS2` is useful for directional translational masses
- but `CMASS2` does not carry rotational inertia naturally
- `CONM2` can represent the translational mass plus `I33` in one card

So the natural mapping is:
- `mass = 1.24224`
- `I33 = 174907.4`
- `I11 = I22 = 0`
- offsets = 0

Important caveat:
- plain `CONM2` in MYSTRAN is isotropic in translational mass (`Mx = My = Mz = mass`)
- SAP uses only X and Y mass, with `Mass3 = 0`

Therefore Stage 1 should accept one of these two positions explicitly:
1. pragmatic baseline: use `CONM2` and accept surplus Z translational mass as a temporary approximation
2. stricter future step: extend concentrated-mass handling so X/Y translational mass and Rz inertia can coexist without adding Z mass

## Constraint model

SAP diaphragm assignments:
- `DIAPH2` includes story-2 nodes and master `49`
- `DIAPH3` includes story-3 nodes and master `50`
- `DIAPHROOF` includes roof nodes and master `51`

For MYSTRAN Stage 1, use rigid diaphragms with master nodes:
- 49 for level 2
- 50 for level 3
- 51 for roof

Dependent DOF should be:
- `1,2,6` = `Ux, Uy, Rz`

This matches the diaphragm behavior used in Problem 1-024 work.

## Stage 1 acceptance scope

Stage 1 for Problem 1-025 should validate:
- modal foundation with diaphragm masters `49/50/51`
- concentrated mass source at those masters
- torsional coupling appears in modal results
- RS methods `SRSS`, `CQC`, `ABS` run under `U1`

Stage 1 does **not** need to claim:
- true two-direction RS combination
- configurable orthogonal percentage combination
- exact X/Y-only translational mass without any Z leakage, unless solver-side concentrated-mass refinement is added

## Practical next build target

Build a first MYSTRAN deck for Problem 1-025 using:
- frame stiffness from SAP `FSEC1`
- rigid diaphragms with masters `49/50/51`
- `CONM2` masses at `49/50/51`
- one-direction response spectrum (`U1` only)
- `SRSS`, `CQC`, `ABS`

Treat that as the Stage 1 baseline before moving on to true directional combination work.
