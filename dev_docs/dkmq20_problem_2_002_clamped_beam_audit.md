# DKMQ20 Problem 2-002 Clamped Beam Audit

Date: 2026-08-25

Scope: CQUAD4 `QUAD4TYP=DKMQ20` compared with `D:\18a\python\linear\DKMQ20_ShellElement_RHR_6dof.py` using `problem_2_002_clamped_beam.py`.

## Finding

The deck/model translation is not the primary issue for DKMQ20 in Problem 2-002. After orienting the DKMQ20 drilling penalty to the element normal and skipping the generic `K6ROT` add-on, LC1, LC2, LC3, LC5, and LC6 track the Python DKMQ20 reference closely. LC4 remains the outlier and still follows the older AU-like behavior.

For `nx=12`:

| Case | Python DKMQ20 | MYSTRAN DKMQ20 | MYSTRAN vs Python |
| --- | ---: | ---: | ---: |
| LC1 | 3.000000E-05 | 3.000000E-05 | 0.000% |
| LC2 | 3.085736E-02 | 3.085692E-02 | -0.001% |
| LC3 | 4.320887E-01 | 4.318321E-01 | -0.059% |
| LC4 | 2.423612E-03 | 2.289612E-02 | 844.711% |
| LC5 | 2.569408E-04 | 2.569370E-04 | -0.001% |
| LC6 | 3.600000E-02 | 3.600000E-02 | 0.000% |

## Implemented Conservative Fix

- `Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90`
  - The drilling penalty now acts on rotation about the element normal, not hard-coded global/local component 6.
- `Source/EMG/EMG1/EMG.f90`
  - DKMQ20 now skips the generic `K6ROT`, like DSQK and CQUADR variants with their own drill treatment.

This removed the vertical-strip AutoSPC issue where the plate normal rotation could be constrained away, and it fixed the severe LC3/LC6 displacement behavior without touching GPSTRESS.

## Tested But Not Kept

The Python DKMQ20 source uses raw `An`, DKMQ20 local `V1/V2` bending rows, and a transform that leaves translations global while mapping rotations to `(alpha, beta, theta_normal)`. A direct transplant into MYSTRAN's standalone DKMQ20 branch made LC4 softer, but introduced mechanisms or large errors in LC2, LC3, and LC6. That path needs a deeper 5-DOF-to-6-DOF embedding audit before it is safe.

## Remaining Work

LC4 still needs a DKMQ20-specific kernel fix. The evidence points to the standalone CQUAD4 DKMQ20 formulation retaining AU-style twist/shear behavior, not to a deck translation/load parsing error.
