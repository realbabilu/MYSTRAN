# Banded Path Audit V1

## Purpose

This note captures the current MYSTRAN banded path as it exists today, with emphasis on:
- where the banded path is selected
- where sparse-to-band conversion happens
- where extra storage is created
- where memory and time are likely being lost

This audit is intentionally limited to the **banded subsystem**.

The sparse path based on:
- `CSR/CSC`
- `SuperLU`
- `SuperLU_MT`

is treated as out of scope and untouched.

## High-Level Summary

The current banded path is already workable, but it is not banded-native.

The main pattern is:
1. build or keep matrices in the default sparse world
2. compute bandwidth
3. allocate dense band arrays
4. convert sparse CRS storage into LAPACK/ARPACK band layouts
5. factor and solve

This means the biggest current risks are:
- peak memory during conversion
- repeated conversion cost
- duplicated operator/storage forms for the same matrix
- unnecessary zero-initialization of large band arrays

## Current Decision Points

## 1. Grid resequencing / ordering

The ordering decision is now centered in:
- [SEQ_PROC.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK1/L1B/SEQ_PROC.f90)

Current behavior:
- `RCM` can be explicitly requested through `GRIDSEQ`
- `BANDEDOPT=Y` plus `SOLLIB='BANDED  '` also activates in-core `RCM`
- `BANDIT` remains only as explicit legacy fallback

Practical result:
- the preferred banded ordering path is already no longer legacy `BANDIT`
- this is a strong base for all future banded optimization work

## 2. Solver-family dispatch

Banded behavior then branches by analysis/solver path, not from one single banded core.

Important current banded users include:
- [SYM_MAT_DECOMP_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/SYM_MAT_DECOMP_LAPACK.f90)
- [FBS_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/FBS_LAPACK.f90)
- [EIG_LANCZOS_ARPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK4/EIG_LANCZOS_ARPACK.f90)
- [EIG_INV_PWR.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK4/EIG_INV_PWR.f90)
- [SOLVE_DLR.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK6/SOLVE_DLR.f90)
- [SOLVE_GOA.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK2/SOLVE_GOA.f90)
- [SOLVE_PHIZL1.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK6/SOLVE_PHIZL1.f90)
- [SOLVE_UO0.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK2/SOLVE_UO0.f90)

So the codebase already has a partial banded subsystem, but it is spread across multiple workflows.

## Current Conversion Path

## 1. Symmetric half-band conversion

The main symmetric band conversion path is:
- sparse CRS input
- bandwidth from `BANDSIZ`
- conversion by [BANDGEN_LAPACK_DPB.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/BANDGEN_LAPACK_DPB.f90)

Used by:
- [SYM_MAT_DECOMP_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/SYM_MAT_DECOMP_LAPACK.f90)

Important detail:
- `BANDGEN_LAPACK_DPB` only stores upper-triangular terms into `KD+1` rows
- this is already the correct storage shape for `DPBTRF/DPBTRS`

This is the best existing starting point for banded memory optimization.

## 2. General band conversion for ARPACK

The Lanczos path in [EIG_LANCZOS_ARPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK4/EIG_LANCZOS_ARPACK.f90) does more work:

1. build sparse `KMSM = KLL - sigma*MLL` or `KLL + sigma*KLLD`
2. compute its bandwidth
3. allocate `RFAC`
4. convert `KMSM` to:
   - `DPB` band form if requested, or
   - `DGB` band form if requested
5. also create `KMSMn`, a nonsymmetric sparse operator form for ARPACK multiply loops

This means Lanczos banded path can temporarily hold:
- `KMSM` sparse
- `RFAC` banded
- `KMSMn` sparse nonsymmetric

That is one of the most important memory hotspots in the current system.

## Current Memory Hotspots

## 1. Full band allocation size is fixed-width

In [SYM_MAT_DECOMP_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/SYM_MAT_DECOMP_LAPACK.f90), memory is estimated as:

- `ABAND ~ (KD + 1) * NROWS`

That is already better than full dense, but still wastes all zero padding inside the band.

## 2. ARPACK `DGB` path is much larger than `DPB`

In [EIG_LANCZOS_ARPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/LK4/EIG_LANCZOS_ARPACK.f90):

- `DPB` factor storage uses `KD + 1`
- `DGB` factor storage uses `3*KD + 1`

So for the same bandwidth, `DGB` is roughly 3x the row storage of `DPB`.

This makes SPD-capable banded problems especially attractive optimization targets.

## 3. Sparse operator duplication in Lanczos

Lanczos currently duplicates sparse operator forms:
- `KMSM`
- `KMSMn`

and may also keep source matrices around until later deallocation.

This is a major peak-memory contributor even before factorization.

## 4. Large arrays are zero-initialized row-by-row

In [ALLOCATE_LAPACK_MAT.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/ALLOCATE_LAPACK_MAT.f90):
- `ABAND`
- `BBAND`

are allocated and then explicitly zero-filled with nested loops.

That is simple and safe, but on large banded arrays it adds:
- noticeable setup time
- full-memory-touch cost before useful work even starts

This is a likely quick win area.

## Current Factor/Solve Path

## 1. SPD half-band path

The main symmetric solve path is:
- decompose in [SYM_MAT_DECOMP_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/SYM_MAT_DECOMP_LAPACK.f90)
- solve in [FBS_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/FBS_LAPACK.f90)

This uses:
- `DPBTRF`
- `DPBTRS`

This is exactly the path we should strengthen for banded optimization.

## 2. Lanczos band operator path

Lanczos still mixes:
- sparse operator storage for matrix-vector work
- band factor storage for shift-invert work

So the banded path there is not a pure band-native workflow.

That is not wrong, but it explains why conversion overhead remains significant.

## SPD Readiness Status

The current code already has a useful starting point.

In [SYM_MAT_DECOMP_LAPACK.f90](D:/fortran/mystran/MYSTRANSolver-18.0.0/Source/UTIL/SYM_MAT_DECOMP_LAPACK.f90), there is already:
- a diagonal positivity diagnostic
- reporting of nonpositive diagonal counts
- a warning path when the matrix is not clearly SPD-ready

This means a future `SPD candidate -> validate -> fallback` design fits naturally into the existing code.

That is better than inventing a new SPD classification system from scratch.

## Main Inefficiencies Found

## 1. Banded conversion is repeated and decentralized

Different workflows build banded forms independently.

This means:
- duplicated code paths
- repeated bandwidth queries
- repeated allocation logic
- repeated sparse-to-band conversion logic

There is no obvious shared banded descriptor/cache layer yet.

## 2. Banded path starts too late

The current path is mostly:
- sparse first
- convert later

That is safe, but it means conversion cost is unavoidable in the current design.

## 3. Lanczos creates too many matrix representations

For `SOL 103` banded Lanczos:
- sparse source matrices exist
- `KMSM` is assembled
- band `RFAC` is created
- nonsymmetric sparse `KMSMn` is created

This is likely the single biggest structural memory issue in the current banded ecosystem.

## 4. Fixed-width band wastes padding

Even after `RCM`, if the bandwidth is still wide, half-band/full-band storage may still waste large memory on zeros.

That does not automatically justify skyline yet, but it confirms why a feasibility gate is necessary.

## Immediate Quick Wins

Based on this audit, the best low-risk early targets are:

1. Add explicit instrumentation around:
- `BANDSIZ`
- `ALLOCATE_LAPACK_MAT`
- `BANDGEN_LAPACK_DPB`
- `BANDGEN_LAPACK_DGB`
- factor time
- solve time

2. Add estimated memory printouts for:
- `DPB` storage
- `DGB` storage
- current sparse operator duplication in Lanczos

3. Centralize banded metadata:
- bandwidth
- storage class
- SPD candidate status
- estimated memory

4. Reduce unnecessary repeated conversion where the same matrix/form is rebuilt.

5. Prefer SPD half-band path whenever safely available.

6. Review whether zero-initialization of full band arrays can be reduced or localized.

## Medium-Risk Next Step

The first more invasive but still reasonable next step is:

- direct-to-band assembly **after** final DOF ordering is known

but only for the banded branch.

This should not happen until instrumentation confirms conversion cost is a major contributor in real runs.

## Recommendation After Audit

The codebase is not ready for skyline-first development yet.

The highest-value next move is:
- instrument current banded path
- strengthen SPD half-band path
- reduce conversion duplication
- quantify memory peaks in Lanczos banded workflows

Only after those measurements should we decide whether:
- optimized fixed-width band is sufficient, or
- skyline/profile is worth the complexity.

## Best Candidate for Phase 1 Implementation

If we start coding immediately after this audit, the best first implementation target is:

### Banded Metrics and Feasibility Layer

Add one shared banded diagnostic/reporting layer that computes and exposes:
- original bandwidth
- post-RCM bandwidth
- `DPB` memory estimate
- `DGB` memory estimate
- optional future skyline estimate
- SPD candidate flag

Then use that layer to:
- inform the user in `.F06`
- guide future auto-dispatch
- provide baseline numbers for later optimization decisions

That is the lowest-risk entry point and will make every later phase more evidence-based.
