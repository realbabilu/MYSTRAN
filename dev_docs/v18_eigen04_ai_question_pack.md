# v18 Eigen04 Failure Pack (for external AI review)

## Scope

This pack is for diagnosing why the Eigen04 axis-debug decks fail in modal factorization on both CBAR and CBEAM, while a related simply-supported shaft deck passes.

## Environment

- Repository/worktree: `E:\mystran17\mystran`
- Validation binary used: `E:\mystran17\mystran\Binaries\mystran.exe`
- Validation date: 2026-04-22

## Failing decks

- `midas_eigen04_cbar_axis_debug.dat`
- `midas_eigen04_cbeam_axis_debug.dat`

## Passing comparator

- `midas_eigen04_simply_supported_shaft_axis_debug.dat`

## Failure symptoms from F06

### Banded/LAPACK path

- `midas_eigen04_cbar_axis_debug.F06`
  - `*ERROR 989: ... KLL ... LAPACK DPBSTF`
  - `LEADING MINOR OF ORDER 62 IS NOT POSITIVE DEFINITE`
- `midas_eigen04_cbeam_axis_debug.F06`
  - `*ERROR 989: ... KLL ... LAPACK DPBSTF`
  - `LEADING MINOR OF ORDER 62 IS NOT POSITIVE DEFINITE`

### INV/SuperLU diagnostic path

- `midas_eigen04_cbar_axis_debug_inv.F06`
  - `*ERROR 981: ... KMSM ... SUPERLU ... INFO = 21`
- `midas_eigen04_cbeam_axis_debug_inv.F06`
  - `*ERROR 981: ... KMSM ... SUPERLU ... INFO = 21`

Interpretation: this is likely a model/constraint indefiniteness or singularity issue, not a single-solver implementation issue.

## What is already known

- Both element formulations (CBAR and CBEAM) fail in the same case.
- Alternative deck (`simply_supported_shaft_axis_debug`) passes in same codebase.
- Broad static/eigen pack mostly passes; failure appears localized to eigen04 axis-debug modeling path.

## Candidate root-cause hypotheses

1. Constraint mapping mismatch for axis-debug conversion (over/under constraint on one rotational/translational DOF).
2. Release semantics mismatch (e.g., member end release interpreted differently between source model and Nastran-like deck).
3. Rigid-link or local-axis coupling introduces a mechanism in modal assembly.
4. Support spring semantics differ from source software assumptions.

## Focused questions for external AI

1. Given same failure for CBAR and CBEAM, where should we prioritize debugging first: converter support semantics or element stiffness/mass?
2. What minimal pre-factorization diagnostics should be added in LINK4 to identify the problematic physical DOF before DPBSTF/SuperLU fails?
3. For mapped beam models with axis-debug orientation, what are common pitfalls that create SPD loss only in modal runs (but not always static)?
4. Is there a robust temporary stabilization strategy (besides global ART_MASS) suitable for isolating mechanism DOFs without masking real mapping bugs?

## Requested output from reviewer AI

- Ranked likely root causes (top 3).
- Concrete, low-risk instrumentation points (file/subroutine level if possible).
- A minimal experiment matrix (3–5 runs) to disambiguate support/release/local-axis causes quickly.
