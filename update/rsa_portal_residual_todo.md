Date: July 14, 2026

# RSA Portal Residual-Structure TODO

## Purpose

This note is for design discussion before changing the current response-spectrum
implementation for portal/base-motion decks that use constrained `SUPORT`
definitions.

The immediate goal is:

- improve commercial-Nastran compatibility for portal/base-excitation RSA decks
- without regressing the beam/frame RSA decks that currently compare reasonably
  well to SAP2000 and MSC

## Problem statement

Current portal-style decks such as:

- `D:\18a\femap_RSA\portal_2d_RS.dat`
- `D:\18a\femap_RSA\portal_3d_RS.dat`

fall into:

- `SCRSPEC COMPATIBILITY SUMMARY (MINIMAL MODES PATH)`

and print:

- `No matching R-set row was found for the primary SUPORT DOF; fallback global-component weighting was used.`

This strongly suggests the current path is not yet reproducing the commercial
residual-structure/base-support excitation treatment.

Observed consequence:

- beam large-mass RSA decks can compare well
- portal/base-motion decks can differ by orders of magnitude in combined
  displacement/velocity/acceleration

## Non-goal

Do not globally replace the current RSA combination path for all decks.

That would risk breaking:

- `problem_1_020_rsa_srss_new.dat`
- `problem_1_024_rsa_srss_new.dat`
- `problem_1_025_rsa_srss_new.dat`
- `beam_RS.dat`

These decks should remain on the currently working path unless a more specific
portal/base-excitation branch is selected.

## Safe interpretation

At the moment, the portal decks should be treated as:

- parser-compatible
- solver-runnable
- numerically limited for commercial portal/base-excitation comparison

until the residual/support path is implemented and validated.

## Working hypothesis

The mismatch is likely not from:

- `TABLED1`
- `DTI,SPECSEL`
- `PARAM,GRAV`
- simple SRSS combination

because those same basic ingredients already behave acceptably for the beam/frame
RSA cases.

The mismatch is more likely from one or more of:

- missing residual-structure support handling
- wrong excitation basis for constrained `SUPORT`
- wrong use of `MPFACTOR_N6` fallback when an `R-set` row is not available
- incomplete treatment of rigid-body/large-mass support behavior

## Current source hotspots

Primary source to review:

- `D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90`

Relevant current logic areas:

- `RSA_COMP_ROW`
- `MPFACTOR_N6`
- `MPFACTOR_NR`
- `TR6_MEFM`
- `RBGLOBAL_GSET`
- `RSA_GAMMA`
- `RSA_MODE_GAMMA_ALT`
- `RSA_MODE_SCALE`

Useful historical snapshot:

- `E:\python\fem_manual\mystran3_project\codex_mod\response_spectra_add_v2\source_snapshot\LINK9.f90`

## Proposed branch rule

Add a more specific RSA branch only when conditions indicate portal/base-motion
support handling is needed.

Candidate trigger:

1. `SCRSPEC == 'Y'`
2. `SOL_NAME(1:5) == 'MODES'` or equivalent SEMODES path
3. primary `SUPORT` DOF exists
4. primary `SUPORT` DOF is constrained or behaves as constrained support excitation
5. no matching `R-set` row is found for the primary `SUPORT` DOF

If these conditions are not met, keep the current working path.

## Implementation options

### Option A: guarded residual-support branch

Add a new branch in `LINK9` that:

- builds/supports a residual-structure excitation basis for the selected support
  direction
- computes modal participation against that basis
- avoids the current plain `MPFACTOR_N6(comp)` fallback for these portal cases

Pros:

- lowest regression risk
- easiest to A/B against current output

Cons:

- more branching in `LINK9`

### Option B: upgrade the existing fallback path

Replace the `MPFACTOR_N6` fallback logic globally when no `R-set` row is found.

Pros:

- simpler surface API

Cons:

- high regression risk for `1-020`, `1-024`, `1-025`, and `beam_RS`
- harder to isolate numerical changes

Recommendation:

- prefer Option A first

## Required validation matrix

Before merging any residual-support implementation, rerun:

- `D:\18a\femap_RSA\beam_RS.dat`
- `D:\18a\femap_RSA\beam_RS_plot.dat`
- `D:\18a\femap_RSA\problem_1_020_rsa_srss_new.dat`
- `D:\18a\femap_RSA\problem_1_024_rsa_srss_new.dat`
- `D:\18a\femap_RSA\problem_1_025_rsa_srss_new.dat`
- `D:\18a\femap_RSA\portal_2d_RS.dat`
- `D:\18a\femap_RSA\portal_3d_RS.dat`

Minimum acceptance intent:

- portal decks move materially closer to MSC
- `1-024` and `1-025` do not move away from SAP-style references
- `beam_RS` does not regress materially

## Output/diagnostic requirements

Keep strong diagnostics in `F06`/`ERR` while this is under development.

Needed diagnostics:

- which RSA branch was selected
- whether the deck used the old minimal-modes path or the new residual-support path
- primary `SUPORT` grid/component
- whether an `R-set` row was found
- whether `MPFACTOR_N6`, `MPFACTOR_NR`, or a residual-support basis drove gamma

## Open technical questions

1. Should constrained `SUPORT` portal decks use a support basis derived from
   `RBGLOBAL_GSET`/`TR6_MEFM`, or do they need a separate residual operator?

2. Is the commercial result closer to:
   - support-direction participation from `MPFACTOR_NR`
   - large-mass rigid-body residual vectors
   - a mixed residual + flexible modal combination

3. Can the branch be identified purely from missing `R-set` row, or should it
   also check:
   - `CONM2` mass at the support grid
   - constrained support dof pattern
   - presence/absence of rigid zero mode

## Recommended next coding step

Do this first:

1. add an explicit internal branch label for
   - `RSA_PATH = MINIMAL_MODES`
   - `RSA_PATH = RESIDUAL_SUPPORT`
2. keep all current numerics unchanged
3. route only portal/base-motion candidates to a new experimental code path
4. log both old and new gamma ingredients side-by-side for comparison

That gives a safe staging point before changing the final combined RSA numbers.
