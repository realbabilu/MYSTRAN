## RSA hinge investigation

Date: 2026-07-14

Decks:
- `D:\18a\femap_RSA\nx_hinge_nx.dat`
- `D:\18a\femap_RSA\nx_hinge_nx_plot_disp.f06`
- `D:\18a\femap_RSA\nx_hinge_msc_plot_disp.f06`

Current state:
- Eigenvalues from MYSTRAN and commercial Nastran are already close for the hinge deck.
- Combined RSA displacement is no longer zero after the latest `LINK9` patch.
- The previous all-zero RSA result was caused by `LINK9` using `MPFACTOR_N6` before `OFP2` had filled it.
- The hinge deck is now under the commercial target by about `0.59x` on the dominant combined displacement components.
- The remaining discrepancy is isolated to modal participation / RSA gamma handling, not to general Lanczos extraction.

What was verified:
- `OFP2` now prints `MPF6_BASIC` diagnostics for `SCRSPEC` runs.
- `MPFACTOR_N6` is not zero. The summary previously showed `gamma(active)=0` because `LINK9` ran the RSA combine before `OFP2` had populated the participation arrays.
- `LINK9` now precomputes `MPFACTOR_N6(JVEC,1:6)` for `SOL MODES` RSA passes using the same rigid-body/basic path as `OFP2`.
- For this hinge deck, the commercial modal participation factors and MYSTRAN `MPF6_BASIC` agree well for many higher flexible modes.
- The main disagreement is concentrated in the first flexible mode.

Observed pattern:
- MYSTRAN mode 2 and above line up closely with the commercial participation-factor table.
- MYSTRAN mode 1 at about `3.634 Hz` does not match the commercial first flexible mode pattern.
- `MPFACTOR_NR(row)` remains zero for the current `SUPORT` row path, even after restoring the old `QGs_COL/DEN` fill for `SCRSPEC`.
- `TR6_MEFM` for the primary support row is:
  - `[0, 1, 0, -1.000006, 0, 5.866710]`

Commercial comparison snapshot:
- NX/MSC first flexible mode participation factors for the hinge deck are approximately:
  - `T2 ~= 0`
  - `R1 ~= +0.942146`
  - `R3 ~= +2.01974`
- MYSTRAN first flexible mode diagnostics currently give:
  - `T2 = -0.318336`
  - `R1 = +1.260539`
  - `R3 = +0.152255`
- This is not a small scaling drift. It indicates a different first-mode participation pattern.

Commercial reference consistency:
- NX and MSC combined response output for the hinge deck is already very close.
- Example near the first grids:
  - grid 2 combined `T2`
    - NX: `5.713637E-03`
    - MSC: `5.713657E-03`
  - grid 2 combined `R1`
    - NX: `7.970944E-04`
    - MSC: `7.970971E-04`
  - grid 2 combined `R3`
    - NX: `1.623207E-03`
    - MSC: `1.623213E-03`
- Therefore the commercial target is stable enough to use as the hinge reference.

Latest hinge result after precomputing `MPFACTOR_N6` in `LINK9`:
- MYSTRAN combined grid 2 values are now:
  - `T2 = 3.365727E-03`
  - `R1 = 4.743716E-04`
  - `R3 = 9.505823E-04`
- Commercial grid 2 values remain:
  - `T2 = 5.713637E-03` (NX) / `5.713657E-03` (MSC)
  - `R1 = 7.970944E-04` (NX) / `7.970971E-04` (MSC)
  - `R3 = 1.623207E-03` (NX) / `1.623213E-03` (MSC)
- This means the hinge deck moved from `0.0` response to a physically active response, but it is still low by about `41%`.

Latest gamma diagnostics after the precompute fix:
- Mode 1 now reports:
  - `freq = 3.634265 Hz`
  - `Sd = 5.702061E-03`
  - `gamma(active) = -3.183358E-01`
  - `gamma(alt-mpf) = -3.183358E-01`
- `max(|UG|)` for the retained hinge modes is now nonzero (`8.487892E-03`), confirming the RSA pass is no longer dead.

Follow-up check on `MPFACTOR_NR(row)`:
- A second isolated `LINK9` patch precomputed the `QS -> QGs -> MPFACTOR_NR` path before the RSA combine, using the same core algebra as `OFP2`.
- This did not change the hinge answer.
- `MPFACTOR_NR(row)` remains `0.0` for the primary hinge `SUPORT` row even after the precompute.
- `beam_RS.dat` remained stable under the same patch.

Interpretation of the follow-up:
- The hinge mismatch is not caused by `MPFACTOR_NR(row)` merely being computed too late.
- For this deck, the active commercial participation path is likely not the current MYSTRAN `R-set row` interpretation at all.
- The residual issue is now narrower:
  - either the commercial hinge deck effectively uses a different support-direction normalization,
  - or MYSTRAN should derive the RSA support participation from a different quantity than the present `TR6_MEFM` row logic and `MPFACTOR_NR(row)` fallback.

Note on available diagnostics:
- The commercial F06 files used here expose:
  - modal participation factors
  - UHVR response matrix
  - combined displacement output
- They do not expose a straightforward printed real-eigenvector table for the flexible mode in this deck, so direct mode-shape comparison will likely need OP2 or another deck/output setup.

What this means:
- The residual hinge mismatch is not caused by `TABLED1`, `DTI,SPECSEL`, or combination-method parsing.
- It is also not caused by all-mode SRSS accumulation in general, because the commercial and MYSTRAN higher-mode patterns are already reasonably aligned.
- The remaining issue is localized to how the hinge RSA path should obtain the support-direction participation scale after `MPFACTOR_N6` becomes available.
- The next likely missing quantity is `MPFACTOR_NR(row)` from the `QGs_COL` path, which is still zero during the RSA combine even though the deck is no longer dead.

Large-mass support finding:
- The commercial hinge deck includes a `0 Hz` rigid support mode in the participation table before the first flexible mode.
- For the current hinge deck that rigid mode is consistent with:
  - `sqrt(CONM2 mass at support grid)` times the rigid support row pattern.
- Here the deck has:
  - `CONM2 244` on `GRID 290` with mass `1000.0`
  - `sqrt(1000.0) = 31.622776...`
  - `TR6_MEFM(row=1,:) = [0, 1, 0, -1.000006, 0, 5.866710]`
- The commercial `0 Hz` participation row matches that pattern closely:
  - `T2 ~= 31.6244`
  - `R1 ~= -31.6341`
  - `R3 ~= 185.5108`
- MYSTRAN currently starts directly at the first flexible mode near `3.634 Hz`, so the present compatibility path is missing that large-mass rigid support mode altogether.

Revised interpretation:
- The hinge mismatch is not just a bad gamma choice inside the first flexible mode.
- A commercial large-mass `SUPORT + CONM2` deck can contribute through a separate rigid `0 Hz` support mode that MYSTRAN does not currently synthesize in the minimal `SCRSPEC` path.
- Because of that, using only flexible-mode `MPFACTOR_N6(component)` underestimates the combined hinge response even after the earlier all-zero bug was fixed.

Safe source note:
- `LINK9` now detects when the primary `SUPORT` grid carries `CONM2` mass but no zero-frequency rigid mode was extracted.
- In that case it emits an explicit warning to `ERR/F06` that the current path is continuing with flexible-mode-only RSA output.

Interpretation:
- The commercial hinge mismatch is not a broad parser problem anymore.
- It is a specific physics/normalization issue in the first flexible hinge mode, or in how MYSTRAN converts that mode into the single-`SUPORT` RSA gamma.
- `MPFACTOR_N6` is now a trustworthy diagnostic.
- `MPFACTOR_NR` is still not the active correction path for this deck.

Safe source changes already made:
- `D:\18a\MYSTRAN\Source\LK9\L92\OFP2.f90`
  - added `SCRSPEC` diagnostics for `MPF6_BASIC`
  - restored `MPFACTOR_NR` fill for `SCRSPEC` in the `MODES` branch as a comparison path
- `D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90`
  - improved `gamma(alt-mpf)` fallback so it uses `MPFACTOR_N6(component)` when the support-scale accumulator is empty or ineffective
  - added summary print of `RS_SUPORT_COMP_SCALE`
  - precomputes `MPFACTOR_N6` during the RSA loop before `OFP2` runs, so combined RSA output is no longer identically zero

Next recommended step:
- Compare the hinge support-direction participation formula directly against the commercial combined result, using the now-bracketed range:
  - old rigid-body mixed gamma path: about `1.27x` commercial
  - current pure `MPFACTOR_N6(component)` path: about `0.59x` commercial
- Do not spend more time trying to revive `MPFACTOR_NR(row)` for this deck unless a commercial reference explicitly shows that quantity is the intended driver.
