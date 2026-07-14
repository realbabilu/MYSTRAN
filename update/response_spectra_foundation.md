# Response Spectrum Foundation

Date: July 13, 2026

## Scope landed in this pass

- Added parser/state foundation for Nastran-style response spectrum bulk data:
  - `PARAM,RSPECTRA,*`
  - `PARAM,SCRSPEC,*`
  - `PARAM,RSTYPE,*`
  - `DLOAD`
  - `RLOAD1`
  - `FREQ1`
  - `TABLED1`
  - `DTI,SPECSEL`
- Reused the earlier `response_spectra_add_v2` direction only where it was narrow and low-risk.
- Kept current solver behavior unchanged for non-RSA runs.
- Added `SUPORT` capture for a minimal Nastran-style excitation direction path.
- Split `TABLED1` storage by table ID so multiple damping-dependent spectra no longer collapse into one global list.
- Added a first solve-side `LINK9` compatibility pass for `SOL 103/SEMODES`:
  - forces internal `MPFACTOR` calculation when `SCRSPEC` needs it
  - forms a minimal combined displacement result in G-set space
  - accepts `PARAM,OPTION,SRSS|CQC|ABS|NRL`
  - currently executes `SRSS`, `CQC`, and `ABS`
  - currently accepts `NRL` as compatibility syntax with warning/fallback
  - writes a summary line to `ERR/F06`
  - writes combined displacement to `F06`
  - writes combined displacement to `OP2` when `DISPLACEMENT(PLOT)` or other OP2-enabled routing is active
  - writes combined displacement to `NEU` when neutral displacement output is active
- For `MODES + SCRSPEC`, `OFP2` now uses the legacy SPC-force / `QGs` participation-factor route instead of the normal rigid-body `MGG*RBG` modal-effective-mass route.
- Relaxed the legacy `R-SET must have at least 6 DOF` fatal only for the narrow case:
  - `SOL_NAME(1:5) == 'MODES'`
  - `PARAM,SCRSPEC,*` active
  - still emits a warning

## Current compatibility intent

- `PARAM,RSPECTRA,0` is accepted as a generation-mode compatibility trigger.
- `PARAM,SCRSPEC,0` is accepted as an application-mode compatibility trigger.
- `PARAM,RSTYPE,FRQG|FRQA|PERG|PERA` controls table interpretation.
- `DTI,SPECSEL` rows are stored so later `SOL 103/SEMODES` RSA application logic can use Nastran-style damping/table selection.

## Not completed in this pass

- Full response spectrum solve flow in `LINK4/LINK5/LINK9`
- Nastran-style `DTI,SPSEL` generation workflow
- Combined RSA stress/strain/force postprocessing beyond grid vectors
- Full modal combination options (`10PCT`, directional 100/30)
- Stress/strain/force response-spectrum postprocessing
- Full validation against `D:\18a\femap_RSA`

## Notes

- This is still a compatibility-transition pass, not a complete RSA implementation.
- The current `SCRSPEC` solve path is intentionally narrow:
  - `SOL_NAME(1:5) == 'MODES'`
  - first `DTI,SPECSEL` pair only
  - `SUPORT`-derived excitation direction
  - acceleration-like interpretation for spectrum ordinates
  - no directional 100/30 combination yet
- Partial `SUPORT`/R-set definitions are accepted only inside that same narrow compatibility path.
- Requested modal stress/strain/force outputs still remain on the normal modal path; they are not yet converted into combined RSA result tables.
- Combined RSA grid-vector outputs now cover:
  - displacement
  - velocity
  - acceleration
- Constraint-force and element-result families are still on the normal modal path in this compatibility pass.
- Combined RSA velocity/acceleration currently use the narrow compatibility interpretation:
  - `V = omega * U`
  - `A = omega^2 * U`
  - no separate absolute/relative spectrum-family split beyond the active `DTI,SPECSEL` kind yet
- Verification against `D:\18a\femap_RSA\beam_rs_msc.f06` on `beam_RS.dat` now gives:
  - displacement at grid 2:
    - `TY`: MYSTRAN `4.172221E-01`, MSC `4.173022E-01`
    - `RZ`: MYSTRAN `1.041884E-01`, MSC `1.041883E-01`
  - velocity at grid 2:
    - `TY`: MYSTRAN `1.330986E+01`, MSC `1.331114E+01`
    - `RZ`: MYSTRAN `3.323730E+00`, MSC `3.323408E+00`
  - acceleration at grid 2:
    - `TY`: MYSTRAN `4.246000E+02`, MSC `4.246000E+02`
    - `RZ`: MYSTRAN `1.060308E+02`, MSC `1.060104E+02`
- Plot-route verification on July 14, 2026:
  - `D:\18a\femap_RSA\beam_RS_plot.dat`
  - `D:\18a\femap_RSA\beam_RS_cqc_plot.dat`
  - both decks now finish normally and write:
    - combined displacement vector
    - combined velocity vector
    - combined acceleration vector
  - current artifact sizes after rerun:
    - `beam_RS_plot.OP2`: `7032` bytes
    - `beam_RS_cqc_plot.OP2`: `7032` bytes
    - `beam_RS_plot.NEU`: `31523` bytes
    - `beam_RS_cqc_plot.NEU`: `31523` bytes
- Ordinary `SOL 103` modal effective mass / participation output remains on the existing rigid-body branch; the `QGs` branch is isolated to `SCRSPEC`.
- Follow-up on July 14, 2026:
  - a zero-response bug remained when `SCRSPEC` used a single `SUPORT` DOF
  - root cause: the `QGs`-derived `MPFACTOR_N6` values are not a reliable basis for a support-row excitation direction
  - fix: `LINK9` now rebuilds the active support influence vector from `TR6_MEFM(row,:)` and `RBGLOBAL_GSET`, then evaluates
    `gamma = phi^T * M * r / genmass` directly for the minimal SRSS accumulation path
  - this change is still isolated to the `MODES + SCRSPEC` compatibility summary path
  - verification deck: `D:\18a\femap_RSA\nx_hinge_msc.dat`
    now reports a nonzero `max(|UG|)` in the SCRSPEC summary
- Follow-up on July 14, 2026, later pass:
  - `TABLED1` continuation parsing now stops cleanly on `ENDT`
  - tiny `SCRSPEC` modal models no longer die in adaptive ARPACK setup; `LINK4` falls back to the condensed dense reference solve
  - combined displacement now routes to:
    - `F06` summary/detail
    - `OP2` as a static-style combined `OUGV1` result when displacement plot output is enabled
    - `NEU` through the normal grid-vector writer when neutral displacement output is enabled
  - verification decks:
    - `D:\18a\femap_RSA\beam_RS.dat`
    - `D:\18a\femap_RSA\beam_RS_cqc.dat`
    - `D:\18a\femap_RSA\beam_RS_plot.dat`
    - `D:\18a\femap_RSA\beam_RS_cqc_plot.dat`
  - numeric verification against `D:\18a\femap_RSA\beam_rs_msc.f06` shows the current combined displacement gap is not from eigen extraction:
    - mode frequencies and modal participation factors already match MSC closely
    - remaining displacement underprediction is explained by two still-missing RSA application pieces:
      - `DLOAD` scale-factor application (`386.` in `beam_RS.dat`)
      - damping interpolation across multiple `DTI,SPECSEL` pairs using the active structural damping (`SDAMP/TABDMP1`)
    - for this deck the missing global factor is:
      - `386.0 * 1.1 = 424.6`
    - this matches the observed ratio between MSC and current MYSTRAN combined displacements almost exactly
- Follow-up on July 14, 2026, RSA scaling pass:
  - `DLOAD` parsing now keeps the actual `S0*Si` scale factors instead of only remembering referenced IDs
  - `SCRSPEC` application in `LINK9` now applies that `DLOAD` scale to the interpolated spectral ordinate
  - `DLOAD` lookup now supports both:
    - `DLOAD -> RLOAD1 -> TABLED1`
    - `DLOAD -> DTI,SPECSEL line-id -> TABLED1`
  - added minimal `SDAMP` and `TABDMP1` support for the narrow RSA compatibility path:
    - `SDAMP = sid` in Case Control is captured
    - `TABDMP1 sid` stores the first valid damping ordinate as the active structural damping target
    - `LINK9` now blends two `DTI,SPECSEL` spectra on the same line by damping when the target lies between them
  - verification against `D:\18a\femap_RSA\beam_rs_msc.f06`:
    - MSC combined displacement at grid 2: `TY=4.173022E-01`, `RZ=1.041883E-01`
    - MYSTRAN combined displacement at grid 2 after the fix: `TY=4.172221E-01`, `RZ=1.041884E-01`
  - verification decks:
    - `D:\18a\femap_RSA\beam_RS.dat`
    - `D:\18a\femap_RSA\beam_RS_plot.dat`
    - `D:\18a\femap_RSA\beam_RS_cqc.dat`
    - `D:\18a\femap_RSA\beam_RS_cqc_plot.dat`
  - current limit of the damping path:
    - `TABDMP1` handling is still intentionally minimal and only captures the first valid damping ordinate for the active SID
  - reporting cleanup:
    - when `SCRSPEC` uses damping interpolation between two `DTI,SPECSEL` tables, the `F06/ERR` summary now reports
      `TABLED1 blend=t1/t2` and the weights `w1/w2` instead of pretending the result came from one table only
  - expected single-table behavior remains for decks that do not define an active structural damping target:
    - `D:\18a\femap_RSA\beam_RS_easy.dat` stays on the first stored `SPECSEL/TABLED1` path and therefore still reports
      the narrow legacy-style `TABLED1=1, damping=5%` summary
- Follow-up on July 14, 2026, `TABDMP1` frequency-dependent pass:
  - `TABDMP1` storage no longer collapses to one damping value per SID
  - `BD_TABDMP1` now keeps all `(frequency,damping)` pairs and ignores trailing `ENDT` tokens safely
  - `LINK9` now evaluates structural damping per retained mode frequency with linear interpolation across the active
    `TABDMP1` pairs
  - mode-by-mode `SPECSEL/TABLED1` selection now follows that interpolated damping instead of a single global damping
    value
  - `CQC` uses the average of the two participating modal damping values in the current compatibility pass
  - summary/reporting change:
    - for frequency-dependent damping, `F06/ERR` now report `TABDMP1 SID`, `damping_min`, and `damping_max`
      instead of one fake constant damping value
  - verification decks:
    - `D:\18a\femap_RSA\beam_RS.dat`
    - `D:\18a\femap_RSA\nx_hinge_nx.dat`
  - observed result on `nx_hinge_nx.dat`:
    - summary now reports `TABDMP1 SID=10, damping_min=2.003634E-02, damping_max=4.000000E-02`
    - previous narrow path incorrectly stayed fixed at `0.02`
- Follow-up on July 14, 2026, RSA SUPORT isolation pass:
  - TSET_PROC no longer converts SUPORT into the structural R-set when SOL MODES + PARAM,SCRSPEC is active
  - intent: keep SUPORT as response-spectrum excitation metadata without altering the validated modal basis
  - LINK9 now guards all MPFACTOR_NR(:,row) accesses with an explicit column-size check so zero-column RSA support cases do not crash
  - status: rebuild completed; SAP-style RSA decks are being rechecked against their modal-only baselines
