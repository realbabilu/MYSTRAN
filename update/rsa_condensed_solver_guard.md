## RSA condensed solver guard

Date: 2026-07-14

Scope:
- `D:\18a\MYSTRAN\Source\LK4\EIGRL_EXTRACT_SOLVERS.F90`

Reason:
- Commercial `nx_hinge_nx` / `nx_hinge_msc` use the low flexible mode at about `3.634265 Hz`.
- MYSTRAN default ARPACK path on the original hinge deck matches that target.
- MYSTRAN condensed extractors (`LANCMETH=DENSE`, `FEAST`, `SUBSP`) currently do not.
- On the same hinge deck they jump to about `17.16302 Hz` as the first mode.

Safe action taken:
- Added an isolated fallback guard:
  - when `SCRSPEC='Y'` and `RS_NUM_SUPORT == 1`
  - `DENSE`, `FEAST`, and `SUBSP` now fall back to ARPACK with a warning

Why this is safe:
- It only affects the current response-spectrum compatibility path with a single `SUPORT` direction.
- It does not change normal modal decks outside that narrow case.
- It avoids silently producing the known-wrong condensed modal basis for hinge-style RSA decks.

Intent:
- This is a compatibility guard, not the final physics fix.
- A later patch can re-enable condensed extractors after their zero/positive-mass condensation matches the commercial low-mode behavior for `SCRSPEC + SUPORT + CONM2`.

Wrapper note:
- `D:\18a\femap_RSA\nx_hinge_nx_bulk.inc` was also cleaned so the helper wrapper decks no longer duplicate `DTI,SPECSEL`.
- The wrapper now starts at the first `TABLED1` block, which matches the intent of the custom `LANCMETH=*` comparison decks.
