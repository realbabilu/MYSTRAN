# CELAS1 Grounded Compatibility

## Scope

This note records the compatibility work for MSC-style grounded `CELAS1` entries
used in the shell validation decks, especially `shell/prob_2_002h_thin.dat`.

## What changed

- `BD_CELAS1` now accepts the short MSC form:
  - `CELAS1, EID, PID, G1, C1`
- The parser no longer forces the deck to provide a second grid/component pair.
- Grounded shorthand is normalized internally so the solver path can treat it as
  a spring-to-ground case instead of a malformed `CELAS1`.
- `ELAS1`, `GET_ELGP`, and the element-data path were hardened so a grounded
  `CELAS1` does not dereference a synthetic `GRID 0` entry.

## Current behavior

- Regular two-grid `CELAS1` behavior is unchanged.
- The grounded shorthand now runs through the MYSTRAN solver path instead of
  failing during bulk-data validation.
- The change is intentionally narrow and only targets the MSC shorthand used in
  the validation deck.
- Smoke test on `shell/prob_2_002h_thin.dat` now reaches normal termination in
  MYSTRAN instead of aborting on `CELAS1`/`GRID 0` validation.

## Limits

- This is not a full generalization of every possible CELAS variant.
- The ground shorthand is handled as a compatibility path, not as a redesign
  of spring element theory.
