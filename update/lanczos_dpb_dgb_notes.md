# Lanczos DPB/DGB Notes

This note records the intended direction for the `DPB` / `DGB` choice used by the Lanczos-family eigen solvers in `mystran18a`.

## What DPB and DGB Are

These are LAPACK matrix-storage / factorization flavors used internally by the banded eigen solve path.

- `DPB`
  - symmetric positive-definite band
  - lower storage cost
  - usually the preferred default when assumptions hold

- `DGB`
  - general band
  - more permissive / more robust for non-ideal cases
  - more expensive than `DPB`

These are implementation-level solver details, not classic user-facing Nastran result-request concepts.

## Compatibility Goal

Move the preferred user control away from `EIGRL` continuation details and toward a global `PARAM` fallback.

Desired direction:

- `PARAM,LANCMETH,<backend>,<DPB|DGB>`

Examples:

- `PARAM,LANCMETH,ARPACK,DPB`
- `PARAM,LANCMETH,ARPACK,DGB`

## Default Direction

Preferred default:

- backend: `ARPACK`
- matrix type: `DPB`

Equivalent explicit form:

- `PARAM,LANCMETH,ARPACK,DPB`

Reason:

- this is the cleaner default for a Nastran-like user experience
- `DGB` remains available as an advanced fallback when needed

## Backward Compatibility

Existing decks may still place `DPB` / `DGB` on `EIGRL` continuation data.

Compatibility policy:

1. keep reading legacy `EIGRL` placement for now
2. emit a warning that the preferred placement is now the `PARAM` fallback
3. keep the old syntax working so legacy decks do not break

## Why This Split Helps

`EIGR` / `EIGRL` should primarily describe:

- search range
- number of roots
- normalization
- extraction setup

while `DPB` / `DGB` is better treated as a solver-policy fallback choice.

That makes decks easier to understand and closer to how users expect Nastran-style input to be organized.

## Recommended User Rule

For ordinary use:

- do not put `DPB` / `DGB` on `EIGRL`
- set it once with `PARAM,LANCMETH,...`

For advanced or legacy decks:

- legacy `EIGRL` placement may remain temporarily supported
- but it should be considered deprecated placement
