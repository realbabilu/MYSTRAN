Date: July 14, 2026

## Scope

This note records the current commercial-comparison status for the portal-style
response-spectrum decks in:

- `D:\18a\femap_RSA\portal_2d_RS.dat`
- `D:\18a\femap_RSA\portal_3d_RS.dat`

Reference commercial decks/results:

- `D:\18a\femap_RSA\portal_2d_RS_msc.dat`
- `D:\18a\femap_RSA\portal_3d_RS_msc.dat`
- `D:\18a\femap_RSA\portal_2d_rs_msc.f06`
- `D:\18a\femap_RSA\portal_3d_rs_msc.f06`

## Current finding

The current MYSTRAN `SCRSPEC` implementation is still a:

- `SCRSPEC COMPATIBILITY SUMMARY (MINIMAL MODES PATH)`

for these portal/base-excitation cases.

The key diagnostic now visible in F06/ERR is:

- `No matching R-set row was found for the primary SUPORT DOF; fallback global-component weighting was used.`

That case occurs when the chosen `SUPORT` excitation dof does not map cleanly into
the rigid-support/residual-structure path used by commercial Nastran.

## Practical interpretation

Beam large-mass decks can still compare well under the current path, but portal-style
base-motion decks can differ substantially because the current implementation does not
yet reproduce the commercial residual-structure support treatment for constrained
support excitation.

This means the following deck family should currently be treated as:

- parser-compatible
- solver-runnable
- numerically limited for commercial portal/base-excitation comparison

until a residual-structure support path is implemented.

## What was changed

`LINK9` now emits a stronger warning in both `ERR` and `F06` when this condition is
detected, so the output is less likely to be mistaken for a fully commercial-equivalent
response-spectrum result.

## Not yet fixed

Not fixed in this pass:

- commercial-style residual-structure support handling for constrained `SUPORT`
- portal/base-motion parity against MSC/NX for combined displacement/velocity/acceleration
