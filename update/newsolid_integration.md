# NEWSOLID Integration

Date: 2026-07-13

This branch now includes the `solid_mystran_add` integration needed to benchmark the newer linear solid paths against the legacy MYSTRAN behavior.

## New Parameter

Use:

```nastran
PARAM,SOLIDTYP,NEWSOLID
```

Accepted values:

- `LEGACY` : keep the original branch behavior
- `NEWSOLID` : enable the new solid branches
- `EAS` : temporary compatibility alias, routed to `NEWSOLID`

If `PARAM,SOLIDTYP` is omitted, the default remains `LEGACY`.

## Bulk Data Support

`CPYRAM` is now accepted with standard Nastran spelling.

Compatibility alias kept for now:

- `CPYRA`

## NEWSOLID Coverage

With `PARAM,SOLIDTYP,NEWSOLID`, the branch now enables:

- `CHEXA8` : EAS9 condensed stiffness path
- `CPENTA6` : EAS9 condensed stiffness path
- `CTETRA4` : smoothed assembly correction path
- `CPYRAM5` : EAS54 condensed stiffness path
- `CPYRAM14` : quadratic pyramid path from the integrated pyramid implementation

The restart persistence path was also updated, so `SOLIDTYP` survives restart I/O through `L1A`.

## Build Status

Integrated files compile successfully in:

- `D:\18a\MYSTRAN\build`
- output binary: `D:\18a\MYSTRAN\Binaries\mystran.exe`

## Notes

- `LEGACY` remains the safe default for existing decks.
- The intent of this integration is side-by-side benchmarking, not silent behavior replacement.
- Benchmarking should compare the same deck under both:

```nastran
PARAM,SOLIDTYP,LEGACY
```

and

```nastran
PARAM,SOLIDTYP,NEWSOLID
```
