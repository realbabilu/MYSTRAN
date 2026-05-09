# NEWSOLID User Guide

Use `PARAM,SOLIDTYP,NEWSOLID` in the bulk data section to enable the new solid formulation branches.

```nastran
PARAM,SOLIDTYP,NEWSOLID
```

To return to the original MYSTRAN behavior, omit the parameter or use:

```nastran
PARAM,SOLIDTYP,LEGACY
```

`PARAM,SOLIDTYP,EAS` is also accepted for now and is routed to `NEWSOLID`.

## What Changes

With `NEWSOLID`, the following linear solid branches are changed:

| Card/order | NEWSOLID behavior |
| --- | --- |
| `CHEXA` 8-node | EAS9 condensed stiffness |
| `CPENTA` 6-node | EAS9 condensed stiffness |
| `CTETRA` 4-node | smooth nodal-patch alpha `0.9` assembly path |
| `CPYRA` 5-node | EAS54 condensed stiffness |
| `CPYRA` 14-node | Liu composite quadratic pyramid stiffness |

The current quadratic baseline elements `CTETRA10`, `CPENTA15`, and `CHEXA20` remain matched to the existing baseline under the guard.

## Example

```nastran
SOL 101
CEND
BEGIN BULK
PARAM,SOLIDTYP,NEWSOLID
MAT1,1,3.0+7,,0.3
PSOLID,1,1
CHEXA,1,1,1,2,3,4,5,6,7,8
ENDDATA
```

For comparison studies, run the same deck twice: once without `PARAM,SOLIDTYP,NEWSOLID`, then once with it. The `.F06` displacement table can be compared directly for static cases.
