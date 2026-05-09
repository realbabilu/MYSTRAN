# solid_mystran_add Development Note

This bundle records the MYSTRAN `NEWSOLID` integration work.

## Scope

- Adds `PARAM,SOLIDTYP,NEWSOLID` as the production switch for the new solid formulation family.
- Keeps legacy solid behavior available with the default `PARAM,SOLIDTYP,LEGACY`.
- Accepts `PARAM,SOLIDTYP,EAS` as a temporary alias for `NEWSOLID`, but `NEWSOLID` is the intended name because not every branch is purely EAS.
- Adds CPYRA card/data plumbing and a new `PYRA.f90` element kernel.
- Keeps the linear tetra legacy element file isolated; the new CTETRA4 path is assembled separately through `CTETRA4S_SMOOTH_ASSEMBLY`.

## Active NEWSOLID Map

| Element | Active branch |
| --- | --- |
| `CHEXA8` | condensed `CHEXA8_EAS9_FROZEN` port |
| `CPENTA6` | condensed `CPENTA6_EAS9_TRIAL` port with EAS modes rotated into MYSTRAN element-local strain axes |
| `CTETRA4` | nodal-patch smooth blend, alpha `0.9`, assembled through `ESP0_FINAL`/`ESP` |
| `CPYRA5` | condensed `CPYRA5_EAS54` stiffness branch |
| `CPYRA14` | Liu composite quadratic pyramid trial |
| `CTETRA10`, `CPENTA15`, `CHEXA20` | quadratic baseline guard remains equivalent to legacy/current baseline |

## Notes

- The older B-bar scaffold is intentionally not the `NEWSOLID` implementation. The live paths are ported from the local Python element candidates.
- CPYRA5 uses the EAS54 branch for stiffness under `NEWSOLID`; mass, thermal load, center stress/strain recovery, and differential stiffness use the standard CPYRA shape-function path.
- CPYRA14 has stiffness, mass, recovery, thermal, and differential stiffness smoke coverage.
- The bundle includes copied source files under `files/` and validation scripts/reports copied from `codex_mod/solids`.
