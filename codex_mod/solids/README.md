# solids_add

Adds selectable 3-D solid formulation options for MYSTRAN solid elements.

## User option

```text
PARAM,SOLIDTYP,NEWSOLID
```

Accepted values:

- `LEGACY` keeps the existing MYSTRAN solid formulation. This is the default.
- `NEWSOLID` is the umbrella switch for the best new formulation by solid shape/order.
- `EAS` is accepted as a temporary alias for `NEWSOLID`, but it is not the best name because not every new solid branch is EAS.

## NEWSOLID Map

Based on the Python solid package in `D:/fortran/mystran2/codex_mod/solid-python`:

- `CHEXA8`: `CHEXA8_EAS9_FROZEN`
- `CHEXA20`: `CHEXA20_BASELINE`
- `CTETRA4`: `CTETRA4_SMOOTH_BLEND_ALPHA09_TRIAL` until a more principled production T4 branch is chosen
- `CTETRA10`: `CTETRA10_BASELINE`
- `CPENTA6`: `CPENTA6_EAS9_TRIAL`
- `CPENTA15`: `CPENTA15_BASELINE`
- `CPYRA5`: `CPYRA5_EAS54_TRIAL` now active for stiffness; mass, thermal, recovery, and KED use the CPYRA5 shape-function path
- `CPYRA14`: `CPYRA14_LIU_COMPOSITE_TRIAL` now active for static/eigen mass using the Python finite-difference derivative prototype

Alternate quadratic custom-next branches:

- `CHEXA20`: `USHEXA20_OOI_TRIAL` for distortion-tolerance study, with an explicit unsymmetric stiffness/eigen policy before production use
- `CTETRA10`: `BT2_KADAPA_DISPLACEMENT_TRIAL` as the Kadapa-style displacement branch; mixed-pressure/B-bar Tet10 remains future work

## Scope

The Fortran `NEWSOLID` path should be ported directly from the element-specific Python candidates, starting with `CHEXA8_EAS9_FROZEN`.

Current status:
- `CHEXA8` has been ported to condensed EAS9 and matches the Python cantilever/MacNeal gates closely. Current rerun: `8x2x2` cantilever `uz/ref ~= 0.978`.
- `CPENTA6` has been ported to condensed EAS9 and now matches the Python reference after rotating Python/basic-basis EAS modes into MYSTRAN element-local strain axes with `TE`. Current reruns: `8x2x2` cantilever `uz/ref ~= 0.985`, `12x2x2` cantilever `uz/ref ~= 0.990`.
- `CTETRA10`, `CPENTA15`, and `CHEXA20` have been checked as quadratic baseline elements against the Python package with `codex_mod/solids/validate_quadratic_mystran.py`; no new Fortran formula branch was needed for these baseline quadratic checks.
- `CTETRA4` has been ported to the mesh-level smooth-blend alpha `0.9` candidate. Legacy `TETRA.f90` remains unchanged; the new path lives in `Source/LK1/L1E/CTETRA4S_SMOOTH_ASSEMBLY.f90` and is routed through `ESP0_FINAL`/`ESP`.
- `CPYRA5/CPYRA14` card/data plumbing and the `PYRA.f90` kernel are now present. CPYRA5 stiffness matches the Python `CPYRA5_EAS54` reference; CPYRA14 stiffness matches the Python Liu reference. Consistent mass, center-point stress/strain recovery matrices, differential stiffness `KED`, and thermal load `PTE`/`STE1` are implemented and smoke-tested.
- Other solid families remain on legacy behavior under `NEWSOLID` until their matching Python candidate is ported.

Current quadratic baseline check, n4 cantilever with total tip load `-1000`:

| case | Python T3 x1000 | MYSTRAN T3 | delta |
|---|---:|---:|---:|
| `CTETRA10` n4 | `-3.851211043e-01` | `-3.851211080e-01` | `-3.740e-09` |
| `CPENTA15` n4 | `-3.872448904e-01` | `-3.872448960e-01` | `-5.631e-09` |
| `CHEXA20` n4 | `-3.890574498e-01` | `-3.890574429e-01` | `+6.931e-09` |

Deck-order notes:
- `CPENTA15` and `CHEXA20` Python/Nastran helper midsides are ordered `bottom, top, vertical`; MYSTRAN shape routines expect `bottom, vertical, top`.
- `CPENTA15` and `CHEXA20` baseline comparison decks use full integration with `PSOLID,1,1,,3,,FULL`.
- `CTETRA10` uses the standard tetra `PSOLID,1,1` form.

Current `CTETRA4` smooth target check with `codex_mod/solids/validate_ctetra4_smooth_target.py`:

| mesh | Python standard T3 x1000 | Python smooth a09 T3 x1000 | MYSTRAN NEWSOLID T3 | MYSTRAN/standard | smooth/standard |
|---|---:|---:|---:|---:|---:|
| n4 | `-2.628385984e-02` | `-1.548880597e-01` | `-1.548880556e-01` | `5.892896` | `5.892896` |
| n8 | `-7.775117808e-02` | `-3.201114590e-01` | `-3.201114333e-01` | `4.117126` | `4.117127` |
| n12 | `-1.233556751e-01` | `-3.950453866e-01` | `-3.950453556e-01` | `3.202490` | `3.202491` |

This confirms the live Fortran `CTETRA4` `NEWSOLID` branch now matches the Python smooth-blend alpha `0.9` target to numerical roundoff. The candidate is not a local element bubble/EAS change; it replaces the assembled elementwise constant strain contribution by nodal-patch averaged strain rows blended with alpha `0.9`.

Linear tetra policy: keep `Source/EMG/EMG5/TETRA.f90` as the legacy/local element. The smooth `CTETRA4` path uses the separate assembly-level entry point `Source/LK1/L1E/CTETRA4S_SMOOTH_ASSEMBLY.f90`, so normal CTETRA behavior stays isolated.

Current CPYRA Python-vs-MYSTRAN check with `codex_mod/solids/validate_cpyra_mystran.py`:

| case | Python T3 x1000 | MYSTRAN T3 | delta |
|---|---:|---:|---:|
| `CPYRA5_EAS54` n4 | `-1.275452610e-01` | `-1.275452667e-01` | `-5.648e-09` |
| `CPYRA14` n2 | `-4.673392309e-02` | `-4.673392080e-02` | `+2.294e-09` |

Current CPYRA SOL 5 smoke check with `codex_mod/solids/validate_cpyra_buckling_smoke.py`:

| case | run | eigen table |
|---|---:|---:|
| `CPYRA5` n2 | ok | true |
| `CPYRA14` n1 | ok | true |

Current CPYRA thermal smoke check with `codex_mod/solids/validate_cpyra_thermal_smoke.py`:

| case | run | fatal |
|---|---:|---:|
| `CPYRA5` | ok | false |
| `CPYRA14` | ok | false |

The working scaffold is wired into:

- `HEXA`
- `PENTA`
- `TETRA`
- `PYRA`

The `NEWSOLID` target is intentionally broader and should be ported element-by-element from the Python candidates above.

## Patch marker convention

Patch-bundle snippets for this package should use:

```fortran
! --- solids_add begin --- !
...
! --- solids_add end --- !
```
