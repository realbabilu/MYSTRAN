# SOLIDTYP NEWSOLID Map

Use this global deck option for the new solid formulation family:

```text
PARAM,SOLIDTYP,NEWSOLID
```

`PARAM,SOLIDTYP,EAS` is accepted as a temporary alias, but `NEWSOLID` is the clearer name because the selected formulation differs by element shape and order.

## Python Source Map

| MYSTRAN family | order | Python candidate | status |
|---|---:|---|---|
| `CHEXA` | 8 | `solid3d_chexa8_eas9_frozen.py` / `CHEXA8_EAS9_FROZEN` | main linear brick candidate |
| `CHEXA` | 20 | `solid3d_chexa20.py` / `CHEXA20_BASELINE` | validated quadratic brick baseline |
| `CTETRA` | 4 | `solid3d_ctetra4_smooth.py` / `CTETRA4_SMOOTH_BLEND_ALPHA09_TRIAL` | best current linear T4 trial |
| `CTETRA` | 10 | `solid3d_ctetra10.py` / `CTETRA10_BASELINE` | production tetra baseline |
| `CPENTA` | 6 | `solid3d_cpenta6_eas.py` / `CPENTA6_EAS9_TRIAL` | main linear wedge candidate |
| `CPENTA` | 15 | `solid3d_cpenta15.py` / `CPENTA15_BASELINE` | production wedge baseline |
| `CPYRA` | 5 | `solid3d_cpyra5_eas.py` / `CPYRA5_EAS54_TRIAL` | active stiffness branch; mass/thermal/recovery/KED use the standard CPYRA5 shape-function path |
| `CPYRA` | 14 | `solid3d_cpyra14_liu.py` / `CPYRA14_LIU_COMPOSITE_TRIAL` | active static/eigen-mass prototype; still uses finite-difference derivatives |

## Alternate Quadratic Custom Next

These are not the default `NEWSOLID` production map yet, but should be kept as custom next branches:

| MYSTRAN family | order | Python candidate | custom-next role |
|---|---:|---|---|
| `CHEXA` | 20 | `solid3d_ushexa20_ooi.py` / `USHEXA20_OOI_TRIAL` | alternate quadratic hex for distortion-tolerance study; unsymmetric on distorted meshes, so eigen/buckling policy must be deliberate |
| `CTETRA` | 10 | `solid3d_bt2_kadapa.py` / `BT2_KADAPA_DISPLACEMENT_TRIAL` | alternate quadratic tetra/Tet10-style Kadapa displacement branch; mixed-pressure/B-bar variant remains future work |

## Important Distinction

Do not use the earlier B-bar volumetric-projection scaffold as the `NEWSOLID` implementation. `NEWSOLID` should be ported directly from the Python candidates listed above, starting with `CHEXA8_EAS9_FROZEN`.

Current Fortran status:
- `CHEXA8`: ported to condensed EAS9 and accepted on the current static checks; current `8x2x2` cantilever rerun is `uz/ref ~= 0.978`.
- `CPENTA6`: ported to condensed EAS9 and accepted on the current static checks after fixing the EAS basis mismatch. Python defines the EAS modes in the basic Cartesian strain basis, while MYSTRAN solid stiffness is assembled in element-local axes, so the Fortran port rotates each enhanced strain tensor with `TE` before condensation. Current reruns: `8x2x2` cantilever `uz/ref ~= 0.985`, `12x2x2` cantilever `uz/ref ~= 0.990`.
- `CTETRA4`: ported to the Python smooth-blend alpha `0.9` nodal-patch formulation in `Source/LK1/L1E/CTETRA4S_SMOOTH_ASSEMBLY.f90`. This is an assembly-level correction on top of legacy Tet4 local stiffness, not a `TETRA.f90` local element branch.
- `CTETRA10`, `CPENTA15`, `CHEXA20`: validated against the Python baseline package with `codex_mod/solids/validate_quadratic_mystran.py`. Current n4 cantilever deltas are at numerical-noise level:

| case | Python T3 x1000 | MYSTRAN T3 | delta |
|---|---:|---:|---:|
| `CTETRA10` n4 | `-3.851211043e-01` | `-3.851211080e-01` | `-3.740e-09` |
| `CPENTA15` n4 | `-3.872448904e-01` | `-3.872448960e-01` | `-5.631e-09` |
| `CHEXA20` n4 | `-3.890574498e-01` | `-3.890574429e-01` | `+6.931e-09` |

Quadratic deck notes:
- `CPENTA15` and `CHEXA20` Python/Nastran helper midsides are ordered `bottom, top, vertical`; MYSTRAN shape routines expect `bottom, vertical, top`.
- `CPENTA15` and `CHEXA20` baseline comparison decks use `PSOLID,1,1,,3,,FULL`; `CTETRA10` uses `PSOLID,1,1`.

- `CTETRA4`: local EAS/bubble diagnostics are intentionally not used; the active candidate is mesh-level smooth-blend alpha `0.9`.
- `CPYRA5/CPYRA14`: active through new `CPYRA` card/data plumbing and `Source/EMG/EMG5/PYRA.f90`. CPYRA5 stiffness matches `CPYRA5_EAS54`; CPYRA14 stiffness matches `CPYRA14_LIU_COMPOSITE_TRIAL`. Center-point stress/strain recovery and differential stiffness are present for SOL 5 buckling smoke. Thermal `PTE`/`STE1` is implemented and smoke-tested.
- Other shapes/orders should stay legacy under `NEWSOLID` until their matching Python candidate is ported.

Current CPYRA check with `codex_mod/solids/validate_cpyra_mystran.py`:

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

Current `CTETRA4` smooth target check with `codex_mod/solids/validate_ctetra4_smooth_target.py`:

| mesh | Python standard T3 x1000 | Python smooth a09 T3 x1000 | MYSTRAN NEWSOLID T3 | MYSTRAN/standard | smooth/standard |
|---|---:|---:|---:|---:|---:|
| n4 | `-2.628385984e-02` | `-1.548880597e-01` | `-1.548880556e-01` | `5.892896` | `5.892896` |
| n8 | `-7.775117808e-02` | `-3.201114590e-01` | `-3.201114333e-01` | `4.117126` | `4.117127` |
| n12 | `-1.233556751e-01` | `-3.950453866e-01` | `-3.950453556e-01` | `3.202490` | `3.202491` |

Fortran implication: the live `CTETRA4` `NEWSOLID` path now follows `CTETRA4_SMOOTH_BLEND_ALPHA09_TRIAL`. The implementation constructs each nodal patch from adjacent tetra elements and assembles the alpha correction through `ESP0_FINAL` and `ESP`.

Isolation rule: legacy `Source/EMG/EMG5/TETRA.f90` stays unchanged. The new smooth Tet4 implementation starts in `Source/LK1/L1E/CTETRA4S_SMOOTH_ASSEMBLY.f90` and should only be routed for linear `CTETRA4` when `PARAM,SOLIDTYP,NEWSOLID` is active.
