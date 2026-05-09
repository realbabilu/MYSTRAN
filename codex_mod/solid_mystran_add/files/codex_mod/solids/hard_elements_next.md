# Hard Solid Elements Next

## CTETRA4 Smooth-Blend

Python source:

- `D:/fortran/mystran2/codex_mod/solid-python/solid3d_ctetra4_smooth.py`
- default target: `CTETRA4_SMOOTH_BLEND_ALPHA09_TRIAL`

This cannot be ported as a normal local `TETRA.f90` EAS/bubble branch. The Python candidate builds a global nodal-patch stiffness:

```text
K = (1 - alpha) * K_standard + alpha * sum_a V_a * Bbar_a.T * D * Bbar_a
```

where each patch `a` is one grid and `Bbar_a` is the volume-weighted average of all adjacent CTETRA4 constant-strain `B` rows. That means one patch contribution couples every node belonging to all CTETRA4 elements connected to grid `a`. The matrix can have a larger stencil than a single tetra element.

Required Fortran direction:

1. Keep standard `TETRA.f90` unchanged for local `KE` generation.
2. Add the new assembly-level path in `Source/LK1/L1E/CTETRA4S_SMOOTH_ASSEMBLY.f90`; this file is the dedicated home for the smooth Tet4 implementation.
3. Reuse the existing `GRID_ELEM_CONN_ARRAY` connectivity table if it is available early enough; otherwise build a CTETRA4-only grid-to-element patch map.
4. For each CTETRA4 element, compute and store:
   - global DOF ids for its four nodes,
   - element volume,
   - basic/global constant-strain `B` matrix.
5. For each grid patch, accumulate `sum(volume * B_global)` over adjacent CTETRA4 elements, divide by total patch volume, then assemble `total_volume / 4 * Bbar.T * D * Bbar`.
6. Blend with alpha `0.9` against the standard CTETRA4 matrix.

Implementation policy:

- Do not modify the numerical behavior of legacy `TETRA.f90`.
- Do not silently change normal `CTETRA` decks.
- Route this only when `PARAM,SOLIDTYP,NEWSOLID` is active and the element is linear `CTETRA4`.
- Keep stress recovery on the legacy Tet4 path until a matching smoothed recovery rule is added.

Current Fortran implementation:

- `CTETRA4S_SMOOTH_ASSEMBLY` has an explicit interface and an isolated entry point for `COUNT` and `ADD` actions.
- `ESP0_FINAL` calls the `COUNT` action after normal element topology counting so the larger nodal-patch stencil is included in `LTERM_KGG`.
- `ESP` calls the `ADD` action on the normal KGG pass only; buckling differential stiffness KGGD remains on the existing element path.
- Internal `ADD_STF_TERM` counts or adds terms into the LINK1 `STF3/STFKEY` linked-list stiffness storage without changing legacy `TETRA.f90`.
- Internal `CTETRA4_B_GLOBAL` reads grid coordinates/translation G-set DOFs and builds the constant-strain basic/global `B` matrix for linear Tet4.
- The current production guard is `PARAM,SOLIDTYP,NEWSOLID`, linear `CTETRA4`, MAT1/MAT9 material stiffness, and single-material nodal patches. Mixed-material patches are left on the standard Tet4 contribution until a material-grouped patch policy is added.

Validation target from `validate_ctetra4_smooth_target.py`:

| mesh | Python standard T3 x1000 | Python smooth a09 T3 x1000 | MYSTRAN NEWSOLID T3 | MYSTRAN/standard | smooth/standard |
|---|---:|---:|---:|---:|---:|
| n4 | `-2.628385984e-02` | `-1.548880597e-01` | `-1.548880556e-01` | `5.892896` | `5.892896` |
| n8 | `-7.775117808e-02` | `-3.201114590e-01` | `-3.201114333e-01` | `4.117126` | `4.117127` |
| n12 | `-1.233556751e-01` | `-3.950453866e-01` | `-3.950453556e-01` | `3.202490` | `3.202491` |

## CPYRA5 / CPYRA14

Python sources:

- `solid3d_cpyra5.py`
- `solid3d_cpyra5_eas.py`
- `solid3d_cpyra14_liu.py`

This MYSTRAN branch now has CPYRA card/data plumbing and the first `PYRA.f90` kernel. The first-pass goal was to make CPYRA usable without touching legacy hex/penta/tetra kernels.

Implemented Fortran plumbing:

1. Add counters and EDAT sizes in `Source/Modules/SCONTR.f90`:
   - `NCPYRA5`, `NCPYRA14`
   - `MEDAT_CPYRA5 = 7`
   - `MEDAT_CPYRA14 = 16`
2. Add element types and node counts in `Source/Modules/MODEL_STUF.f90`:
   - `PYRA5`, 5 nodes
   - `PYRA14`, 14 nodes
3. Add bulk card readers:
   - `Source/LK1/L1A-BD/BD_CPYRA0.f90`
   - `Source/LK1/L1A-BD/BD_CPYRA.f90`
   - matching interface and `USE_IFs` modules.
4. Wire `CPYRA` into `LOADB0`, `LOADB`, element-property/material mapping, and EMG solid handling.
5. Add EMG support:
   - include `PYRA5/PYRA14` in solid-type checks,
   - add geometry orientation handling,
   - add `PYRA.f90` and `PYRA_Interface.f90`.
6. For `PYRA5`, port the Python 2x2x2 baseline first, then promote the EAS54 stiffness branch.
7. For `PYRA14`, port the Liu composite branch directly from Python, including finite-difference derivatives for this prototype stage.

Current Fortran limits:

- Static stiffness and consistent mass are implemented.
- CPYRA5 stiffness uses the condensed `CPYRA5_EAS54` branch. Mass, thermal, recovery, and KED continue to use the standard CPYRA5 shape-function path.
- Center-point `BE1/BE2` and `SE1/SE2` are implemented for stress/strain helper use.
- Differential stiffness `KED` is implemented with the same `CBAR^T KWW CBAR` pattern used by HEXA/PENTA/TETRA and passes SOL 5 smoke checks.
- Thermal `PTE` and `STE1` are implemented using shape-function temperature interpolation and pass smoke tests.
- `PYRA14` should later replace finite-difference derivatives with analytic derivatives before production hardening.

Current CPYRA check with `validate_cpyra_mystran.py`:

| case | Python T3 x1000 | MYSTRAN T3 | delta |
|---|---:|---:|---:|
| `CPYRA5_EAS54` n4 | `-1.275452610e-01` | `-1.275452667e-01` | `-5.648e-09` |
| `CPYRA14` n2 | `-4.673392309e-02` | `-4.673392080e-02` | `+2.294e-09` |

Current CPYRA SOL 5 smoke check with `validate_cpyra_buckling_smoke.py`:

| case | run | eigen table |
|---|---:|---:|
| `CPYRA5` n2 | ok | true |
| `CPYRA14` n1 | ok | true |

Current CPYRA thermal smoke check with `validate_cpyra_thermal_smoke.py`:

| case | run | fatal |
|---|---:|---:|
| `CPYRA5` | ok | false |
| `CPYRA14` | ok | false |

Current Python reference status:

| element | status | bending note |
|---|---|---|
| `CPYRA5_EAS54_TRIAL` | active Fortran stiffness branch | much better static bending than baseline; still stiff in eigen/buckling |
| `CPYRA14_LIU_COMPOSITE_TRIAL` | validated formula prototype | good convergence, but Python uses finite-difference derivatives |
