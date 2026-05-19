## Composite Buckling Benchmark Stage 1

This note captures the first clean benchmark bridge from the Python reference
package to MYSTRAN3 for laminated shell buckling.

Reference Python source:

- `D:\fortran\mystran3\codex_mod\composite_python_ref\composite_buckling_v5_operator.py`

Element mapping:

- Python `CompDKMQ24` <-> MYSTRAN `CQUADR + PCOMP`
- Python `CompDKMT18` <-> MYSTRAN `CTRIAR + PCOMP`

## What Python is actually solving

The Python operator is not a plain "apply edge force and buckling solve" case.
It builds laminate prebuckling resultants from prescribed generalized strain:

- `N = A eps0 + B kappa`
- `M = B eps0 + D kappa`

and then solves the reciprocal operator form of:

- `K phi = lambda KG phi`

with `KG` assembled from the in-plane resultants `N`.

That means the fairest MYSTRAN parity path is:

1. match the laminate `ABD/T` basis,
2. match the intended prebuckling state,
3. then compare buckling factors.

## First benchmark target

Use the default Python benchmark shape first:

- plate size `a = 1.0`
- total thickness `h = 0.01`
- mesh `n = 4`
- `nmodes = 3`
- preload generalized strain:
  - `eps = [-1e-5, 0, 0]`
- `kap = [0, 0, 0]`

## Degeneracy check before composite buckling

Before comparing laminated buckling, the Python reference now supports both:

- `isotropic`
- `oneply_iso`

where `oneply_iso` is a literal one-ply laminate built through
`Laminate.from_plies(...)`, not the direct `Laminate.isotropic(...)` helper.

Observed results:

- `test_comp_elements.py`
  - `DKMQ24 oneply_iso reproduction relerr = 3.44e-16`
  - `DKMT18 oneply_iso reproduction relerr = 3.06e-16`
- `composite_patch_tests.py`
  - `23 tests, 0 failed`
- `composite_buckling_v5_operator.py`
  - `oneply_iso` gives exactly the same buckling factors as `isotropic`

This means the composite wrappers already degenerate cleanly to the isotropic
DKMQ24/DKMT18 paths at benchmark level.

## Python anchor values

### Symmetric cross-ply

Command:

```text
D:\python312\python.exe D:\fortran\mystran3\codex_mod\composite_python_ref\composite_buckling_v5_operator.py --all --laminate sym_crossply --n 4 --eps=-1e-5,0,0 --kap=0,0,0 --nmodes 3
```

Laminate resultants:

- `Nxx = -7.292349645348e+03`
- `Nyy = -2.816355725100e+02`
- `Nxy = -8.443763395802e-15`
- `Mxx = 0`
- `Myy = 0`
- `Mxy = 0`

Buckling factors:

- `CompDKMQ24`: `18.59959718065`, `67.00449388219`, `72.47469327626`
- `CompDKMT18`: `18.29520954378`, `58.06311182480`, `62.15739066448`

### Isotropic / one-ply isotropic

Commands:

```text
D:\python312\python.exe D:\fortran\mystran3\codex_mod\composite_python_ref\composite_buckling_v5_operator.py --all --laminate isotropic --n 4 --eps=-1e-5,0,0 --kap=0,0,0 --nmodes 3
D:\python312\python.exe D:\fortran\mystran3\codex_mod\composite_python_ref\composite_buckling_v5_operator.py --all --laminate oneply_iso --n 4 --eps=-1e-5,0,0 --kap=0,0,0 --nmodes 3
```

Both produce identical resultants:

- `Nxx = -7.692307692308e+03`
- `Nyy = -2.307692307692e+03`
- `Nxy = 0`
- `Mxx = 0`
- `Myy = 0`
- `Mxy = 0`

Both produce identical buckling factors:

- `CompDKMQ24`: `20.50494044634`, `44.40949932444`, `73.03413019448`
- `CompDKMT18`: `25.21133911545`, `51.91792677963`, `95.23584846675`

### Unsymmetric cross-ply

Command:

```text
D:\python312\python.exe D:\fortran\mystran3\codex_mod\composite_python_ref\composite_buckling_v5_operator.py --all --laminate unsym_crossply --n 4 --eps=-1e-5,0,0 --kap=0,0,0 --nmodes 3
```

Laminate resultants:

- `Nxx = -7.292349645348e+03`
- `Nyy = -2.816355725100e+02`
- `Nxy = -8.443763395802e-15`
- `Mxx = 1.571627078739e+01`
- `Myy = 0`
- `Mxy = -2.110940848950e-17`

Buckling factors:

- `CompDKMQ24`: `9.650901686627`, `22.95096430089`, `46.13741023666`
- `CompDKMT18`: `10.58824015653`, `24.57931337449`, `53.05078124380`

## Why Stage 1 should start with sym_crossply

`sym_crossply` keeps `B = 0`, so the first parity run isolates:

- composite stiffness path,
- shell differential stiffness path,
- buckling extraction path,

without also mixing in bending-membrane coupling from `B`.

## MYSTRAN3 implementation direction

Use the existing shell buckling pattern from:

- `D:\fortran\mystran3\codex_mod\shell_renovation_fix\buckling06_out\buckling06_dkmq24.dat`
- `D:\fortran\mystran3\codex_mod\shell_renovation_fix\buckling06_out\buckling06_dkmt18.dat`

Those files confirm the working MYSTRAN pattern:

- `SOL 5`
- static preload in `SUBCASE 1`
- buckling extraction in `SUBCASE 2`
- `METHOD = 20`
- `EIGRL,20,...`

## Practical caution

Python v5 operator is a prescribed-prestrain benchmark, not a generic load
case. If MYSTRAN stage-1 uses only edge forces, we should expect some mismatch
unless the force pattern reproduces the same membrane resultant field.

So the immediate goal of Stage 1 is not "perfect one-shot equality", but:

1. confirm `CQUADR + PCOMP` and `CTRIAR + PCOMP` both run a clean composite
   buckling case,
2. confirm mode ordering and magnitude family are reasonable against Python,
3. then refine the preload specification if needed.

## Recommended next deck family

- `comp_buckling_cquadr_pcomp_sym_n4.dat`
- `comp_buckling_ctriar_pcomp_sym_n4.dat`

First comparison metrics:

- first 3 buckling load factors,
- whether the first mode family matches between quad and tri paths,
- whether `CQUADR` and `CTRIAR` stay near the Python pair:
  - around `18.6 / 18.3` for mode 1.

## Preload audit update (2026-05-19)

The first MYSTRAN3 `sym_crossply n=4` buckling decks were rerun with
`SPCFORCE/GPFORCE/OLOAD/STRESS/STRAIN` output enabled.

Key result:

- The static preload state is not the source of the remaining mismatch.
- In the `CTRIAR + PCOMP` preload stress table, the ply stresses integrate back
  to the Python target laminate resultants to within normal roundoff:
  - Python target: `Nxx = -7.292349645348e+03`, `Nyy = -2.816355725100e+02`
  - MYSTRAN stress table shows approximately:
    - 0 deg plies: `sigma_1 ~ -1.35789E+06`
    - 90 deg plies x-equivalent: `sigma_2 ~ -1.00584E+05`
    - with `t_ply = 2.5E-03`, the integrated membrane resultants recover the
      same `Nxx/Nyy` family as the Python benchmark.

Therefore the remaining buckling gap for `sym_crossply n=4`:

- Python:
  - `CompDKMQ24`: `18.5996`, `67.0045`, `72.4747`
  - `CompDKMT18`: `18.2952`, `58.0631`, `62.1574`
- MYSTRAN:
  - `CQUADR + PCOMP`: `10.1652`, `36.1459`, `41.8727`
  - `CTRIAR + PCOMP`: `9.99146`, `31.3156`, `35.2067`

is now attributed primarily to the laminated shell buckling / differential
stiffness path (`KG`) rather than to:

- shell routing,
- laminate `ABD/T` property construction,
- isotropic/one-ply composite degeneration,
- or static preload/resultant reproduction.

## Unsym prescribed-state bridge update (2026-05-19)

We also tested a more literal "Python-like" bridge for `unsym_crossply n=8`
using MYSTRAN Case Control:

- `ENFORCED = <file.enf>`

with an `.enf` file that prescribes:

- `u_x = eps_x * x`
- `u_y = 0`
- `u_z = 0`
- `r_x = r_y = r_z = 0`

for every grid, intending to mimic:

- `eps0 = [-1e-5, 0, 0]`
- `kappa0 = [0, 0, 0]`

This route **does pass through the static preload chain**, but it is not yet a
usable parity path for buckling because of a current MYSTRAN limitation:

1. `ENFORCED` in the present implementation places **all DOFs** in the `SE`
   set, not only the in-plane DOFs we would want for this benchmark.
2. That collapses the remaining buckling problem to an ultra-small reduced
   system (`KLL/KLLD` effectively `1x1` in the test run).
3. LINK4 then crashes in the progress formatter while assembling the tiny
   buckling operator:
   - runtime failure in `COUNTER_PROGRESS.f90`
   - seen during `PUT RFAC MATRIX IN ARPACK BAND FORM`

Practical conclusion:

- `ENFORCED` is useful as an experiment because it confirms MYSTRAN can ingest
  a prescribed displacement field,
- but it is **not yet the right production bridge** for the unsymmetric
  composite buckling benchmark,
- because we need a way to prescribe the membrane state without globally
  enslaving all rotational/out-of-plane DOFs.

So the clean next-step for `unsym` is now better defined:

1. add a **partial enforced-state bridge** (in-plane only), or
2. add a **direct benchmark preload injection path** for shell buckling that
   feeds the intended `(eps0, kappa0)` state into the differential stiffness
   construction without using a full static edge-load solve.

## Unsym direct state injection update (2026-05-19)

We then implemented the second option narrowly for the composite shell kernels
only:

- `CQUADR_DKMQ24.f90`
- `CTRIAR_DKMT18.f90`

using:

- `DEBUG,239,1`

as a benchmark-only toggle in the composite buckling `KED` path.

When `PCOMP_PROPS = 'Y'` and `DEBUG(239) > 0`, the shell buckling differential
stiffness now injects the Python-v5 prescribed laminate state directly:

- `eps0   = [-1.0e-5, 0, 0]`
- `kappa0 = [ 0,      0, 0]`

without relying on the static edge-load solution for `STRAIN(1:6)`.

This was rerun for `unsym_crossply n=8`:

- Python anchors:
  - `CompDKMQ24`: `9.4040`, `18.7418`, `39.0515`
  - `CompDKMT18`: `9.7572`, `19.4696`, `41.6094`
- MYSTRAN with direct state injection:
  - `CQUADR + PCOMP`: `16.2808`, `33.7357`, `66.4137`
  - `CTRIAR + PCOMP`: `16.4446`, `33.7769`, `67.0471`

Interpretation:

- this removes the old ambiguity about whether the `unsym` gap is caused only
  by the static edge-load preload bridge,
- because even with the Python-like prescribed laminate state injected
  directly into `KED`, the MYSTRAN `unsym` benchmark remains much higher than
  Python,
- while the `sym_crossply n=8` case still matches closely.

So the remaining `unsym` gap is now better localized:

- not shell routing,
- not `PCOMP` parsing,
- not static composite stiffness parity,
- not the crude edge-load state bridge alone,
- but most likely the **laminated shell buckling operator under unsymmetric
  `B`-coupled behavior** for `CQUADR/CTRIAR`.

For clarity, the benchmark decks are now split into:

- ordinary edge-load unsym benchmark:
  - `comp_buckling_cquadr_pcomp_unsym_n8.dat`
  - `comp_buckling_ctriar_pcomp_unsym_n8.dat`
- Python-v5 prescribed-state benchmark:
  - `comp_buckling_cquadr_pcomp_unsym_n8_pyv5.dat`
  - `comp_buckling_ctriar_pcomp_unsym_n8_pyv5.dat`

## Static K cross-check update (2026-05-19)

Before going deeper into the unsymmetric buckling operator, we checked whether
the composite static stiffness `K` itself was already healthy in MYSTRAN for
`CQUADR` and `CTRIAR`.

Added benchmark helpers:

- `calc_field_patch_energy.py`
- `composite_python_ref/make_mystran_composite_field_patch_decks.py`
- `composite_python_ref/probe_composite_buckling_model.py`

The field-patch generator now writes three exact enforced-field families:

- `membrane`
- `bending`
- `coupled`

Python reference energies from `composite_field_patch_tests.py` for `n=4`:

- `sym_crossply, membrane`: `4.432525742087e+00`
- `unsym_crossply, membrane`: `4.432525742087e+00`
- `sym_crossply, bending`: `2.066868735223e-02`
- `unsym_crossply, bending`: `1.477508580696e-02`
- `sym_crossply, coupled`: `4.453194429440e+00`
- `unsym_crossply, coupled`: `4.211556766084e+00`

MYSTRAN energies reconstructed from `0.5 u^T f` using the enforced field and
printed `SPC FORCE` reactions:

### CQUADR

- `sym_crossply, coupled`: `4.453194488075e+00`
- `unsym_crossply, coupled`: `4.211557338275e+00`
- `unsym_crossply, membrane`: `4.432525800000e+00`
- `unsym_crossply, bending`: `1.477508565000e-02`

### CTRIAR

- `sym_crossply, coupled`: `4.453194488105e+00`
- `unsym_crossply, coupled`: `4.211557339030e+00`
- `unsym_crossply, membrane`: `4.432525800000e+00`
- `unsym_crossply, bending`: `1.477508550510e-02`

Small-model Python probe (`n=2`, prescribed `eps0=[-1e-5,0,0]`, `kappa0=0`)
showed:

- `KG_norm` is identical for `sym_crossply` and `unsym_crossply`
- `K_norm` changes only slightly between `sym` and `unsym`
- the `n=2` buckling eigenvalue is identical in Python for `sym` and `unsym`
  for both `CompDKMQ24` and `CompDKMT18`

Takeaway:

- static composite `K` for `CQUADR` and `CTRIAR` is already matching Python to
  very tight tolerance for membrane, bending, and coupled fields,
- therefore the remaining `unsym_crossply` discrepancy is not a general
  `PCOMP` or static stiffness problem,
- the unresolved part is now isolated more tightly to the **buckling operator
  path** for unsymmetric `B`-coupled laminates.
