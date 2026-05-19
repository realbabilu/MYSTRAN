## Composite Python vs MYSTRAN Plan

Reference package copied to:

- `D:\fortran\mystran3\codex_mod\composite_python_ref`

Primary Python references:

- `composite_patch_tests.py`
- `composite_field_patch_tests.py`
- `composite_plate_bending.py`
- `composite_buckling_v5_operator.py`

Element mapping:

- Python `CompDKMQ24` <-> MYSTRAN `CQUADR + PCOMP`
- Python `CompDKMT18` <-> MYSTRAN `CTRIAR + PCOMP`

Current MYSTRAN routing status:

- `CQUADR` has explicit isotropic/composite path labeling in F06/BUG
- `CTRIAR` has explicit isotropic/composite path labeling in F06/BUG
- legacy `PCOMP` input behavior is preserved

## Benchmark ladder

### Stage 0: Composite element sanity

Python:

- `test_comp_elements.py`
- `composite_patch_tests.py`

Purpose:

- isotropic laminate reproduces isotropic base element
- symmetric laminate gives `B ~= 0`
- unsymmetric laminate gives `B != 0`
- stiffness remains symmetric
- laminate `ABD/As` algebra is healthy

Status:

- copied package runs clean in `mystran3`
- `composite_patch_tests.py`: `19 tests, 0 failed`

Interpretation:

- this is the right first gate before any MYSTRAN parity claim

### Stage 1: Field patch parity

Python:

- `composite_field_patch_tests.py`

Purpose:

- assembled mesh
- prescribed exact displacement field
- compare `Ufe = 0.5 u^T K u`
- compare `Uref = 0.5 area [eps0,kappa]^T ABD [eps0,kappa]`

Cases:

- `membrane`
- `bending`
- `coupled`

Recommended first MYSTRAN parity target:

- `coupled`
- laminate = `sym_crossply`
- then laminate = `unsym_crossply`

Why:

- this is the cleanest bridge from laminate theory to assembled shell behavior
- it checks the `B` coupling branch explicitly

MYSTRAN note:

- MYSTRAN does not naturally expose `0.5 u^T K u` from a prescribed field in the same way as the Python harness
- for production parity, use this stage first as a qualitative/structural anchor and then decide whether to:
  - add a small debug energy print path, or
  - compare resulting nodal force reactions under prescribed displacements

### Stage 2: Structural static bending

Python:

- `composite_plate_bending.py`

Problem:

- square plate
- simply supported edges with `w = 0`
- minimal in-plane anchors
- uniform transverse pressure

Recommended first comparison:

- laminate = `sym_crossply`
- `n = 4, 8`
- compare:
  - center displacement `center_w`
  - max transverse displacement `max_abs_w`
  - corner in-plane drift `corner_ux`, `corner_uy`

Current Python anchors for `sym_crossply` with `pressure = 1000.0`, `a = 1.0`, `h = 0.01`:

- `CompDKMQ24`
  - `n=4`: `center_w = -1.202206841986e-03`
  - `n=8`: `center_w = -1.159953254336e-03`
- `CompDKMT18`
  - `n=4`: `center_w = -1.102182097757e-03`
  - `n=8`: `center_w = -1.127396064442e-03`

Why this stage matters:

- easier to mirror in MYSTRAN than field-energy patch
- gives a structural benchmark before buckling

### Stage 3: Composite buckling operator

Python:

- `composite_buckling_v5_operator.py`

Problem:

- benchmark-level buckling operator
- laminate prebuckling resultants from:
  - `N = A eps0 + B kappa`
  - `M = B eps0 + D kappa`
- reciprocal operator route:
  - `K phi = lambda KG phi`

Recommended first comparison:

- laminate = `sym_crossply`
- `eps = [-1e-5, 0, 0]`
- `kap = [0, 0, 0]`
- `n = 4, 8`
- `nmodes = 3`

Current Python anchors:

- `CompDKMQ24`
  - `n=4`: `18.5996, 67.0045, 72.4747`
  - `n=8`: `18.4700, 59.2993, 63.8817`
- `CompDKMT18`
  - `n=4`: `18.2952, 58.0631, 62.1574`
  - `n=8`: `18.5403, 57.4626, 62.5623`

Why this is not first:

- harder to mirror in MYSTRAN
- better after static composite path is already trusted

## Recommended MYSTRAN implementation order

1. Confirm composite shell routing with tiny smoke decks
   - already done for `CQUADR + PCOMP`
   - already done for `CTRIAR + PCOMP`

2. Build first MYSTRAN structural comparison from Stage 2
   - composite plate bending
   - easiest practical parity case

3. Add Stage 1 style field patch support only if needed
   - especially if we need to isolate `B` coupling numerically

4. Move to Stage 3
   - composite buckling operator style benchmark

## First MYSTRAN parity deck set

Suggested first deck family:

- `comp_plate_bending_cquadr_pcomp_n4.dat`
- `comp_plate_bending_cquadr_pcomp_n8.dat`
- `comp_plate_bending_ctriar_pcomp_n4.dat`
- `comp_plate_bending_ctriar_pcomp_n8.dat`

Shared settings:

- laminate = `sym_crossply`
- geometry = square plate `a = 1.0`
- thickness = `0.01`
- pressure = `1000.0`
- simply supported `w = 0` on all edges
- minimal in-plane anchors matching Python intent

Metrics to compare:

- center `UZ`
- max `|UZ|`
- corner `UX`
- corner `UY`

## First MYSTRAN run status

Generated decks:

- `D:\fortran\mystran3\codex_mod\composite_plate_bending_mystran\comp_plate_bending_cquadr_pcomp_n4.dat`
- `D:\fortran\mystran3\codex_mod\composite_plate_bending_mystran\comp_plate_bending_cquadr_pcomp_n8.dat`
- `D:\fortran\mystran3\codex_mod\composite_plate_bending_mystran\comp_plate_bending_ctriar_pcomp_n4.dat`
- `D:\fortran\mystran3\codex_mod\composite_plate_bending_mystran\comp_plate_bending_ctriar_pcomp_n8.dat`

Load implementation note:

- switched from `PLOAD2` to equivalent nodal `FORCE`
- this matches the current Python benchmark intent more closely for both quad and tri paths

Observed F06 path confirmation:

- `CQUADR`: `DKMQ24 composite shell path`
- `CTRIAR`: `DKMT18 composite shell path`

Observed center `UZ` from MYSTRAN:

- `CQUADR + PCOMP`
  - `n=4`: node `13` -> `UZ = -2.737718E-03`
  - `n=8`: node `41` -> `UZ = -2.613753E-03`
- `CTRIAR + PCOMP`
  - `n=4`: node `13` -> `UZ = -6.381260E-04`
  - `n=8`: node `41` -> `UZ = -6.601814E-04`

Python anchors for same benchmark:

- `CompDKMQ24`
  - `n=4`: `center_w = -1.202206841986e-03`
  - `n=8`: `center_w = -1.159953254336e-03`
- `CompDKMT18`
  - `n=4`: `center_w = -1.102182097757e-03`
  - `n=8`: `center_w = -1.127396064442e-03`

Immediate reading:

- composite routing works
- benchmark deck family is alive
- parity is not yet there
- `CQUADR + PCOMP` is currently softer than Python by about `2.3x`
- `CTRIAR + PCOMP` is currently stiffer than Python by about `0.58x` displacement ratio

This is a useful first parity case because it shows:

- the laminated composite paths are active
- the mismatch is now in shell implementation behavior, not in routing or parser selection

## Practical conclusion

Use this ladder:

1. `composite_patch_tests.py`
2. `composite_field_patch_tests.py`
3. `composite_plate_bending.py`
4. `composite_buckling_v5_operator.py`

For MYSTRAN parity, start at Stage 2 first. It is the best balance of:

- structural meaning
- implementation effort
- diagnostic value
