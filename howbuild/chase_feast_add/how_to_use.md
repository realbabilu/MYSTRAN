# CHASE/FEAST - How To Use

This note is for runtime usage after build integration.

## 1) Pick one example deck
Pass-focused set (Eigen8 only):
- `examples/midas_eigen08_cbar_axis_debug.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbar_axis_debug_feast.dat`
- `examples/midas_eigen08_cbar_axis_debug_chase.dat`
- `examples/midas_eigen08_cbeam_axis_debug.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbeam_axis_debug_feast.dat`
- `examples/midas_eigen08_cbeam_axis_debug_chase.dat`
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli.dat` (ARPACK baseline)
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli_feast.dat`
- `examples/midas_eigen08_cbeam_axis_debug_bernoulli_chase.dat`

## 2) How to switch method in deck
Use `PARAM,LANCMETH,...`:
- `PARAM,LANCMETH,ARPACK` (default/stable path)
- `PARAM,LANCMETH,FEAST`
- `PARAM,LANCMETH,CHASE`

Recommended solver pairing for modal runs:
- `PARAM,SOLLIB,SPARSE`
- `PARAM,SPARSEFLAVOR,SUPERLU`

## 3) Run command
Example:

```powershell
E:\mystran17\mystran\Binaries\mystran.exe E:\mystran17\mystran\howbuild\chase_feast_add\examples\midas_eigen08_cbeam_axis_debug_feast.dat
```

## 4) What to check
- `*.ERR`: warnings/fallback info
- `*.F06`: eigenvalue table / frequencies
- Compare FEAST/CHASE results against ARPACK baseline deck with same model.

## 5) Current behavior note
`LANCMETH=FEAST/CHASE` is integrated as experimental path with guarded behavior in LINK4 flow.  
This package intentionally uses Eigen8 decks because they are the most stable/pass-focused validation set for quick checks.
