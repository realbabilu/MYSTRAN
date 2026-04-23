# CHASE/FEAST - How To Use

This note is for runtime usage after build integration.

## 1) Pick one example deck
From:
- `examples/midas_eigen09_pyramid_modal_axis_debug_nev6.dat` (ARPACK baseline)
- `examples/midas_eigen09_pyramid_modal_axis_debug_nev6_feast.dat`
- `examples/midas_eigen09_pyramid_modal_axis_debug_nev6_chase.dat`
- `examples/midas_eigen09_pyramid_modal_cbeam_axis_debug_nev6.dat` (ARPACK baseline)
- `examples/midas_eigen09_pyramid_modal_cbeam_axis_debug_nev6_feast.dat`
- `examples/midas_eigen09_pyramid_modal_cbeam_axis_debug_nev6_chase.dat`
- `examples/feast_modal_chain_fullrank.dat` (small surrogate check)

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
E:\mystran17\mystran\Binaries\mystran.exe E:\mystran17\mystran\howbuild\chase_feast_add\examples\midas_eigen09_pyramid_modal_axis_debug_nev6_feast.dat
```

## 4) What to check
- `*.ERR`: warnings/fallback info
- `*.F06`: eigenvalue table / frequencies
- Compare FEAST/CHASE results against ARPACK baseline deck with same model.

## 5) Current behavior note
`LANCMETH=FEAST/CHASE` is integrated as experimental path with guarded behavior in LINK4 flow.  
Keep ARPACK as production default unless validation deck set is explicitly passing for your case.
