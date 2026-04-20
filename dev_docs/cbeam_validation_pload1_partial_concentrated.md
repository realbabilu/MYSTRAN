# CBEAM Validation: PLOAD1 Partial-Span And Concentrated-In-Element

Validation input:

- [cbeam_pload1_partial_concentrated_validation.dat](E:/mystran17/mystran/Binaries/cbeam_pload1_partial_concentrated_validation.dat)

Run output:

- [cbeam_pload1_partial_concentrated_validation.F06](E:/mystran17/mystran/Binaries/cbeam_pload1_partial_concentrated_validation.F06)

## Model

Single `CBEAM` cantilever (`L=1.0`) with local beam `+y` aligned to global `+y`.

## Subcase 1: partial-span uniform local-y load

Input:

- `TYPE=FYE`, `P1=P2=1.0`
- `X1=0.2`, `X2=0.7`

Result checks from `F06`:

- fixed-end reaction shear `Fy = -5.000000E-01` (matches total applied `q*(X2-X1)*L = 0.5`)
- fixed-end reaction moment `Mz = -2.250000E-01`
- element engineering forces at end A: shear `+5.000000E-01`, moment `+2.250000E-01`

This confirms partial-span distributed load is converted and recovered consistently.

## Subcase 2: concentrated local-y load in element

Input:

- `TYPE=FYE`, `P1=P2=1.0`
- `X1=X2=0.5` (point load at midspan)

Result checks from `F06`:

- fixed-end reaction shear `Fy = -1.000000E+00`
- fixed-end reaction moment `Mz = -5.000000E-01`
- element engineering forces at end A: shear `+1.000000E+00`, moment `+5.000000E-01`
- element engineering moment at end B is approximately zero (`-5.551115E-17`, roundoff)

This confirms the concentrated-in-element path (`X1=X2`) is working and matches cantilever statics.

## Conclusion

The `PLOAD1` implementation for `CBAR/BART/CBEAM` now supports:

- full-span uniform loads
- partial-span distributed loads
- concentrated-in-element loads (`X1=X2`)

for local component types currently supported (`FYE/FZE/FXE/MXE` plus `FY/FZ/Y/Z` aliases).
