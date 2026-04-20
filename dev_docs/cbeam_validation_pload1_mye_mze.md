# CBEAM Validation: PLOAD1 MYE/MZE

Validation input:

- [cbeam_pload1_mye_mze_validation.dat](E:/mystran17/mystran/Binaries/cbeam_pload1_mye_mze_validation.dat)

Run output:

- [cbeam_pload1_mye_mze_validation.F06](E:/mystran17/mystran/Binaries/cbeam_pload1_mye_mze_validation.F06)

## What was validated

Single-element cantilever `CBEAM` (`L=1`) with two subcases:

- Subcase 1: `PLOAD1,MYE,FR,0,1,1`
- Subcase 2: `PLOAD1,MZE,FR,0,1,1`

## Results from F06

Subcase 1 (`MYE`):

- root reaction `R2 = -1.000000E+00`
- element engineering force at end A, plane 2 moment = `-1.000000E+00`

Subcase 2 (`MZE`):

- root reaction `R3 = -1.000000E+00`
- element engineering force at end A, plane 1 moment = `+1.000000E+00`

Both runs terminate normally and pass parser -> load processing -> assembly -> solve -> recovery.
