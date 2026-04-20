# CBEAM PLOAD1 Validation

Validation deck:
- [cbeam_pload1_validation.dat](E:/mystran17/mystran/Binaries/cbeam_pload1_validation.dat)

Run output:
- [cbeam_pload1_validation.F06](E:/mystran17/mystran/Binaries/cbeam_pload1_validation.F06)

## What this checks

This case validates the new limited `PLOAD1` beam-load path:

- `PLOAD1 -> LINK1Q`
- `PRESSURE_DATA_PROC -> PPNT/PDATA/PTYPE`
- `ELMDAT2 -> PRESS`
- `BAR1/BEAM -> PPE`
- `EPTL -> SYS_LOAD -> PG`

The model is a one-element cantilever `CBEAM` with:

- local `y = global y`
- local `z = global z`
- `q_y = +1.5`
- `q_z = -2.0`
- `L = 1.0`

## Expected result

For uniform full-span distributed loads on a cantilever:

- total shear in local `y`: `q_y L = 1.5`
- total shear in local `z`: `q_z L = -2.0`
- fixed-end moment from `q_y`: `q_y L^2 / 2 = 0.75`
- fixed-end moment from `q_z`: `q_z L^2 / 2 = -1.0`

In global coordinates for this orientation:

- reaction `T2 = -1.5`
- reaction `T3 = +2.0`
- reaction `R3 = -0.75`
- reaction `R2 = -1.0`

## MYSTRAN result

From the SPC force table:

- `T2 = -1.500000E+00`
- `T3 =  2.000000E+00`
- `R2 = -1.000000E+00`
- `R3 = -7.500000E-01`

These match the expected fixed-end reactions exactly within printed precision.

From the displacement table at the free end:

- `T2 =  1.526400E-04`
- `T3 = -2.035200E-04`
- `R2 =  2.560000E-04`
- `R3 =  1.920000E-04`

This confirms the distributed beam load is entering the solve and producing bending in both transverse directions.

## Engineering-force recovery

After the recovery-side fix in LINK9, the `ELEMENT ENGINEERING FORCES` table now reports:

- end A plane 1 moment: `+7.500000E-01`
- end A plane 2 moment: `-1.000000E+00`
- end B plane 1 moment: `-1.110223E-16`
- end B plane 2 moment: `-2.220446E-16`
- plane 1 shear: `+1.500000E+00`
- plane 2 shear: `-2.000000E+00`

These values are consistent with the cantilever fixed-end actions from the distributed load:

- root bending moments match `qL^2/2`
- free-end moments are zero within roundoff
- root shears match `qL`

So the current status is:

- global load assembly: validated
- reactions: validated
- displacement response: validated
- element engineering-force recovery for `PLOAD1`: validated for the current uniform full-span `FY/FZ` implementation
