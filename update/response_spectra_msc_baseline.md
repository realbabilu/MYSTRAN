# MSC Response-Spectrum Baseline

Date: July 14, 2026

## Scope

Commercial reference runs were executed with:

- `E:\MSC_Nastran\2025.2\bin\nastran.exe`

using normalized copies of the active MYSTRAN RSA decks:

- `D:\18a\femap_RSA\problem_1_020_rsa_srss_msc.dat`
- `D:\18a\femap_RSA\problem_1_024_rsa_srss_msc.dat`

These copies were created only to satisfy MSC bulk-data parsing. The original working
MYSTRAN decks were not overwritten.

## Input-normalization required by MSC

The active MYSTRAN decks did not run directly in MSC because of formatting that MYSTRAN
accepts but MSC rejects:

- `TABDMP1` lines using `.05ENDT` / `.04ENDT`
- `PBAR` trailing fields in `problem_1_020`

After bulk normalization, MSC 2025.2 still rejected:

- `MEFFMASS = ALL`
- `MPFACTOR = ALL`

The working MSC form for these decks is:

- `MEFFMASS = YES`

Also, to force legacy external output generation, these decks should keep:

- `PARAM,POST,-1`
- `PARAM,POSTEXT,YES`

With that combination, MSC produced `.op2` files successfully for both decks. As a
result, the baseline below contains:

- eigenvalues
- modal/eigenvector output
- response-spectrum displacement
- response-spectrum velocity
- response-spectrum acceleration
- response-spectrum bar force output
- printed effective-mass summary

## Problem 1-020 baseline

Source:

- `D:\18a\femap_RSA\problem_1_020_rsa_srss_msc.f06`

Eigenvalues:

- mode 1:
  - eigenvalue = `1.617799E+01`
  - cycles = `6.401511E-01`
- mode 2:
  - eigenvalue = `1.146448E+02`
  - cycles = `1.704109E+00`

Representative response-spectrum displacement:

- `GRID 2`:
  - `T1 = 4.205205E-02`
  - `R2 = 3.712308E-04`
- `GRID 3`:
  - `T1 = 1.068353E-01`
  - `R2 = 3.761756E-04`
- `GRID 7`:
  - `T1 = 4.205222E-02`
  - `R2 = 1.856120E-04`
- `GRID 8`:
  - `T1 = 1.068355E-01`
  - `R2 = 1.880832E-04`

Representative response-spectrum velocity:

- `GRID 2`:
  - `T1 = 1.867535E-01`
  - `R2 = 1.502405E-03`
- `GRID 3`:
  - `T1 = 4.340617E-01`
  - `R2 = 1.842201E-03`
- `GRID 7`:
  - `T1 = 1.867549E-01`
  - `R2 = 7.511893E-04`
- `GRID 8`:
  - `T1 = 4.340625E-01`
  - `R2 = 9.210822E-04`

Representative response-spectrum acceleration:

- `GRID 2`:
  - `T1 = 1.132627E+00`
  - `R2 = 6.300186E-03`
- `GRID 3`:
  - `T1 = 1.865169E+00`
  - `R2 = 1.347277E-02`
- `GRID 7`:
  - `T1 = 1.132646E+00`
  - `R2 = 3.150048E-03`
- `GRID 8`:
  - `T1 = 1.865175E+00`
  - `R2 = 6.736296E-03`

Representative response-spectrum bar forces:

- `CBAR 1`:
  - `M2A = 6.861236E+01`
  - `M2B = 3.345469E+01`
  - `V2 = 8.447395E-01`
  - `AXIAL = 6.906425E-01`
- `CBAR 5`:
  - `M2A = 5.568427E+01`
  - `M2B = 9.965016E-15`
  - `V2 = 4.640356E-01`
  - `AXIAL = 5.869367E-01`

Printed effective-mass summary:

- total effective mass fraction:
  - `T1 = 1.000000E+00`
  - `R2 = 1.000000E+00`
- effective mass matrix key terms:
  - `(T1,T1) = 1.554600E+00`
  - `(T1,R2) = 2.487360E+02`
  - `(R2,R2) = 4.477248E+04`

## Problem 1-024 baseline

Source:

- `D:\18a\femap_RSA\problem_1_024_rsa_srss_msc.f06`

Eigenvalues:

- mode 1:
  - eigenvalue = `7.657227E+02`
  - cycles = `4.404087E+00`
- mode 2:
  - eigenvalue = `8.490395E+02`
  - cycles = `4.637502E+00`
- mode 3:
  - eigenvalue = `7.338652E+03`
  - cycles = `1.363416E+01`
- mode 4:
  - eigenvalue = `7.614361E+03`
  - cycles = `1.388791E+01`

MSC warning:

- modes 3 and 4 use extrapolated damping because `TABDMP1` stops at `10 Hz`
- printed value:
  - damping = `8.000000E-02`

Representative response-spectrum displacement:

- `GRID 10`:
  - `T1 = 2.426231E-04`
  - `T2 = 1.098048E-05`
  - `R2 = 1.478662E-05`
- `GRID 14`:
  - `T1 = 2.529781E-04`
  - `T2 = 1.576009E-05`
  - `R2 = 8.997470E-06`
- `GRID 19`:
  - `T1 = 4.744178E-04`
  - `T2 = 2.002564E-05`
  - `R2 = 8.387294E-06`

Representative response-spectrum velocity:

- `GRID 10`:
  - `T1 = 7.003608E-03`
  - `T2 = 4.210367E-04`
  - `R2 = 4.093171E-04`
- `GRID 14`:
  - `T1 = 7.307261E-03`
  - `T2 = 5.852829E-04`
  - `R2 = 2.497602E-04`
- `GRID 19`:
  - `T1 = 1.367623E-02`
  - `T2 = 8.319996E-04`
  - `R2` should be read from full table if needed for exact deck-to-deck comparison

Representative response-spectrum acceleration:

- `GRID 10`:
  - `T1 = 2.583290E-01`
  - `T2 = 2.692284E-02`
  - `R2 = 1.136507E-02`
- `GRID 14`:
  - `T1 = 2.703952E-01`
  - `T2 = 3.659472E-02`
  - `R2 = 7.116074E-03`
- `GRID 19`:
  - `T1` and companions are available in the same response block and can be scraped in a
    second pass if needed

Representative response-spectrum bar forces:

- `CBAR 1`:
  - `M1A = 2.781558E+00`
  - `M2A = 1.364042E-01`
  - `M1B = 1.803660E+00`
  - `M2B = 1.042826E-01`
  - `V1 = 3.524217E-01`
  - `V2 = 1.845931E-02`
  - `AXIAL = 2.343958E-01`
- `CBAR 14`:
  - `M1A = 2.314137E+00`
  - `M2A = 1.745968E-01`
  - `M1B = 2.606628E+00`
  - `M2B = 1.844232E-01`
  - `V1 = 3.784511E-01`
  - `V2 = 2.760562E-02`

Printed effective-mass summary:

- total effective mass fraction:
  - `T1 = 1.000000E+00`
  - `T2 = 1.000000E+00`
  - `R1 = 1.000000E+00`
  - `R2 = 1.000000E+00`
  - `R3 = 1.000000E+00`
- effective mass matrix key terms:
  - `(T1,T1) = 1.242240E+01`
  - `(T1,R2) = 2.422368E+02`
  - `(T1,R3) = -3.354048E+02`
  - `(T2,R1) = -2.422368E+02`
  - `(R2,R2) = 5.248464E+03`
  - `(R3,R3) = 2.699388E+04`
