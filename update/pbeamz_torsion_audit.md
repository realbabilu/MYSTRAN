## PBEAMZ torsion audit for SAP2000 Example 1-006a

Date: July 16, 2026

### Scope

Compared:

- `D:\18a\pbeamz_test\prob_006a_nonprismatic_pbeamz_clean.dat`
- `D:\18a\pbeamz_test\prob_006a_nonprismatic_pbeam_original.dat`

Main check:

- free-end response at grid 5 for all seven subcases
- special focus on subcase 5 torsion (`Rx`)

### Root cause found

`PBEAMZ` properties that had already been expanded into explicit stations were still being routed through the special 2-end taper stiffness path in `BEAM.f90`.

That meant:

- explicit multi-station `PBEAM` used the generic stationed beam stiffness path
- expanded `PBEAMZ` used a different stiffness path

This was the source of the larger torsion mismatch.

### Fix applied

In `D:\18a\MYSTRAN\Source\EMG\EMG3\BEAM.f90`:

- keep `BUILD_PBEAMZ_TAPERED_BEAM_KE` only for abstract 2-station `PBEAMZ`
- send expanded multi-station `PBEAMZ` through `BUILD_TAPERED_BEAM_KE`, same as explicit `PBEAM`

Condition changed from:

- `CBEAM_ACTIVE_TAPER_MODE > 0 .AND. NSTA >= 2`

to:

- `CBEAM_ACTIVE_TAPER_MODE > 0 .AND. NSTA == 2`

### Result after fix

Free-end response at grid 5:

- SC1: exact match
- SC2: exact match
- SC3: exact match
- SC4: exact match
- SC6: exact match
- SC7: exact match

Torsion subcase SC5:

- explicit `PBEAM`: `7.981836E-02`
- `PBEAMZ clean`: `7.984913E-02`

Residual difference:

- about `3.08E-05`
- about `0.0386%`

### Interpretation

The major torsion-path defect is fixed.

The small remaining SC5 difference is no longer a stiffness-path mismatch. It is most likely from section-property generation differences between:

- explicit hand-entered `PBEAM` station `J` values in the reference deck
- `PBEAMZ` section-property generation from the PBEAML shape database

### Next audit if needed

If torsion must be matched even tighter, compare the generated station `J` values from `PBEAMZ` against the explicit `PBEAM` `J` values for each segment in Example 1-006a.

### Follow-up audit

Additional tracing after the main fix showed:

- generated `PBEAMZ` station `J` values match the explicit `PBEAM` station list for Example 1-006a to the rounding shown in the deck
- the runtime loader in `Source/EMG/EMG1/ELMDAT1.f90` copies:
  - `PBEAM_NSTATIONS`
  - `PBEAM_XL`
  - `PBEAM_RPROPS`
  directly into:
  - `CBEAM_ACTIVE_NSTATIONS`
  - `CBEAM_ACTIVE_XL`
  - `CBEAM_ACTIVE_RPROPS`

That means the remaining small SC5 torsion delta is not explained by:

- the old special `PBEAMZ` stiffness path
- the `BAR` torsion-constant generator itself
- the property-to-active-property loader

Current conclusion:

- the remaining SC5 difference is small enough to classify as a residual numerical mismatch
- it is no longer an obvious logic defect in `PBEAMZ`
