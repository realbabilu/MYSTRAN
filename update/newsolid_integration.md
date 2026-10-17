# NEWSOLID Integration

Date: 2026-10-15

This branch includes the revised Bedrosian CPYRAM13 and frame-invariant EAS9 v2 stiffness paths. See `solid_revision_results_2026-10-15.md` for executed results and independent feature limitations.

## New Parameter

Use:

```nastran
PARAM,SOLIDTYP,NEWSOLID
```

Accepted values:

- `LEGACY` : keep the original branch behavior
- `NEWSOLID` : enable the new solid branches
- `EAS` : temporary compatibility alias, routed to `NEWSOLID`

If `PARAM,SOLIDTYP` is omitted, the default remains `LEGACY`.

## Bulk Data Support

`CPYRAM` is now accepted with standard Nastran spelling.

Compatibility alias kept for now:

- `CPYRA`

## NEWSOLID Coverage

With `PARAM,SOLIDTYP,NEWSOLID`, the branch now enables:

- `CHEXA8` : EAS9 v2 natural-to-Cartesian strain mapping and local condensation
- `CPENTA6` : EAS9 v2 natural-to-Cartesian strain mapping and local condensation
- `CTETRA4` : smoothed assembly correction path
- `CPYRAM5` : standard linear pyramid in both LEGACY and NEWSOLID
- `CPYRAM13` : Bedrosian 13-node, eight-point Jacobi stiffness in both modes

CPYRA14/Liu and pyramid EAS54 were removed. Standard CPYRAM cards accept exactly 5 or 13 nodes. No global enrichment DOFs are introduced.

EAS9 v2 stiffness passes the supplied 5-001 patch, rotation of geometry and loads, and 5-004 Python reproduction. General enhanced thermal RHS and corner recovery remain KNOWN_FAIL; pressure, Kg and dynamic/output parity have independent acceptance gates. PYRAM13 OP2 stress/strain is suppressed with a warning until center-plus-corner recovery is implemented. For modal PYRAM13 mass use `PARAM,COUPMASS,1`; indefinite row-sum modal mass is rejected.

The restart persistence path was also updated, so `SOLIDTYP` survives restart I/O through `L1A`.

## Build Status

Integrated files compile successfully in:

- `H:\Projects\MYSTRAN\source\build`
- output binary: `H:\Projects\MYSTRAN\source\Binaries\mystran.exe`
- tools: `C:\gcc\bin\cmake.exe`, `C:\gcc\bin\make.exe`, `C:\gcc\bin\gfortran.exe`

## Notes

- `LEGACY` remains the safe default for existing decks.
- The intent of this integration is side-by-side benchmarking, not silent behavior replacement.
- Benchmarking should compare the same deck under both:

```nastran
PARAM,SOLIDTYP,LEGACY
```

and

```nastran
PARAM,SOLIDTYP,NEWSOLID
```


### CPYRAM13 v3 thermal/buckling update (2026-10-15)

CPYRAM13 K remains 8-point in both SOLIDTYP modes. Consistent mass uses independent Jacobi27. Physical initial-stress Kg now uses stresses recovered from the displacement field at each integration point, with thermal strain subtracted. PSOLID IN=2/default selects Kg8; IN=3 selects Kg27 while K remains K8 and mass27. Thermal 5-011 passes analytic, GP and full SOL 101 output checks. The 5-010 nominal-uniform-prestress screen matches Python, but full loaded SOL 105/SAP/Nastran certification is pending. See `cpyram13_v3_results_2026-10-15.md` for all gates and retained risks.


### CHEXA8 EAS9 v2 feature update (2026-10-15)

NEWSOLID CHEXA8 uses v2 covariant modes and element-local Schur condensation. Condensed thermal RHS and full enhanced center/corner recovery now include the local thermal parameter contribution. Kg uses enhanced prestress at each Gauss point and the physical initial-stress gradient operator. SOL 101/103/105 manufactured checks pass independently; no nonlinear finite-strain EAS tangent is claimed. PSOLID IN=2 is required; incompatible order now errors instead of silently changing formulation. Solid PLOAD4 G1/G3 and stress-output OP2 readback remain KNOWN_FAIL, with matched Nastran feature parity pending. See `hexa8_eas9_v2_results_2026-10-15.md`.


### Confirmed user selection (2026-10-16)

NEWSOLID CPENTA6 remains EAS9 v2; CTETRA4 remains nodal-smoothed alpha=0.9. The elastic smoothing correction is now explicitly excluded from buckling Kg assembly. Full static CTETRA4 agreement with v4 Python passes, while CTETRA4/CPENTA6/CPYRAM5 buckling remains NOT FINAL PASS. V4/v5 execute standard-element screens and do not establish EAS9 v2 or smoothed buckling parity. See `newsolid_selection_2026-10-16.md`.


### Explicit experimental CPENTA6 SOL 105 selector (2026-10-16)

PENTABUC defaults to EAS9. With NEWSOLID and PENTABUC=ANSEAS, both phases of SOL 105 use Claude-original ANS/EAS K and recovery, and Kg uses its recovered GP prestress with the physical displacement-gradient operator. Standalone static and other solutions remain EAS9 v2. LEGACY ignores this selector. Distorted affine patch failure is reproduced and warned; this is NOT FINAL PASS or certified nonlinear tangent/SAP parity. Details: `cpenta6_buckling_anseas_2026-10-16.md`.
