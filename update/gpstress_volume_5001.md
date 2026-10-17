> Superseded output coverage: see `op2_solid_5001_audit.md`. Both direct and
> principal OP2 volumes are now implemented, with reference corner-grid rows
> and metadata. Historical limitations/tables below describe earlier builds.

# Solid GPSTRESS VOLUME: attached 5-001, BASIC

**Subsequent quadratic audit:** CPENTA15 and CHEXA20 selective-integration
implementation errors have now been corrected. Both pass the FIXED-input
Python/Fortran/full-solver audit below 1e-7. See
`quadratic_patch_5001_audit.md`; quadratic failures shown below describe the
earlier implementation, before that repair.

## Scope and status

User scope was narrowed to `5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH.dat`, BASIC
coordinates, with Simcenter Nastran **2512** as the reference target. No 5-002
GPSTRESS work was performed. This is an initial, explicitly identified recovery
mode, not a claim of Simcenter numerical equivalence.

The unchanged attached input produces all eight F06/OGS1 volumes, containing
590 grid rows. Five volumes pass the analytical constant-stress comparison;
three fail. The acceptance criterion that all eight pass is **not met**.
Independent Python recovery agrees with the three failing volumes using the
displacements actually solved by MYSTRAN. No stiffness formulation, global
DOF, mesh coordinate, connectivity, material, or load was changed to force a pass.

## Implementation map

| Component | Reused infrastructure and change |
|---|---|
| `LOADC`, `CC_STRE` | Parse POST VOLUME cards and retain separate subcase GP requests. Restore original ordinary stress requests and location for solid-only volume decks. |
| `CC_OUTPUT_DESCRIBERS` | Dynamic volume definitions and POST SET storage; no fixed volume/SET-count ceiling. Retain volume ID, element SET ID and CID. |
| `GPSTRESS_SURFACE_UTILS` | Resolve element ranges, validate membership, collect connected grid IDs, recover and average per volume, write F06 and OGS1. Existing shell routines remain available. |
| `SHP3DT`, `SHP3DP`, `SHP3DH`, `JAC3D`, `B3D_ISOPARAMETRIC` | Reuse existing shape, Jacobian and strain-displacement operations at interior recovery samples. |
| `ELEM_STRE_STRN_ARRAYS` | Optional recovery coordinate for enhanced solids; existing callers keep their previous behavior. Reuse constitutive law, displacement reconstruction and stress transformation. |
| `HEX8_EAS9_STATE`, `PENTA6_ANSEAS_STATE` | Reconstruct condensed enhanced variables at the actual recovery location; no standard-B replacement. |
| OP2 helpers | Reuse `WRITE_TABLE_HEADER`, `WRITE_ITABLE`, `END_OP2_TABLE`, and existing Fortran record conventions. |

## Syntax verification

The installed `C:/Program Files/Siemens/Femap 2512/nastran/bin/nastran.ini`
identifies version 2512.0000. The bundled Quick Reference Guide, VOLUME entry
(PDF page 744, printed page 4-538), defines element SET membership and BASIC,
ELEMENT or CORD output systems; VOLUME follows OUTPUT(POST). The attached
comma-separated `VOLUME id, SET sid, SYSTEM BASIC` form is accepted here.
The parser also accepts `DIRECT`. PRINCIPAL output is not implemented and is
rejected explicitly; this implementation writes the six direct components.

The bundled guide's default requests both DIRECT and PRINCIPAL when neither
is specified. This initial implementation intentionally implements the direct
tensor requirement only; it does not claim complete VOLUME-command parity.

Reference command:

```powershell
& 'C:/Program Files/Siemens/Femap 2512/nastran/bin/nastran.exe' `
  simcenter_2512_attached.dat scr=yes sdir=H:/Projects/MYSTRAN/source/tests/gpstress
```

The first automated attempts reported license error 3060. The user's subsequent
manual run succeeded in `nastran_test/5-001_GPSTRESS_BASIC/nastran`, with no
USER FATAL or SYSTEM FATAL messages and all eight OGS1 volumes present. The
input SHA256 is unchanged. The completed reference supersedes the earlier
license-blocked status.

### Completed Simcenter 2512 comparison

Both the standalone reader and pyNastran decode all eight reference volumes.
Simcenter stores several volumes in a single OGS1 table, each with separate
table-3 metadata; the standalone reader now enumerates all of them.
pyNastran additionally skips PVT0 parameter metadata because its installed
version rejects the `SOLIDTYP` string parameter. OGS results are not skipped.

| Volume | Element | Simcenter grids | MYSTRAN grids | Simcenter analytical max relative error | MYSTRAN vs Simcenter at common grids |
|---|---|---:|---:|---:|---:|
| 1 | CTETRA4 | 47 | 47 | 0 | 0 |
| 2 | CTETRA10 | 47 | 273 | 0 | 0 |
| 3 | CPENTA6 | 16 | 16 | 0.1854621 | 0.3894738 |
| 4 | CPENTA15 | 16 | 56 | 0 | 0.5632124 |
| 5 | CHEXA8 | 16 | 16 | 0 | 0 |
| 6 | CHEXA20 | 16 | 48 | 0 | 0.6505704 |
| 7 | CPYRAM5 | 23 | 23 | 0 | 0 |
| 8 | CPYRAM13 | 23 | 111 | 0 | 0 |

Relative differences use the corresponding analytical component as scale.
Simcenter passes 7/8 analytical volume checks; CPENTA6 differs by up to 18.55%.
CPENTA15 and CHEXA20 pass in Simcenter but fail in MYSTRAN. Their MYSTRAN
recovery agrees with independent Python recovery from the actual MYSTRAN
displacements; that does not establish parity of the underlying formulations.

Simcenter outputs corner grids only for quadratic elements in this deck.
MYSTRAN currently also outputs midsides. Comparisons therefore join GRID IDs,
not array positions. Extra midside rows are a documented output coverage
difference, not exact Simcenter row-layout parity.

`tests/gpstress/compare_simcenter.py` exports `simcenter_2512_comparison.csv`
and `simcenter_2512_summary.csv`. Reference metadata uses approach/device 13,
S-code 1 and output-axis code 3; exact metadata parity has not yet been implemented.

## Recovery mode

Recovery samples use positive-Jacobian interior locations. Tetrahedra use four
symmetric tetrahedral locations; wedges use three triangular locations on each
of two thickness planes; hexahedra use 2x2x2 locations. The existing element
operators evaluate stress there. A linear tetrahedral, wedge, or trilinear
hexahedral tensor field is extrapolated to every connected grid, including
midsides. PYRAM5/13 currently project their existing center tensor to connected
grids. This pyramid behavior preserves constant fields but is not a general
nodal stress extrapolator.

Each contribution is transformed to BASIC before equal-element tensor
averaging. No principal or von Mises values are averaged. Constitutive shear
stresses are used directly; shear strains are not substituted for stress.
Pressure and von Mises are derived from the averaged tensor for OGS1.

Element IDs and grid IDs are sorted, and contributions are accumulated in that
order. Volume identity remains separate even when volumes share the same
coordinates, grid IDs, or elements. Storage scales with model counts and
connected grids rather than maximum external IDs. Invalid Jacobians, missing
IDs/SETs, unsupported families, empty domains and nonfinite stresses produce
diagnostics. Invalid samples are not replaced by zeros.

## OP2 layout and independent decoding

OGS1 real SORT1 volume-direct uses table code **27**, format **1**, width **9**:

| Word | Meaning |
|---|---|
| 1 | `10*GRID + device_code` (integer) |
| 2â€“7 | SXX, SYY, SZZ, SXY, SYZ, SZX (float32) |
| 8 | Pressure, minus one third of the normal-stress trace |
| 9 | von Mises computed from the averaged tensor |

Table-3 word 3 is volume ID; word 4 is subcase ID; word 8 is reference CID.
Table-3 has 50 numeric words and three 128-byte text fields. Existing OP2
helpers frame the records with 32-bit little-endian words on this build.
This layout was checked against the installed pyNastran OGS table-27 decoder.

`tests/gpstress/read_ogs.py` independently checks all Fortran record lengths
and enumerates each volume. pyNastran independently reads the same eight
volumes through `grid_point_stresses_volume_direct`; its result keys include
volume identity. It is not accessed through a generic stress attribute.

The installed pyNastran requires local NumPy compatibility aliases for `in1d`
and `chararray`; no installed package was modified. Existing OES1X1/OSTR1X
tables are skipped in that reader because its pyramid type-255 decoder is
unsupported. OGS1 is decoded, not skipped.

## Attached input results

The fixture has SHA256
`B4EFE70CD735D82BD037C6F93307E4117137ECF084A746D52203CDCDDEC04074`, identical
to the file supplied in Downloads. Analytical BASIC tensor:
`[2000, 2000, 2000, 400, 400, 400] psi`.

| Volume | Element | Grids | Maximum relative component error | Patch |
|---|---|---:|---:|---|
| 1 | CTETRA4 | 47 | 0 at OP2 precision | PASS |
| 2 | CTETRA10 | 273 | 0 at OP2 precision | PASS |
| 3 | CPENTA6 | 16 | 0.5125284 | FAIL |
| 4 | CPENTA15 | 56 | 0.5632124 | FAIL |
| 5 | CHEXA8 | 16 | 0 at OP2 precision | PASS |
| 6 | CHEXA20 | 48 | 0.8311227 | FAIL |
| 7 | CPYRAM5 | 23 | 0 at OP2 precision | PASS |
| 8 | CPYRAM13 | 111 | 0 at OP2 precision | PASS |

CPENTA6/15 and CHEXA20 solved displacements depart from the affine field.
After removing a rigid translation per independent patch, the maximum range
of affine displacement residuals is respectively 2.8930e-4, 1.7397e-4 and
1.5879e-4 inches. Other volumes are within approximately 1.54e-10 inches at
OP2 precision. This is evidence of an existing solution/formulation issue,
not a reason to replace recovered stresses by the analytical target.

Independent Python element kernels evaluate stresses from the actual solved
displacements, then use a separate polynomial linear solve and grid averaging.
Maximum absolute discrepancies are 0.000192 psi (CPENTA6), 0.001243 psi
(CPENTA15), and 0.001618 psi (CHEXA20). These calculations use float32 OP2
displacements, so differentiation amplifies displacement rounding. The audit
tolerance is 1e-5 relative plus 0.002 psi absolute.

The existing ordinary OES/OSTR corner pipeline reports negative corner
Jacobians for CHEXA8 500002 and CHEXA20 600002 (about -1.6e-3), plus associated
EAS recovery messages. The baseline executable already reports these errors.
They remain in F06/ERR; the baseline's final "terminated normally" message is
not sufficient to establish an error-free run. Interior GPSTRESS samples are
valid and all eight GPSTRESS tables remain readable. This does not certify
the mesh as valid throughout its entire natural-coordinate domain.

## Verification and reproduction

```powershell
cmake --build source/build -j 8
Set-Location source/tests/gpstress
H:/Projects/MYSTRAN/.venv/Scripts/python.exe -m unittest test_gpstress -v
H:/Projects/MYSTRAN/source/Binaries/mystran.exe 5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH.dat
H:/Projects/MYSTRAN/.venv/Scripts/python.exe audit_attached.py
```

All **10** regression tests pass. Regular affine regressions pass all eight element types in both LEGACY and
NEWSOLID. Tests also cover F06/OGS consistency, selected domains, subcase
scope, one volume, more than 256 simultaneous volumes, coordinate rotation,
missing elements/ranges/SETs/CIDs, duplicate volume IDs, empty domains and
unsupported families. The regular decks are supplementary regressions, not
substitutes for the attached mesh.

For both the regular fixture and the unchanged attached mesh, enabling
GPSTRESS leaves all existing OUG/OES/OSTR binary payloads byte-identical to
the run without GPSTRESS. No shell, beam, stiffness, or eigenvalue algorithm
was edited. No commits were pushed.

Deliverables are under `tests/gpstress`: the attached fixture, generator,
regression suite, independent OGS reader, Python recovery audit, full grid
tensor CSV, per-component analytical/Python comparisons, and volume summary.
The primary executable is `Binaries/mystran.exe`.

## Remaining limits

All-eight analytical acceptance and Simcenter 2512 numerical parity remain
unachieved. Pyramid recovery remains center projection, and averaging remains
an explicitly identified equal-element initial method. Arbitrary curvilinear
coordinates, nonlinear/modal volume recovery and full PRINCIPAL output are
unsupported. Mixed shell SURFACE/solid VOLUME request side effects have not
been validated; the current accepted scope is the supplied solid-only BASIC
deck. No stiffness changes are proposed as part of this GPSTRESS work.
