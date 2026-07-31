# MYSTRAN 18a July 2026 Release Notes

Branch:

- `v18.00.a`

This release note summarizes the main backports, fixes, and feature additions
carried into the July 2026 `18a` line.

It is intentionally release-facing and focuses on behavior and usage rather
than the full internal patch history.

---

## Summary

| Area | Status |
| :--- | :--- |
| External MUMPS / FEAST / SuperLU / OpenBLAS build support | Done |
| Multi-subcase buckling & STATSUB (`EIGR` / `EIGRL`) | Fixed |
| RFORCE load processing | Fixed |
| K6ROT / element assembly / CELAS1 | Fixed |
| Shell / 1D stress, strain, principal values, curvature output | Fixed |
| Grid-level AUTOSPC stabilization | Fixed |
| Case-control / output descriptor updates | Fixed |
| Response spectrum support (beam path) | Added |
| New solid formulations (`EAS9`, `EAS54`) & restart I/O | Added |
| MITC3+ option for `CTRIA3` | Added |
| `CQUADR` (`DKMQ24`) & `CTRIAR` (`DKMT18`) | Added |
| `PARAM,STR_CID,0` basic/global shell stress output | Added |
| `GPSTRESS`, `GSTRESS`, `STRFIELD`, `SET`, `SURFACE`, `VOLUME` | Added |
| Thin-shell `PSHELL/CQUAD4` blank `MID3` compatibility | Fixed |
| `PBEAMZ` deck support | Added |

Notes:

- targets MinGW builds against a local GCC/OpenBLAS/MUMPS/FEAST stack
- solver-dispatch reporting is clearer at runtime
- this line also includes focused writer work for `F06`, `OP2`, and `NEU`

---

## GPSTRESS / GSTRESS update

### What was added

- `GPSTRESS` and `GSTRESS` are accepted in Case Control
- `STRFIELD` is accepted for the same workflow
- MSC-style `OUTPUT(POST)` shell decks with `SET` / `SURFACE` are now parseable
- the OP2 writer now emits a readable `OGS1` table for the investigated shell
  patch-test workflow

This is an initial compatibility implementation for grid-point surface stress
work. It is separate from ordinary shell element stress.

### Important distinction

- ordinary shell stress in OP2 is the element-stress family (`OES1X1`)
- `GPSTRESS` / `GSTRESS` in OP2 is the grid-point surface-stress family (`OGS1`)

`pyNastran` reads `OGS1` through:

```python
op2.grid_point_surface_stresses
```

### Example usage

```nastran
SET 10 = 100
GPSTRESS = 10
STRFIELD = 10

SUBCASE 1
  SPC = 1
  LOAD = 1
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL

SUBCASE 2
  SPC = 2
  LOAD = 2
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL

OUTPUT(POST)
SET 1 ALL
SURFACE 100 SET 1 NORMAL Z

BEGIN BULK
$PARAM,POST,-1
$PARAM,POSTEXT,YES
```

### Current notes

- in the MSC-style shell patch deck used during this work, `GPSTRESS` and `STRFIELD`
  are tied to the surface request through `SET 10 = 100`
- when `OUTPUT(POST)` is used, `PARAM,POST,-1` and `PARAM,POSTEXT,YES` should remain
  commented out or omitted; enabling both paths together causes a conflict
- for current compatibility, when the deck does not explicitly set
  `PARAM,STR_CID`, the `GPSTRESS/GSTRESS` path defaults the stress-coordinate
  request toward the basic/global style used by the patch-test comparisons
- validation of invariants is safest when derived from `NX/NY/TXY`
  rather than assuming all later raw `OGS1` payload slots have identical
  solver-independent semantics

Useful local references:

- [gpstress_patch_test_example.md](D:/18a/MYSTRAN/update/gpstress_patch_test_example.md)
- [gpstress_recovery_design.md](D:/18a/MYSTRAN/update/gpstress_recovery_design.md)
- [pynastran_gpstress_notes.md](D:/18a/MYSTRAN/update/pynastran_gpstress_notes.md)
- [ogs1_writer_audit_2026-07-19.md](D:/18a/MYSTRAN/update/ogs1_writer_audit_2026-07-19.md)

---

## Thin-shell `PSHELL/CQUAD4` blank `MID3` fix

### What changed

The `PSHELL` parser now restores the classic compatibility fallback:

- if `MID3` is blank, `MID3` inherits `MID2`

This matters for thin-shell style decks that specify membrane and bending
material but intentionally leave the transverse-shear material field blank.

### Why it matters

Without this compatibility path, older thin-shell decks can lose the intended
material carry-over when `MID3` is omitted.

For July 2026 this was especially relevant to shell verification and patch-test
style `CQUAD4` workflows.

### Example usage

```nastran
PSHELL,1,1,0.001,1
```

Meaning in the current `18a` branch:

- `MID1 = 1` supplies membrane material
- `MID2 = 1` supplies bending material
- blank `MID3` now falls back to `MID2` again for compatibility

Field 6 (`12I/TM^3`) still defaults to `1.0` when blank, so ordinary bending
inertia remains active.

### Practical outcome

- old thin-shell decks do not need to be rewritten just to force an explicit
  `MID3`
- this is a compatibility fix, not a new shell formulation

---

## `PARAM,STR_CID,0` shell stress output

MYSTRAN now supports:

```nastran
PARAM,STR_CID,0
```

to request shell stress and strain output in the basic/global coordinate
system.

This is useful for verification decks whose reference values are expressed as
global engineering components, such as SAP2000 patch-test tables.

Supported shell families in this pass:

- `CTRIA3`
- `CQUAD4`
- `CQUADR`
- `CQUAD8`

Default local/material-style shell stress behavior remains unchanged when
`PARAM,STR_CID,-2` is used or implied.

---

## `PBEAMZ` deck support

### Purpose

`PBEAMZ` was added as a beam-property input path built around `PBEAML`-style
section definitions plus a compact wizard-style metadata tail.

It is intended to separate:

- physical section geometry
- stiffness-only modifiers
- optional taper metadata
- optional rigid-offset metadata
- optional station-control metadata

### Recognized metadata

- `STIFFMOD`
- `NSM`
- `TAPER`
- `STATIONS`
- `RIOFFSET` / `ROFSET`

### Minimal non-tapered example

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,END
```

### Modifier example

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,STIFFMOD,1000.,1.0,1.0,0.0,0.0
+,END
```

Interpretation of that example:

- axial-area stiffness modifier = `1000.0`
- major/minor bending modifiers remain `1.0`
- shear modifiers are driven to `0.0`
- omitted trailing modifier fields keep the default `1.0`

### Station-control example

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,STATIONS,10,0.33,0.66
+,END
```

This yields the default base segmentation plus inserted extra breakpoints,
conceptually like:

- `0.0, 0.1, 0.2, 0.3, 0.33, 0.4, 0.5, 0.6, 0.66, 0.7, 0.8, 0.9, 1.0`

### Taper note

A tapered member adds `TAPER` and an end-section block such as `DIM1A`.

Current taper modes:

- `LINEAR` or `1`
- `PARABOLIC` or `2`
- `CUBIC` or `3`

The inertia interpolation follows the SAP-style exponent rule:

```text
I(x) = [ (I1^(1/n)) * (1 - x/L) + (I2^(1/n)) * (x/L) ]^n
```

with:

- `n = 1` linear
- `n = 2` parabolic
- `n = 3` cubic

Torsion `J` is kept linear in the current `PBEAMZ` path.

### Current validation direction

Main local comparison decks:

- `D:\18a\pbeamz_test\prob_001_inclined_frame_pbeaml.dat`
- `D:\18a\pbeamz_test\prob_001_inclined_frame_pbeamz_nomod.dat`
- `D:\18a\pbeamz_test\prob_001_inclined_frame_pbeamz_mod.dat`

Broader notes:

- [pbeamz_phase1.md](D:/18a/MYSTRAN/update/pbeamz_phase1.md)
- [pbeamz_torsion_audit.md](D:/18a/MYSTRAN/update/pbeamz_torsion_audit.md)

---

## Other main July 2026 items

- multi-subcase buckling and `STATSUB` flow backported from the newer path
- RFORCE sign/unit/origin handling corrected
- K6ROT helper introduced and wired into the shell/assembly path
- shell stress/strain/principal/curvature writers updated
- response spectrum support added for the current beam path
- new solid formulation path added through `PARAM,SOLIDTYP,NEWSOLID`
- `MITC3` option added for `CTRIA3`
- `CQUADR` / `CTRIAR` support integrated
- `CELAS1` grounded shorthand compatibility added

---

## Practical outcome

The July 2026 `18a` line is no longer just a small build refresh. It is a
rolling maintenance branch that now includes:

- external solver / BLAS integration work
- solver-dispatch cleanup
- shell output and compatibility fixes
- GPSTRESS / OGS1 groundwork
- beam-property extensions through `PBEAMZ`
- response spectrum beam support
- new solid formulations and restart-state plumbing

For deeper implementation details, use:

- [v18_backport_summary.md](D:/18a/MYSTRAN/update/v18_backport_summary.md)
