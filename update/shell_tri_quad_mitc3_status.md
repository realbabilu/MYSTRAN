# Shell Triangle/Quadrilateral Element Status

## Scope

This note records the current `mystran18a` shell-element output status after the `CTRIAR`, `CQUADR`, `MITC3+`, global stress-coordinate, and baseline `OGS1` work.

The purpose is to avoid mixing three different topics:

- Bulk Data card compatibility.
- Element formulation used for stiffness/recovery.
- F06/OP2/NEU output support.

## Current element mapping

| Bulk Data card | Internal `TYPE` | Current formulation path | Current output path |
| --- | --- | --- | --- |
| `CTRIA3` | `TRIA3` | MIN3/MITC3+ path, depending branch settings | TRIA3 stress/strain/force writers |
| `CTRIAR` | `TRIA3` with DKMT sentinel | `CTRIAR_DKMT18` for DKMT-style shell stiffness | Reuses TRIA3 stress/strain/force writers |
| `CQUAD4` | `QUAD4` | MIN4/MITC4+ path, depending branch settings | QUAD4 stress/strain/force writers |
| `CQUADR` | `QUADR` | `CQUADR_DKMQ24` | QUAD4/CQUADR stress/strain/force writers |
| `CQUAD8` | `QUAD8` | MITC8 path | QUAD8 shell stress/strain/force writers |
| `CTRIA6` | not present | not implemented | not implemented |

## `CTRIAR`

`CTRIAR` is parsed by the existing `BD_CTRIA` path and stored internally as `TYPE='TRIA3'`.

The distinguishing marker is stored in `EDAT` as a DKMT18 sentinel:

```text
EDAT(... thickness key / CTRIAR sentinel) = -18
```

This lets the element-generation side dispatch to `CTRIAR_DKMT18` while still reusing the existing triangular output infrastructure.

Practical consequence:

- `CTRIAR` benefits from the same F06/OP2/NEU triangular stress/strain/force output path as `CTRIA3`.
- The new baseline `OGS1` writer includes `CTRIAR` because it checks internal `TYPE(1:5) == 'TRIA3'`.
- Output is currently corner-derived, not full MSC grid-point stress averaging.

## `CQUADR`

`CQUADR` is parsed as a separate internal shell type:

```text
TYPE = 'QUADR   '
```

Its stiffness path dispatches to:

```text
CQUADR_DKMQ24
```

The LINK9 output paths have been expanded so `QUADR` is treated with the QUAD4 family where appropriate:

- Element engineering forces.
- Element stresses.
- Element strains.
- FEMAP neutral element result vectors.
- Baseline `OGS1` grid-point surface stress output.

Practical consequence:

- `CQUADR` is not hidden as `QUAD4`; it has its own formulation path.
- Output compatibility intentionally reuses the QUAD4-family output layout where the result topology is four-corner shell data.
- The baseline `OGS1` writer includes `CQUADR` because it checks `TYPE == 'QUADR   '`.

## MITC3+

`MITC3+` is an alternate triangular shell formulation path for `CTRIA3`-family elements.

Current status:

- The active output type remains `TRIA3`.
- F06/OP2/NEU result writers do not need a separate `MITC3+` table name.
- Stress/strain/force values depend on the formulation/recovery arrays filled upstream, then reuse the standard TRIA3 output path.

Important distinction:

- `MITC3+` is a formulation/recovery choice.
- `CTRIA3`/`CTRIAR` are Bulk Data element names.
- OP2/F06 still see these as triangular shell result tables, not a new Nastran element table family.

## Baseline `OGS1` support

The current baseline `OGS1` OP2 writer supports:

- `CTRIA3`
- `CTRIAR`
- `CQUAD4`
- `CQUADR`
- `CQUAD8`

Current pyNastran readback on the local `GPSTRESS` patch-test deck confirms:

```text
['OUGV1', 'OES1X1', 'OGS1', 'OUGV1', 'OES1X1', 'OGS1']
grid_point_surface_stresses = 2
```

Limit:

This is corner-derived `OGS1` output. It is structurally OP2-compatible, but it is not yet the full MSC/Nastran `GPSTRESS` surface/volume averaging and grid-point stress recovery algorithm.

## `CTRIA6` status

`CTRIA6` is not currently implemented in this branch.

No active source support was found for:

- `CTRIA6`
- `TRIA6`
- `NCTRIA6`
- `MEDAT_CTRIA6`
- `BD_CTRIA6`

Therefore, adding only an OP2/F06 writer would be incorrect. There is no element stiffness, recovery, or stress array source for `CTRIA6` yet.

Minimum safe implementation path for real `CTRIA6` support:

1. Add parser/counting support:

```text
BD_CTRIA6 or extended BD_CTRIA
NCTRIA6
MEDAT_CTRIA6
ELMTYP entry 'TRIA6   '
NELGP = 6
```

2. Add EMG formulation:

```text
TRIA6 stiffness/mass/recovery path
material angle/offets/thickness handling
pressure/thermal/differential stiffness where applicable
```

3. Add LINK9 recovery sizing:

```text
NUM_SEi for TRIA6
MAXREQ_OGEL handling
ELEM_STRE_STRN_ARRAYS handling
CALC_ELEM_STRESSES / CALC_ELEM_STRAINS handling
```

4. Add output support:

```text
F06 stress/strain/force rows
OP2 OES/OEF/OSTR table mapping
baseline OGS1 rows
NEU element result vectors
```

5. Validate against MSC/NX:

```text
displacement
element force
center stress/strain
corner/grid stress/strain
OP2 readback through pyNastran
```

Until those steps exist, `CTRIA6` should stay documented as unsupported rather than partially accepted.

## Practical validation note

For current shell patch tests:

- Use `CTRIA3`/`CTRIAR` for 3-node triangle comparisons.
- Use `CQUAD4`/`CQUADR` for 4-node quadrilateral comparisons.
- Use `CQUAD8` only where the branch limitations are acceptable.
- Do not use `CTRIA6` as a validation target yet.

## SAP2000 shell problem 2-004

`prob_2_004_thin.dat` and `prob_2_004_thick.dat` are the SAP2000 twisted
cantilever shell benchmarks.  The source `.s2k` deck uses a regular shell
section with drilling DOF enabled, unit modifiers, and a 12 x 2 curved/twisted
quadrilateral mesh.

Validation targets in `cases_shell.txt` are the SAP2000/independent published
values from Problem 2-004.  The MSC reference deck is already close to those
published values, so MSC remains a useful reference for this case.

| Quantity | SAP2000/independent target | MSC/reference behavior | MYSTRAN `CQUADR` |
| --- | ---: | ---: | ---: |
| Thin `OUT`, tip `TZ` | `0.001749` | within about `1%` | `0.0016041`, about `8.3%` low |
| Thin `IN`, tip `TY` | `0.005429` | within about `1%` | `0.0053329`, about `1.8%` low |
| Thick `OUT`, tip `TZ` | `0.001749` | within about `1%` | `0.0016110`, about `7.9%` low |
| Thick `IN`, tip `TY` | `0.005429` | within about `1%` | `0.0053585`, about `1.3%` low |

Local MYSTRAN element-path checks showed the regular `CQUAD4` path is not usable
for this specific geometry:

| Variant | Maximum displacement error vs SAP/PDF target |
| --- | ---: |
| `CQUAD4` | about `3.5e4%` |
| `CQUAD4` + `PARAM,K6ROT` | about `4.9e3%` |
| `CQUAD4` + MITC4+ path | about `1.8e5%` |
| `CQUADR` / DKMQ24 | about `8.3%` thin, `7.9%` thick |

Conclusion:

- `K6ROT` reduces the singular drilling-mode symptom but does not solve the
  twisted-shell benchmark.
- `CQUADR`/DKMQ24 is the correct current MYSTRAN validation choice for this
  curved/twisted quadrilateral benchmark.
- The remaining roughly `8%` error is a real formulation/recovery limitation
  relative to SAP2000/MSC targets, not an input conversion issue.
- The validation deck was updated to use `CQUADR` so the case measures the
  best available MYSTRAN shell path instead of the known-bad regular `CQUAD4`
  path.

## SAP2000 shell problem 2-008c

`prob_2_008c_thick.dat` and `prob_2_008c_thin.dat` are modal shell checks
using `SOL 103` and an MSC-reference comparison.

The validator summary can look severe because the current `msc;` bulk compare
prints absolute eigenvalue differences for `REALEIGENVALUES`.  Since the
eigenvalues are order `1e5` to `1e7`, a visually large absolute difference can
still be a small relative error.

Current local F06 comparison:

| Deck | Worst mode | MYSTRAN | MSC/reference | Absolute diff | Relative diff |
| --- | ---: | ---: | ---: | ---: | ---: |
| `prob_2_008c_thick.dat` | 5 | `1.004925E+07` | `1.003714E+07` | `1.2110E+04` | about `0.12%` |
| `prob_2_008c_thin.dat` | 2 | `7.847476E+05` | `7.859006E+05` | `1.1530E+03` | about `0.15%` |

Conclusion:

- The large `Error = 1.2e+04` / `8.6e+03` lines in `failed_mystran18a.txt`
  are absolute-difference artifacts, not large physical/modal errors.
- This case should not outrank true large relative displacement/stress errors.
- A future validation cleanup should either report relative eigenvalue error for
  modal `msc;` checks or use explicit `pth;` eigenvalue checks with percent
  tolerance.

## Current large-error triage notes

The current shell validation list mixes true large relative errors, intentional
formulation-tracking failures, parser/deck compatibility failures, and absolute
modal differences.  Do not treat the largest printed number as the highest
engineering priority without checking what the validator is reporting.

| Problem/deck family | Current observation | Priority interpretation |
| --- | --- | --- |
| `prob_2_004_thin/thick` | `CQUAD4` is unusable on the twisted shell; `CQUADR` reduces the error to about `8%`. | Real shell formulation gap. Battle `MIN4`, `MITC4+`, `DKMQ24`, MSC, and SAP2000 per problem. |
| `prob_2_008c_thin/thick` | Printed failure is an absolute eigenvalue difference; relative error is only about `0.15%` or less. | Validation-reporting cleanup, not a top formulation bug. |
| `prob_2_001_thin/thick` | Irregular elements 1/2/5 are intentionally tracked; MSC CQUAD4 also misses SAP/independent patch-test values on those elements. | Formulation benchmark case. Compare `MIN4`, `MITC4+`, `DKMQ24`, MSC, and SAP2000 in a table. |
| `prob_2_002h_thin` | Grounded `CELAS1` short form is rejected by MYSTRAN as grid/component zero. | Parser/deck compatibility issue for scalar springs, separate from shell formulation. |
| `prob_2_002h_thin` after `CELAS1` patch | Short grounded `CELAS1` form is being made MSC-compatible in `BD_CELAS1`. | Input-compatibility fix; reclassify after rebuild/test. |

Recommended next battle-table columns:

```text
Problem | Quantity | SAP2000/independent | MSC/Nastran | MIN4/CQUAD4 |
MITC4+/CQUAD4 | DKMQ24/CQUADR | K6ROT variant | Decision
```
