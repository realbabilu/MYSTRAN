# MYSTRAN Response Spectrum Load-to-Mass Spec

## Scope
This note defines the intended direction for response spectrum support in MYSTRAN, starting from simple planar-frame verification problems such as SAP2000 Problem 1-020 and then generalizing to broader models.

## 1. Active DOF handling for planar-frame RS problems
For planar frame verification problems like Problem 1-020, only the in-plane structural DOF should be active:

- `Ux`
- `Uz`
- `Ry`

All other DOF should be inactive in the effective analysis model.

### Practical interpretation in MYSTRAN
There are two acceptable ways to achieve this:

1. Explicitly restrain unused DOF in the deck.
2. Let unused/out-of-plane DOF be removed naturally through `AUTOSPC`, provided this does not distort the intended physics.

For verification work, the preferred reading is:
- if the benchmark clearly intends a plane frame, the model should behave as a plane frame;
- extra DOF should not carry mass or spectrum participation accidentally.

## 2. Mass only in the intended excitation direction
For Problem 1-020, the building mass is only active in the global X direction.
That means:

- no vertical excitation mass contribution;
- no out-of-plane excitation mass contribution;
- participation should be controlled by the intended translational excitation component only.

### If mass comes from explicit lumped mass
If the benchmark mass is defined from nodal/applied mass, this is straightforward:
- assign only the intended mass component;
- leave other translational directions at zero.

### If mass comes from density or self weight
If mass is generated from density/self-weight style inputs, the solver should have a clean conversion rule from weight to mass.

## 3. Mass model structure
The correct structure is that `LOAD2MASS` is only one additional contributor to the total analysis mass. It is not a replacement for the existing mass sources.

### Intended total mass model
For modal and response spectrum analysis:

```text
M_total =
  M_density
+ M_NSM / M_NSML
+ M_CONM2 / M_CMASS
+ M_LOAD2MASS
```

This means:
- density-derived structural mass still contributes;
- explicit nonstructural mass still contributes;
- lumped/equipment mass still contributes;
- selected load combinations may contribute additional equivalent mass through `LOAD2MASS`.

## 4. Existing mass sources remain active
### 4.1 `CONM2` / `CMASS`
These remain the correct mechanism for explicit lumped mass such as:
- equipment,
- machinery,
- bearing mass,
- discrete attached hardware.

Example:

```text
CONM2,1001,501,0,2.5
```

This must still enter the assembled mass matrix directly.

### 4.2 `NSM` / `NSML`
These remain the right mechanism for nonstructural mass that is naturally attached to elements as distributed mass.

Examples:
- asphalt,
- cladding,
- parapet,
- permanent utilities,
- finishes that are better represented as mass per area or per length.

Example:

```text
NSM,10,PSHELL,100,0.025
```

and if/when supported similarly:

```text
NSML,20,100,0.025
```

These must still remain valid and still contribute to the mass matrix.

## 5. `LOAD2MASS` is an additional mass source from load combinations
`LOAD2MASS` is intended for situations where mass source data already exists naturally in load form instead of direct mass form.

Typical use cases:
- asphalt already modeled as load,
- parapet/barrier load already modeled,
- utility/service load already modeled,
- some fraction of live load (for example 30 percent) treated as seismic mass.

### Example concept
```text
LOAD,990,1.0,1.0,800,0.30,900
PARAM,LOAD2MASS,YES
PARAM,LOAD2MASSID,990
```

Meaning:
- the selected `LOAD` combination `990` is converted into additional equivalent mass;
- that equivalent mass becomes `M_LOAD2MASS`;
- it is then added to the normal mass assembly, not substituted for it.

So the result is:

```text
M_total = existing mass + mass from LOAD 990
```

not:

```text
M_total = mass from LOAD 990 only
```

## 6. Gravity and weight-to-mass conversion
The clean rule is to let the user declare gravity directly:

```text
PARAM,GRAV,9810.0
```

and then let MYSTRAN derive internally:

```text
WTMASS = 1.0 / GRAV
```

### Example: N-mm-tonne style
```text
PARAM,GRAV,9810.0
```

then:

```text
WTMASS = 1 / 9810.0
       = 1.01936799185E-4
```

### Example: N-m-kg style
```text
PARAM,GRAV,9.80665
```

then:

```text
WTMASS = 1 / 9.80665
       = 0.101971621
```

So the exact `WTMASS` depends on the user’s declared unit system through `GRAV`.

## 7. Proposed precedence rule for `GRAV` and `WTMASS`
There are two practical user styles:

1. user provides `GRAV` only;
2. user provides explicit `WTMASS` only;
3. user provides both.

### Recommended behavior
#### Case A: only `PARAM,GRAV,...` exists
- derive `WTMASS = 1.0 / GRAV` internally;
- use this for weight-to-mass conversion in response spectrum / auto-mass workflows.

#### Case B: only `PARAM,WTMASS,...` exists
- use the supplied `WTMASS` directly;
- no inferred gravity is required for conversion.

#### Case C: both `WTMASS` and `GRAV` exist
- use `WTMASS` as the controlling factor for weight-to-mass conversion;
- use `GRAV` for self-weight/static load direction and magnitude;
- issue a warning if `WTMASS` and `1.0/GRAV` are not numerically consistent within tolerance.

## 8. Proposed auto-mass control
Suggested parameter family:

```text
PARAM,AUTOMASS,YES
PARAM,GRAV,9810.0
PARAM,GRAVDIR,0.0,0.0,-1.0
```

### Intent
- `AUTOMASS,YES`
  - enable automatic conversion from weight-style loading or density/self-weight interpretation into analysis mass for response spectrum use.
- `GRAV`
  - gravity magnitude in current unit system.
- `GRAVDIR`
  - gravity direction vector.

## 9. Overhauled response spectrum concept
This is the key conceptual change.

Response spectrum processing should be staged as:

1. build modal basis from the structural model and total mass;
2. interpret the response spectrum input according to its declared format;
3. compute directional response results separately:
   - `RSX`
   - `RSY`
   - optionally `RSZ`
4. only after directional results exist, perform the requested combination.

So the architecture becomes:

```text
modal solve
   -> read response spectrum (`PARAM,RSTYPE,PERG`)
   -> RSX
   -> RSY
   -> RSZ (optional)
   -> combination stage
```

This is cleaner than embedding a fixed combination rule too early in the solver path.

### Why this is better
- easier to validate each direction independently;
- easier to compare against commercial programs and benchmark tables;
- easier to expose writer-side output for raw directional sets;
- easier to support multiple combo rules without rebuilding the whole RS pipeline.

## 10. Response spectrum input format
For civil and structural workflows, the response spectrum input should be controlled by one compact parameter family:

```text
PARAM,RSTYPE,PERG
```

with the intended meanings:

```text
PARAM,RSTYPE,PERG   = period vs Sa/g
PARAM,RSTYPE,PERA   = period vs absolute acceleration
PARAM,RSTYPE,FRQG   = frequency vs Sa/g
PARAM,RSTYPE,FRQA   = frequency vs absolute acceleration
```

This keeps the user deck compact and readable while still covering both the civil-facing and solver-facing spectrum conventions.

### Preferred civil default
The most natural first-class civil input style remains:

```text
PARAM,RSTYPE,PERG
```

Meaning:
- the response spectrum table is defined as `Period` vs `Sa/g`;
- the abscissa is structural period;
- the ordinate is spectral acceleration normalized by gravity.

This is the style commonly used in civil software and benchmark documents, and it should be supported directly instead of forcing users to pre-convert the spectrum into frequency-based absolute acceleration form.

### Explicitly rejected alternative
The response spectrum format should not be split into separate parameters such as:

```text
PARAM,RSABSC,PERIOD
PARAM,RSORD,SG
```

That style is more verbose, easier to mis-pair, and less natural for ordinary deck authoring than a single `RSTYPE` code.

### Intended solver behavior for `RSTYPE=PERG`
When `PARAM,RSTYPE,PERG` is present:

1. the solver reads `TABLED1` as `T` vs `Sa/g`;
2. the solver converts `T` to the internal working frequency form as needed;
3. the solver converts `Sa/g` to absolute acceleration using the active gravity convention;
4. the converted spectrum is then used in the modal response spectrum pipeline.

### General interpretation rules
- `PERG`
  - read `TABLED1` as `Period` vs `Sa/g`
- `PERA`
  - read `TABLED1` as `Period` vs absolute acceleration
- `FRQG`
  - read `TABLED1` as `Frequency` vs `Sa/g`
- `FRQA`
  - read `TABLED1` as `Frequency` vs absolute acceleration

In all four cases, the solver may convert the input into one internal working representation, but the user-facing deck should preserve the stated meaning of the original spectrum data.

### Gravity coupling
The acceleration conversion should follow the gravity rule in this note:
- if only `GRAV` exists, derive `WTMASS = 1/GRAV`;
- if `WTMASS` also exists, `WTMASS` remains authoritative for weight-to-mass conversion;
- `GRAV` still provides the acceleration scale for `Sa/g -> absolute acceleration`.

This keeps the user-facing deck civil-friendly while preserving a consistent internal representation.

## 11. Combination stage
The combination stage should be explicit and selectable.

### Conceptual user controls
For now the intent can be expressed as:

```text
PARAM,OPTION,SRSS
```

or

```text
PARAM,OPTION,CQC
```

Meaning:
- `SRSS` means combine directional or modal response through SRSS logic;
- `CQC` means combine through CQC logic when correlation data is available.

### Important design reading
The combination selection is post-directional:
- first produce `RSX`, `RSY`, `RSZ`;
- then apply `SRSS` or `CQC` as requested.

This is the conceptual overhaul.

## 12. Directional combination behavior
At present, the local MYSTRAN RS path uses fixed orthogonal combination behavior equivalent to 100/30 through an internal constant.
That is too rigid for the intended next stage.

### Desired behavior
Directional combination should be configurable rather than hard-coded.

Instead of only fixed:
- `100/30`

we want generic percentage-style control, such as:
- `100/40`
- `100/20`
- `100/30/30`
- or any explicit orthogonal percentage rule the user requests.

### Better conceptual model
There are two distinct layers:

1. Directional result generation
   - produce `RSX`, `RSY`, `RSZ`
2. Directional/result combination
   - `SRSS`
   - `CQC`
   - orthogonal percentage rule
   - signed linear combinations if needed for workflow compatibility

These should not be tangled together more than necessary.

## 13. Proposed combo model
A future response spectrum combo specification should support at least:

1. directional result generation
   - `RSX`
   - `RSY`
   - `RSZ`
2. combination method
   - `SRSS`
   - `CQC`
   - orthogonal percentage combination
3. configurable orthogonal percentage factor(s)

### Example conceptual parameters
Illustrative only:

```text
PARAM,OPTION,SRSS
PARAM,RSDIRCOMB,ORTHO
PARAM,RSORTHO,30.0
```

or:

```text
PARAM,OPTION,CQC
PARAM,RSDIRCOMB,ORTHO
PARAM,RSPCT,100.0,30.0,30.0
```

The exact syntax can evolve later. The important part is:
- combination happens after `RSX/RSY/RSZ` are known;
- the solver should not permanently bake in one fixed 30 percent rule.

## 14. Recommendation for verification sequence
The clean order is:

1. Problem 1-020
   - planar frame
   - only `Ux, Uz, Ry`
   - only X-direction mass and excitation
   - use this to debug participation and directional logic first
2. Problem 1-024
   - larger packaging and workflow case
   - use this after the simpler directional logic is trustworthy
3. then extend to writer-side full NEU packaging and combo validation

## 15. Immediate implementation guidance
Short-term:
- keep current solver-side RS path running;
- debug directional participation using Problem 1-020;
- do not assume `AUTOSPC` alone is enough unless the resulting active mass and participation are checked;
- make `PARAM,RSTYPE,PERG` the first-class civil input path before adding broader spectrum-format families.

Next-stage:
- add mass derivation controls (`AUTOMASS`, `GRAV`, `GRAVDIR` or equivalent);
- respect explicit `WTMASS` when present;
- warn on `WTMASS` vs `GRAV` inconsistency;
- add `LOAD2MASS` as an additional mass contributor, not a replacement;
- support direct `Period vs Sa/g` input through `PARAM,RSTYPE,PERG`;
- split directional result generation from final combination;
- replace fixed `100/30` logic with configurable percentage-based orthogonal combination;
- document the combo rules in the response spectrum notes and validation package.

## 16. Summary
The target behavior is:
- plane-frame RS benchmarks behave as true plane-frame models;
- mass participation is only present in intended directions;
- automatic load-to-mass conversion is available when density or self-weight style input is used;
- `GRAV` can drive derived `WTMASS` when `WTMASS` is absent;
- explicit `WTMASS` remains authoritative if both are present;
- `LOAD2MASS` contributes additional equivalent mass on top of normal mass sources;
- `PARAM,RSTYPE,PERG` lets users input spectra directly as `Period vs Sa/g`;
- `RSX`, `RSY`, `RSZ` are generated first;
- `SRSS` or `CQC` is applied after directional results are available;
- directional response spectrum combination is configurable and not permanently fixed to `100/30`.


