## RSA OP2 writer audit - 2026-07-15

### Scope

Audit focused on MYSTRAN response-spectrum (`PARAM,SCRSPEC`) OP2 output behavior versus MSC-style reference behavior observed from:

- `D:\18a\femap_RSA\msc\beam_rs_msc2.bdf.op2`

### Findings

1. RSA combined grid results in `LINK9`

`D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90`

The current SCRSPEC compatibility path writes combined RSA:

- displacement via `CALL WRITE_GRD_OP2_OUTPUTS ( 0, RSA_NUM_OUT, 'DISP', ITABLE, NEW_RESULT )`
- velocity via `CALL WRITE_GRD_OP2_OUTPUTS ( 0, RSA_NUM_OUT, 'VELO', ITABLE, NEW_RESULT )`
- acceleration via `CALL WRITE_GRD_OP2_OUTPUTS ( 0, RSA_NUM_OUT, 'ACCE', ITABLE, NEW_RESULT )`

2. Header writer choke point

`D:\18a\MYSTRAN\Source\UTIL\OUTPUT2_WRITE_TABLE.f90`

Grid-result OP2 metadata is controlled by:

- `WRITE_OUG3_STATIC`
- `WRITE_OUG3_EIGN`
- `WRITE_OUG3`

3. Mismatch identified in existing MYSTRAN RSA OP2

Old `beam_RS.op2` readback showed a velocity table header with:

- `table_code = 10`
- `analysis_code = 6`
- `thermal = 0`

That combination is not MSC-style response-spectrum metadata and causes pyNastran to reject or misroute the result.

4. Safe writer patch applied

`WRITE_OUG3` now sets:

- `THERMAL = 4`

when all of the following are true:

- `SOL_NAME(1:5) == 'MODES'`
- `SCRSPEC == 'Y'`
- `TABLE_CODE` is one of:
  - `1` displacement
  - `10` velocity
  - `11` acceleration
  - `3` SPC force

This is intended to align MYSTRAN SCRSPEC grid/SPCF OP2 output with the MSC/NX response-spectrum convention used by downstream readers.

### Current limitation

Fresh verification against `beam_RS.dat` could not be completed in this pass because that deck currently stops earlier in LINK0 on unrelated deck parsing issues (`SPC1` formatting / prior deck compatibility noise), so the new OP2 header path still needs rerun confirmation on a clean RSA deck.

### Follow-up verification on clean deck

A MYSTRAN-safe reference-aligned deck was created and run successfully:

- `D:\18a\femap_RSA\beam_RS_msc2_mystran.dat`

Observed behavior from the resulting OP2:

- file produced: `D:\18a\femap_RSA\beam_RS_msc2_mystran.OP2`
- readback showed normal modal content:
  - eigenvectors present
  - beam stress present
- but no combined RSA grid-result containers were observed by the current OP2Query readback:
  - no `SC/1/DISPLACEMENTS/GID/...`
  - no `SC/1/VELOCITY/GID/...`
  - no `SC/1/ACCELERATION/GID/...`

This means the writer patch to tag RSA-style grid/SPCF result headers is necessary but not sufficient: the combined RSA OUG/OQG payloads are still either:

1. not being written at all,
2. being written with malformed table/subtable structure, or
3. colliding with the existing modal OUG content in a way that prevents downstream readers from materializing them as separate RSA results.

### Most likely next choke points

Primary source files to inspect next:

- `D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90`
- `D:\18a\MYSTRAN\Source\LK9\L91\WRITE_GRD_OP2_OUTPUTS.f90`
- `D:\18a\MYSTRAN\Source\UTIL\OUTPUT2_WRITE_TABLE.f90`

Specific suspicion:

- combined RSA sets currently call `WRITE_GRD_OP2_OUTPUTS(JSUB=0, ...)`
- that path currently routes them through static-style `WRITE_OUG3_STATIC`
- modal OUG content in the same file is also present
- result separation may require a distinct table/header convention closer to MSC response-spectrum OUPV1 handling rather than plain static-style OUG reuse

### Next recommended check

Run a clean SCRSPEC deck that reaches LINK9 and then verify:

- displacement table lands in OP2 as RSA-compatible
- velocity table no longer appears as `analysis_code = 6 / thermal = 0`
- acceleration table no longer appears as `analysis_code = 6 / thermal = 0`
- if written, SPC force table is also tagged consistently
