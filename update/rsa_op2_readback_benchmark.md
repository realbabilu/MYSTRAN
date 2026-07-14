## RSA OP2 readback benchmark

Date: 2026-07-14

### Scope

Readback and compare RSA-related OP2 files using:
- `D:\18a\MYSTRAN_Validation-main\op2_query.py`
- `D:\18a\MYSTRAN_Validation-main\compare_op2_vs_op2.py`
- `D:\18a\MYSTRAN_Validation-main\compare_f06_vs_op2.py`

Python runtime used:
- `C:\Users\arypr\AppData\Local\Programs\Python\Python312\python.exe`
- `pyNastran 1.4.1`

### Files checked

- `D:\18a\femap_RSA\beam_RS.OP2`
- `D:\18a\femap_RSA\beam_RS_abs.OP2`
- `D:\18a\femap_RSA\beam_RS_cqc.OP2`
- `D:\18a\femap_RSA\nx_hinge_nx_dense.OP2`
- `D:\18a\femap_RSA\nx_hinge_nx.OP2`
- `D:\18a\femap_RSA\nx_hinge_msc.OP2`

### Result 1: beam RSA OP2 readback was initially empty

Observed with:
- `compare_f06_vs_op2.py D:\18a\femap_RSA\beam_RS --deck D:\18a\femap_RSA\beam_RS.dat`
- `probe_op2_query_tree.py`

Findings:
- `beam_RS.OP2` is a valid binary file and pyNastran sees table names like:
  - `ROUGV1`
  - `srss.displacements`
  - `abs.displacements`
  - `spc_forces`
- However, `OP2Query` currently builds an empty tree for `beam_RS.OP2`:
  - `TOP_KEYS=[]`
  - `LEAF_COUNT=0`
- `compare_f06_vs_op2.py` therefore reports all requested values as missing in OP2:
  - combined `DISPLACEMENTS`
  - modal `SPCFORCES`
  - `REALEIGENVALUES`

Interpretation:
- The current RSA beam OP2 writer was producing OP2 tables that pyNastran identified by name.
- But the nodal/eigenvalue payload was not ending up in the standard pyNastran containers that `OP2Query` currently reads.
- So for RSA beam cases, the initial bottleneck was `OP2 writer and/or pyNastran container population`, not `F06`.

### Result 1A: combined RSA displacement OP2 readback now works after LINK9 routing fix

Rebuilt on 2026-07-14 after updating `Source/LK9/LINK9/LINK9.f90` so the combined RSA
`DISP/VELO/ACCE` write path emits OP2 whenever the OP2 unit is already open and the family is
requested for either `PRINT` or `PLOT`.

Observed with:
- `probe_op2_query_tree.py D:\18a\femap_RSA\beam_RS.OP2`
- `compare_f06_vs_op2.py D:\18a\femap_RSA\beam_RS --deck D:\18a\femap_RSA\beam_RS.dat`
- same for `beam_RS_abs` and `beam_RS_cqc`

Results:
- `beam_RS.OP2`
  - `TOP_KEYS=['SC']`
  - `SUBCASES=['1']`
  - `LEAF_COUNT=12`
- `beam_RS_abs.OP2`
  - `Matching leaves = 12`
  - `Differing leaves = 0`
- `beam_RS_cqc.OP2`
  - `Matching leaves = 12`
  - `Differing leaves = 0`

What now matches:
- the combined RSA displacement vector values requested through `F06`

What still does not appear in OP2 for the standard non-`PLOT` RSA beam decks:
- modal `SPCFORCES`
- `REALEIGENVALUES`

Interpretation:
- The empty-tree failure is fixed for the combined nodal displacement payload.
- The remaining gap is narrower: modal supporting tables still rely on the classic plot-oriented path
  and are not yet mirrored into the standard RSA OP2 path for `beam_RS*.dat`.

### Result 1B: standard `beam_RS.dat` now carries modal OP2 content without `PLOT`

Additional fixes applied on 2026-07-14:
- `Source/LK9/LINK9/LINK9.f90`
  - for `SOL MODES + PARAM,SCRSPEC`, if the OP2 unit is already open, `PRINT` requests for
    `DISP/VELO/ACCE/SPCF/MPCF` are promoted to OP2-family output before `OFP1/OFP2` run
- `D:\18a\MYSTRAN_Validation-main\op2_query.py`
  - modal nodal blocks now map to `SC/#/MODE/#/...` when pyNastran exposes `modes`
  - RSA eigenvector `mode_cycles` are normalized from rad/s to cycles/s when needed

Observed with:
- `compare_f06_vs_op2.py D:\18a\femap_RSA\beam_RS --deck D:\18a\femap_RSA\beam_RS.dat`

Latest result:
- `Matching leaves = 56`
- `Missing in OP2 = 0`
- `Missing in F06 = 4`
- `Differing leaves = 4`

Interpretation:
- From the OP2 side, all requested benchmark leaves for `beam_RS.dat` are now present.
- The remaining 4 `Differing leaves` are only floating-point precision drifts in
  `REALEIGENVALUES/{EIGENVALUE,CYCLES}`.
- The remaining 4 `Missing in F06` are negative-sign modal eigenvector components that the current
  F06 parser does not retain, not missing OP2 content.

### Result 1C: `beam_RS_plot.dat` remains a separate solver issue

After rerunning `beam_RS_plot.dat` with the updated binary:
- LINK4 stops with `*ERROR 7112` from `dsaupd`, `INFO = -9`
- message text: initial Arnoldi residual vector is zero

Additional note on 2026-07-14:
- `Source/LK4/LINK4.f90` was adjusted so the small-model warning path no longer claims
  `LINK4 WILL USE THE CONDENSED DENSE REFERENCE SOLVER` when the active RSA compatibility rule
  is `SCRSPEC + single SUPORT`, because that path is immediately forced back to ARPACK in
  `EIGRL_EXTRACT_SOLVERS.F90`.
- This is a logging/dispatch-consistency fix only. It does not yet solve the underlying
  `INFO = -9` ARPACK startup failure for `beam_RS_plot.dat`.

Interpretation:
- this is not an OP2 readback issue
- it is the small-model `SCRSPEC + single SUPORT + plot-path` eigen extraction problem already
  known in the ARPACK compatibility path
- do not use `beam_RS_plot.dat` as the acceptance benchmark for the OP2 readback work until that
  solver-side issue is stabilized

### Result 2: hinge OP2 element-result readback matches NX/MSC closely

Observed with:
- `compare_op2_vs_op2.py D:\18a\femap_RSA\nx_hinge_nx_dense.OP2 D:\18a\femap_RSA\nx_hinge_nx.OP2 --path-filter SHELLSTRESSES`
- `compare_op2_vs_op2.py D:\18a\femap_RSA\nx_hinge_nx_dense.OP2 D:\18a\femap_RSA\nx_hinge_msc.OP2 --path-filter SHELLSTRESSES`

Results:
- MYSTRAN vs NX:
  - `MATCHING_LEAVES=17756`
  - `DIFFERING_LEAVES=4`
  - `MISSING_IN_TEST=0`
  - `MISSING_IN_REF=0`
- MYSTRAN vs MSC:
  - `MATCHING_LEAVES=17756`
  - `DIFFERING_LEAVES=4`
  - `MISSING_IN_TEST=0`
  - `MISSING_IN_REF=0`

The only reported differences were tiny shell-stress drifts:
- `SC/1/SHELLSTRESSES/EID/140/CORNER/3/Z1/YY`
  - test `8.497752999999999`
  - ref `8.497790999999999`
- `SC/1/SHELLSTRESSES/EID/140/CORNER/3/Z2/YY`
  - test `-8.497752999999999`
  - ref `-8.497790999999999`
- `SC/1/SHELLSTRESSES/EID/144/CORNER/4/Z1/MAJOR`
  - test `4.51924`
  - ref `4.519252`
- `SC/1/SHELLSTRESSES/EID/144/CORNER/4/Z2/MINOR`
  - test `-4.51924`
  - ref `-4.519252`

Interpretation:
- For hinge RSA OP2, the element-result portion that `OP2Query` currently maps is already benchmark-quality against both NX and MSC.
- The residual mismatch is at floating-point noise level for shell stresses.

### Result 3: nodal RSA OP2 readback is not currently available even for the references

Observed with direct `pyNastran` inspection:
- `len(op2.displacements) == 0`
- `len(op2.spc_forces) == 0`
- `len(op2.eigenvectors) == 0`
- `len(op2.op2_results.srss.displacements) == 0`
- same outcome for:
  - `beam_RS.OP2`
  - `nx_hinge_nx_dense.OP2`
  - `nx_hinge_nx.OP2`

Interpretation:
- For these RSA OP2 files, missing nodal readback is not specific to MYSTRAN.
- The NX reference OP2 in this benchmark folder also does not populate the standard pyNastran nodal namespaces for the RSA combined result.
- Therefore:
  - comparing RSA nodal combined displacement through `OP2Query` is not yet a valid benchmark path
  - `F06` remains the reliable benchmark source for RSA nodal combined vectors right now

### Practical conclusion

Current status is split:
- `F06`:
  - usable for RSA nodal combined output benchmarking
- `OP2`:
  - usable for hinge element-result benchmarking
  - not yet usable for RSA nodal combined vector benchmarking through current `pyNastran + OP2Query`

### Next recommended work

1. Extend the OP2 audit around the actual RSA tables present in the file:
   - `ROUGV1`
   - `srss.*`
   - `abs.*`
2. Determine whether MYSTRAN's RSA nodal OP2 payload layout is nonstandard, or whether `pyNastran 1.4.1` simply does not map these tables for this result family.
3. Keep using `F06` as the source of truth for RSA combined nodal displacement/velocity/acceleration until that binary readback path is solved.
