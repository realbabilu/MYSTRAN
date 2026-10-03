# Native T6 OP2 validation — 2026-10-03

This supersedes the earlier OP2 limitation in `T6_UPGRADE_VALIDATION_2026-10-03.md` for the tested T6 stress/force/GPSTRESS results.

## Changes

- Register TRIA6 in the OES and OEF table-name dispatchers. Previously stress payloads could be written without a valid OES header, and force dispatch entered `ERROR STATE TRIA6`.
- Emit native CTRIA6 type 75, rather than labeling T6 as CTRIA3 type 74.
- OES uses 70 words per element: EID, CEN/, center plus three corner locations, two fiber layers, eight stress values per layer.
- OEF uses 38 words per element: EID, CEN/, center plus three corner locations, eight force/moment values per location.
- Recover all seven T6 samples internally even for CENTER-only requests, so the quadratic OP2 payload contains real corner values. F06 selectors still control visible rows; midside F06 and GPSTRESS output remain available.
- Use the correct per-element OGEL stride for CENTER force output.
- Reopen OES when its ITABLE state is zero after an OGS table closes.
- Write recovered OGS data separately for each requested SURFACE using its actual ID. The old hardcoded ID 100 combined different surfaces and caused pyNastran array-size failures in mixed decks.

No element stiffness formula changed in this OP2 work. Q8 formulation was not changed; shared OES/OGS state handling now supports the mixed deck.

## Siemens reference

Ran `C:/Program Files/Siemens/Femap 2512/nastran/bin/nastran.exe` on `test/t6_upgrade/op2_reference/nastran_t6.dat`. Job completed normally. The counterpart `mystran_t6.dat` has the same T6 patch geometry/material/boundary conditions, requests PLOT/PRINT results and POST=-1, and selects the final MYSTRAN MITC6 formulation. The Siemens deck removes MYSTRAN-specific formulation/debug parameters.

pyNastran reads both files. For both subcases, their native T6 stress tables have type 75, width 70, shape `(1,80,8)` and element IDs 51–60. Their subcase-2 T6 force tables have type 75, width 38, shape `(1,40,8)` and the same element IDs. This comparison checks the table contract and identity, not equality between different solver element formulations.

## Numerical and regression checks

- All four T6 families: OES stress has 80 F06-compared rows per subcase (10 elements × 4 locations × 2 fibers); OEF force has 40 compared rows in bending; OGS surface 600 has 50 compared rows per subcase (25 grids × 2 fibers). All pass.
- `test/duel3a.OP2`: complete pyNastran read succeeds, including both displacement subcases, linear/quadratic result tables and six surface IDs. T6 CENTER stresses, CENTER moments, and surface-600 GPSTRESS match F06. The original mixed deck selects CENTER, so local corner values are additionally checked using the four-family CENTER/CORNER decks.
- F06 stress A/B/C regression: 8/8 pass. Standalone CENTER/CORNER/GPSTRESS selectors: 32/32 pass.
- Build successful: `MYSTRAN/Binaries/mystran.exe`.

Artifacts:

- `test/t6_upgrade/op2_validation.txt`
- `test/duel3a_pynastran.txt` and `test/duel3a.pynastran.json`
- `test/t6_upgrade/stress_recovery/<FAMILY>.pynastran.txt` and `.pynastran.json`
- `test/t6_upgrade/op2_reference/{nastran_t6,mystran_t6}.pynastran.json`

Linear/Q8 numerical OP2 completeness, strain/thermal/modal OP2, and direct Femap GUI import are outside this validation. Reading a mixed file successfully does not establish numerical completeness for every other family.

## Example pyNastran reader

From `C:/PROJECTAI/18a`:

```powershell
python test/read_t6_op2.py test/duel3a.OP2
python test/check_t6_op2_reference.py
```

Basic API:

```python
from pyNastran.op2.op2 import OP2

model = OP2(debug=False)
model.read_op2(r'C:\PROJECTAI\18a\test\duel3a.OP2')
print(model.get_op2_stats(short=True))
stress = model.op2_results.stress.ctria6_stress
force = model.op2_results.force.ctria6_force
for key, result in stress.items():
    print(key, result.element_node, result.data.shape)
```

pyNastran may return tuple result keys containing subcase, analysis and surface information; do not assume every key is integer 1 or 2. The supplied validation reader uses `result.isubcase` and `result.ogs_id`.
