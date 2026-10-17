# NEWSOLID default and CPENTA6 validation — 2026-10-17

Global default is now `SOLIDTYP=NEWSOLID`. `PARAM,SOLIDTYP,LEGACY` explicitly selects LEGACY. Within NEWSOLID, CPENTA6 defaults to ABE; `PARAM,PENTATYP,EAS9` selects the frozen EAS9_V2 formulation. Both selections were built and run separately. No empirical stiffness factors were introduced.

## Output corrections

Corrected solid stress recovery-row/EID indexing and OP2 record sizes, retained resolved recovery CID, and corrected principal/equivalent strain calculation from engineering shear. Linear and quadratic solids sharing an OP2 element type now append to one result record rather than appearing as duplicate time steps and losing EIDs.

`check_newsolid_default.py` passes bitwise omitted-parameter versus explicit NEWSOLID/ABE equality for displacement, SPC forces, stress and strain across three subcases. All eight solid groups are present and identical in F06. Seven groups are covered in OP2: PYRAM13 stress/strain OP2 remains intentionally suppressed by the existing writer. The reader uses a temporary in-memory element-name mapping for code 255; binary output is unchanged.

`check_cpenta6_op2.py` passes ABE and EAS9 recovery parity against float64 Python for regular and rotated two-element fixtures, 5-004, and 5-011 thermal. EID, point GRID, CID, six components, principal and equivalent values are checked. Large-field input preserves coordinates consistently between solvers. Nastran 2512 outputs were generated with the supplied executable and successful analysis logs, then read back. Uniform thermal results agree; ABE/EAS9 versus Nastran standard CPENTA are explicitly classified as formulation differences. The LEGACY warped 5-004 discrepancy (about 0.115–0.116% displacement, 0.125–0.129% stress, 0.133–0.136% strain) is reported only, per user instruction.

## Optimization

An exact-key, thread-local one-entry CPENTA6 state cache avoids repeated condensation at recovery points. All state dependencies invalidate the cache; point-specific recovery is evaluated each time. Cache-bypassed comparisons give maximum difference exactly zero. Latest 20,000-call microbenchmark: 1.921875 s uncached, 0.171875 s cached, 11.18x. An earlier run measured 13x. These are routine microbenchmarks, not full-model speedup claims.

## Known formulation failures and open coverage

ABE 5-001 patch remains KNOWN_FAIL: worst frozen displacement error 0.0793928 and stress error 0.462125. Native displacement agrees with Python; corner stress parity is INVALID where the supplied distorted elements have nonpositive corner Jacobians. EAS9_V2 buckling remains KNOWN_FAIL. These labels preserve observed FAIL status.

For 5-010 recovered-prestress strong-axis buckling, ABE is 269.7684 kip and EAS9 is 269.8246 kip: difference 0.0562 kip (0.0208%). Independent Reference is 127.91 kip and SAP2000 is 127.12 kip; both variants are about 110.9% above Independent and FAIL. Weak-axis ABE is 0.8912347 kip (PASS), EAS9 about 9.680284 kip (FAIL), versus Independent 0.888 kip. Native/Python parity passes; that does not imply accuracy.

5-013 native pore-pressure/effective-stress validation remains missing. Pressure 5-008 acceptance remains INVALID pending original restraint and stress-convention resolution. Plate reference provenance and PYRAM13 OP2 coverage remain limited as recorded in the retained reports. Overall completion remains INCOMPLETE_NATIVE_VALIDATION.

## Reproduction and artifacts

Use the repository `.venv` Python with scripts in `tests/solids`: `check_penta6_abe.py`, `check_penta6_buckling.py`, `check_cpenta6_op2.py` (add `--nx` to regenerate Nastran outputs), `check_newsolid_default.py`, `check_penta6_cache.py`, and `report_penta6_abe.py`. Nastran execution requires the working host license context and a writable scratch directory.

Final build and regressions passed: kernel operators, patch/frame, EAS feature limits, CHEXA8_V2, CPYRAM13_V3, and tetra smoothing selection. Existing reported pressure/parser and coverage limits are retained.

Retained benchmark deliverables: `tests/solids/cpenta6_abe_validation` (37 identical-model ABE/EAS9 DAT pairs, parity CSV, Independent/SAP comparison, known_failures.csv, hashes and manifests). Current OP2/Nastran/default artifacts are retained under its `op2_validation` folder; full generated solver outputs remain in `tests/solids/generated/cpenta6_op2_nastran`.
