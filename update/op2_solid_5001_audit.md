# MYSTRAN solid GPSTRESS OP2 — 5-001

Implemented the two result types actually present in the supplied successful
Simcenter 2512 5-001 reference: OGS1 direct (27, 9 words) and principal
(28, 15 words). Default VOLUME writes both; explicit DIRECT or PRINCIPAL
selects the corresponding result. Output uses corner grids for all eight
solid types. Metadata words 1–23 and ordered GRID IDs match the reference.
The file remains named for the MYSTRAN input model, with the MYSTRAN writer
and executable version 18.0.A from MYSTRAN_Version. The standard OP2 tape
signature and XXXXXXXX format marker are retained for reader compatibility;
these are not a claim to a Simcenter solver version. No Simcenter product/version
label is added.

The observed table-28 payload differs from the installed pyNastran docstring:
three principal values precede the row-wise direction-cosine matrix, pressure
and von Mises. Reference principal order is maximum, minimum, intermediate.
Repeated eigenvalues have zero direction columns in this reference. General
eigenvector signs and degenerate bases are not unique. LAPACK DSYEV computes
principal results from the recovered average tensor, not averaged eigenvalues.

Validation: build successful; 11 GPSTRESS regression tests pass. Existing
OUG/OES/OSTR payloads remain unchanged. The independent record parser checks
Fortran framing, table widths, metadata and ordered grid IDs. The principal
test checks eigenvalues, the eigen-equation for emitted nonzero axes, pressure,
von Mises, and DIRECT/PRINCIPAL selection. Installed pyNastran reads direct
volumes; its principal reader does not populate data, so raw principal records
are independently checked instead.

## Full 5-001 comparison

Actual output: `nastran_test/mystran/quadratic_patch_audit/5-001_PENTATYP_EAS9.OP2`.
Reference input retains the original topology convention; MYSTRAN input is the
previously audited equivalent mapped FIXED_QUADRATIC deck, NEWSOLID/EAS9.
No coordinates or boundary values changed for this OP2 work.

Both result types contain eight volumes, with row counts
47, 47, 16, 16, 16, 16, 23, 23 (204 rows each).
All raw float payload slots, including direction cosines, match exactly at
OP2 single precision for volumes 1,2,4,5,6,7,8.
Volume 3 CPENTA6 differs: maximum relative direct/invariant difference 0.227690;
principal-value/invariant difference 0.151397. The reference CPENTA6 fails the
analytical constant-strain patch, whereas MYSTRAN EAS9 passes. Do not alter
correct EAS9 stresses merely to reproduce the reference's distorted-patch error.
This establishes output coverage/layout and 7/8 numerical reference matches,
not full solver numerical parity. Existing ordinary CHEXA8 corner-Jacobian
diagnostic remains separate from GPSTRESS results in this combined deck.

Reproduce: run `tests/gpstress/audit_op2_5001.py`; results are exported to
`tests/gpstress/op2_5001_comparison.csv`. No commit or push performed.

5-001.dat 
---------


