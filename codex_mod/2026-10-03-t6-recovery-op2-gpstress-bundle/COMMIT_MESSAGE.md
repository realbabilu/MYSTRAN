Align T6 shell formulations, recovery and native OP2 output

Align SIMOT6 v2, MITC6 v4, MH6T v3 and REZAIEE v3 with the final
Python references. Recover CENTER and all six T6 nodes with consistent
fiber-stress and moment signs, and include midsides in GPSTRESS recovery.

Write native CTRIA6 OES/OEF type-75 tables, reopen OES after OGS, and use
separate OGS records with the actual surface IDs. Suppress empty surface
force headers and remove temporary matrix/recovery printing.

Include the existing local recovery-buffer and output-selector changes,
plus a dated codex_mod bundle containing the baseline patch, changed-file
snapshots, Python references, validation inputs and recorded evidence.

Validation:
- Displacement 2-001 through 2-004: 88/88 refined cases and 88/88 recovery
  regression cases pass across the four T6 families.
- Patch-test local CENTER/CORNER and GPSTRESS: 8/8 cases pass.
- Standalone CENTER/CORNER/GPSTRESS selectors: 32/32 checks pass.
- Four-family native OES/OEF/OGS reads match F06 values with pyNastran.
- Siemens Nastran 2512 reference matches T6 table layouts and element IDs.
- Mixed duel3a OP2 reads successfully; empty-header removal preserves
  force, GPSTRESS and displacement values.

Curved/twisted stress recovery, numerical completeness of other shell
families, thermal/strain/modal OP2 and direct Femap GUI import are outside
the recorded validation scope. Q8 formulation was not upgraded here.
