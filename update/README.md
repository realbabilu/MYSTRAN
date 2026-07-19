# Update Notes

This folder collects concise notes for the `v18.00.a` source backport work done in this repository.

Files:

- `v18_backport_summary.md`
  Summary of the MYSTRAN 18a source-level changes in `CMakeLists.txt` and the Fortran source tree under `Source/`.
  This now also includes short notes on:
  - the July 2026 `NEU` writer architecture cleanup
  - the July 2026 `MEFFMASS/MPFACTOR` compatibility bridge and backend corrections

- `output_routing_guide.md`
  Practical guidance for current `F06`/`OP2`/`NEU` routing in `v18.00.a`, including how to request:
  - `OP2` without detailed result tables in `F06`
  - `NEU` without detailed result tables in `F06`
  - `OP2` + `NEU`
  - classic full output

- `neu_writer_map.md`
  A code-level map of the current `.neu` writer path in `LINK9`, including:
  - file/block structure
  - nodal vector writing
  - element force/stress/strain writing
  - `MODES` and `BUCKLING` set mapping
  - current implementation boundaries

- `neu_reference_v9_audit.md`
  Short audit of the FEMAP v9 text neutral files in `reference_msc`, including:
  - block `100/450/451` observations
  - modal and buckling set-header examples
  - what the current MYSTRAN writer already matches
  - what remains custom or lighter than the reference exports

- `min4_basic_stress_patch_test.md`
  Patch-test note for ordinary shell element stress validation using `PARAM,STR_CID,0`.
  This records the `OES1X1`/`SHELLSTRESSES/EID` path for SAP2000 `2-001` style global element stress checks.

- `gpstress_patch_test_example.md`
  Separate GPSTRESS/GSTRESS example note for MSC-style grid-point stress output.
  This records the `OUTPUT(POST)`/`SURFACE` shape and clarifies that GPSTRESS maps to `OGS1`, not ordinary element stress.

- `gpstress_recovery_design.md`
  Phase-1 implementation design for moving GPSTRESS from the current corner-copy
  compatibility writer toward surface-driven grid-point recovery.
  The current source tree now includes the first metadata foundation for this:
  - structured `SURFACE` parsing in `LOADC`
  - stored `surface_id/set_id/normal` descriptors
  - a helper path to resolve shell element/grid membership from the referenced set
