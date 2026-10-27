# Linear shell benchmark — reference-anchored comparisons

Current MYSTRAN results: 27 October 2026 rebuild. Nastran results retained from prior unchanged-deck runs. CQUAD4, CQUADR, CTRIA3 and CTRIAR only.

[All linear comparisons](comparison_all_linear.md) · [Full ZIP](comparison_all_linear.zip) · [Tested DAT](DAT/linear_shell_benchmark_DAT.zip) · [Linear publication QA](PUBLICATION_QA_linear.md)

Case IDs retain source-model filenames; actual element and formulation are given by their columns and deck selectors. Two-axis 2-002 tip observables use the two physical corner nodes.

NATIVE means the Nastran result for the corresponding element card: CQUAD4, CQUADR, CTRIA3 or CTRIAR. It does not mean the same proprietary interpolation as a MYSTRAN selector. Curves use independent/analytical references; published SAP2000 values are separate table columns. A completed solve is not a physical-accuracy pass. Zero references use absolute residuals; missing and failed output is never replaced by zero.

The tested source snapshot includes DKMT18 SF as an explicitly labelled candidate and corrected DSG3 local-z Kg. It is documented in [source provenance](data_linear/provenance_linear.json), with a reproducible element-local patch. This benchmark snapshot is not a declaration that SF passed every production-promotion gate. No shared eigensolver prototype is included.

## Case index

- [2-001: Membrane, bending, shear and twist patches](benchmark_2_001_linear.md)
- [2-002: Straight beam: six load cases](benchmark_2_002_linear.md)
- [2-003: Curved beam: Model C](benchmark_2_003_linear.md)
- [2-004: Twisted cantilever beam](benchmark_2_004_linear.md)
- [2-005: 2 x 10 plate: simple and clamped](benchmark_2_005_linear.md)
- [2-006: Scordelis-Lo roof](benchmark_2_006_linear.md)
- [2-007: Hemispherical shell with opening](benchmark_2_007_linear.md)
- [2-008: Plate vibration: Model C](benchmark_2_008_linear.md)
- [2-009: Plate on elastic foundation](benchmark_2_009_linear.md)
- [2-010: Cylinder under internal pressure](benchmark_2_010_linear.md)
- [2-012: Thick annular plate](benchmark_2_012_linear.md)
- [2-013: Uniform temperature](benchmark_2_013_linear.md)
- [2-014: Through-thickness thermal gradient](benchmark_2_014_linear.md)
- [2-016: Plate out-of-plane buckling](benchmark_2_016_linear.md)
- [2-017: Column in-plane buckling](benchmark_2_017_linear.md)
- [2-050: Pinched cylinder](benchmark_2_050_linear.md)
- [T: Rank and supplementary analytical tests](rank_and_patch_tests_linear.md)

## Formulation matrix

| Element card | PARAM selector | Tested values |
|---|---|---|
| CQUAD4 | QUAD4TYP | `DKMQ20`, `DSQK`, `MIN4`, `MIN4T`, `MITC4`, `MITC4+`, `SIMO` |
| CQUADR | QUADRTYP | `DKM24AU`, `DKM24EA`, `DKMQ24`, `MBP1C0`, `MITC4PD`, `Q4EASANS`, `Q4RS`, `SIMO` |
| CTRIA3 | TRIA3TYP | `DSG3`, `MIN3`, `MITC3+`, `T3FF` |
| CTRIAR | TRIARTYP | `DKMT18`, `KMITC3+`, `T3FFD` |

DKMT18 is SF in this tested candidate. KMITC3+ replaces retired MITC3+HB. The current formulation kernel tags are also columns in the CSV.

## Solver completion

| Solver | Element / formulation | Completed | Failed |
|---|---|---:|---:|
| MYSTRAN | CQUAD4_DKMQ20 | 83 | 0 |
| MYSTRAN | CQUAD4_DSQK | 83 | 0 |
| MYSTRAN | CQUAD4_MIN4 | 83 | 0 |
| MYSTRAN | CQUAD4_MIN4T | 83 | 0 |
| MYSTRAN | CQUAD4_MITC4 | 83 | 0 |
| MYSTRAN | CQUAD4_MITC4plus | 83 | 0 |
| MYSTRAN | CQUAD4_SIMO | 83 | 0 |
| MYSTRAN | CQUADR_DKM24AU | 83 | 0 |
| MYSTRAN | CQUADR_DKM24EA | 83 | 0 |
| MYSTRAN | CQUADR_DKMQ24 | 83 | 0 |
| MYSTRAN | CQUADR_MBP1C0 | 83 | 0 |
| MYSTRAN | CQUADR_MITC4PD | 83 | 0 |
| MYSTRAN | CQUADR_Q4EASANS | 83 | 0 |
| MYSTRAN | CQUADR_Q4RS | 83 | 0 |
| MYSTRAN | CQUADR_SIMO | 83 | 0 |
| MYSTRAN | CTRIA3_DSG3 | 82 | 1 |
| MYSTRAN | CTRIA3_MIN3 | 83 | 0 |
| MYSTRAN | CTRIA3_MITC3plus | 83 | 0 |
| MYSTRAN | CTRIA3_T3FF | 83 | 0 |
| MYSTRAN | CTRIAR_DKMT18 | 83 | 0 |
| MYSTRAN | CTRIAR_KMITC3plus | 83 | 0 |
| MYSTRAN | CTRIAR_T3FFD | 83 | 0 |
| NASTRAN | CQUAD4_NATIVE | 83 | 0 |
| NASTRAN | CQUADR_NATIVE | 83 | 0 |
| NASTRAN | CTRIA3_NATIVE | 83 | 0 |
| NASTRAN | CTRIAR_NATIVE | 83 | 0 |

## Data and reproduction

The CSV contains 5,564 metric rows from 2,158 solver configurations and is 3.28 MB. It is suitable for direct analysis in Python or ChatGPT; preserve status, units, references and hashes when filtering.

[Run status](data_linear/run_status_linear.csv) · [Model audit](data_linear/model_audit.csv) · [Source snapshot](data_linear/provenance_linear.json) · [DAT archive](DAT/linear_shell_benchmark_DAT.zip)

The DAT archive contains every tested input, including failed cases, with separate MYSTRAN and NASTRAN directories and portable run/extract scripts. Licensed executables, raw licensed reference files and full solver outputs are not redistributed.

## Additional convergence and preload verification

2-003c now includes six meshes (6x1 through 180x6),132 additional MYSTRAN and24 Nastran completed runs. Totals:1826 MYSTRAN and332 Nastran configurations; the CSV contains5564 metrics. One additional SOL101 confirms DSG3 BEND LC1 preload (99.8876ksi/100ksi); the SOL105 LC2 six-mode residual failure remains open. Earlier counts above describe the original publication snapshot.

## Rebenchmark - 27 October 2026

All 1,826 MYSTRAN configurations rerun using rebuilt DSQK mixed hybrid; 1,825 completed and one DSG3 BEND six-mode known failure remains. The 332 Nastran reference configurations retain prior unchanged-deck outputs. Default residual acceptance is 1e-6; BKRTOL is absent. See [current report](../Linear_Shell_Rebenchmark_20261027.md), [run status](data_linear/run_status_linear.csv), and [source/build provenance](data_linear/provenance_linear.json). Earlier sections describe historical snapshots.
