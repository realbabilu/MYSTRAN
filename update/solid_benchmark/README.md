> **18 October CHEXA8 update:** [TIED vs EAS9v2](comparison_all.md#chexa8) uses the existing beam/dome benchmark meshes and previously published Nastran results. [Comparison data (ZIP)](comparison_all.zip).

# Solid benchmarks: MYSTRAN vs Nastran

Comparison reports for SAP2000-derived solid benchmarks, with independent reference values, published SAP2000 values, mesh counts and refinement curves. All text is in English.

| Report | Coverage | Read on GitHub | Standalone HTML |
|---|---|---|---|
| MYSTRAN solid comparison | Five cases: 5-001, 5-004, 5-007, 5-010 and 5-011; 71 response comparisons | [Report](mystran_solid_comparison.md) | [Download HTML](mystran_solid_comparison.html) |
| All solid elements - SAP2000 benchmarks | Eight solid types on 5-004 and 5-007; 56 matched meshes, 63 MYSTRAN selections and 126 responses | [Report](comparison_all.md) | [Download HTML](comparison_all.html) |

The all-solid report covers CTETRA4/10, CPENTA6/15, CHEXA8/20 and CPYRAM5/13. CPENTA6 ABE and EAS9 are shown separately. Beam curves use three mesh levels; dome curves use four. The bold ratio reference is 1.0; the percentage-error plot uses a separate 1% reference line.

## 5-001 GPSTRESS VOLUME and solid-element patch-test example

- [Validated MYSTRAN input](5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH_FIXED_QUADRATIC_MYSTRAN_EAS9.dat), with normalized BDF formatting, reconciled CPENTA15/CHEXA20 midside ordering and `PARAM,PENTATYP,EAS9`.

The example contains eight independent solid groups: CTETRA4/10, CPENTA6/15, CHEXA8/20 and CPYRAM5/13. It requests GPSTRESS VOLUME direct and principal OP2 results in the basic coordinate system. Coordinates, material and prescribed displacement values retain the supplied 5-001 patch definition.

Use the validated MYSTRAN input to reproduce the audited run. Seven volumes match the tested Nastran OP2 float payloads; CPENTA6 EAS9 passes the analytical patch while differing from the Nastran result. The [18 October recovery correction](../solid_element_recovery_5001_2026-10-18.md) restores ordinary element stress/strain for all 462 elements, including CHEXA8 and CPYRAM13, while preserving this patch's GPSTRESS payloads. See the [OP2 audit](../op2_solid_5001_audit.md) and [quadratic patch audit](../quadratic_patch_5001_audit.md) for the earlier numerical results and scope.

## Benchmark input archive

[Download input decks only](solid_benchmark_inputs.zip). The ZIP contains input DAT files, instructions and an input manifest. **No Nastran or MYSTRAN output files are included.** See input checksums (research data retained locally; see embedded tables).

## Numerical summaries

- Five-case comparison CSV (research data retained locally; see embedded tables)
- All-solid comparison CSV (research data retained locally; see embedded tables)
- Mesh inventory (research data retained locally; see embedded tables)
- Geometry audit (research data retained locally; see embedded tables)
- Constraint audit (research data retained locally; see embedded tables)

## Reading the comparison

Nastran and MYSTRAN are compared on matched inputs. Independent values and published SAP2000 values are separate reference columns. SAP2000 values retain their original mesh provenance; no SAP2000 result is fabricated at an unreported mesh size. The all-solid analyses were run with both solvers; the original five-case report also uses previously supplied Nastran reference outputs as documented.

Accuracy depends on the case, mesh and formulation. Agreement between solvers does not establish agreement with independent. The reports retain actual numerical differences and documented geometry-check settings. No production stiffness code or empirical scaling factor was changed for these comparison runs.

GitHub renders the Markdown reports and PNG figures directly. Download the HTML reports to use the local element selector in a browser.

## Reference Hex Benchmarks
https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Solids/Problem%205-001.pdf to https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Solids/Problem%205-013.pdf
Computers & Structures, Inc. (2024). SAP2000 is a proven FEA well-known product for civil engineering by Computers & Structures, Inc.


[All-solid comparison data (ZIP)](comparison_all.zip) · [Five-case comparison data (ZIP)](mystran_solid_comparison.zip)
