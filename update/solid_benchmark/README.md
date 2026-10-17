# Solid benchmarks: MYSTRAN vs Nastran

Comparison reports for SAP2000-derived solid benchmarks, with independent reference values, published SAP2000 values, mesh counts and refinement curves. All text is in English.

| Report | Coverage | Read on GitHub | Standalone HTML |
|---|---|---|---|
| MYSTRAN solid comparison | Five cases: 5-001, 5-004, 5-007, 5-010 and 5-011; 71 response comparisons | [Report](mystran_solid_comparison.md) | [Download HTML](mystran_solid_comparison.html) |
| All solid elements - SAP2000 benchmarks | Eight solid types on 5-004 and 5-007; 56 matched meshes, 63 MYSTRAN selections and 126 responses | [Report](comparison_all.md) | [Download HTML](comparison_all.html) |

The all-solid report covers CTETRA4/10, CPENTA6/15, CHEXA8/20 and CPYRAM5/13. CPENTA6 ABE and EAS9 are shown separately. Beam curves use three mesh levels; dome curves use four. The bold ratio reference is 1.0; the percentage-error plot uses a separate 1% reference line.

## Input download

[Download input decks only](solid_benchmark_inputs.zip). The ZIP contains input DAT files, instructions and an input manifest. **No Nastran or MYSTRAN output files are included.** See [input checksums](input_manifest.json).

## Numerical summaries

- [Five-case comparison CSV](data/mystran_solid_comparison.csv)
- [All-solid comparison CSV](data/comparison_all.csv)
- [Mesh inventory](data/mesh_inventory.csv)
- [Geometry audit](data/geometry_audit.csv)
- [Constraint audit](data/constraint_audit.csv)

## Reading the comparison

Nastran and MYSTRAN are compared on matched inputs. Independent values and published SAP2000 values are separate reference columns. SAP2000 values retain their original mesh provenance; no SAP2000 result is fabricated at an unreported mesh size. The all-solid analyses were run with both solvers; the original five-case report also uses previously supplied Nastran reference outputs as documented.

Accuracy depends on the case, mesh and formulation. Agreement between solvers does not establish agreement with independent. The reports retain actual numerical differences and documented geometry-check settings. No production stiffness code or empirical scaling factor was changed for these comparison runs.

GitHub renders the Markdown reports and PNG figures directly. Download the HTML reports to use the local element selector in a browser.

## Reference Hex Benchmarks
https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Solids/Problem%205-001.pdf to https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Solids/Problem%205-013.pdf
Computers & Structures, Inc. (2024). SAP2000 is a proven FEA well-known product for civil engineering by Computers & Structures, Inc.
