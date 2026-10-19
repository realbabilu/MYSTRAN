# Publication QA — 19 October 2026

The supplied CQUAD8/CTRIA6 publication draft is reconciled with locally audited data. Nine individual theory pages, their combined manuscript, and16 benchmark case pages are retained in Markdown and HTML. No solver runs or formulation code changes were made during publication.

## Gates

- Current publication trees contain Markdown/HTML, SVG/PNG, license text and input archives; no tracked CSV/Python/notebook/raw solver outputs in Solid or Shell publication paths.
- Research source files removed from Git tracking remain byte-identical on the local research computer. This is current-tree cleanup, not rewriting Git history.
- Local file/image targets and Markdown anchors pass static verification. Regenerated numerical tables share the same private input records as their plots; private CSV is not copied into public folders.
- Plot sample inspection confirms readable labels/markers; multimesh studies use distinct markers/dashes and independent ratio/error references; single-mesh thermal cases use bars.
- HTML TeX is preserved and uses MathJax from its public CDN. Browser screenshot automation could not connect to the local preview server in this environment; interactive rendering is not represented as screenshot-verified. PNG numerical figures were inspected directly.
- DAT archive CRC and every manifest SHA256 are verified. The archive contains593 actual input decks across14 verified-run cases plus Markdown README/hash manifest, with no results, Python, CSV or private filesystem paths.

## Reconciliation

2-005 uses the correctedC/D supplied independent targets, audited quarter-only4/8/16/32 observations, and the complete PLOAD4-to-nodal audit. No old zero probes or CW result is substituted. Independent-error and normalized-ratio graphs replace parity-based accuracy plots.

2-013B uses the latest native GPSTRESS for both Q8/T6, target geometry, Z1/Z2 and basic-aligned SURFACE100. Nastran-only MID rows and MYSTRAN-only midsides are not silently paired. 2-014 retains the actual Python28×32 mesh and declared radial/tangential moment recovery; older cross-case thermal figures are omitted.

2-008 uses diagonal CONM2 and frequency references1/T. CSI mode shapes were unavailable, so ordered-frequency comparisons do not claim MAC identity. 2-016/017 normalization and physical mode matching remain PENDING; provisional stress/force overlays have been removed. Raw signed mode-index1 diagnostics and the complete eigenvalue tables are labeled accordingly. Nastran t=.01 fatal7340 is a failure/missing result, notzero. 2-009 remains EXCLUDED. Case2-015 is PENDING and has no input archive entry without a verified-run receipt.

The theory pages preserve the supplied mathematical content, primary-reference list and limitations. Bibliographic items explicitly labeled for later title/DOI confirmation in the supplied manuscript remain provisional; publication here is research documentation, not a claim those bibliographies were independently verified for journal submission.

## Solid privacy cleanup

No Solid formulation math or result tables were rewritten. Only research-data links/prose were adjusted to retain embedded public tables while removing links to private source files. The following28 tracked research files were untracked with local bytes preserved:

- `update/Formulation/Solid/data/comparison_all.csv`
- `update/Formulation/Solid/data/constraint_audit.csv`
- `update/Formulation/Solid/data/documentation_checks.json`
- `update/Formulation/Solid/data/finest_mesh_summary.csv`
- `update/Formulation/Solid/data/geometry_audit.csv`
- `update/Formulation/Solid/data/mesh_inventory.csv`
- `update/Formulation/Solid/data/provenance.json`
- `update/Formulation/Solid/data/sap2000_reference.csv`
- `update/solid_benchmark/data/comparison_all.csv`
- `update/solid_benchmark/data/constraint_audit.csv`
- `update/solid_benchmark/data/element_recovery_5001_20261018.csv`
- `update/solid_benchmark/data/finest_mesh_summary.csv`
- `update/solid_benchmark/data/five_case_mesh_inventory.csv`
- `update/solid_benchmark/data/geometry_audit.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/artifact_hashes.json`
- `update/solid_benchmark/data/hexa8_tied_20261018/benchmark_matrix.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/eas9_baseline_parity.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/feature_matrix.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/manifest.json`
- `update/solid_benchmark/data/hexa8_tied_20261018/manual_nastran_import.json`
- `update/solid_benchmark/data/hexa8_tied_20261018/operator_matrix.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/orientation_matrix.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/output_coverage.csv`
- `update/solid_benchmark/data/hexa8_tied_20261018/runs.csv`
- `update/solid_benchmark/data/hexa8_tied_vs_eas9v2_archive.csv`
- `update/solid_benchmark/data/mesh_inventory.csv`
- `update/solid_benchmark/data/mystran_solid_comparison.csv`
- `update/solid_benchmark/data/sap2000_reference.csv`

## Input archive coverage

|Case|DAT count|
|---|---:|
|2-001|6|
|2-002|97|
|2-003|42|
|2-004|61|
|2-005|88|
|2-006|22|
|2-007|30|
|2-008|33|
|2-010|18|
|2-012|6|
|2-013|24|
|2-014|13|
|2-016|120|
|2-017|33|

Existing Solid input archives were inspected for forbidden result/Python/CSV entries; none were found. One retained Solid archive includes its existing Nastran batch launcher; it contains no solver output. No production Python scripts elsewhere in the repository were blanket-ignored or removed.
