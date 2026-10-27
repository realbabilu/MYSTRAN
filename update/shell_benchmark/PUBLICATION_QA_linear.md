# Linear shell publication QA

25 October 2026. 1826 MYSTRAN and 332 Nastran configurations, plus eight Nastran unit-displacement rank decks. All results are fresh native runs, not substituted Python values. All Nastran comparisons use the label Nastran.

## Recorded corrections

- KMITC3+ replaces retired MITC3+HB in CTRIAR; all 22 current selectors are explicit.
- MYSTRAN snapshot is source base21540039 plus the supplied element-local SF adaptation and DSG3 transverse-Kg correction. Exact patch, source hashes and executable SHA256 are retained; candidate status is explicit. No shared solver defaults or iterative-refinement/Jacobi code is included.
- SOL105 uses explicit per-deck ARP_TOL=1e-10, unchanged residual acceptance1e-6 and default DARPACK. Nastran varying-load STATSUB selection is kept separate from MYSTRAN PRELOAD syntax. First Pcr/stress comparisons are not six-mode mode-shape certification.
- SPC1 cards are split to six IDs. 2-017 inactive246 covers every node. Root/support/load definitions and the common physical models are preserved.
- Four centroid triangles per parent quadrilateral remain the linear triangle convention. 2-003c and2-008c are the selected variants. No invented results for absent problems2-011/2-015.
- Two-axis 2-002 observables use the two physical tip corners; averaging interior tip nodes would understate torsion and is excluded. The 192x16 SF ratio is 0.99610044.
- Added LC4-only two-axis beam refinements and retained hemisphere n8/16/32 convergence. Hemisphere GRID-CD translations are converted to basic XYZ; the negative-Y response is no longer mistaken for a local tangential component.
- Independent references remain primary. SAP2000 printed values and source meshes are separate secondary columns. The 2-010 principal-stress envelopes are not relabelled as pointwise hoop stress. Zero references use absolute errors; missing/failed data are blank and explicitly flagged.

## Limits carried forward

Corrected DSG3 BEND six-mode acceptance still fails in this snapshot and remains SOLVER_ERROR. DKMT18 SF has confirmed pipe-stress trade-offs, inherited rigid-rotation stabilization and open shear-patch/output gates. Uniform-shear operator consistency is not a blanket physical-patch pass. Thermal displacement/recovery and rank figures are only the actual recorded observables; this report does not infer all directional tensors or physical mode matching from them. Completion is not accuracy acceptance.

## Publication controls

Quadratic Markdown, HTML, figures and DAT archives are preserved. New linear reports use `_linear.md` and new figures use `figures/linear_...png`; no linear HTML is generated. The existing README gains linear navigation and omits its unsupplied2-015 index entry. Each metric carries units, source/deck/output hashes, candidate tag, archive entry and native run status. CSV size is below10,000,000 bytes. Raw outputs and licensed binaries remain local.

[Publication validation](data_linear/publication_validation.json) records CSV arithmetic, unique keys, links, archive CRC, deck hashes and quadratic preservation checks.

The current rebuilt binary and checksums are distributed through the [October release](https://github.com/realbabilu/MYSTRAN/releases/tag/mystran18a-october2026). Use the 27 October executable hash in current provenance; earlier dated snapshot descriptions above are historical. ZIP archives are marked binary in Git attributes.

## Additional convergence and preload verification

2-003c now includes six meshes (6x1 through 180x6),132 additional MYSTRAN and24 Nastran completed runs. Totals:1826 MYSTRAN and332 Nastran configurations; the CSV contains5564 metrics. One additional SOL101 confirms DSG3 BEND LC1 preload (99.8876ksi/100ksi); the SOL105 LC2 six-mode residual failure remains open. Earlier counts above describe the original publication snapshot.

## Rebenchmark - 27 October 2026

All 1,826 MYSTRAN configurations rerun using rebuilt DSQK mixed hybrid; 1,825 completed and one DSG3 BEND six-mode known failure remains. The 332 Nastran reference configurations retain prior unchanged-deck outputs. Default residual acceptance is 1e-6; BKRTOL is absent. See [current report](../Linear_Shell_Rebenchmark_20261027.md), [run status](data_linear/run_status_linear.csv), and [source/build provenance](data_linear/provenance_linear.json). Earlier sections describe historical snapshots.
