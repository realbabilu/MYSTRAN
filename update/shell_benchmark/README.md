# Shell benchmark library

## Quadratic shells
[All comparisons](comparison_all.md)

- [2-001: Membrane/bending patch; pointwise recovery](benchmark_2_001.md) — Joint and element stresses require correct frame; Q8 transverse shear recovery LC2 under review.
- [2-002: Fixed A/Q8 and E/T6 with midpoint BC](benchmark_2_002.md) — Single source mesh each; midpointSPC126, six load cases, bar errors and retained pivot warnings.
- [2-003: Plate/shell validation](benchmark_2_003.md) — Q8 7-mesh full-clamp convergence and 6×1 values; T6 support parity pending.
- [2-004: Fixed-source twisted beam Q8/T6](benchmark_2_004.md) — Source12×2; root midsides clamped, T6 two triangles per quad; tip-load sensitivity bars and explicit MACQ8D SPCF-output failure. No convergence study.
- [2-005: 2×10 in C/D plate quarter convergence](benchmark_2_005.md) — Corrected reference: C uniform 12.97 in; C point 16.95 in; D uniform 2.60 in; D point 7.25 in. T6 point-load under investigation.
- [2-006: Scordelis-Lo roof](benchmark_2_006.md) — Free-edge error bars and seven-point displacement profiles re-exported; moment axes pending.
- [2-007: S2K-based hemispherical shell thick convergence](benchmark_2_007.md) — Complete Q8/T6 native OP2, nx4/8/16/32; independent-reference ratios and errors, 44 matching runs.
- [2-008: Thick plate modal / diagonal CONM2](benchmark_2_008.md) — CSI five independent periods: 0.01723, 0.00704, 0.00281, 0.00219, 0.00193 s; plot f/fref.
- [2-009: Plate on elastic foundation](benchmark_2_009.md) — EXCLUDED invalid older model; foundation A/B stiffness must be corrected and rerun.
- [2-010B: Fixed24×48 pressure cylinder](benchmark_2_010.md) — TargetsUr+0.12414in, topUZ−0.12414in, hoop60ksi; native OLOAD verified; fiber/corner stress nonuniformity explicit.
- [2-012A: Fixed6×96 annulus, thick shear/bending](benchmark_2_012.md) — One source mesh; midpoint inner support and consistent line load; targetUZ−0.00534in. Nastran9137 failures and BAILOUT warnings retained.
- [2-013: Thermal constant through thickness](benchmark_2_013.md) — 2-013A Ux=0.000132, Uy=0.000066 in at joint8; 2-013B sigma_xx=sigma_yy=-733.33 psi, with measured GPSTRESS.
- [2-014: Annular thermal gradient](benchmark_2_014.md) — 9×16 CSI thick reference: Uz=0.01931 in, R2=0.00352 rad, M11=-2.181 and M22=-4.406 kip-in/in; Python-generated 28×32 Q8 / 1792 T6 triangles is not the CSI mesh.
- [2-016: Out-of-plane buckling](benchmark_2_016.md) — CSI shear/bend references t=1: 425/1165 ksi; t=0.01: 0.0425/0.1165 ksi. Buckling lambda to critical stress requires load normalization, including t=1 provisional.
- [2-017: In-plane buckling](benchmark_2_017.md) — CSI reference Pcr=127.91 kip; lambda requires load basis and physical buckling mode matching.

[Publication and Solid privacy audit](PUBLICATION_QA.md) · [Actual tested DAT archive](DAT/quadratic_shell_benchmark_DAT.zip)

## Linear shells

[All comparisons](comparison_all_linear.md)

Case IDs retain source-model filenames; actual element and formulation are given by their columns and deck selectors. Two-axis 2-002 tip observables use the two physical corner nodes.

NATIVE means the Nastran result for the corresponding element card: CQUAD4, CQUADR, CTRIA3 or CTRIAR. It does not mean the same proprietary interpolation as a MYSTRAN selector. Curves use independent/analytical references; published SAP2000 values are separate table columns. A completed solve is not a physical-accuracy pass. Zero references use absolute residuals; missing and failed output is never replaced by zero.

The tested source snapshot includes DKMT18 SF as an explicitly labelled candidate and corrected DSG3 local-z Kg. It is documented in [source provenance](data_linear/provenance_linear.json), with a reproducible element-local patch. This benchmark snapshot is not a declaration that SF passed every production-promotion gate. No shared eigensolver prototype is included.

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

[Linear comparisons](comparison_all_linear.md) · [Linear ZIP](comparison_all_linear.zip) · [Linear publication QA](PUBLICATION_QA_linear.md) · [Linear DAT archive](DAT/linear_shell_benchmark_DAT.zip)

CQUAD4, CQUADR, CTRIA3 and CTRIAR use separate `_linear.md` reports. NATIVE denotes Nastran for the corresponding element card.


[Quadratic comparison tables (ZIP)](comparison_all.zip) — tables at published precision, with source reports and qualifications.
