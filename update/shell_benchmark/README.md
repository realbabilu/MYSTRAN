# Shell benchmark library

[All comparisons](comparison_all.md)

- [2-001: Membrane/bending patch; pointwise recovery](benchmark_2_001.md) — Joint and element stresses require correct frame; Q8 transverse shear recovery LC2 under review.
- [2-002: Cantilever six loading branches](benchmark_2_002.md) — Support/rotational DOF parity and low-mesh pivot warnings under review.
- [2-003: Plate/shell validation](benchmark_2_003.md) — Q8 123456 clamped fix audited; T6 support parity pending.
- [2-004: Thin and thick shell response](benchmark_2_004.md) — Comparison only when location and same independent reference fixed.
- [2-005: 2×10 in C/D plate quarter convergence](benchmark_2_005.md) — Corrected reference: C uniform 12.97 in; C point 16.95 in; D uniform 2.60 in; D point 7.25 in. T6 point-load under investigation.
- [2-006: Scordelis-Lo roof](benchmark_2_006.md) — Several newer quadratic variants require regression checks.
- [2-007: Shell benchmark](benchmark_2_007.md) — Do not use old exported-probe envelope as reference PASS.
- [2-008: Thick plate modal / diagonal CONM2](benchmark_2_008.md) — CSI five independent periods: 0.01723, 0.00704, 0.00281, 0.00219, 0.00193 s; plot f/fref.
- [2-009: Plate on elastic foundation](benchmark_2_009.md) — EXCLUDED invalid older model; foundation A/B stiffness must be corrected and rerun.
- [2-010: Shell benchmark](benchmark_2_010.md) — No verified independent-reference figure in consolidated dataset.
- [2-012: Shell benchmark](benchmark_2_012.md) — Do not infer validated reference from solver parity.
- [2-013: Thermal constant through thickness](benchmark_2_013.md) — 2-013A Ux=0.000132, Uy=0.000066 in at joint8; 2-013B sigma_xx=sigma_yy=-733.33 psi, with measured GPSTRESS.
- [2-014: Annular thermal gradient](benchmark_2_014.md) — 9×16 CSI thick reference: Uz=0.01931 in, R2=0.00352 rad, M11=-2.181 and M22=-4.406 kip-in/in; Python-generated 28×32 Q8 / 1792 T6 triangles is not the CSI mesh.
- [2-015: Shell benchmark](benchmark_2_015.md) — Pending independent-reference and run metadata verification.
- [2-016: Out-of-plane buckling](benchmark_2_016.md) — CSI shear/bend references t=1: 425/1165 ksi; t=0.01: 0.0425/0.1165 ksi. Buckling lambda to critical stress requires load normalization, including t=1 provisional.
- [2-017: In-plane buckling](benchmark_2_017.md) — CSI reference Pcr=127.91 kip; lambda requires load basis and physical buckling mode matching.

[Publication and Solid privacy audit](PUBLICATION_QA.md) · [Actual tested DAT archive](DAT/quadratic_shell_benchmark_DAT.zip)
