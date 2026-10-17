# CPYRAM13 v3 implementation plan, 2026-10-15

Reference: `python-solids/MYSTRAN_CPYRAM13_SAP5010_5011_Thermal_Buckling_v3.zip`.
README.md, CODEX_TASK.md, both runners and cpyram13_operators.py were read completely.

1. Preserve Bedrosian interpolation, standard 13-node connectivity, 39 mechanical DOFs and K8. Preserve CPYRAM5 and other elements.
2. Change independent consistent mass integration to 27-point Jacobi; compare actual Fortran mass with v3 Python.
3. Add integration-point recovery of sigma = D (B u - alpha (N T - Tref)).
4. Replace only CPYRAM13 geometric stiffness with physical initial-stress integration: grad(Ni)^T sigma grad(Nj) times I3. Compression is negative. Support independent 8/27-point Kg without changing K8.
5. Reproduce thermal 5-011 displacement, stress and residual checks using actual Fortran operators, and compare independent Python results.
6. Reproduce 5-010 uniform-prestress buckling screening for both Kg orders. This is not acceptance of a full SAP/Nastran prebuckling solution: supplied S2K lacks the applied compression load.
7. Run existing stiffness/patch/frame regressions. Record independent operator gates and all unresolved failures/risks. Nastran acceptance requires separate results; do not infer it from Python agreement.
