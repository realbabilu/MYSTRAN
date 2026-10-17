# CTETRA4 NEWSOLID — Python buckling parity fix

Only `Source/EMG/EMG5/TETRA.f90` production source was changed in this fix.
The new branch applies only to `TYPE=TETRA4` and `SOLIDTYP=NEWSOLID` at OPT(6),
LOAD_ISTEP>1. LEGACY and TETRA10 execute the historical code unchanged.

Root cause: the old path calls ELEM_STRE_STRN_ARRAYS, which rotates STRESS to
the requested output CID, then combines it with local CBAR/shape gradients.
Changing PSOLID CORDM from 0 to -1 changed the 4x1x24 buckling factor from
1.494327 to 1.438734 without changing physical geometry, load or restraint.
The historical rotational operator also differs from the initial-stress
operator used by the chosen Python reference.

The NEWSOLID T4 path now recovers local prestress directly from ES*B*UEL,
subtracts local thermal stress where applicable, and integrates
Kg_ij = integral(grad Ni transpose sigma grad Nj) I3, matching the Python
reference. Elastic nodal smoothing remains alpha=.9; no scalar calibration.

Build: cmake --build source/build -j 12, successful.
Executable SHA256: 0cdab20974b032d1429f86e9fc4fc10068c2e75062f8102607d42140a331e713.

Two meshes, 14 eigenpairs each, CORDM=0/-1: eigenvalue relative error <4e-7,
minimum MAC >.99999994. Output-CID eigenvalues are bitwise equal and mode MAC
exceeds 1-1e-10. LEGACY factors and eigenvectors are bitwise unchanged.
Patch 5-001, twisted beam 5-004 and uniform thermal 5-011 controls have
bitwise unchanged displacement and stress.

Reproduce: `.venv/Scripts/python.exe source/tests/solids/check_ctetra4_python_kg.py`.
Evidence: `CTETRA4_NEWSOLID_CODEX_HANDOFF_v3_20261017/CASE_5010/parity_fixed`.
This matches the Python reference's SMOOTH elastic K plus standard recovered
initial-stress Kg. It is not a claim of a newly derived fully nonlocal
SMOOTH geometric tangent.
