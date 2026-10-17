## October 17 update: solids and GPSTRESS - 17 October 2026

The [October release](https://github.com/realbabilu/MYSTRAN/releases/tag/mystran18a-october2026) includes the latest MYSTRAN 18.0.A Windows binary from `v18.00.a`.

- `NEWSOLID` is now the default; select `PARAM,SOLIDTYP,LEGACY` for historical solid paths.
- Enhanced CHEXA8 EAS9, CPENTA6 ABE/EAS9, CTETRA4 smoothing and Bedrosian CPYRAM13 are integrated. CPENTA6 defaults to ABE; `PARAM,PENTATYP,EAS9` selects EAS9. `EAS6` is not supported. CPYRAM accepts 5 or 13 nodes.
- CPENTA15 and CHEXA20 selective shear integration is corrected. Both pass the 5-001 six constant-strain patches and rigid-body rotation audit below 1e-7 relative error, with no extra global DOFs.
- Solid GPSTRESS VOLUME writes OP2 direct tensors, principal stresses and directions, pressure and von Mises for all eight 5-001 solid groups. Corner-grid rows and metadata match the tested Nastran reference. In the tested 5-001 case, seven volumes produce OP2 float payloads identical to Nastran. CPENTA6 EAS9 MYSTRAN is closer to the analytical solution than the Nastran reference. Thus, the tested GPSTRESS results are equal to or better than that reference for this case; this comparison does not extend to all solid benchmarks.
- Solid OP2 recovery indexing, mixed linear/quadratic records and engineering-shear principal strains are corrected. CTETRA4 NEWSOLID geometric stiffness matches the audited Python buckling operator and is independent of output coordinates.

Validation remains scoped: ABE on the distorted 5-001 patch and some CPENTA6 buckling benchmarks remain outside the established analytical acceptance tolerance. Ordinary CHEXA8 corner recovery can report a nonpositive Jacobian on the supplied distorted mesh. GPSTRESS VOLUME coverage does not imply complete ordinary element-stress OP2 coverage: CPYRAM13 OES/OSTR output remains limited.

See the [5-001 OP2 audit](update/op2_solid_5001_audit.md), [quadratic patch audit](update/quadratic_patch_5001_audit.md), [solid defaults and limitations](update/newsolid_op2_default_2026-10-17.md), and [CTETRA4 buckling audit](update/ctetra4_python_kg_parity_2026-10-17.md).

