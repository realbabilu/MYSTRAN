# Quadratic shell force and SIMOQ8 recovery update — 2026-10-19

The generic GPFORCE recovery treated basic-coordinate quadratic shell stiffness as a local-coordinate matrix. It multiplied KE_basic by U_local and then applied a further local-to-basic force rotation. Rotated element axes produced large, nonphysical force-balance residuals even though assembled applied loads and solved displacements were unchanged.

GP_FORCE_BALANCE_PROC now explicitly uses UEB for the five active Q8 and four active T6 basic-DOF kernels. TRANSFORM_NODE_FORCES has an optional BASIC_INPUT argument so these forces pass directly to the existing basic-to-GRID-CD output conversion. Omitted arguments retain the legacy local-input path. Assembly stiffness, quadratic pressure integration and constraints are unchanged. The existing thermal-vector subtraction is retained; the mechanical tests below do not certify thermal/dynamic/nonlinear recovery. Other nodal-force output paths are outside this change.

The SIMOQ8 output path now shares SHELL_POINT_TENSOR_TO_FRAME, using the basis at the same evaluation point as the recovered membrane/bending tensor. It preserves the native center-frame OP2 convention. This makes the previously equivalent point-to-center operation explicit; it does not claim to fix a production displacement or stiffness error. Shear recovery and curved-shell analytical accuracy are not certified by the pointwise tests.

## Validation

- Build succeeds; a final smoke run completes normally and verifies banner `18.0.A OCT 19 2026` in terminal and F06.
- GPFORCE: 38 subcase checks cover all nine formulations on native 4x4 quarter meshes, rotated geometry with GRID CD10, and a legacy CQUAD4 guard. All displacement, SPCFORCES and OLOAD arrays are bitwise identical before/after. The maximum corrected force residual is below1.50e-11 lb. Legacy CQUAD4 GPFORCE also remains bitwise identical. The rotated diagnostic uses rounded coordinates and is not an exact rotational-invariance assertion about displacement.
- SIMOQ8: existing validated source snapshots are byte-identical to the three included pointwise-recovery source changes. The analytical/Python-v13/Fortran diagnostic covers100 flat CW/CCW samples and1800 tensor components. Seven earlier before/after solver regressions preserve displacement/eigenvector/stress/force arrays wherever produced. These CSVs are retained as historical validation of the identical recovery source; they are not represented as freshly rerun after the banner change.
- Separate PLOAD4 audit of72 MYSTRAN and16 researcher-run Nastran outputs verifies uniform pressure loads on flat regular quarter meshes4/8/16/32. Applied nodal vectors agree to native OP2 precision, and displacement/SPCF arrays remain bitwise identical. This supports a recovery-only diagnosis; existing formulation convergence differences remain.

The GPFORCE CSV executable hash precedes the date-only banner rebuild. That rebuild changes executable provenance, not the tested recovery source. Native OP2/F06, diagnostic DAT, executable backups and internal Python stay local. Only numeric validation CSVs are committed here. The local codex_mod archive records source, patch and build hashes.

The earlier Nastran CCW LC2 analytical discrepancy in 2-001 remains unresolved. No full patch-test, shear or curved-shell accuracy PASS is inferred. No kernel change was made to force agreement.
