# DKMQ24 corrections and validation

24 October 2026. Applies to the plain CQUADR/DKMQ24 branch.

- Project basic-coordinate nodal rotations onto each shell normal in the consistent drilling operator. A fixed basic RZ term artificially stiffened rotated shells.
- Assemble equivalent membrane loads for mean temperature change, including TREF. Through-thickness temperature-gradient bending remains unsupported by this branch.
- Apply the scalar initial-stress kernel to all three translations, without a direct rotational term. The former transverse-only default suppressed in-plane buckling modes.

## Validation

Isolated incremental builds retained existing configured objects/libraries; this was not a full clean rebuild. Exact-DAT Python comparisons covered all six load cases on six 2-002 meshes (36 observables). Mean-temperature sign/zero/TREF checks passed 5/5. Corrected 2-013a free expansion: UX=0.000132 in, UY=0.000066 in.

Hemisphere 2-007 refinement n8/n12/n16/n24/n32 gives response/reference 1.012039, 0.998687, 0.994732, 0.993590, 0.993875. The corrected buckling elastic/prestress/geometric matrices and first mode match Python; 2-017 critical loads are 153.6169 and 133.2703 kip, against 127.91 kip (fine-mesh error 4.19%). All 71 retained DKMQ24 configurations completed in the geometric-fix build.

Closure applies to these documented defects, not all shell features. The original benchmark measurements and default executable were preserved. No shared solver correction is included. An earlier thermal-only incremental build hit a recovery bounds error that did not recur in the later 71-case build; its cause remains untraced. Pure temperature-gradient response in 2-014b is still zero and remains open.

Local retained evidence lives under shell_linear_benchmark_20261023/investigations: dkmq24_python_parity, case_06_dkmq24_2_007_lc1, case_11_dkmq24_2_013a_lc1, and case_15_dkmq24_2_017_lc2. These development artifacts are not bundled in this source-only publication. Each local build records source/executable/archive hashes and commands.
