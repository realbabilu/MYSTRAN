# Linear shell thermal recovery and DSG3 output correction

Linear triangular shells (including CTRIAR stored internally as TRIA3) store TEMPP1 gradient at DT(4). The shared default recovery previously read DT(5), omitting thermal bending stress. Read DT(4) for TRIA3 only; retain quad, quadratic and solid behavior.

DSG3 engineering-force output already rotates local tensors into basic axes. Remove the subsequent special reversal of its three moment components: output moments now satisfy M=integral(z*sigma(z) dz), consistently with recovered fiber stresses. Preserve tensor rotation, stress recovery and every stiffness/load/solver calculation.

## Verification

Controlled incremental build: 220 restrained runs (before/after all22 formulations, five thermal loads), 110 candidate fiber-stress checks, 44 annulus runs with all22 displacement arrays unchanged. The load support baseline contains separately developed all22 gradient equivalent-load support; this recovery commit does not add those load kernels.

Independent NumPy reconstruction of DSG3 BB recovery reproduces all3584 annulus element fiber stresses within7.50e-5 kip/in2 and original engineering-force moments within1.20e-5 kip-in/in, including printed displacement rounding. The output correction passes8 runs and7461 element/subcase fiber-integrated moment checks; max discrepancy8.34e-6 kip-in/in (printed stress rounding). POST on/off agree. A mechanical 2-005 baseline rerun confirms both load cases retain identical fiber stresses and displacements; the original moment-sign inconsistency is corrected for mechanical loads too.

DSG3 common-sign annulus outer-adjacent center Mrr=-2.179963, inner-adjacent center Mtt=-4.240787 kip-in/in. Published edge references -2.181,-4.406 refer to edge locations, not these adjacent-center samples. No edge-error certification is asserted.

These are flat isotropic default-recovery checks. Free thermal curvature physical-sign certification, exact edge recovery, DEBUG176 SE/STE paths, anisotropy, laminates and arbitrary normals remain open. No eigensolver tolerance/default, other formulation or default executable is changed.

Local evidence: shell_linear_benchmark_20261023/investigations/case_12_tempgradient_2_014b_lc1/thermal_recovery_20261024; portable Python/patch audit: codex_mod/research_exports_20261024/DSG3_Thermal_Recovery_Audit_20261024.zip.
