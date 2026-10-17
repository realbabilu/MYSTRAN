# CHEXA8 TIED integration — 18 October 2026

Status: **IMPLEMENTED / PROVISIONAL. Results from 20 Nastran comparisons are included; accuracy acceptance remains scoped and provisional.** No surrogate run is counted as official acceptance.

## Selection

For CHEXA8 on NEWSOLID, omitted HEXATYP now selects TIED. `PARAM,HEXATYP,TIED` selects the corrected Stage 6 operator. `PARAM,HEXATYP,EAS9` selects the original formulation, with its original production stiffness path and unchanged legacy state algorithm. CHEXA20 and other solid formulations are not changed. There are still only 24 translational element DOFs.

`PARAM,HEXATAX,ZETA` is the default orientation policy. XI and ETA explicitly select another natural transverse axis. AUTO is a separate **opt-in**, not an assurance of universal orientation safety. HEXAXIS retains its pre-existing SPLITD/SIDE12 meaning. Selectors accept lowercase values and reject unknown values. F06 records formulation, policy, selected natural axis and fallback, and explicitly labels TIED provisional.

AUTO uses the mean length of four corresponding edges for each natural direction. It accepts the shortest direction only when its length is at most half the next-shortest mean length, the maximum edge-length ratio within any family is at most 1.5, and the maximum absolute cosine between mean edge directions is at most 0.3. Otherwise it explicitly falls back to ZETA and logs AUTO_FALLBACK=T. Near-cubic and strongly distorted bricks therefore do not silently select a minimum-length direction. Explicit XI/ETA/ZETA remains available. These criteria are reproducible policy gates, not a proof for arbitrary distorted geometries.

Cyclic axes are implemented by cycling the sampling directions and centroid frame in the original nodal coordinates. This is algebraically equivalent to right-handed natural-node permutation followed by inverse DOF scattering. The matrix tests explicitly compare against that independently permuted Python construction.

Restart files now store both selectors after PENTATYP. Older restarts without HEXATYP retain EAS9, rather than reinterpreting saved stiffness as TIED; invalid stored selectors are rejected.

## Operators and integration

TIED uses full 2x2x2 physical integration, positive determinant checks at evaluated points, physical engineering strains xx/yy/zz/xy/yz/zx, centroid physical shear projectors and covariant EAS9 mapping. The volume-weighted mean of the tying modification is subtracted from every B operator before stiffness, thermal and recovery evaluation. Kaa is solved by LAPACK Cholesky with DPOCON reciprocal-condition rejection below 1e-12; there is no explicit inverse.

The production HEXA stiffness assembly has a distinct TIED branch. EAS9 retains the previous assembly. Thermal RHS includes the internal enhanced load and condensed state; recovery uses corrected B and enhanced strains. Existing corner extrapolation and GPSTRESS interior extrapolation remain separate from direct interior evaluation. KG uses the physical gradients and the stress recovered from the selected formulation and actual static preload. Compatible mass integration is unchanged.

## Evidence

- 36 independent operator cases: regular, distorted, thin and spatially rotated bricks; XI/ETA/ZETA; nu=0.3 and 0.4999. Largest relative error across K, thermal RHS, effective B, strain, stress and KG: **5.53e-13**. Exactly six rigid modes and no resolved negative/hourglass modes, using spectral tolerance 1e-11 of the largest eigenvalue. EAS9 operators are bitwise equal to the frozen original source.
- 20 full-model baseline comparisons: current explicit EAS9 and executable saved before implementation have identical displacement/eigenvalue and solid-stress results at output precision.
- 80 full MYSTRAN solves completed normally. All 80 OP2 files were read and solid stress coverage checked. This is solver execution evidence, not automatic benchmark acceptance.
- 5-001 uses the exact seven-brick distorted patch with eight free interior grids and the original outer grids. Six separate affine engineering strains and a rigid rotation were run. Displacement, all 63 center/corner stress rows and all 16 GPSTRESS grids pass analytical checks. Maximum center/corner stress absolute error over tested formulations/states is 2.15e-11. BASIC output coordinates are explicit; comparing local-frame stress components directly against BASIC targets would be incorrect.
- Default/explicit/lowercase selector OP2 results agree byte-for-byte; invalid formulation/axis values are rejected. Both COUPMASS=0 and 1 modal checks agree with independently assembled K/M eigenvalues within 2e-7 relative error.
- AUTO policy tests cover cubic ambiguity, clear thin XI/ETA/ZETA directions and skew fallback. Full 8x8 dome cyclic node relabeling changes output displacements by at most 7.46e-9 (single-precision OP2 output).
- The existing multi-solid kernel regression passes with historical CHEXA8 references explicitly selecting EAS9.

The compact table below lists MYSTRAN and published/analytical values. The full benchmark CSV now also includes all 20 newly completed Nastran results, their explicit IN=2 selection, matched differences and buckling mode MAC.

| Case | Response | TIED default ZETA | TIED + AUTO | Reference |
|---|---|---:|---:|---:|
| 5-004 | IN mean loaded tip displacement | see CSV | 0.00532654 | 0.005429 |
| 5-004 | OUT mean loaded tip displacement | see CSV | 0.00171882 | 0.001749 |
| 5-007 | 8x8x1 point A displacement | 0.0152698 | 0.0691112 | 0.094 |
| 5-007 | 16x16x1 point A displacement | 0.0700826 | 0.0923827 | 0.094 |
| 5-007 | 24x24x1 point A displacement | 0.0871244 | 0.0932866 | 0.094 |
| 5-010 | weak critical factor, actual preload | see CSV | 0.8909257 | 0.888 |
| 5-010 | strong critical factor, actual preload | see CSV | 127.1197 | 127.91 |
| 5-011 | restrained uniform SXX/SYY | about -733.333 | about -733.333 | -733.333333 |

The default ZETA policy is **not** equivalent to AUTO on the dome: its 24x24 error remains **7.314%**; AUTO gives **0.759%**. The coarse 8x8 AUTO dome also has a **26.48%** error. These accuracy limitations are retained and reported, not tuned away or hidden by fallback.

5-004 uses the native CSI S2K geometry, material, supports and nodal forces. Dome 16x16 retains the supplied A model; 8x8 is obtained by retaining every second coordinate of A, and 24x24 retains every second coordinate of the supplied 48x48 B model. No retained joint is moved. This is documented source-mesh coarsening; it is not the bundled reconstructed Python beam/dome benchmark.

5-008 uses CSI S2K meshes C (5 radial elements) and E (20 radial elements). The verification document explicitly states that Uz is inactive: every axial translation is restrained to enforce plane strain. Local tangential restraints are represented by MPCs. Unit pressure on the actual planar inner quad is integrated exactly into consistent nodal forces in both solver decks, with the resultant recorded. There is no arc-area surrogate. Cases cover nu=0.3, 0.49, 0.499 and 0.4999 on the five-element mesh, and nu=0.49 on the twenty-element mesh. Inner displacement and radial/hoop/axial corner stresses are tabulated. Default TIED inner displacement errors are approximately 0.95–1.39% on the five-element mesh and 0.0735% on the twenty-element mesh; these are analytical comparisons only. Native PLOAD4 parser functionality is not certified by an equivalent-nodal-load test.

5-010 uses the existing CHEXA8 SOL105 source model with one-kip static preload, not constant prescribed prestress or Euler-beam-only validation. MYSTRAN has explicit static and buckling subcases. New Nastran mode matches and eigenvalue differences are now reported in the CSV; archived-reference fields, where present, remain separately labeled. 5-011 covers source restrained uniform temperature, free uniform expansion and an additional through-thickness gradient; the gradient response is measured, not certified against an independent gradient solution.

At 24x24, elapsed times in a two-job run were approximately 11.50 s for EAS9, 16.01 s for TIED and 16.10 s for TIED AUTO. These are indicative wall times with concurrent jobs, not isolated performance measurements. Peak memory and larger-mesh scaling remain unmeasured.

## Nastran execution and acceptance scope

Twenty Nastran cases produced F06/OP2 results. Each imported DAT is byte-identical to its prepared regression deck. `manual_nastran_import.json` records input/output hashes.

The new Nastran inputs explicitly use PSOLID IN=2. They are a different integration selection from the archived source/default Nastran inputs: on the 24x24 dome, new Nastran A=0.00211390248; in buckling the weak/strong factors are 26.7418/155.2266, versus TIED/AUTO 0.8909257/127.1197. The corresponding translation mode MAC values remain about 0.999989/1.000000. Factor agreement and mode-shape agreement are distinct checks. Historical default-integration Nastran numbers are retained only in the historical archive, not copied into the new table.

All six benchmark categories now have actual Nastran execution evidence. The new CSV includes Nastran scalar responses, differences and mode MAC. Patch GPSTRESS is read from the OP2 volume-direct table rather than a page-dependent F06 parse. Analytical accuracy and the large new full-integration Nastran differences remain explicit; this is not a blanket release acceptance. General native pressure-card compatibility, gravity/body-force/global initial-strain regression and large-mesh memory/scaling remain uncertified. TIED stays provisional with EAS9 available unchanged.

## Reproduce from H:/Projects/MYSTRAN/source

PowerShell (the repository already has its configured CMake build and GCC runtime):

```powershell
$env:PATH='C:\gcc\bin;C:\gcc\openblas32\bin;'+$env:PATH
C:/gcc/bin/make.exe -C build -j8
../.venv/Scripts/python.exe tests/solids/hexa8_tied/build_probes.py
../.venv/Scripts/python.exe tests/solids/hexa8_tied/check_operators.py
../.venv/Scripts/python.exe tests/solids/hexa8_tied/run_regression.py
../.venv/Scripts/python.exe tests/solids/hexa8_tied/collect_regression.py
../.venv/Scripts/python.exe tests/solids/hexa8_tied/check_features.py
../.venv/Scripts/python.exe tests/solids/hexa8_tied/check_orientation.py
../.venv/Scripts/python.exe tests/solids/check_kernels.py
```

To rerun the Nastran comparisons, run `run_regression.py --nastran-only` and then `collect_regression.py`. This preserves MYSTRAN results and reruns the Nastran jobs. The collector accepts a new Nastran result only after a completed solve, and then records matched displacement differences. New results still require benchmark accuracy review before release acceptance.

The original historical `check_hexa8_v2.py` probe segfaults in its production recovery harness on the saved baseline executable; this pre-existing harness failure was preserved in `update/hexa8_eas9_baseline_20261018.log`. It is not listed as a passing test. The new direct operator probe, full solver patch/recovery checks and legacy kernel checks provide separate evidence.

## Files

Changed production files: `Source/Modules/PARAMS.f90`, `Source/Modules/EAS9_FRAME.f90`, `Source/LK1/L1A-BD/BD_PARAM.F90`, `Source/UTIL/READ_L1A.f90`, `Source/UTIL/WRITE_L1A.f90`, `Source/EMG/EMG5/HEXA.f90`.

Reproducible scripts and frozen operator baseline: `tests/solids/hexa8_tied/`. Every executed model has a DAT deck, F06, OP2 (when the solve completed) and stdout under `regression/`. `manifest.json` records model sources, mesh derivation and pressure conventions. `benchmark_matrix.csv`, `runs.csv`, `eas9_baseline_parity.csv`, `output_coverage.csv`, `operator_matrix.csv`, `feature_checks/feature_matrix.csv` and `feature_checks/orientation_matrix.csv` contain the evidence.

The final rebuild adds only the provisional F06 notice after the numerical regression build. Final selector/mass/policy checks and rebuilt direct operator probes passed against that rebuild; no formulation changed in that last rebuild. Executable/source hashes are recorded in `artifact_hashes.json`; individual numerical runs retain their own executable hashes in `runs.csv`.
