# Solid element stress and strain recovery: 5-001 correction

18 October 2026. Scope: the supplied eight-solid 5-001 distorted patch and the regression cases listed below.

## Result

The corrected MYSTRAN executable writes ordinary OP2 **element stress and strain for all 462 elements**, with **2506 center/corner rows in each result**. The earlier executable wrote 413 elements, omitting all seven CHEXA8 and 42 CPYRAM13 elements. CPYRAM5's 210 zero corner-stress rows are also corrected.

On the unchanged 5-001 input, all six tensor components for CTETRA4/10, CPENTA15, CHEXA8/20 and CPYRAM5/13 match the freshly computed Nastran results at OP2 single precision. CPENTA6 EAS9 retains the analytical constant-strain solution; the Nastran prism differs. Matching tensor values is not a claim of byte-identical entire OP2 files or identical eigenvector direction conventions.

[Before/after comparison CSV](solid_benchmark/data/element_recovery_5001_20261018.csv).

| Element | Before | Corrected | Nastran | Corrected stress point rows |
|---|---:|---:|---:|---:|
| CTETRA4 | 168 | 168 | 168 | 840 |
| CTETRA10 | 168 | 168 | 168 | 840 |
| CPENTA6 | 14 | 14 | 14 | 98 |
| CPENTA15 | 14 | 14 | 14 | 98 |
| CHEXA8 | 0 | 7 | 7 | 63 |
| CHEXA20 | 7 | 7 | 7 | 63 |
| CPYRAM5 | 42 | 42 | 42 | 252 |
| CPYRAM13 | 0 | 42 | 42 | 252 |
| **Total** | **413** | **462** | **462** | **2506** |

## Recovery changes

- **CHEXA8 NEWSOLID:** ordinary corner tensors and recovery operators are extrapolated from the eight valid enhanced states at 2x2x2 Gauss points. Center/interior requests retain direct enhanced recovery. Output no longer inverts the distorted corner mapping. Condensed enhancement and thermal contributions are included in the sampled states.
- **CPYRAM5/13:** populate six recovery locations: center, four basal corners and apex. The center retains direct recovery at collapsed coordinates (0,0,0.25). Corners use tensor-product extrapolation from the existing eight integration points; the apex uses the averaged collapsed direction (a=b=0, z=1), avoiding a singular apex Jacobian. Actual output-point shape values provide thermal stress contributions.
- **CPYRAM13 OES/OSTR:** six populated locations now support the standard center-plus-five-corner record, so the center-only suppression guards are removed. This is ordinary element strain output, **not GPSTRAIN**.
- **Batch writing:** flush successfully recovered solid records after processing the type, independently of the original requested row count. A failed element no longer hides the other valid elements in its batch.
- **Error reporting:** report recovery failures using the failing TYPE/EID immediately. The error summary no longer attributes earlier CHEXA8 errors to the final PYRAM13 element visited.
- **Thermal GPSTRESS:** request thermal data loading before interior recovery in thermal subcases. With every EMG option disabled, its early return previously left temperatures unloaded. A free uniform-expansion regression now checks zero thermal stress for both ordinary output and GPSTRESS.

No stiffness formulation, mesh, prescribed displacement, global DOF or coordinate was changed. The existing kernel comparisons still pass for stiffness, mass and thermal loading. The original 5-001 displacements remain bitwise identical. Its GPSTRESS direct and principal payloads also remain bitwise identical, retaining eight volumes and 204 grid rows.

## Geometry remains deliberately distorted

For CHEXA8 500002, local corner 8 (GRID 500005) has detJ=-0.001581269375; the minimum sampled 2x2x2 and 3x3x3 integration determinants are +0.006335311428 and +0.002483659322. These values are unchanged. Extrapolated output does not make the mapping orientation-preserving everywhere or excuse an invalid integration point. The LEGACY direct-corner path still diagnoses this corner; its valid neighboring elements are now retained in output. Mesh improvement remains a separate production-model decision.

## Verification

Build: `cmake --build build -j 6`.

1. `python -m unittest discover -s tests/gpstress -p test_solid_element_recovery.py -v`: five tests pass, covering the supplied full patch, six independent constant-strain states and rigid rotation, independent nonconstant pyramid strain/temperature fields, uniform thermal expansion, and a real LEGACY recovery failure with successful neighboring records preserved. Constant-field stress/strain checks use the established 1e-7 relative criterion; solved displacement absolute errors are below 1e-8 in these tests.
2. `python -m unittest discover -s tests/gpstress -p test_gpstress.py -v`: all 11 existing volume regressions pass. The old test expectation that the attached deck must emit a corner-Jacobian error is updated to reflect corrected recovery.
3. `python tests/solids/check_kernels.py`: Python/Fortran kernel comparisons pass for regular, warped and rotated geometries, both NEWSOLID and LEGACY selections, including matrix rank, mass and thermal checks.
4. `python tests/solids/check_newsolid_default.py`: omitted selectors and explicit NEWSOLID/ABE produce bitwise-identical results for displacement, SPC forces and all eight solid groups' element stress/strain over three subcases.
5. Independent OP2 readback against the fresh Nastran run verifies all 462 IDs and the same center/corner keys for stress and strain. Parameter-table skipping and a temporary in-memory element-name map accommodate the installed reader; binary result records are not patched.

The full matched-model audit and raw output files remain in the local `battle/solid_sap_comparison/element_stress_5001_20261018/` directory, with final results under `corrected/`. Only numerical summary data accompany this note.

## Scope

These results establish recovery/output correctness for the stated cases and field checks. They do not establish universal distorted-mesh quality, full finite-strain behavior, GPSTRAIN, or all-element/all-analysis numerical parity with Nastran. The 17 October benchmark curves remain historical results; no new beam/dome analyses are claimed here.
