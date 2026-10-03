# T6 final Python alignment — 2026-10-03

## Scope and reference classes

T6 only. Existing unrelated changes in the MYSTRAN working tree were preserved.
Q8 formulation and its Python registry were not upgraded in this work.

| Family | Final Python reference |
|---|---|
| SIMOT6 | `Simo1993_Tri6_ShellElement_v2.py` |
| MITC6 | `MITC6_Tri_v4.py` |
| MH6T | `MacNeal_MH6T_Tri_v3.py` |
| REZAIEE | `Rezaiee2017_Tri6_v3.py` |

## Formula alignment

- Coherent geometric directors and pointwise covariant-to-local mapping follow final Python defaults. Explicit SNORM cards remain supported; automatic averaged normals no longer replace geometric defaults.
- SIMO/MITC drilling rotations follow the point normal. MH6T fixed frame follows the centroid geometry, including flat elements; its drilling uses the physical metric projection.
- MITC6/REZAIEE bending includes director derivatives with the final cross-product/sign convention.
- REZAIEE membrane and shear tying locations, affine reconstruction, shear centroid and unnormalized interpolated shear director now follow v3.

## Recovery alignment

- All four kernels recover directly at CENTER and all six nodes in order: `(1/3,1/3)`, `(0,0)`, `(1,0)`, `(0,1)`, `(1/2,0)`, `(1/2,1/2)`, `(0,1/2)`. Stiffness quadrature remains unchanged.
- MYSTRAN uses membrane minus z times recovered curvature and a negative bending-force conversion. Therefore recovery BE2 is minus Python Bb; stiffness still uses Python Bb. Fiber stresses and moments now match Python signs.
- CENTER-only uses one point; CORNER/GPSTRESS internally reserve seven. F06 CORNER includes the three midside nodes.
- Engineering force recovery copies all T6 samples; subsequent nodal rows are no longer left uninitialized.
- GPSTRESS includes all six element nodes and uses the native signed centroid basis to transform local stresses to the requested surface basis.
- `compare_f06_vs_python.py` and `dump_python_debug.py` use final T6 references and their natural-coordinate node order. Historical aliases such as `MHT6v3` resolve to MH6T.
- Dump labels use the model element ID (`E51`), and each dump identifies the actual reference module/class.

## Validation

MYSTRAN built successfully with `C:/gcc/bin/make.exe -j4 mystran` in `MYSTRAN/build`; executable is `MYSTRAN/Binaries/mystran.exe`.

| Check | Result | Artifact |
|---|---|---|
| Displacement 2-001 through 2-004, four T6 families, mesh 8/24 for 2-002/3/4 | 88/88 pass | `test/t6_upgrade/displacement_refined/results.json` |
| Displacement regression after final recovery build, mesh 2/4 | 88/88 pass | `test/t6_upgrade/recovery_final/results.json` |
| 2-001 STEP A/B/C, membrane + bending, four families | 8/8 pass | `test/t6_upgrade/stress_recovery/results.json` |
| CENTER/CORNER/hidden GPSTRESS selections separately | 32/32 pass | `test/t6_upgrade/output_selection/results.json` |
| Final dump reference, ID, nodal interpolation, BM/BB/BS presence | All pass | `test/t6_upgrade/debug_validation.txt` |

Each stress case checks both fibers at 10 centers, 60 element-node samples, and 25 global grid averages, including midsides. Bending also checks 70 local moment samples (center plus six nodes per element). Missing outputs fail the gate. Maximum relative discrepancy is approximately 2.86e-6 for stress and 3.78e-7 for moments, consistent with F06 formatting.

Displacement compares all six DOFs at every node using absolute + component-scaled relative tolerances: `1e-10 + 1e-4*scale`. Stress/moment tolerance is `1e-4`. Models exported for displacement retain the same geometry, loads, prescribed displacements, material and final class as Python.

Human-readable comparator reports: `test/t6_upgrade/stress_recovery/{SIMOT6,MITC6,MH6T,REZAIEE}.compare.txt`. Eight element-51 dumps with attributes: `test/t6_upgrade/debug_<FAMILY>_<membrane|bending>.txt`.

This establishes displacement matching on 2-001/2/3/4 and F06 stress/moment A/B/C matching on the 2-001 patch. Curved/twisted stress recovery, strain output, thermal recovery, and OP2 quadratic table layouts have not been validated by these gates. Do not infer their validation from the displacement results.

## Reproducible checks

Run from `C:/PROJECTAI/18a`:

```powershell
python test/t6_displacement_gate.py --phase recheck --mesh 2 4
python test/t6_stress_gate.py
python test/t6_output_selection_gate.py
python test/verify_final_t6_debug.py
```

Source backups before formula/recovery changes are in `test/t6_upgrade/source_before_formula` and `source_before_recovery`. The mutation scripts `upgrade_t6_formula.py`, `refine_t6_formula.py`, `upgrade_t6_recovery.py`, and `finalize_t6_recovery.py` are one-time migration records, not rerunnable validation tools.
