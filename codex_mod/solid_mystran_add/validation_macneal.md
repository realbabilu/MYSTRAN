# MacNeal Solid Validation

Validation used MacNeal Static-37 style solid cantilever cases. Reports and scripts are copied in `files/codex_mod/solids/`.

## MYSTRAN Legacy vs NEWSOLID

| Element | Mesh | Case | Legacy/ref | NEWSOLID/ref |
| --- | --- | --- | ---: | ---: |
| `CHEXA8` | n12 | in-plane | `0.931` | `0.981` |
| `CHEXA8` | n12 | out-of-plane | `0.938` | `0.977` |
| `CPENTA6` | n12 | in-plane | `0.212` | `0.352` |
| `CPENTA6` | n12 | out-of-plane | `0.255` | `0.403` |
| `CPYRA5` | n12 | in-plane | `0.160` | `0.249` |
| `CPYRA5` | n12 | out-of-plane | `0.268` | `0.385` |
| `CTETRA4` | n12 | in-plane | `0.088` | `0.628` |
| `CTETRA4` | n12 | out-of-plane | `0.145` | `0.748` |
| `CTETRA10` | n8 | in-plane | `0.712` | `0.712` |
| `CTETRA10` | n8 | out-of-plane | `0.733` | `0.733` |

## Python vs MYSTRAN

- CPENTA6, CPYRA5, CPYRA14, CTETRA4 smooth, and quadratic guard checks match the local Python targets to near roundoff in the dedicated validation scripts.
- CHEXA8 `NEWSOLID` is close on the MacNeal gate; the local Python standard/legacy CHEXA comparison is not treated as the oracle for MYSTRAN legacy.
- CPYRA buckling, thermal, and differential stiffness smoke checks completed without fatal errors.

## Commands Used

```powershell
cmake --build build_mingw --target mystran --config Release -- -j8
python codex_mod/solids/validate_cpyra_mystran.py
python codex_mod/solids/validate_cpyra_buckling_smoke.py
python codex_mod/solids/validate_cpyra_thermal_smoke.py
python codex_mod/solids/validate_quadratic_mystran.py
python codex_mod/solids/compare_macneal_hexa8_legacy_newsolid.py
python codex_mod/solids/compare_macneal_cpenta6_legacy_newsolid.py
python codex_mod/solids/compare_macneal_cpyra5_legacy_newsolid.py
python codex_mod/solids/compare_macneal_ctetra_legacy_newsolid.py
python codex_mod/solids/compare_macneal_python_vs_mystran.py
```
