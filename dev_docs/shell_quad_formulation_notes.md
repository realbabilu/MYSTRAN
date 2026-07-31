# Shell QUAD formulation notes: DKMQ20, DKMQ24, and MITC4+

Context: these notes document the current comparison between the Python prototype elements under `D:\18a\bending_only\Shell\gemini2\shit` and the MYSTRAN Fortran ports under `D:\18a\MYSTRAN\Source\EMG\EMG4`.

The main benchmark referenced here is Problem 2-004 MacNeal twisted beam, thick case, after rerunning the stale MYSTRAN F06 files for `MY-DKMQ20`.

## Current MacNeal 2-004 thick results

Reference target used by the script:

- `Uy_ref = 0.005429`
- `Uz_ref = 0.001749`

Last available mesh result:

| Label | Mesh N | Uy | Uz | Notes |
|---|---:|---:|---:|---|
| Python `DKMQ20-AU` | 24 | `5.4019394e-03` | `1.7053021e-03` | Python `DKMQ20_ShellElement_RHR_6dof_AU_K6ROT.py` |
| MYSTRAN `MY-DKMQ20` | 24 | `5.4019393e-03` | `1.7053023e-03` | `CQUAD4_DKMQ20_RHR.f90`, `PARAM,QUAD4TYP,DKMQ20` |
| MYSTRAN `MY-DKMQ24` | 24 | `5.4019393e-03` | `1.7053023e-03` | Standard CQUADR branch currently equivalent to AU-style kernel in this benchmark |
| Python original `DKMQ24` | 24 | `7.5919117e-03` | `2.3553501e-03` | Python `DKMQ24_ShellElement_RHR.py`, non-AU/original path |
| Python `MITC4+ 6DOF` | 24 | `6.8564643e-03` | `2.3885624e-03` | Python `MITC4p_ShellElement_6dof.py` |
| MYSTRAN `MY-MITC4P` | 24 | `5.4071923e-03` | `1.7076490e-03` | MYSTRAN CQUAD4 MITC4+ |

The cleaned `MY-DKMQ20` convergence sequence is:

| N | MY-DKMQ20 Uy | MY-DKMQ20 Uz |
|---:|---:|---:|
| 2 | `4.7161867e-03` | `7.6440647e-04` |
| 4 | `5.0876373e-03` | `1.2303757e-03` |
| 6 | `5.2309257e-03` | `1.4046257e-03` |
| 8 | `5.2981517e-03` | `1.5043113e-03` |
| 12 | `5.3585167e-03` | `1.6110387e-03` |
| 18 | `5.3904300e-03` | `1.6775640e-03` |
| 24 | `5.4019393e-03` | `1.7053023e-03` |

Earlier non-monotonic `MY-DKMQ20` behavior at N=8 was stale F06 output from before the AU/K6ROT kernel replacement. After rerun, the N=8 value matches the Python AU value.

## 1. DKMQ20 MYSTRAN vs Python DKMQ20 5DOF non-AU

These are not the same formulation.

The old/simple Python `DKMQ20_ShellElement_RHR.py` should be treated as the 5-DOF/non-AU baseline. It does not include the current AU projection machinery that is now used in MYSTRAN `CQUAD4_DKMQ20_RHR.f90`.

Current MYSTRAN `CQUAD4_DKMQ20_RHR.f90` is closer to:

```text
DKMQ20_ShellElement_RHR_6dof_AU_K6ROT.py
```

The important differences are:

- MYSTRAN DKMQ20 uses a 24x24 shell slot, i.e. 6 DOF per node at the solver interface.
- The physical drilling stiffness is not part of the DKMQ20 bending/membrane/shear kernel; the sixth DOF is handled by the external MYSTRAN K6ROT stabilization path.
- It uses the AU/ADELTA projection path:

```text
AU = BUILD_AU(XYZ, NORMALS)
ADELTA = BUILD_ADELTA(XYZ, thickness)
AINV_AU = inv(ADELTA) * AU
```

- `Bb` has the additional translational curvature contribution caused by nodal normal variation.
- `Bs` is built through the AU shear projection path, not the old 5-DOF shear expression.

Result: `MY-DKMQ20` should be compared against Python `DKMQ20-AU`, not against the old Python `DKMQ20` 5-DOF/non-AU file.

MacNeal confirms this:

```text
Python DKMQ20-AU N=24 : Uy=5.4019394e-03, Uz=1.7053021e-03
MYSTRAN MY-DKMQ20 N=24: Uy=5.4019393e-03, Uz=1.7053023e-03
```

The difference is only F06/parser rounding.

## 2. MYSTRAN DKMQ24/AU vs original Python DKMQ24

These are also not the same formulation.

The current MYSTRAN DKMQ24-like result in MacNeal is numerically identical to the AU-family DKMQ20 result:

```text
MY-DKMQ24 N=24 : Uy=5.4019393e-03, Uz=1.7053023e-03
MY-DKMQ20 N=24 : Uy=5.4019393e-03, Uz=1.7053023e-03
DKMQ20-AU N=24 : Uy=5.4019394e-03, Uz=1.7053021e-03
```

The original Python `DKMQ24_ShellElement_RHR.py` gives:

```text
Python DKMQ24 original N=24: Uy=7.5919117e-03, Uz=2.3553501e-03
```

That original Python DKMQ24 path is softer on MacNeal 2-004 than the MYSTRAN/AU family. It should not be used as the equivalence target for the current MYSTRAN CQUADR/DKMQ24-style result unless the MYSTRAN branch is intentionally reverted to the old non-AU formulation.

Practical interpretation:

- `DKMQ24_ShellElement_RHR.py` = original Python DKMQ24 reference path.
- `DKMQ24_MystranCQUADR_ShellElement_RHR.py` / MYSTRAN `MY-DKMQ24` = MYSTRAN-like AU-compatible path.
- If the purpose is to reproduce MYSTRAN CQUADR, compare against `DKMQ24_MystranCQUADR_ShellElement_RHR.py`, not the original Python `DKMQ24_ShellElement_RHR.py`.

This is why `MY-DKMQ24` and Python original `DKMQ24` differ substantially in MacNeal even when both are named DKMQ24.

## 3. Why MYSTRAN MITC4+ is better on MacNeal 2-004 than Python MITC4+ 5/6DOF

Current MacNeal N=24:

```text
MY-MITC4P         Uy=5.4071923e-03, Uz=1.7076490e-03
Python MITC4+6    Uy=6.8564643e-03, Uz=2.3885624e-03
Reference         Uy=5.4290000e-03, Uz=1.7490000e-03
```

MYSTRAN MITC4+ is closer to the benchmark reference in this problem. The main reason is not just "5 DOF vs 6 DOF"; the implementation details differ.

Key MYSTRAN MITC4+ traits:

- MYSTRAN uses the full shell infrastructure: 6-DOF CQUAD4 assembly slot, element coordinate transforms, and the established K6ROT handling for the drilling DOF.
- The MYSTRAN MITC4+ path uses the Ko-Lee-Bathe MITC4+ covariant strain direct interpolation for membrane/bending terms.
- The MYSTRAN shear terms go through the existing MITC tying-point/covariant-to-local transform path in `MITC4_B.f90`.
- The MYSTRAN implementation applies the standard MYSTRAN element orientation and shell-local basis machinery consistently through the CQUAD4 element pipeline.

By contrast, the Python MITC4+ prototypes are useful for isolated formulation experiments but do not necessarily reproduce MYSTRAN's complete production path:

- depending on the file, the Python prototype may be 5-DOF or manually expanded to 6-DOF;
- drilling handling may be fictitious or absent rather than MYSTRAN's K6ROT-style stabilization;
- local/global transform conventions are not guaranteed to match MYSTRAN's production element transform exactly;
- the shear/director convention can be materially different in twisted geometry.

MacNeal twisted beam is sensitive to director interpolation, local basis handling, and shear/bending coupling. Therefore a seemingly small difference in transform or drilling treatment can produce a large response difference.

Current evidence suggests:

```text
MYSTRAN MITC4+ is not simply "Python MITC4+ 5DOF plus one drilling DOF".
```

It is a more complete MYSTRAN shell pipeline implementation. That is why it can be better on MacNeal 2-004 than the standalone Python MITC4+ prototype.

## Practical comparison map

Use this mapping for future head-to-head tests:

| MYSTRAN label | PARAM selector | Fortran kernel | Correct Python comparison |
|---|---|---|---|
| `MY-DKMQ20` | `PARAM,QUAD4TYP,DKMQ20` | `CQUAD4_DKMQ20_RHR.f90` | `DKMQ20_ShellElement_RHR_6dof_AU_K6ROT.py` |
| `MY-DKMQ24` | default CQUADR DKMQ24 path | CQUADR DKMQ24/MYSTRAN path | `DKMQ24_MystranCQUADR_ShellElement_RHR.py` |
| Python original DKMQ24 | n/a | n/a | `DKMQ24_ShellElement_RHR.py`; do not assume it equals MYSTRAN DKMQ24 |
| `MY-MITC4P` | `PARAM,QUAD4TYP,MITC4+` or default CQUAD4 MITC4+ path | MYSTRAN MITC4/MITC4+ files | `MITC4p_MystranCQUAD4_ShellElement_6dof.py` if comparing MYSTRAN-like behavior; not the older plain MITC4+ prototype |

## Recommendation

For validation plots, avoid label names that hide formulation changes.

Recommended labels:

- `DKMQ20-old5` for old Python non-AU DKMQ20.
- `DKMQ20-AU` for Python AU/K6ROT DKMQ20.
- `MY-DKMQ20` for MYSTRAN CQUAD4 DKMQ20.
- `DKMQ24-orig` for original Python DKMQ24.
- `DKMQ24-MY` or `MY-DKMQ24` for the MYSTRAN-like DKMQ24/AU path.
- `MITC4+PY` for standalone Python MITC4+.
- `MITC4+MY` or `MY-MITC4P` for MYSTRAN CQUAD4 MITC4+.

This avoids treating different kernels as failed ports when they are actually different formulations.
