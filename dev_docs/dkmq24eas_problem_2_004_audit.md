# DKMQ24EAS Problem 2-004 Audit

## Scope

Problem 2-004 MacNeal twisted beam showed `MY-DKMQ24EA` tracking the poor
unfixed Python `DKMQ24` behavior instead of the current Python
`DKMQ24_EAS4_ShellElement_RHR_standalone.py` behavior.

The target behavior is the standalone DKMQ24EAS implementation:

- flat elements use the EAS4 membrane enhancement
- warped elements disable the EAS4 correction and keep the native DKMQ24
  membrane path
- this is not intended to make DKMQ24EAS follow DKMQ24AU as its formulation;
  it follows the Python standalone's flatness guard

## MYSTRAN Change

File changed:

- `Source/EMG/EMG4/CQUADR_DKMQ24N.f90`

Change made:

- `DKMQ20_MODE` is no longer forced on for every `DKM24EA` element.
- `EAS4_ACTIVE` is now the explicit flat-element gate:
  `DKM24EA_MODE .AND. IS_FLAT_QUAD(XYZ)`.
- Only when `EAS4_ACTIVE` is true does `DKMQ20_MODE` switch on for the EAS4
  branch.
- Warped `DKM24EA` elements therefore retain the native DKMQ24 membrane path,
  matching the Python standalone guard.

## Verification

Build:

- `make -j4` in `D:\18a\MYSTRAN\build`
- result: passed, produced `D:\18a\MYSTRAN\Binaries\mystran.exe`

Problem 2-004 rerun:

- decks: `D:\18a\python\linear\working_mystran_2_004\prob_2_004_nx*_cquadr_dkm24ea.dat`
- all `nx = 2, 4, 8, 12, 16, 24` completed successfully

Final comparison against Python `DKMQ24EAS` standalone:

| nx | qty | Python DKMQ24EAS | MYSTRAN DKM24EA | error vs Python |
| --- | --- | ---: | ---: | ---: |
| 2 | UY | 0.00547601671153 | 0.00471618666667 | -13.8756% |
| 2 | UZ | 0.000892452540881 | 0.000764406466667 | -14.3477% |
| 4 | UY | 0.00530238322771 | 0.00508763733333 | -4.04999% |
| 4 | UZ | 0.00128810215724 | 0.00123037566667 | -4.48151% |
| 8 | UY | 0.00534996803711 | 0.00529815166667 | -0.968536% |
| 8 | UZ | 0.00152216534988 | 0.00150431133333 | -1.17294% |
| 12 | UY | 0.00537713281261 | 0.00535851666667 | -0.34621% |
| 12 | UZ | 0.00161921611514 | 0.00161103866667 | -0.505025% |
| 16 | UY | 0.00539102875629 | 0.00538342 | -0.141137% |
| 16 | UZ | 0.00166649230711 | 0.00166198633333 | -0.270387% |
| 24 | UY | 0.00540363165485 | 0.00540193933333 | -0.0313182% |
| 24 | UZ | 0.00170717044684 | 0.00170530233333 | -0.109427% |

Flat problem 2-003 spot check was also rerun for `DKM24EA` at
`nx = 4, 8, 12, 16`; UY stayed matched to Python while UZ stayed about
1.5% high, consistent with the previous 2-003 audit.

## MITC3+ 2-004 Notes

MITC3+ was investigated but no change was kept.

Tested and reverted experiments in `Source/EMG/EMG4/TPLT_MITC3P.f90`:

- switching MITC3+ bending/shear quadrature from current 7-point integration
  to the Python 3-point rule
- bypassing the extra post-condensation rotational sign/transform layer
- replacing the current MYSTRAN ANS shear construction with the direct Python
  `j_mat @ N_gamma @ B_mid` shear operator
- pairing the direct Python shear operator with direct `SHELL_T`

Results:

- 3-point quadrature made the 2-004 MITC3+ result worse.
- The direct Python shear substitutions collapsed the displacements to roughly
  5% or less of the Python target, so the Python shear operator is not a safe
  local substitution into MYSTRAN's current MITC3+ kernel.
- The original MITC3+ source was restored. MITC3+ still needs deeper frame and
  shear-operator derivation work before applying a fix.

