# Daftar file berubah dan snapshot pendukung

Baseline GitHub: `931dde4df56d9c39b43d884bda941f7fd27c6036`. Branch: `v18.00.a`.

## 38 file tracked yang berubah

| File di repo | Ringkasan |
|---|---|
| `Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUAD4_SIMO1989.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUAD8_MACQ8D.f90` | Existing local change: move RV/SV parameters to host scope; no new Q8 formula upgrade in this session. |
| `Source/EMG/EMG4/CQUADR_DKM24EA.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_DKMQ24.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_DKMQ24R.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_MBP1C0.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_Q4EASANS.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_Q4RS.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CQUADR_SIMO1993.f90` | Existing local changes to MELDOF-sized recovery buffers / bounded BE assignment; included in the initial-to-current delta. |
| `Source/EMG/EMG4/CTRIA6_MH6T.f90` | Align formulation, signed directors, metric/drilling and seven-point recovery with final Python; preserve explicit SNORM. |
| `Source/EMG/EMG4/CTRIA6_MITC6.f90` | Align formulation, signed directors, metric/drilling and seven-point recovery with final Python; preserve explicit SNORM. |
| `Source/EMG/EMG4/CTRIA6_REZAIEE.f90` | Align formulation, signed directors, metric/drilling and seven-point recovery with final Python; preserve explicit SNORM. |
| `Source/EMG/EMG4/CTRIA6_SIMO1993.f90` | Align formulation, signed directors, metric/drilling and seven-point recovery with final Python; preserve explicit SNORM. |
| `Source/LK1/L1A-CC/CHK_CC_CMD_DESCRIBERS.f90` | Track CENTER/CORNER request flags independently for stress/strain/force. |
| `Source/LK1/L1A/LOADB.f90` | Reserve recovery storage for seven T6 samples (MAX_STRESS_POINTS plus center slot). |
| `Source/LK9/L91/WRITE_ELEM_ENGR_FORCE.f90` | Native CTRIA6 type-75 OEF width 38 and corrected T6 CENTER stride/selectors. |
| `Source/LK9/L91/WRITE_ELEM_STRESSES.f90` | T6 F06 CENTER/six nodes; native OES width 70; signed basis GP conversion; actual OGS surface IDs and per-surface tables. |
| `Source/LK9/L92/ELEM_STRE_STRN_ARRAYS.f90` | Use basic-coordinate element displacement UEB for T6 recovery; remove temporary printing. |
| `Source/LK9/L92/OFP1.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP2.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP3.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP3_ELFE_1D.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP3_ELFE_2D.f90` | Recover/copy all seven T6 force samples; suppress empty surface force headers. |
| `Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90` | T6 direct seven-point output; per-element metadata; native signed centroid basis. |
| `Source/LK9/L92/OFP3_STRE_PCOMP.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP3_STRN_NO_PCOMP.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/OFP3_STRN_PCOMP.f90` | Preserve local MAXREQ*5 output-buffer capacity/initialization change. |
| `Source/LK9/L92/POLYNOM_FIT_STRE_STRN.f90` | Formatting/indentation delta; no behavior change claimed. |
| `Source/LK9/L92/SHELL_STRESS_OUTPUTS.f90` | Formatting/indentation delta; no behavior change claimed. |
| `Source/LK9/LINK9/ALLOCATE_LINK9_STUF.f90` | Preserve local MAXREQ*5 allocation of OGEL and shell transform arrays. |
| `Source/LK9/LINK9/MAXREQ_OGEL.f90` | T6 row counts 7 force / 14 stress; preserve existing local output capacity policy. |
| `Source/Modules/CC_OUTPUT_DESCRIBERS.f90` | Add CENTER/CORNER request flags. |
| `Source/Modules/GPSTRESS_SURFACE_UTILS.f90` | Include TRIA6/QUAD8 in shell patches and collect all six/eight nodes. |
| `Source/Modules/MODEL_STUF.f90` | Existing local TRIA6 NUM_SEi change 4 to 3; T6 stress/force output now explicitly reserves seven. |
| `Source/UTIL/OUTPUT2_WRITE_ELFORCE.f90` | Accept TRIA6 in OEF dispatch. |
| `Source/UTIL/OUTPUT2_WRITE_STRESS.f90` | Accept TRIA6 in OES dispatch; reopen OES after OGS closes it. |

## Python dengan bukti perubahan terhadap backup lokal

Folder sibling `python` dan `test` bukan bagian clone GitHub. Diff-nya tidak disamakan dengan diff GitHub.

- `battle_2_002_quadratic_clamped.py`
- `battle_2_003_quadratic_curved.py`
- `battle_2_004_twistedv2.py`
- `battle_2_006_scordelisv2.py`
- `compare_disp_v4c.py`
- `compare_f06_vs_python.py`
- `test_patch2001_q8_t6_v4c.py`
- `test_quadraticv2a.py`
- `test_result_2_001_t6_disp.py`

`dump_python_debug.py` juga diperbarui pada sesi ini, tetapi tidak memiliki backup pra-edit.
`duel3a_mht6_1.dat` diganti nama menjadi `duel3a_mh6t_1.dat`; helper perbandingan aktif disesuaikan.

## Snapshot tambahan

Referensi T6 final, core/adaptor/dependensi Python yang diimpor, skrip validasi, deck, laporan, dan bukti hasil ikut disalin agar dapat direview. Snapshot dependensi tidak berarti file tersebut diedit pada sesi ini.

Lokasi dan checksum setiap salinan tercantum dalam `FILE_MANIFEST.csv` dan `manifest.json`.
