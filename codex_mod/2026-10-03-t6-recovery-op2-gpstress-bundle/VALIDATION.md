# Validasi dan reproducibility

## Evidence yang telah direkam

| Pemeriksaan | Hasil | Evidence dalam bundle |
|---|---|---|
| Displacement 2-001–2-004, empat T6, mesh 8/24 untuk 2-002/3/4 | 88/88 | `evidence/t6_upgrade/displacement_refined/results.json` |
| Regression displacement setelah recovery, mesh 2/4 | 88/88 | `evidence/t6_upgrade/recovery_final/results.json` |
| STEP A/B/C patch 2-001, membrane + bending | 8/8 | `evidence/t6_upgrade/stress_recovery/results.json` |
| Selector CENTER/CORNER/hidden GPSTRESS | 32/32 | `evidence/t6_upgrade/output_selection/results.json` |
| Dump referensi, ID, shape nodal, BM/BB/BS | PASS | `evidence/t6_upgrade/debug_validation.txt` |
| OP2 OES/OEF/OGS empat T6 vs F06 | PASS | `evidence/t6_upgrade/op2_validation.txt` |
| Layout OES/OEF dan ID T6 vs Siemens 2512 | MATCH | `evidence/t6_upgrade/op2_reference/` |
| Mixed duel3a OP2 read + nilai T6 | PASS | `evidence/duel3a.pynastran.json` |
| Empty force header cleanup | PASS, nilai hasil tetap | `evidence/duel3a_force_header.stdout.txt`; pemeriksaan readback dilakukan pada working tree |

Displacement memakai absolute + component-scale tolerance `1e-10 + 1e-4*scale`; stress/moment gates memakai relative tolerance `1e-4`. Selisih maksimum patch stress sekitar `2.86e-6`, sesuai pembulatan F06. Missing output gagal gate. Bukti build dipilih pada folder `evidence/t6_upgrade/`.

Pada CENTER/CORNER patch, setiap kasus memeriksa 10 center, 60 node-elemen, 25 grid global dan kedua fiber; bending juga memeriksa 70 momen lokal. OP2 T6 mempunyai 80 stress rows/subcase, 40 force rows pada bending, serta 50 OGS fiber rows/surface/subcase.

Regression mesh 2/4 dilakukan sesudah formula/recovery final. Perubahan setelah itu berada di serialization dan pencetakan output; gates F06/OP2 serta readback deck mixed dijalankan kembali sesuai evidence. Tidak ada klaim seluruh daftar tes dijalankan ulang setelah setiap perubahan whitespace/header terakhir.

## Layout workspace untuk rerun

Skrip validasi menggunakan layout awal:

```text
workspace/
  MYSTRAN/     # clone + executable hasil build
  python/      # core, final references, benchmark/helper
  test/        # gates, decks, output sementara
```

Pada workspace asli layout ini sudah ada. Untuk clone lain yang bernama `MYSTRAN`, dari root repo:

```powershell
python codex_mod/2026-10-03-t6-recovery-op2-gpstress-bundle/tools/restore_validation.py --workspace C:/work/18a
```

Tool menyalin `validation/python` dan `validation/test`; existing file yang berbeda ditolak kecuali opsi `--overwrite` diberikan. `--evidence` ikut mengembalikan evidence ke folder test untuk pemeriksaan layout referensi. Paket tidak menyalin ulang archival Fortran ke Source; source produksi sudah merupakan bagian commit yang direncanakan.

## Build dan runtime

Build source dengan workflow CMake/Fortran/OpenBLAS yang tersedia di repo. Pada mesin validasi digunakan gfortran/make dari `C:/gcc/bin`, target `mystran`, executable `MYSTRAN/Binaries/mystran.exe`.

```powershell
C:/gcc/bin/make.exe -C MYSTRAN/build -j4 mystran
```

Perintah di atas dijalankan dari parent workspace setelah build directory dikonfigurasi. Runtime compiler dan OpenBLAS DLL harus tersedia melalui PATH; mereka tidak disalin ke bundle. File `validation/requirements-tested.txt` merekam versi Python packages terpasang saat packaging, bukan jaminan kompatibilitas setiap versi lain.

## Gates yang dapat dijalankan ulang

Dari parent workspace:

```powershell
python test/t6_displacement_gate.py --phase recheck --mesh 2 4
python test/t6_displacement_gate.py --phase recheck_refined --mesh 8 24
python test/t6_stress_gate.py
python test/t6_output_selection_gate.py
python test/verify_final_t6_debug.py
python test/read_t6_op2.py test/t6_upgrade/stress_recovery/MITC6.OP2
python test/check_t6_op2_reference.py
```

OP2 checks membutuhkan hasil OP2 yang dibentuk oleh stress gate. Layout reference comparison juga membutuhkan recorded JSON referensi pada `test/t6_upgrade/op2_reference`; gunakan restore `--evidence`, atau jalankan solver referensi kembali melalui deck yang disediakan. Siemens executable memerlukan instalasi/lisensi pengguna dan tidak masuk bundle.

Benchmark 2-006/full legacy convergence bukan gate final di atas. Modul historis `Simo1993_Q8_ShellElement_V2_MacNealPatched`, `Simo1993_Q8_ShellElement_v3`, dan `solver_f06_convergence` tidak tersedia di workspace ini. Snapshot legacy tidak dianggap telah memenuhi dependency tersebut.
