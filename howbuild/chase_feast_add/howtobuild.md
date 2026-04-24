# CHASE/FEAST Add - How To Build

## Commit Group
`chase_feast_add`

Related:
- `how_to_use.md` for run-time usage.
- `examples/` for ready-made FEAST/CHASE decks.

## CMakeLists.txt Snippet
Tambahkan snippet ini di `CMakeLists.txt`:

```cmake
# !--- CHASE and FEAST --- begin!
# FEAST/ChASE wrapper sources use preprocessor guards for optional native hooks.
set_source_files_properties(
  "${CMAKE_SOURCE_DIR}/LK4/EIG_LANCZOS_FEAST.f90"
  "${CMAKE_SOURCE_DIR}/LK4/EIG_LANCZOS_CHASE.f90"
  PROPERTIES COMPILE_OPTIONS "-cpp"
)
# !--- CHASE and FEAST --- end!
```

## Build Command

```powershell
cmake ..\mystran `
  -DMYSTRAN_USE_EXTERNAL_FEAST=ON `
  -DMYSTRAN_USE_EXTERNAL_CHASE=ON `
  -DMYSTRAN_FEAST_LIBRARY="...\libfeast.a" `
  -DMYSTRAN_CHASE_LIBRARY="...\libchase_f.a" `
  -DMYSTRAN_CHASE_EXTRA_LIBS="...\libchase_c.a;stdc++;gomp"
cmake --build . -j 8
```

Reference from `E:\mystran17\eigenbenchmark\build4.bat`:

```powershell
cmake -S E:\mystran17\mystran -B E:\mystran17\mystran\build-codex -G "MinGW Makefiles" `
  -DMYSTRAN_USE_EXTERNAL_CHASE=ON `
  -DMYSTRAN_CHASE_LIBRARY="C:/gcc/chase32/nompi/libchase_f.a" `
  -DMYSTRAN_CHASE_EXTRA_LIBS="C:/gcc/chase32/nompi/libchase_c.a;stdc++;gomp" `
  -DMYSTRAN_USE_EXTERNAL_FEAST=ON `
  -DMYSTRAN_FEAST_LIBRARY="C:/gcc/feast32/libfeast.a"
cmake --build E:\mystran17\mystran\build-codex -j 8
```

Notes:
- Current CMake wiring for CHASE requires `MYSTRAN_USE_EXTERNAL_CHASE=ON` + `MYSTRAN_CHASE_LIBRARY`.
- If CHASE backend is not linked, runtime uses ARPACK fallback with warning `4912`.
- FEAST include/library variables in this package are legacy placeholders from earlier notes; only library path is used in current wiring.

## Info
- v18: ARPACK sudah bisa jalan lewat SuperLU saat `SOLLIB=SPARSE`, jadi bukan LAPACK-only.
- Jalur ini terlihat di `ARPACK_LANCZOS_EIG` yang factor/solve `KMSM` via `SYM_MAT_DECOMP_SUPRLU` + `FBS_SUPRLU`.

## Roadmap Integrasi FEAST/CHASE
Referensi:
- `dev_docs/eigensolver_feast_chase_integration_plan.md`

Rekomendasi praktis:
- Pertahankan ARPACK+SuperLU sebagai default dulu (sudah stabil dan cepat di benchmark).
- Tambah `LANCMETH=FEAST/CHASE` sebagai opsi eksperimental bertahap.
- Implement lewat adapter di `LINK4` (bukan patch langsung di jalur ARPACK) supaya fallback ke ARPACK tetap aman.

## Modified Files (CHASE/FEAST add)
Dengan snippet `# !--- CHASE and FEAST --- begin!` s.d. `# !--- CHASE and FEAST --- end!`

- `mystran/Source/USE_IFs/LINK4_USE_IFs.f90`
- `mystran/Source/LK4/LINK4.f90`
- `mystran/Source/LK1/L1A-BD/BD_PARAM.f90`
- `mystran/Source/LK4/EIG_LANCZOS_CHASE.f90`
- `mystran/Source/LK4/EIG_LANCZOS_FEAST.f90`
- `mystran/Source/USE_IFs/EIG_LANCZOS_FEAST_USE_IFs.f90`
- `mystran/Source/USE_IFs/EIG_LANCZOS_CHASE_USE_IFs.f90`
- `mystran/Source/Modules/PARAMS.f90`
- `mystran/CMakeLists.txt`

## Supporting Files (CHASE/FEAST add)
- `mystran/Binaries/midas_eigen08_cbar_axis_debug.dat`
- `mystran/Binaries/midas_eigen08_cbeam_axis_debug.dat`
- `mystran/Binaries/midas_eigen08_cbeam_axis_debug_bernoulli.dat`
- `mystran/dev_docs/eigensolver_feast_chase_surrogate_validation.md`
- `dev_docs/v18_static_eigen_validation_report.md`
- `dev_docs/eigensolver_feast_chase_integration_plan.md`
