! --- response_spectrum_mystran_add begin --- !
# User Guide (Response Spectrum Writer)

## Inputs
- Geometry NEU template (importable in FEMAP)
- `SOL103` modal NEU + F06
- `SOL112` RSX NEU
- `SOL112` RSY NEU

## Print/Post Syntax (Nastran)
- `PRINT` -> F06 text
- `PLOT` -> OP2 binary
- `PRINT,PLOT` -> keduanya
- `NEU` -> post-step (writer), bukan syntax langsung di deck

Referensi profil siap pakai:
- `docs/femap_output_profiles.md`

## Solver Recommendation
- `DENSE`: quick baseline/smoke
- `SUBSPACE`: larger practical models
- `FEAST`: when FEAST path is available and validated in your environment

## Run Command
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\scripts\femap_rs_neu_writer.py \
  --geometry E:\mystran17\dynamics\verification\nxnastra\model-028-1-024-full-real-vectors.neu \
  --modal-neu E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.NEU \
  --modal-f06 E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.F06 \
  --rsx-neu E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.NEU \
  --rsy-neu E:\mystran17\dynamics\verification\Example-1-024_sol112_rsy.NEU \
  --out E:\mystran17\dynamics\verification\nxnastra\model-033-1-024-rs-writer.neu
```

## One-Command FEMAP Mode
```powershell
& E:\mystran17\codex_mod\response_spectrum_mystran_add\scripts\run_femap_mode.ps1 `
  -Bdf E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.bdf `
  -Mode FEMAP
```

## Ready-to-Run Deck Templates (Example 1-024)
- Fast (`PLOT`-heavy):
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_srss_femap_fast.bdf`
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_cqc_femap_fast.bdf`
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_rsy_femap_fast.bdf`
- Check (`PRINT,PLOT`):
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_srss_femap_check.bdf`
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_cqc_femap_check.bdf`
  - `E:\mystran17\dynamics\verification\Example-1-024_sol112_rsy_femap_check.bdf`

Dengan generate NEU:
```powershell
& E:\mystran17\codex_mod\response_spectrum_mystran_add\scripts\run_femap_mode.ps1 `
  -Bdf E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.bdf `
  -Mode FEMAP `
  -WriteNeu `
  -GeometryNeu E:\mystran17\dynamics\verification\nxnastra\model-028-1-024-full-real-vectors.neu `
  -ModalNeu E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.NEU `
  -ModalF06 E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.F06 `
  -RsxNeu E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.NEU `
  -RsyNeu E:\mystran17\dynamics\verification\Example-1-024_sol112_rsy.NEU `
  -OutNeu E:\mystran17\dynamics\verification\nxnastra\model-femap-mode-output.neu
```

## Output Sets
- Modal 1..4
- RS X Direction
- RS Y Direction
- RSX+RSY, RSX-RSY, -RSX-RSY, -RSX+RSY
- 1.0RSX+0.3RSY, 0.3RSX+1.0RSY, and sign variants

## Custom Combo (No Code Edit)
Gunakan `--combo "NAME:a:b"` (repeatable), dengan rumus `a*RSX + b*RSY`.

Contoh:
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\scripts\femap_rs_neu_writer.py \
  --geometry E:\mystran17\dynamics\verification\nxnastra\model-028-1-024-full-real-vectors.neu \
  --modal-neu E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.NEU \
  --modal-f06 E:\mystran17\dynamics\verification\Example-1-024_sol103_modal.F06 \
  --rsx-neu E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.NEU \
  --rsy-neu E:\mystran17\dynamics\verification\Example-1-024_sol112_rsy.NEU \
  --no-default-combos \
  --combo "100X+30Y:1.0:0.3" \
  --combo "30X+100Y:0.3:1.0" \
  --combo "100X-30Y:1.0:-0.3" \
  --combo "-30X+100Y:-0.3:1.0" \
  --out E:\mystran17\dynamics\verification\nxnastra\model-034-1-024-custom-combo.neu
```

## Validate Quickly
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\tools\neutral_parser.py \
  E:\mystran17\dynamics\verification\nxnastra\model-033-1-024-rs-writer.neu
```

## Compare SOL111 vs SOL112 (NEU Value Check)
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\tools\compare_neu_sets.py \
  E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.NEU \
  E:\mystran17\dynamics\verification\Example-1-024_sol111_from112_compat_clean.NEU \
  --abs-tol 0 --rel-tol 0 \
  --json-out E:\mystran17\codex_mod\response_spectrum_mystran_add\validation\sol111_vs_sol112_example_1_024.json
```

## SOL111 Compatibility Run (from SOL112-style deck)
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\tools\sol112_to_sol111_compat.py \
  E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.bdf \
  E:\mystran17\dynamics\verification\Example-1-024_sol111_from112_compat_clean.bdf \
  --clean-sol111
```
```powershell
& E:\mystran17\mystran\Binaries\mystran.exe E:\mystran17\dynamics\verification\Example-1-024_sol111_from112_compat_clean.bdf
```

! --- response_spectrum_mystran_add end --- !
