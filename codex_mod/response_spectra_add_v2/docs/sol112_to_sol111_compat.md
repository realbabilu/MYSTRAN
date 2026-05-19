! --- response_spectrum_mystran_add begin --- !
# SOL112 -> SOL111 Compatibility Layer (External)

Tool:
- `tools/sol112_to_sol111_compat.py`

Purpose:
- Accept input deck style currently used for `SOL112` response spectrum setups.
- Convert to `SOL111` deck so direct/modal FRF path can process it.

What it changes:
1. `SOL 112` -> `SOL 111`
2. Ensures each `SUBCASE` has `FREQ = <sid>` (auto from first `FREQ1`, fallback `40`)
3. Updates title text containing `SOL112` -> `SOL111` when found
4. Optional clean mode: comments `ACCELERATION(...)` request in case control for cleaner SOL111 logs
5. Keeps bulk cards unchanged (DAREA/RLOAD1/TABLED1/DLOAD etc.)
6. Adds trace marker comment in bulk data

Usage:
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\tools\sol112_to_sol111_compat.py \
  E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.bdf \
  E:\mystran17\dynamics\verification\Example-1-024_sol111_from112_compat.bdf
```

Usage (clean SOL111 output):
```powershell
& 'C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' \
  E:\mystran17\codex_mod\response_spectrum_mystran_add\tools\sol112_to_sol111_compat.py \
  E:\mystran17\dynamics\verification\Example-1-024_sol112_srss.bdf \
  E:\mystran17\dynamics\verification\Example-1-024_sol111_from112_compat_clean.bdf \
  --clean-sol111
```

Then run MYSTRAN normally on the generated SOL111 deck.

! --- response_spectrum_mystran_add end --- !
