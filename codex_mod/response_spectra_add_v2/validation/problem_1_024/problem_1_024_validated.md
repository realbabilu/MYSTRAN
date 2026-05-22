# Problem 1-024 Validated

Date:
- `2026-05-22`

Scope:
- MYSTRAN validation package for `SAP2000 Problem 1-024`
- package target: `codex_mod\1-24validated`

Validated artifacts:
- [`problem_1_024_mystran_dense_modal_beam_swap.dat`](D:/user/doc/New%20project/1-24validated/problem_1_024_mystran_dense_modal_beam_swap.dat)
- [`problem_1_024_static_unit_load.dat`](D:/user/doc/New%20project/1-24validated/problem_1_024_static_unit_load.dat)

Executable used:
- `D:\mystran2\MYSTRANSolver-18.0.0\Binaries\mystran.exe`

Validation basis:
- `Problem 1-024` is considered validated in MYSTRAN after the `RBE2` and `SPC1` continuation fixes
- modal and static checks both match the established OpenSees/reference baseline closely

Reference modal periods:
- `T1 = 0.227100 s`
- `T2 = 0.215600 s`
- `T3 = 0.073300 s`
- `T4 = 0.072000 s`

Validated MYSTRAN modal result:
- F06 frequencies from [`problem_1_024_mystran_dense_modal_beam_swap.F06`](D:/user/doc/New%20project/problem_1_024_mystran_dense_modal_beam_swap.F06):
  - `f1 = 4.404087 Hz`
  - `f2 = 4.637502 Hz`
  - `f3 = 13.634160 Hz`
  - `f4 = 13.887910 Hz`
- corresponding periods:
  - `T1 = 0.227062 s`
  - `T2 = 0.215633 s`
  - `T3 = 0.073345 s`
  - `T4 = 0.072005 s`

Absolute modal deltas:
- `|dT1| = 0.000038 s`
- `|dT2| = 0.000033 s`
- `|dT3| = 0.000045 s`
- `|dT4| = 0.000005 s`

Validated static check:
- from [`problem_1_024_static_unit_load.F06`](D:/user/doc/New%20project/problem_1_024_static_unit_load/problem_1_024_static_unit_load.F06)
- node `29`:
  - OpenSees `1.714145015e-04`
  - MYSTRAN `1.714145000e-04`
- node `19`:
  - OpenSees `1.645314940e-04`
  - MYSTRAN `1.645315000e-04`

Deck notes:
- modal deck:
  - `SOL 103`
  - `EIGRL,1,,,4`
  - continuation `,DENSE,,64`
  - validated configuration uses `CBAR` plus beam inertia swap interpretation
- static deck:
  - `SOL 101`
  - unit load at node `29`
  - same diaphragm/master-mass semantics as the validated benchmark model

Status:
- `Problem 1-024 modal = validated`
- `Problem 1-024 static unit-load check = validated`
- this package freezes the validated MYSTRAN deck set and the validation record
