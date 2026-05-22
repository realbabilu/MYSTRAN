**Problem 1-020 Validated**

Scope:
- SAP2000 `Example 1-020`
- frame response-spectrum analysis of a two-dimensional rigid frame
- MYSTRAN package frozen from a dense modal deck plus deterministic SRSS post-processing

Validated artifacts:
- [problem_1_020_modal.dat](D:/user/doc/New%20project/1-20validated/problem_1_020_modal.dat)
- [benchmark_1_020_modal.py](D:/user/doc/New%20project/1-20validated/benchmark_1_020_modal.py)
- [benchmark_1_020_validate.py](D:/user/doc/New%20project/1-20validated/benchmark_1_020_validate.py)

Executable used:
- `D:\mystran2\MYSTRANSolver-18.0.0\Binaries\mystran.exe`

Primary evidence files:
- [problem_1_020_modal.F06](D:/user/doc/New%20project/benchmark_1_020_modal/problem_1_020_modal.F06)
- [Problem 1-020.pdf](D:/mystran2/verification/beam/Problem%201-020.pdf)
- [Example 1-020.s2k](D:/mystran2/verification/beam/Example%201-020.s2k)

Reference results from SAP2000 / Chopra:
- `T1 = 1.5620 s`
- `T2 = 0.5868 s`
- `Ux@joint 2 = 7.576 in`
- `Ux@joint 3 = 18.84 in`
- `M33`:
  - `elem 1 @ joint 1 = 12636 k-in`
  - `elem 1 @ joint 2 = 6793 k-in`
  - `elem 2 @ joint 2 = 6023 k-in`
  - `elem 2 @ joint 3 = 5222 k-in`
  - `elem 5 @ joint 2 = 9810 k-in`
  - `elem 6 @ joint 5 = 9810 k-in`
  - `elem 7 @ joint 3 = 5222 k-in`
  - `elem 8 @ joint 6 = 5222 k-in`

MYSTRAN validated results:
- modal:
  - `T1 = 1.562131 s`
  - `T2 = 0.586817 s`
- response spectrum SRSS:
  - `Ux@joint 2 = 7.576084 in`
  - `Ux@joint 3 = 18.838518 in`
- `M33` SRSS:
  - `elem 1 @ joint 1 = 12635.826 k-in`
  - `elem 1 @ joint 2 = 6793.146 k-in`
  - `elem 2 @ joint 2 = 6023.038 k-in`
  - `elem 2 @ joint 3 = 5222.469 k-in`
  - `elem 5 @ joint 2 = 9810.580 k-in`
  - `elem 6 @ joint 5 = 9810.580 k-in`
  - `elem 7 @ joint 3 = 5222.619 k-in`
  - `elem 8 @ joint 6 = 5222.619 k-in`

Method note:
- The MYSTRAN deck is a dense `SOL 103` modal run with:
  - bending-only frame semantics
  - `CMASS2` masses attached only to `U1` at joints `7` and `8`
  - full 3D beam elements, with inactive SAP DOFs restrained
- SRSS response-spectrum quantities are reconstructed from:
  - modal periods
  - MYSTRAN eigenvectors
  - MYSTRAN modal participation factors
  - the explicit SAP response-spectrum ordinates from the benchmark PDF

Validation status:
- `Problem 1-020 modal = validated`
- `Problem 1-020 SRSS displacement = validated`
- `Problem 1-020 SRSS beam/column M33 = validated`
