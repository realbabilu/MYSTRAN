# CSI 2-013 — Thermal constant through thickness

Independent accuracy compares each solver with the same supplied CSI/analytical observable. Missing references and unverified normalization remain PENDING. Nastran–MYSTRAN differences are not accuracy metrics. Zero references use absolute residuals. Numerical data are retained locally; public tables embed the compared values.

2-013B uses current native OGS1 surface100. Joint8=(.24,.12,0)in found by geometry; Z1/Z2=-.0005/+.0005in. Q8:5elements/20GRID; T6:10elements/25GRID. E=1e6psi,nu=.25,alpha=5.5e-6/F,DeltaT=100F,no gradient,MID3=MID2=1. Both solvers share physical inputs.

Nastran outputs8cornerGRID andMID/Z1/Z2; MYSTRAN outputsall20Q8/all25T6GRID atZ1/Z2. Only common native fibers/locations are paired. Z2 normal stresses equal Z1; table below showsZ1 once. MYSTRAN averaging uses area/4 orarea/3; uniform target OES contributors agree withOGS. Nastran nonuniform averaging identity is not established by this uniform field. Surface axes align basicXY,NORMALZ,CID0. No principal-angle comparison. Rounded -733.33psi errors are dominated by reference rounding; exact analytical target=-733.3333333333psi.

2-013A supplied free-expansion targets Ux=.000132in,Uy=.000066in; current public quantitative re-export here focuses onB. No legacy parity chart establishes A accuracy.

| Family | MYSTRAN | Output | Fiber | Reference | Nastran | N error%/abs | MYSTRAN | M error%/abs |
|---|---|---|---|---|---|---|---|---|
| Q8 | SIMOQ8 | SXX | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | SIMOQ8 | SYY | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | SIMOQ8 | TXY | Z1 | 0 | 1.7114463e-13 | abs 1.7114463e-13 | 2.9338538e-14 | abs 2.9338538e-14 |
| Q8 | SIMOQ8 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | SIMOQ8 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | SIMOQ8 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MITC8 | SXX | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | MITC8 | SYY | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | MITC8 | TXY | Z1 | 0 | 1.7114463e-13 | abs 1.7114463e-13 | 2.9338538e-14 | abs 2.9338538e-14 |
| Q8 | MITC8 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MITC8 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MITC8 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | ANS8BDG6 | SXX | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | ANS8BDG6 | SYY | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | ANS8BDG6 | TXY | Z1 | 0 | 1.7114463e-13 | abs 1.7114463e-13 | 2.9338538e-14 | abs 2.9338538e-14 |
| Q8 | ANS8BDG6 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | ANS8BDG6 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | ANS8BDG6 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | HBQ8 | SXX | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | HBQ8 | SYY | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | HBQ8 | TXY | Z1 | 0 | 1.7114463e-13 | abs 1.7114463e-13 | 2.9338538e-14 | abs 2.9338538e-14 |
| Q8 | HBQ8 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | HBQ8 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | HBQ8 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MACQ8D | SXX | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | MACQ8D | SYY | Z1 | -733.33 | -733.33331 | 0.00045177318 | -733.33331 | 0.00045177318 |
| Q8 | MACQ8D | TXY | Z1 | 0 | 1.7114463e-13 | abs 1.7114463e-13 | 0 | abs 0 |
| Q8 | MACQ8D | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MACQ8D | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| Q8 | MACQ8D | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | SIMOT6 | SXX | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | SIMOT6 | SYY | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | SIMOT6 | TXY | Z1 | 0 | 4.9571241e-14 | abs 4.9571241e-14 | 2.5895334e-13 | abs 2.5895334e-13 |
| T6 | SIMOT6 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | SIMOT6 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | SIMOT6 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MITC6 | SXX | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | MITC6 | SYY | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | MITC6 | TXY | Z1 | 0 | 4.9571241e-14 | abs 4.9571241e-14 | 4.8632702e-13 | abs 4.8632702e-13 |
| T6 | MITC6 | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MITC6 | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MITC6 | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MH6T | SXX | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | MH6T | SYY | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | MH6T | TXY | Z1 | 0 | 4.9571241e-14 | abs 4.9571241e-14 | 3.1579679e-14 | abs 3.1579679e-14 |
| T6 | MH6T | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MH6T | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | MH6T | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | REZAIEE | SXX | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | REZAIEE | SYY | Z1 | -733.33 | -733.33344 | 0.00046841921 | -733.33331 | 0.00045177318 |
| T6 | REZAIEE | TXY | Z1 | 0 | 4.9571241e-14 | abs 4.9571241e-14 | 2.6526928e-13 | abs 2.6526928e-13 |
| T6 | REZAIEE | UX | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | REZAIEE | UY | GRID | 0 | 0 | abs 0 | 0 | abs 0 |
| T6 | REZAIEE | UZ | GRID | 0 | 0 | abs 0 | 0 | abs 0 |


## Figures

![2_013B_GPSTRESS_reference_error.png](figures/2_013B_GPSTRESS_reference_error.png)

[Benchmark index](README.md) · [Audited input DAT archive](DAT/quadratic_shell_benchmark_DAT.zip)
