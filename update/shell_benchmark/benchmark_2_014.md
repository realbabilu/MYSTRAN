# CSI 2-014 — Annular thermal gradient

Independent accuracy compares each solver with the same supplied CSI/analytical observable. Missing references and unverified normalization remain PENDING. Nastran–MYSTRAN differences are not accuracy metrics. Zero references use absolute residuals. Numerical data are retained locally; public tables embed the compared values.

One actual Python mesh per family:28radial×32circumferential=896Q8 or1792T6 elements, not the originalCSI9×16mesh. Radius3..30in,h=1in,MID3=MID2=1,outerclamped/innerfree. References UZ=.01931in,R2=.00352rad,M11=-2.181kip-in/in,M22=-4.406kip-in/in.

Displacement/rotation are atinner(3,0,0);M11atouter(30,0,0),M22atinner. Native tensors transformed once to radial/tangential axes, followed by declared area-per-corner averaging. This does not assert identicalSAP averaging. Mesh comes from Python, notS2K; S2K is used to audit physical definitions/orientation. Executed modelB identity is retained although the supplied table labelsModelA.

| Family | MYSTRAN | Output | Reference | Nastran | N error(%) | MYSTRAN | M error(%) |
|---|---|---|---|---|---|---|---|
| Q8 | SIMOQ8 | UZ | 0.01931 | 0.019121716 | 0.9750606 | 0.019121403 | 0.97668113 |
| Q8 | SIMOQ8 | R2 | 0.00352 | 0.0034965987 | 0.66481008 | 0.0034966003 | 0.66476378 |
| Q8 | SIMOQ8 | M11 | -2.181 | -2.1857838 | 0.21933882 | -2.1859088 | 0.22507076 |
| Q8 | SIMOQ8 | M22 | -4.406 | -4.4607646 | 1.2429552 | -4.4600211 | 1.2260797 |
| Q8 | MITC8 | UZ | 0.01931 | 0.019121716 | 0.9750606 | 0.019120781 | 0.97990289 |
| Q8 | MITC8 | R2 | 0.00352 | 0.0034965987 | 0.66481008 | 0.0034967663 | 0.66004764 |
| Q8 | MITC8 | M11 | -2.181 | -2.1857838 | 0.21933882 | -2.1858854 | 0.22399956 |
| Q8 | MITC8 | M22 | -4.406 | -4.4607646 | 1.2429552 | -4.4588968 | 1.2005619 |
| Q8 | ANS8BDG6 | UZ | 0.01931 | 0.019121716 | 0.9750606 | 0.019120781 | 0.97990289 |
| Q8 | ANS8BDG6 | R2 | 0.00352 | 0.0034965987 | 0.66481008 | 0.0034967663 | 0.66004764 |
| Q8 | ANS8BDG6 | M11 | -2.181 | -2.1857838 | 0.21933882 | -2.1858854 | 0.22399956 |
| Q8 | ANS8BDG6 | M22 | -4.406 | -4.4607646 | 1.2429552 | -4.4588968 | 1.2005619 |
| Q8 | HBQ8 | UZ | 0.01931 | 0.019121716 | 0.9750606 | 0.019120781 | 0.97990289 |
| Q8 | HBQ8 | R2 | 0.00352 | 0.0034965987 | 0.66481008 | 0.0034967663 | 0.66004764 |
| Q8 | HBQ8 | M11 | -2.181 | -2.1857838 | 0.21933882 | -2.1858854 | 0.22399956 |
| Q8 | HBQ8 | M22 | -4.406 | -4.4607646 | 1.2429552 | -4.4588968 | 1.2005619 |
| Q8 | MACQ8D | UZ | 0.01931 | 0.019121716 | 0.9750606 | 0.019121073 | 0.97838847 |
| Q8 | MACQ8D | R2 | 0.00352 | 0.0034965987 | 0.66481008 | 0.0034965111 | 0.66729714 |
| Q8 | MACQ8D | M11 | -2.181 | -2.1857838 | 0.21933882 | -2.1855822 | 0.21009447 |
| Q8 | MACQ8D | M22 | -4.406 | -4.4607646 | 1.2429552 | -4.4608814 | 1.2456052 |
| T6 | SIMOT6 | UZ | 0.01931 | 0.019077269 | 1.2052338 | 0.019077269 | 1.2052338 |
| T6 | SIMOT6 | R2 | 0.00352 | 0.0034948536 | 0.71438581 | 0.0034948536 | 0.71438581 |
| T6 | SIMOT6 | M11 | -2.181 | -2.183165 | 0.099266247 | -2.1831648 | 0.099255468 |
| T6 | SIMOT6 | M22 | -4.406 | -4.4112425 | 0.11898608 | -4.4112421 | 0.11897611 |
| T6 | MITC6 | UZ | 0.01931 | 0.019077269 | 1.2052338 | 0.019147793 | 0.84001641 |
| T6 | MITC6 | R2 | 0.00352 | 0.0034948536 | 0.71438581 | 0.003498537 | 0.60974431 |
| T6 | MITC6 | M11 | -2.181 | -2.183165 | 0.099266247 | -2.1815198 | 0.023832765 |
| T6 | MITC6 | M22 | -4.406 | -4.4112425 | 0.11898608 | -4.4139646 | 0.18076631 |
| T6 | MH6T | UZ | 0.01931 | 0.019077269 | 1.2052338 | 0.019087791 | 1.1507434 |
| T6 | MH6T | R2 | 0.00352 | 0.0034948536 | 0.71438581 | 0.0034979682 | 0.62590355 |
| T6 | MH6T | M11 | -2.181 | -2.183165 | 0.099266247 | -2.1814193 | 0.019223032 |
| T6 | MH6T | M22 | -4.406 | -4.4112425 | 0.11898608 | -4.4153917 | 0.21315662 |
| T6 | REZAIEE | UZ | 0.01931 | 0.019077269 | 1.2052338 | 0.019147793 | 0.84001641 |
| T6 | REZAIEE | R2 | 0.00352 | 0.0034948536 | 0.71438581 | 0.003498537 | 0.60974431 |
| T6 | REZAIEE | M11 | -2.181 | -2.183165 | 0.099266247 | -2.1815198 | 0.023832765 |
| T6 | REZAIEE | M22 | -4.406 | -4.4112425 | 0.11898608 | -4.4139646 | 0.18076631 |


## Figures

![CSI_2_014_Q8.png](figures/CSI_2_014_Q8.png)

![CSI_2_014_T6.png](figures/CSI_2_014_T6.png)

[Benchmark index](README.md) · [Audited input DAT archive](DAT/quadratic_shell_benchmark_DAT.zip)
