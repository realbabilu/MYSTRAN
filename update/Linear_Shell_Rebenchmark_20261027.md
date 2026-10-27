# Linear shell rebenchmark — 27 October 2026

Rebuilt with DSQK mixed hybrid from b1da4a6a plus the same DKMT18 SF candidate / corrected DSG3 local-z Kg element deltas used in the prior linear snapshot. BKRTOL is absent; default 1e-6 physical-pencil residual acceptance and solver defaults are unchanged. Quadratic material is unchanged.

## Run completion

| Solver | Status | Configurations |
|---|---|---:|
| MYSTRAN | completed | 1825 |
| MYSTRAN | solver_error | 1 |
| NASTRAN | completed | 332 |

All MYSTRAN configurations above were rerun with one OpenBLAS/OMP thread and four concurrent jobs. Nastran results are retained from prior runs for unchanged decks; no new Nastran solve or license claim is made. One additional fresh SOL101 independently checks DSG3 BEND preload. Successful completion is separate from physical benchmark accuracy.

## Known failures

- prob_2_016c_thick_beamseg_bend / CTRIA3_DSG3: *ERROR: buckling eigenpair residual exceeds 1E-6           2   1.9954529446686861E-006

The DSG3 BEND LC2 six-mode failure remains visible and is not replaced with a one-mode pass. Its LC1 static preload is independently refreshed.

## Regression and trade-offs

Compared 4114 measured non-DSQK MYSTRAN metrics against the prior published dataset; 0 differ at stored output precision. Differences including near-zero roundoff values are retained in the delta CSV.

| DSQK case / observable | New result / reference |
|---|---:|
| prob_2_002_twoaxis_nx48_ny4 / tip_abs_UY | 0.902094721 |
| prob_2_002_twoaxis_nx96_ny8 / tip_abs_UY | 0.970073021 |
| prob_2_002_twoaxis_nx192_ny16 / tip_abs_UY | 0.999486804 |
| prob_2_003c_nx180_ny6 / tip_UZ | 0.990062750 |
| prob_2_003c_nx180_ny6 / tip_UZ_alternative_shell_torsion | 1.089809503 |

2-002 LC4 improves, but the finest tested mesh is not a demonstrated asymptotic limit. The 2-003 beam and shell references remain distinct; mixed hybrid improves the beam comparison and worsens the shell comparison. Twisted 2-004 remains deficient. Higher-mode dynamics and constant-shear enriched bending limitations remain documented; no universal DSQK closure is claimed.

## Reproduction and outputs

Executable SHA256: `c76a0c35f980c579d1de0e803f32163939eb493e50b5d527221788526a775d91`.

[All linear results](shell_benchmark/comparison_all_linear.md) · [CSV](shell_benchmark/comparison_all_linear.csv) · [Run status](shell_benchmark/data_linear/run_status_linear.csv) · [Deltas](shell_benchmark/data_linear/rebenchmark_deltas_20261027.csv) · [Build/source hashes](shell_benchmark/data_linear/provenance_linear.json)

Full raw native runs and build logs are retained locally under H:/Projects/MYSTRAN/codex_mod/linear_rebenchmark_20261027. The source snapshot is preserved there; the dirty main source checkout was not overwritten.
