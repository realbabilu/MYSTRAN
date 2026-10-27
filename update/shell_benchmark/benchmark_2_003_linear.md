# Linear shells — 2-003c: mesh convergence

27 October 2026 rebenchmark: all MYSTRAN entries were rerun using the rebuilt mixed-hybrid DSQK binary. Nastran entries retain previous outputs for unchanged decks. Default residual acceptance remains 1e-6. See [rebenchmark report](../Linear_Shell_Rebenchmark_20261027.md).


MYSTRAN convergence rerun with the 27 October mixed-hybrid DSQK build; Nastran outputs retained.

Meshes: 6x1, 12x1, 30x1, 60x2, 90x3, 180x6. This is a combined angular/radial refinement path, not constant-aspect-ratio refinement. Geometry r=4.12..4.32 in, 90-degree arc, t=.1 in, E=1e7 psi, nu=.25. Root DOFs1234 constrained, rotations5/6 free, as in source Model C. Unit tip forces use radial tributary weights; observe the average of the inner and outer physical tip corners.

All 132 MYSTRAN and 24 Nastran runs completed. Nastran needed a local scratch directory and 1GB sequential retry after shared-scratch and memory-allocation failures; no physical deck changes. The generated90x3 reproduces the imported Model C results within 1.014e-7 relative across all26 formulations and three observables.

References: IN=.0886 in; OUT primary=.5004 in; OUT shell-torsion alternative=.4546 in. The two OUT references are distinct. SAP Model C thick=.0884/.4558 in. A curve need not converge to1 against the primary OUT reference because shell torsion and end restraint differ from the hand calculation.

## Reading the two OUT convergence plots

**Both OUT plots use exactly the same computed tip UZ. Only the normalization denominator changes; they are not different solver runs.**

| Plot | Vertical axis | Meaning of the reference |
|---|---|---|
| `tip_UZ` | UZ / 0.5004 in | Original independent beam calculation, using the rectangular-section torsion constant |
| `tip_UZ_alternative_shell_torsion` | UZ / 0.4546 in | The same independent calculation corrected for the shell-theory torsion constant |

The [CSI SAP2000 Problem 2-003 verification PDF](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-003.pdf), in its out-of-plane discussion, explicitly explains this distinction. The beam calculation uses J=0.000045776 in^4, whereas the shell approximation uses J=b*t^3/3. Their ratio is 0.68664. Only the twisting contribution is corrected:

- Original displacement: 0.35414 bending + 0.00010 shear + 0.14619 twisting = approximately **0.5004 in**.
- Shell-torsion adjustment: 0.35414 + 0.00010 + 0.68664*0.14619 = approximately **0.4546 in**.

For **DSQK mixed hybrid at 180x6**, UZ is **0.49542740 in**. The same displacement gives **0.990063** against 0.5004 and **1.089810** against 0.4546. This is an improvement against the beam target and a regression against the shell alternative; both targets are retained.

The corrected reference is still approximate: CSI also notes a mismatch between the modeled end restraints and the hand calculation. Conversely, closeness to 0.5004 alone does not establish better shell-theory accuracy. The different asymptotes of the other formulations warrant a separate torsional-response investigation; these plots do not establish its cause. Both references remain visible to avoid selecting a target after seeing the results.

### IN displacement: reference 0.0886 in

![](figures/linear_003_convergence_2003_tip_UY.png)

### OUT displacement: original beam reference 0.5004 in

![](figures/linear_003_convergence_2003_tip_UZ.png)

### Same OUT displacement: shell-torsion reference 0.4546 in

![](figures/linear_003_convergence_2003_tip_UZ_alternative_shell_torsion.png)

## Finest mesh:180x6

| Formulation | IN UY (in) | IN/reference | OUT UZ (in) | OUT/0.5004 | OUT/0.4546 |
|---|---:|---:|---:|---:|---:|
| CQUAD4_DKMQ20 | 0.0871887 | 0.984071 | 0.4956924 | 0.9905924 | 1.090393 |
| CQUAD4_DSQK | 0.08718 | 0.9839729 | 0.4954274 | 0.9900627 | 1.08981 |
| CQUAD4_MIN4 | 0.08836968 | 0.9974004 | 0.4978447 | 0.9948935 | 1.095127 |
| CQUAD4_MIN4T | 0.08836968 | 0.9974004 | 0.4974966 | 0.9941979 | 1.094361 |
| CQUAD4_MITC4 | 0.08718872 | 0.9840714 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUAD4_MITC4plus | 0.08718872 | 0.9840714 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUAD4_NATIVE | 0.0885334 | 0.9992483 | 0.4983749 | 0.9959531 | 1.096293 |
| CQUAD4_SIMO | 0.08853374 | 0.9992521 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUADR_DKM24AU | 0.0871887 | 0.984071 | 0.4956924 | 0.9905924 | 1.090393 |
| CQUADR_DKM24EA | 0.0885337 | 0.9992517 | 0.4956924 | 0.9905924 | 1.090393 |
| CQUADR_DKMQ24 | 0.0871887 | 0.984071 | 0.4956924 | 0.9905924 | 1.090393 |
| CQUADR_MBP1C0 | 0.0885337 | 0.9992517 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUADR_MITC4PD | 0.0871887 | 0.984071 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUADR_NATIVE | 0.08866204 | 1.0007 | 0.4980181 | 0.99524 | 1.095508 |
| CQUADR_Q4EASANS | 0.0885337 | 0.9992517 | 0.4969579 | 0.9931214 | 1.093176 |
| CQUADR_Q4RS | 0.08718872 | 0.9840714 | 0.4987362 | 0.9966751 | 1.097088 |
| CQUADR_SIMO | 0.08852773 | 0.9991844 | 0.4969579 | 0.9931214 | 1.093176 |
| CTRIA3_DSG3 | 0.08647893 | 0.9760602 | 0.4953932 | 0.9899943 | 1.089734 |
| CTRIA3_MIN3 | 0.08647932 | 0.9760645 | 0.4974966 | 0.9941979 | 1.094361 |
| CTRIA3_MITC3plus | 0.08647547 | 0.9760212 | 0.4977151 | 0.9946346 | 1.094842 |
| CTRIA3_NATIVE | 0.08647842 | 0.9760545 | 0.4979638 | 0.9951314 | 1.095389 |
| CTRIA3_T3FF | 0.08647932 | 0.9760645 | 0.4960001 | 0.9912071 | 1.091069 |
| CTRIAR_DKMT18 | 0.08647932 | 0.9760645 | 0.4949886 | 0.9891859 | 1.088844 |
| CTRIAR_KMITC3plus | 0.08806428 | 0.9939535 | 0.498324 | 0.9958513 | 1.096181 |
| CTRIAR_NATIVE | 0.08832367 | 0.9968812 | 0.4976903 | 0.9945848 | 1.094787 |
| CTRIAR_T3FFD | 0.08647932 | 0.9760645 | 0.4960001 | 0.9912071 | 1.091069 |
