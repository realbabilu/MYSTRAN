# DSQK MacNeal patch audit â€” 27 October 2026

## Decision

Retain the existing mixed hybrid DSQK implementation. The supplied handoff selects the same mixed An already implemented on `research/dsqk-mixed-20261026`; no new element source change or rebuild is required by this package. MixedShear remains a separate Python research candidate. This audit adds validation of the existing mixed hybrid. Source and reports are published together; release binaries and solver tolerances are unchanged.

The selected formulation keeps covariant nodal shear, consistent-spin drilling alpha=0.03, mixed elimination An, and inherited bending expression Kbu + An^T KbDelta An. It is not certification of the complete published mixed functional. Source inspection confirms the artificial-shear SHELL_T_MOD fallback remains present.

## Reproduction

The supplied Python audit was rerun offline using the local published DAT archive, with immutable hashes recorded below. Four distorted decks contain five physical patch fields: MH membrane LC1 and bending LC2, shear LC3, bending LC4 and twist LC5. Original source decks and supplied results remain under input; rerun results are separate. No boundary conditions were changed.

| Check | Result |
|---|---:|
| Native vs Python assembled stiffness, maximum relative Frobenius difference | 3.77e-16 |
| Native vs Python translation, maximum absolute difference | 2.91e-11 in |
| Native vs Python rotation, maximum absolute difference | 1.41e-16 rad |
| Energy from native KGG/F06 vs Python, maximum relative difference | 9.84e-10 |
| Reconstructed nodal reaction vectors, maximum relative difference | 2.63e-8 |

Native results use the already tested mixed binary SHA256 `131bb8959464757f244ae9014bebcc8f7f165d1c09419ba9c3c47be040658635`, without BKRTOL. Four SOL101 runs completed. Reaction comparison is reconstructed from KGG times the displacement vector; it is not an independent F06 SPC reaction or prescribed boundary-traction audit. F06 rounding and cancellation limit the energy/reaction comparison.

Python hybrid membrane and bending/twist fields reproduce the analytical strains/curvatures and energy to roundoff (MH bending curvature error 5.86e-15). Native operator and response parity support the same kinematics; an independent native recovered-stress table audit has not been added here.

For the constrained shear patch (all nodal bending rotations held):

| Python variant | Translation error, in | Total energy / target |
|---|---:|---:|
| DK baseline | 2.90705e-11 | 1.00000012817 |
| Mixed hybrid | 2.90705e-11 | 1.00000012817 |
| MixedShear | 1.79313e-4 | 1.28029e-7 |

Hybrid and DK retain a nonzero enriched curvature in this shear test (~7.486e-4); do not call this a zero-curvature certification. MixedShear satisfies discrete free equilibrium but nearly suppresses the target shear energy. This constrained test is not proof about all mixed shear equilibrium fields.

## Existing trade-offs remain

The prior 87-run native regression remains the validation basis for the unchanged binary. Two-axis 2-002 LC4 hybrid ratios are 0.902095 / 0.970073 / 0.999487 at 48x4 / 96x8 / 192x16. The handoff's finer Python runs indicate a limit near 1.012: the 192x16 agreement is not the asymptotic limit. Fine-mesh values were not independently rerun in this audit.

2-003 OUT 180x6 remains approximately +8.98% against the shell alternative reference, although close to the beam reference. Both targets must remain visible. Twisted 2-004 remains poor: n24 ratios approximately 0.50840 (Fz) and 0.35214 (Fy). No new universal fix is claimed.

The handoff's additional buckling and dynamics figures are supplied evidence, not newly certified native results. Exact 2-016c candidate-to-candidate prestress/Kg and six-mode audits, independent boundary reaction distributions, multiple distortion seeds and rotated/reversed assembled patch tests remain future work. Existing DSG3 Case 13 known failure is unaffected.

## Reproducibility

Run `audit_offline.py`, then `audit_native.py` with the workspace Python environment. The first script uses a local archive path; the second records the exact binary hash. The ZIP contains input modules/decks, scripts, raw native F06/KGG and JSON outputs. Paths can be adjusted when moving the package. No network refresh is needed.

[Native/Python evidence ZIP](validation_dsqk_macneal_20261027/DSQK_MacNeal_Native_Verified_20261027.zip) · [Native results](validation_dsqk_macneal_20261027/native_results.json) · [Python results](validation_dsqk_macneal_20261027/python_results.json)
