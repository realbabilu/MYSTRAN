# Discussion cases 1–16: published-source retest, 25 October 2026

**Still unresolved in the published implementation: cases 2, 12 and 13. Cases 14 and 16 are retired because MITC3+HB was removed; KMITC3+ is assessed separately as the successor.**

Source: commit `1f0c8f8f51683a6b8fd933d9270fad677612e729`. Fresh CMake build of the publication checkout, using the documented external linear algebra libraries. No Fortran modifications were made in this audit. All 70 requested runs completed normally; accuracy status is assessed separately.

| Case | Formulation / benchmark | Status | Fresh evidence / remaining limit |
|---:|---|---|---|
| 1 | DKMQ24 / 2-002 LC3+LC6 | Closed | n24 UY=0.4320889 vs 0.4321; RZ=0.036 exactly. Six meshes rerun. |
| 2 | DKMT18 / 2-002 LC4 | Known fail | n24 UYabs=0.002352926 vs 0.00341: 31.00% low. Baseline selector; opt-in SF does not close the baseline case. |
| 3 | DKMT18 / 2-005 clamped LC1 | Closed | n24 UZ=-2.614101 vs -2.60: 0.5423% error. Four meshes. |
| 4 | MITC3+ Lee / 2-005 simple LC1 | Closed | n24 UZ=-12.97491 vs -12.97: 0.0379% error. Four meshes; pressure response remains restored. |
| 5 | DKMT18 / 2-005 simple LC1 | Closed | n24 UZ=-12.97928 vs -12.97: 0.0715% error. Four meshes. |
| 6 | DKMQ24 / 2-007 LC1 | Closed | n32 UX=0.09342421 vs 0.094: 0.6125% error; n8/n16/n32 rerun. |
| 7 | DSQK spin / 2-007 LC1 | Closed for hemisphere | n32 UX=0.09326114: 0.7860% error. Twisted beam remains a separate known accuracy failure. |
| 8 | DSG3 spin / 2-007 LC1 | Closed for hemisphere | n32 UX=0.09301597: 1.0468% error. Separate twisted LC1 regression remains documented. |
| 9 | MITC3+ Lee / 2-007 LC1 | Closed for hemisphere | n32 UX=0.09367761: 0.3430% error. Shear-patch trade-off remains separate. |
| 10 | DKMT18 / 2-010 LC1 | Closed for duplicate pressure | 24x48 axial UZ=-0.1239743 vs -0.12414: 0.1335% error. Radial Ur=0.1116052 remains 10.10% low on this mesh; no universal pipe-accuracy closure. |
| 11 | DKMQ24 / 2-013a LC1 | Closed | UX=0.000132 and UY=0.000066 exactly match the free-expansion targets. |
| 12 | Multiple / 2-014b LC1 | Open in published source | 17/22 formulations still zero. Nonzero: MIN4, MIN4T, MIN3, MITC3+ Lee and KMITC3+. Earlier all-formulation load fixes remain isolated/unpublished; recovery is a distinct audit. |
| 13 | DSG3 / 2-016c bending LC2 | Open in published source | Critical edge stress=703.6821 ksi vs 0.1165, for both ND1 and ND6. Published DSG3 lacks the previously tested local-z Kg contribution. Solver completion does not make this an accuracy pass. |
| 14 | MITC3+HB / 2-016c bending LC2 | Retired / replaced | HB no longer exists. KMITC3+ successor first critical stress=0.1157315 ksi (0.6597% error); default ND1/ND6 complete. This is successor evidence, not a repair certification of HB. |
| 15 | DKMQ24 / 2-017 LC2 | Closed for divergent Pcr defect | First Pcr=153.6169 -> 133.2703 kip (2x12 -> 4x24), toward 127.91 kip; finest error 4.19%. First-Pcr acceptance only. |
| 16 | MITC3+HB / 2-017 LC2 | Retired / replaced | HB removed as requested. KMITC3+ successor Pcr=131.5733 -> 127.3982 kip. Fine result agrees with prior exact-mesh Python=127.3981539; default eigensolver, one requested Pcr. |

## Important distinctions

Closure concerns the originally reported defect, not universal formulation accuracy. Case 7/8 hemisphere closure does not close DSQK twisted-beam accuracy or the DSG3 twisted-LC1 trade-off. Case 10 closes doubled pressure, while radial mesh sensitivity remains. Case 15 closes the worsening first-Pcr trend, not every higher eigenpair.

Cases 12 and 13 expose a publication gap: the earlier isolated thermal-gradient and DSG3 Kg-z corrections were not carried into this published source. This audit does not describe those isolated successful tests as current production support. In particular, the previously accepted corrected-DSG3 mode2 residual limitation is distinct from the current published kernel returning the wrong first critical stress.

The 22-formulation gradient test uses KMITC3+ instead of the removed HB selector. KMITC3+ bending first-Pcr is close to the target and both default ND1/ND6 requests complete. Its reported preload maximum fiber principal stress is 147.1633 ksi against the nominal 100 ksi boundary target; this is retained as a separate stress-recovery/prestress diagnostic, not certified by first-Pcr agreement.

## Input and scope controls

SPC1 lists are split into at most six GRID IDs per physical line. 2-017 uses the source SAP active UX/UZ/RY on all nodes, including centroids, with original support restraints, and requests one Pcr. No DARPACK override or tolerance/backend change was introduced. Hemisphere displacements are transformed from each GRID output coordinate system to basic axes and evaluated at the actual force application nodes.

Files: [comparison CSV](validation_discussion_20261025/comparison.csv), [raw extracted results](validation_discussion_20261025/results.json), [source/build provenance](validation_discussion_20261025/provenance.json). Local workspace `H:/Projects/MYSTRAN/codex_mod/discussion_retest_20261025` retains the runner `retest.py`, DAT/F06/OP2, run.json input/executable hashes, CMake cache and logs. Source hash and executable hash are recorded separately from older isolated test binaries.

Counts: 11 cases closed for their documented defects (1,3,4,5,6,7,8,9,10,11,15); 3 unresolved (2,12,13); 2 retired (14,16). No automatic threshold was used to declare every completed solver run an accuracy pass.


## Case 12 follow-up — 25 October 2026

Case 12 is now fixed and verified; owner approved publication: 22/22 gradient responses, 396/396 principal runs and 440/440 assertions pass, with 22 multi-gradient/TREF and 22 OP2 checks. The earlier table remains the historical audit of published commit 1f0c8f8f. Case 13 Kg-z is not included. See [case 12 completion](CASE_12_Temperature_Gradient_20261025.md).
