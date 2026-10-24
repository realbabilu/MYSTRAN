# 2-017: corrected Python / MYSTRAN parity — 2026-10-25

First critical load parity is confirmed for MITC3+ Lee + Spin and KMITC3+ Nguyen on identical meshes and boundary conditions. These meshes use four centroid triangles per parent quadrilateral.

| Formulation | Parent mesh | T3 elements | Python Pcr (kip) | MYSTRAN Pcr (kip) |
|---|---|---:|---:|---:|
| MITC3+ Lee + Spin | 2 x 12 | 96 | 167.9630195 | 167.9630 |
| MITC3+ Lee + Spin | 4 x 24 | 384 | 136.8640141 | 136.8640 |
| KMITC3+ Nguyen | 2 x 12 | 96 | 131.5733113 | 131.5733 |
| KMITC3+ Nguyen | 4 x 24 | 384 | 127.3981539 | 127.3982 |

Maximum relative difference is 3.62e-7 (0.0000362%), within printed F06 precision. All four runs completed using one requested mode and no DARPACK override.

## Correction to the earlier higher-mode diagnosis

The earlier apparent Python / MYSTRAN spectrum discrepancy was caused by different effective constraints. Original free-field SPC1 lines listed eight GRID IDs; MYSTRAN BD_SPC1 reads six GRID fields per physical line (fields 4–9). The Python reader accepted every listed ID. The coarse native system consequently had 27 extra free DOFs: components 2, 4 and 6 at GRID 7, 8, 15, 16, 23, 24, 31, 32 and 39.

Splitting SPC1 cards into at most six GRID IDs made all six diagnostic eigenvalues agree for both formulations on both meshes. The earlier interpretation as a formulation or eigensolver parity defect is superseded by this input-format finding. Historical outputs remain unchanged for traceability.

For the final SAP comparison, inactive UY, RX and RZ are constrained at every node, including added centroids. UX, UZ and RY remain active, with the original base supports. The published SAP verification acceptance concerns the first Pcr: SAP 127.22 kip and independent reference 127.91 kip.

Production Fortran and the shared eigensolver were not changed by this audit. The instrumented executable used to inspect native matrices is isolated from production.

## Reproduction and evidence

- `minimal_Pcr_parity.json`: final four-run results and input hashes.
- `inputs/`: corrected one-mode decks and native outputs.
- `corrected_SPC1_results.json`: six-mode diagnostic agreement after card splitting.
- `boundary_dof_audit.json`: missing-constraint evidence.
Original malformed-card outputs and the parity runner remain in the local investigation workspace; they are not included in this published evidence subset.

Differences from another convergence plot should be assessed using actual triangle connectivity and element count, not parent mesh labels alone.
