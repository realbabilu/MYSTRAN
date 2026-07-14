## RSA output status

Date: 2026-07-14

Decks checked:
- `D:\18a\femap_RSA\beam_RS.dat`
- `D:\18a\femap_RSA\beam_RS_abs.dat`
- `D:\18a\femap_RSA\beam_RS_cqc.dat`
- `D:\18a\femap_RSA\beam_RS_dense.dat`
- `D:\18a\femap_RSA\nx_hinge_nx.dat`
- `D:\18a\femap_RSA\nx_hinge_nx_dense.dat`

Current parser / solver status:
- `TABLED1`
- `DTI,SPECSEL`
- `PARAM,SCRSPEC`
- `PARAM,OPTION`

These remain stable on the current local tests.

F06 output families confirmed present on RSA runs:
- real eigenvalue table
- combined `DISPLACEMENTS`
- combined `VELOCITIES`
- combined `ACCELERATIONS`
- combined `SPC FORCES`
- combined element `STRESSES`
- `SCRSPEC COMPATIBILITY SUMMARY (MINIMAL MODES PATH)`

Observed deck behavior:
- `beam_RS.dat`
  - default ARPACK path runs normally
  - combined `SPCFORCE` and beam stress tables are emitted
- `beam_RS_dense.dat`
  - still runs normally after the latest guard work
  - beam RSA baseline remains unchanged
- `beam_RS.dat` and `beam_RS_dense.dat`
  - RSA combined grid-vector tables now print with standard F06 headers
  - `SOL 103 + SCRSPEC` combined `DISP/VELO/ACCE` no longer reuses the modal `EIGENVECTOR` label
- `beam_RS_abs.dat`
  - `PARAM,OPTION,ABS` was run directly after the parser cleanup
  - `ERR/F06` now correctly report `Combined RSA displacement output was assembled with ABS across retained modes.`
  - the old stale warning claiming `ABS` still fell back to `SRSS` is removed
- `beam_RS_cqc.dat`
  - `PARAM,OPTION,CQC` was run directly on the same beam model
  - `ERR/F06` report `Combined RSA displacement output was assembled with CQC across retained modes.`
  - for this 2-mode beam deck, reported combined displacement is numerically the same to printed precision as `SRSS/ABS`
- `nx_hinge_nx.dat`
  - commercial-compatible modal basis remains the ARPACK path
  - combined `SPCFORCE` and element stress blocks are emitted
- `nx_hinge_nx_dense.dat`
  - now safely falls back to ARPACK because `SCRSPEC + single SUPORT` is guarded
  - wrapper deck was cleaned so it no longer duplicates `DTI,SPECSEL`

Important limitation:
- Output tables now exist beyond displacement-only RSA.
- However, for hinge-style `SUPORT + CONM2` cases, the residual physics mismatch is still in the participation / rigid-support contribution, not in whether the F06 blocks are written.
- So:
  - table presence is acceptable
  - combined values for hinge-like cases are still not fully commercial-equivalent

Practical interpretation:
- It is safe to keep extending the minimal RSA compatibility path in isolated writer / reporting steps.
- It is not yet safe to claim full commercial equivalence for all combined force/stress quantities on hinge-style support-mass problems.

Writer note:
- The latest fix is formatting/reporting only:
  - `LINK9` now requests full grid-table headers for RSA combined `DISP/VELO/ACCE`
  - `WRITE_GRD_PRT_OUTPUTS` suppresses `OUTPUT FOR EIGENVECTOR` and uses `DISPLACEMENTS` instead of `EIGENVECTOR` when `SCRSPEC='Y'`
- This does not change combined-response math; it only restores proper F06 table identity.

Combination-method note:
- `SRSS`, `CQC`, and `ABS` are active in the current `LINK9` RSA combine path.
- `NRL` is still accepted only for compatibility parsing and currently warns that it falls back to `SRSS`.

Current local blocker:
- A quick `pyNastran` readback audit of the RSA `OP2` files could not be completed in the current shell environment because the local Python runtime returned `Access is denied`.
- The RSA `OP2` files are being written (`beam_RS*.OP2` exist), but binary readback still needs to be checked once the local Python execution path is usable again.
