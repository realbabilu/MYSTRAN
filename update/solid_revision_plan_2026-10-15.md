# Solid revision audit and implementation plan

Date: 2026-10-15. Branch: v18.00.a. Initial source working tree clean.
References: python-solids/handoff/CODEX_HANDOFF.md and SOURCE_MANIFEST.md,
read in full. Archived evidence is historical, not a new MYSTRAN validation.

## Scope updates after the initial audit

The user subsequently required removal of CPYRA14. Steps below that originally
proposed preserving it are superseded: type slot 26 now represents PYRAM13,
METYPE remains 27, and only 5/13-node pyramid cards are accepted.

Executed comparisons established that the handoff's older EAS source does not
reproduce its v2 evidence. The unchanged v2 reference was recovered from
nastran_test/6/SAP5004_CPENTA6_CPYRA5_Nastran_Python_Comparison_v1.zip.
Use its solid_5004/solid3d_eas9_frame_v2.py and natural_to_cartesian helper.
No tuning is authorized before v2 patch, rotation and twisted-beam reproduction.
Those three gates now pass; see solid_revision_results_2026-10-15.md for
the executed metrics and remaining independent feature gates.

## Architecture audit

* Parsing: LOADB0/LOADB dispatch CPYRAM and CPYRA to BD_PYRAM0/BD_PYRAM.
  BD_PYRAM currently requires a third card for node 14; standard 13-node
  connectivity is not supported. BD_PYRAM0 conservatively allocates 16 EDAT
  entries for any continuation. Keep that allocation compatible with both.
* Metadata: SCONTR holds element counts and METYPE; MODEL_STUF holds parallel
  ELMTYP, ELEM_ONAME, NELGP, NUM_SEi tables. Add a new type without renumbering
  existing ones. Material/property lookup is ELEM_PROP_MATL_IIDS.
* EMG dispatch and transformation: EMG initializes local coordinates, material,
  matrices and output. PYRAM uses compact 3*n matrices internally and maps to
  existing indices 6*(node-1)+1:3. Global GRID layout and sparse assembly stay.
* PYRAM5 uses collapsed coordinates and standard quadrature, but NEWSOLID
  currently selects EAS54. Remove that selection; retain standard PYRAM5 in
  both modes. Preserve the separate existing PYRAM14 kernel.
* PYRAM14 uses Liu interpolation with finite-difference derivatives. It is not
  a substitute for Bedrosian PYRAM13; do not reuse those formulas.
* PENTA and HEXA already assemble Kuu/Kua/Kaa and condense nine local modes.
  PENTA transforms basic-axis strain modes to local axes and regularizes Kaa
  as in the supplied Python. HEXA uses local-axis modes. Both silently skip
  correction on solve failure. Verify coordinates, pivots and recovery.
* Mass, thermal loads, geometric stiffness and recovery are separate OPT
  branches. EAS correction currently affects stiffness only. Recovery uses
  ordinary B rather than B-M*inv(Kaa)*Kau; nonuniform thermal loads need the
  matching condensed enhanced thermal RHS. These are validation gaps.
* Pressure enters through BD_PLOAD4 and separate face-load routines. Existing
  FORCE-loaded 5-004 evidence does not validate PLOAD4 implementation.
* Output spans LK9 SOLID_*_OUTPUTS, OFP3_*, FEMAP and OUTPUT2_WRITE_STRESS.
  Pyramid recovery is currently populated at center only; node locations and
  OP2 metadata require their own verification.
* Available: gfortran, CMake, configured build, Python NumPy/SciPy/SymPy and
  pyNastran, local Nastran DAT/OP2/F06 evidence. No fresh Nastran solver run
  has yet been established.

## Ordered implementation and regression gates

1. Capture baseline executable and benchmark outputs. Add reproducible local
   kernel tests; retain source references unchanged. Record exact evidence
   paths, tolerances and current run results separately from archived CSVs.
2. Port Bedrosian analytic shape/derivatives in collapsed coordinates, with
   node order corners 1-4, apex 5, base edges 6-9, lateral edges 10-13. Test
   nodal interpolation including apex limit, partition/derivative sums,
   coordinate reproduction, face trace and finite-difference derivatives.
3. Add eight-point stiffness quadrature with correct collapsed Jacobian
   weighting. Compare full 39x39 matrix against Python on affine and warped
   geometries; check symmetry, six rigid motions, affine energy, rank 32,
   reversed orientation and invalid Jacobians. No stabilizer or extra DOFs.
4. Integrate new PYRAM13 metadata, parsing, EMG/material dispatch. Distinguish
   5/13/14 nodes without changing existing type IDs. Test free/fixed/large
   cards, continuation placement and missing/excess midsides. Preserve 14
   node historical path. Standard PYRAM5 must remain selected by NEWSOLID.
5. Compare actual MYSTRAN displacements with Python and existing Nastran OP2
   on 5-004 consistent loads and 5-007 meshes; run 5-001 patch and 5-008
   incompressibility screens. Evaluate scaled global eigenvalues where
   feasible. Reject claimed parity if geometry, constraints or loads differ.
6. Validate existing PENTA6 and HEXA8 EAS9 separately: full matrices, patch,
   rigid modes, rotated/distorted geometries, Schur residual, Kaa condition,
   5-004 IN/OUT. Report singular condensation explicitly. Preserve LEGACY
   default; do not port EAS18 or change working quadratic solids.
7. Independently gate consistent/lumped mass (total mass, moments, dynamics),
   PLOAD4 (all faces, signs, midside forces and moments), thermal (uniform
   expansion, gradients, EAS RHS condensation), geometric stiffness
   (pre-stress and buckling), and recovery (strain/stress, EAS alpha,
   coordinates, center/corner locations, F06/OP2). A stiffness pass grants
   no acceptance to these features. Mark untested cells NOT VALIDATED.
8. Run existing solid build cases and LEGACY/NEWSOLID regressions. Publish
   Python/Nastran/MYSTRAN table and limitations. Do not call the collection
   production-ready while any required gate remains unvalidated.

## Known failures and risks

* KNOWN_FAIL: THIN_DOME_LOCKING, CPYRAM5 5-007. Standard ~0.00077869,
  EAS18 ~0.00080589, NS alpha=.8 ~0.00296729 versus target 1 at 24x24.
  Standard remains the default despite this documented limitation.
* KNOWN_RISK: LOCAL_HOURGLASS, CPYRAM13 eight-point rank 32/39: six rigid
  modes plus one underintegration mode. Passing tested assemblies cannot
  certify arbitrary meshes. Rank-one v4 stabilization degrades bending;
  27-point stiffness over-stiffens the thin dome. Neither is a default fix.
* KNOWN_RISK: EAS stiffness-only recovery/thermal inconsistency; silent
  condensation fallback; coordinate-frame sensitivity. Validate before
  accepting EAS beyond displacement benchmarks.
* Excluded: global enrichment, CPYRAM5 smoothing/EAS/bubble default, PYRA14
  as standard Nastran topology, and experimental stabilization.

## Acceptance ledger (initial state)

| Formulation | Stiffness | Mass | Pressure | Thermal | Kg | Recovery |
|---|---|---|---|---|---|---|
| PYRAM13 | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED |
| PENTA6 EAS9 | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED |
| HEXA8 EAS9 | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED | NOT VALIDATED |

Existing source coverage and historical evidence do not change this ledger
without executed checks against the revised source.
