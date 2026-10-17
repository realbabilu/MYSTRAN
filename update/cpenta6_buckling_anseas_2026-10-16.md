# Optional CPENTA6 Claude ANS/EAS buckling branch, 2026-10-16

Git branch: `cpenta6-buckling-anseas`, created from the current `v18.00.a` checkout. Existing uncommitted solid work was preserved; no commit or push was performed.

## Selection contract

```nastran
PARAM,SOLIDTYP,NEWSOLID
PARAM,PENTABUC,ANSEAS
```

| Solid selector | PENTABUC | Analysis | CPENTA6 behavior |
| --- | --- | --- | --- |
| NEWSOLID | omitted / EAS9 | All | Existing EAS9 v2 |
| NEWSOLID | ANSEAS | SOL 105, both stages | Claude original edge-tied ANS + one locally condensed thickness EAS parameter |
| NEWSOLID | ANSEAS | SOL 101 or other nonbuckling solutions | Existing EAS9 v2 |
| LEGACY | either | All | Existing LEGACY |

The static prestress stage inside SOL 105 belongs to the buckling analysis and therefore also uses ANS/EAS when selected. It is not mixed with EAS9 elastic K. Standalone SOL 101 remains EAS9 v2. Existing CPENTA15 and other element formulations are unaffected. `PENTABUC` defaults to EAS9 in PARAMS.f90; unsupported values are fatal input errors. It is serialized as an optional trailing L1A restart record; older records without it default to EAS9. Full restart migration has not been separately certified.

## Source audit

Sources read: original `CPENTA6_ANS_EAS_original.py` and the patch-projector README/audit, plus the new five-formulation follow-up README, head-to-head README and assembled patch script/results. The original source is ported, **not** the symmetric affine projector and not mean-shear correction. Exact original bytes and SHA256 provenance are in `tests/solids/references/cpenta6_anseas/`.

Claude original: transverse natural shear rows are tied at (.5,0,z), (0,.5,z), (.5,.5,z). They are mapped using the constant center-frame natural/Cartesian strain transformation. One thickness strain mode is (detJ0/detJ)*z in natural epsilon_zz, with the original max(1e-14*abs(Kaa),1e-14) regularization. Its internal parameter is locally eliminated. No global enrichment DOF or connectivity change is introduced.

Thermal RHS and stress recovery use the same local parameter equation (Fa-Kau*u)/Kaa. Thermal support is a derived extension; the original Python file has mechanical recovery only. Tests independently integrate the thermal equations and compare every recovery station. Kg integrates the physical enhanced prestress tensor against standard displacement shape gradients. It does not copy the inherited same-Kg uniform-prestress screening and does not claim a derived finite-strain consistent ANS/EAS tangent.

## Validation

Production-linked tests exercise PENTA, ELMDIS, OPT(6) and ELEM_STRE_STRN_ARRAYS; detached helper success alone is not used as validation.

* K agrees with the unchanged Claude original Python on regular/distorted/rotated wedges to about 4.2e-16 relative error. Mechanical stress/strain, independently derived gradient thermal RHS/recovery and recovered-stress Kg comparisons remain below 6e-16. Kg rotation covariance error is 3.29e-16.
* Regular 3x3x3 assembled affine patch residual: 7.75e-16.
* Distorted 3x3x3 assembled affine patch residual: **0.004709884304**, reproducing the follow-up **KNOWN_FAIL**. Relative interior displacement correction: 0.001799808232.
* Distorted single-element weak affine nodal-force diagnostic: **0.07724339560** relative error for the specific random affine gradient used by this regression. This is not the archive's different ~9.77% diagnostic; both indicate failure, not acceptance.
* Full SOL 105 manufactured square cantilever (20 wedges, unit end compression, E=1e6, nu=.25): selected ANSEAS first factor **3098.513**, independent actual-static/recovered-stress-Kg result **3098.513109**; relative printed difference 3.52e-8. Default EAS9 factors begin **2197.873, 2213.842** on this verification case. These numbers are specific to this mesh/operator and do not establish SAP/Nastran accuracy or superiority of ANSEAS.
* Default SOL 105 and explicit PENTABUC=EAS9 eigenvalue lists are identical.
* Standalone SOL 101 with/without PENTABUC=ANSEAS gives identical OP2 displacement arrays. LEGACY with/without the selector also gives identical arrays. Invalid selectors are rejected.
* Existing default solid kernel comparisons were rerun; EAS9 v2 remains unchanged. Existing patch/frame and CPENTA6 5-004 regressions are retained and rerun separately.

## Retained failures and risks

**KNOWN_FAIL: DISTORTED_AFFINE_PATCH** remains intrinsic to the selected original ANS/EAS formulation. An explicit F06/ERR warning identifies this experimental branch; no projector, scaling or tuning hides it.

**KNOWN_RISK: PRESTRESS_TANGENT_NOT_CERTIFIED**. The new full SOL 105 route is internally verified against the same independently assembled small-strain stress-gradient operator. A variational finite-strain ANS/EAS tangent, follower-load terms and matched Nastran/SAP prestress/recovery/mode-shape comparison remain pending. It is **NOT FINAL PASS**. The archive's ~0.891 kip number uses imposed uniform prestress and inherited standard Kg, not this full loaded analysis; that screening value is not an acceptance target for tuning.

**KNOWN_RISK: THERMAL_EXTENSION_SCOPE**. Thermal condensation/recovery pass independent mathematical comparisons, but no matched Nastran thermal-buckling case or nonlinear thermal tangent is supplied.

Solid PLOAD4 G1/G3 input and stress-request OP2 serialization limitations remain unchanged and unaccepted. No result here certifies mass/modal ANSEAS behavior: the selector is intentionally inactive outside SOL 105. EAS9 v2's existing general thermal/recovery limitations outside the ANSEAS branch remain separately documented.

Run `tests/solids/check_penta6_anseas.py` after building. Results: `tests/solids/penta6_anseas_results.json`. Generated BDFs, output and linked probes are ignored. The additional five-formulation follow-up archive is extracted under `python-solids/CPENTA6_FiveFormulations_MultiPatch_Followup_v2/` for source evidence; no files from that archive were modified.
