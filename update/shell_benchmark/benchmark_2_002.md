# CSI 2-002 — fixed source meshes A/Q8 and E/T6

20 October 2026. **One source mesh per model. No remeshing or convergence study.** Existing S2K vertices/topology are enriched with quadratic edge nodes only. Nine MYSTRAN formulations and two Nastran models are compared for six load cases; the root-midpoint boundary diagnostic is kept separate.

## Source models and reference

A comes from Example 2-002a-thick.s2k: 6 rectangular areas → 6 CQUAD8 / 33 nodes. E comes from Example 2-002e-thick.s2k: 24 triangles, **four per parent panel** → 24 CTRIA6 / 63 nodes. E is not the old two-triangle Model D. Original joint coordinates/IDs, panel centers and area IDs are preserved, including area53. No new mesh level was generated.

The source E triangle windings are mixed (12 normal +Y, 12 normal −Y). DAT order is consistently normalized to −Y like A; per-element changes are recorded. Original S2K files remain unchanged. This reorder changes element-frame/fiber conventions, not geometry or the basic physical loads.

Units lb/in, length6 in, width0.2 in, shell h0.1 in, E=10,000,000 psi, G=3,846,153.84615385 psi, nu0.3. Both source sections are Shell-Thick; PSHELL MID1=MID2=MID3=1, transverse shear factor5/6. Basic XZ is the shell plane; no Y/Z swap is applied to the observables.

Source loads are retained: LC1 FX=+0.5 at tips7/14; LC2 FZ=+0.5 at both tips and **−0.5 at root joint8**; LC3 FY=+0.5 at both tips; LC4 FY=−5/+5 at7/14; LC5 FX=−5/+5; LC6 MZ=+0.5/+0.5. The negative LC2 load acts on a free UZ DOF and is not discarded as a reaction.

| LC | Observable at physical joints7/14 | Independent | Units | CSI SAP2000 A | CSI SAP2000 E |
|---|---|---|---|---|---|
| 1 | signed mean UX | 3e-05 | in | 3e-05 | 3e-05 |
| 2 | signed mean UZ | 0.1081 | in | 0.1072 | 0.0066 |
| 3 | signed mean UY | 0.4321 | in | 0.4321 | 0.4298 |
| 4 | mean absolute UY | 0.00341 | in | 0.00224 | 0.00458 |
| 5 | mean absolute UX | 0.0009 | in | 0.000899 | 5.5e-05 |
| 6 | signed mean RZ | 0.036 | rad | 0.036 | 0.036 |


CSI printed SAP2000 results use its source linear shell elements; they are not native quadratic outputs and are not used as independent denominators. Numerical errors below are 100×abs(solver−independent)/abs(independent). Fixed tip coordinates are (6,0,0) and (6,0,0.2) in, GRID7/14, verified within 1e-12 in. All GRID CD=0. LC4/5 use mean absolute values so opposite tip signs do not cancel; no arbitrary maximum-response probe is substituted.

## Boundary interpretation: source versus quadratic midpoint diagnostic

The source support is **not a full six-DOF clamp**. Joint1 has SPC1236; joint8 SPC126, with upper-root UZ deliberately free for the Poisson treatment. Source-corner baseline retains precisely those assigned joints and adds no midside constraints.

The second, explicitly labelled diagnostic adds **SPC126 only at (0,0,0.1)**: GRID18 in A/Q8 and GRID25 in E/T6. Those are the DOFs fixed at both edge endpoints. Midpoint UZ, RX and RY remain free; no SPC123456 is imposed. This is an additional support assumption for the quadratic boundary trace, not an exact copy of the source joint list. It is not selected merely because a reference error is smaller.

Both boundary variants keep the identical mesh, material, element formulation and loads. Baseline DAT and outputs remain intact; the variant lives separately. No formulation stiffness was edited.

## Single-mesh bar comparisons

![A/Q8 fixed-mesh boundary comparison](figures/2_002_A_Q8_fixed_mesh_BC_errors.png)

![E/T6 fixed-mesh boundary comparison](figures/2_002_E_T6_fixed_mesh_BC_errors.png)

## Nastran: boundary effect

| Model | LC | Observable | Source corner BC | + midpoint SPC126 | Baseline error (%) | Midpoint error (%) |
|---|---|---|---|---|---|---|
| A | 1 | UX | 3.03520756e-05 | 3.01724831e-05 | 1.17359 | 0.574944 |
| A | 2 | UZ | 0.108060114 | 0.108060114 | 0.0368972 | 0.0368972 |
| A | 3 | UY | 0.447650611 | 0.428939611 | 3.59885 | 0.731402 |
| A | 4 | abs_UY | 0.00304914406 | 0.00304914406 | 10.5823 | 10.5823 |
| A | 5 | abs_UX | 0.000899999985 | 0.000899999985 | 1.71777e-06 | 1.71777e-06 |
| A | 6 | RZ | 0.0368562452 | 0.0363403894 | 2.37846 | 0.945526 |
| E | 1 | UX | 3.06351394e-05 | 3.03175711e-05 | 2.11713 | 1.05857 |
| E | 2 | UZ | 0.107494935 | 0.107494935 | 0.559727 | 0.559727 |
| E | 3 | UY | 0.440287083 | 0.427932769 | 1.89472 | 0.964414 |
| E | 4 | abs_UY | 0.00334176095 | 0.00334176095 | 2.00115 | 2.00115 |
| E | 5 | abs_UX | 0.000899999985 | 0.000899999985 | 1.71777e-06 | 1.71777e-06 |
| E | 6 | RZ | 0.0365954377 | 0.036247056 | 1.65399 | 0.686267 |


Adding midpoint SPC126 changes LC1/3/6 but leaves LC2/4/5 unchanged in these runs. For Nastran A, LC3 error falls from3.60% to0.73% and LC6 from2.38% to0.95%; LC4 remains10.58%. For Nastran E, LC3 goes from1.89% to0.96% and LC6 from1.65% to0.69%; LC2 remains0.560%. Midpoint BC therefore matters, but does not account for every discrepancy and does not warrant silently replacing the source baseline.

## All native quadratic results

| Model | Formulation | LC | Units | Baseline value | Baseline error (%) | Midpoint126 value | Midpoint126 error (%) |
|---|---|---|---|---|---|---|---|
| A | NASTRAN | 1 | in | 3.03520756e-05 | 1.17359 | 3.01724831e-05 | 0.574944 |
| A | NASTRAN | 2 | in | 0.108060114 | 0.0368972 | 0.108060114 | 0.0368972 |
| A | NASTRAN | 3 | in | 0.447650611 | 3.59885 | 0.428939611 | 0.731402 |
| A | NASTRAN | 4 | in | 0.00304914406 | 10.5823 | 0.00304914406 | 10.5823 |
| A | NASTRAN | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | NASTRAN | 6 | rad | 0.0368562452 | 2.37846 | 0.0363403894 | 0.945526 |
| A | SIMOQ8 | 1 | in | 3.03065281e-05 | 1.02176 | 3.01533728e-05 | 0.511243 |
| A | SIMOQ8 | 2 | in | 0.107539706 | 0.518311 | 0.107539706 | 0.518311 |
| A | SIMOQ8 | 3 | in | 0.433270723 | 0.270938 | 0.428156883 | 0.912547 |
| A | SIMOQ8 | 4 | in | 0.0030258426 | 11.2656 | 0.0030258426 | 11.2656 |
| A | SIMOQ8 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | SIMOQ8 | 6 | rad | 0.0362165868 | 0.60163 | 0.0360687077 | 0.190855 |
| A | MITC8 | 1 | in | 3.03065281e-05 | 1.02176 | 3.01533746e-05 | 0.511249 |
| A | MITC8 | 2 | in | 0.107539706 | 0.518311 | 0.107539706 | 0.518311 |
| A | MITC8 | 3 | in | 0.436072707 | 0.919395 | 0.428902864 | 0.739906 |
| A | MITC8 | 4 | in | 0.00303439144 | 11.0149 | 0.00303439144 | 11.0149 |
| A | MITC8 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | MITC8 | 6 | rad | 0.0362164974 | 0.601382 | 0.0360195935 | 0.0544263 |
| A | ANS8BDG6 | 1 | in | 3.03065281e-05 | 1.02176 | 3.01533746e-05 | 0.511249 |
| A | ANS8BDG6 | 2 | in | 0.107539706 | 0.518311 | 0.107539706 | 0.518311 |
| A | ANS8BDG6 | 3 | in | 0.436072707 | 0.919395 | 0.428902864 | 0.739906 |
| A | ANS8BDG6 | 4 | in | 0.00303439144 | 11.0149 | 0.00303439144 | 11.0149 |
| A | ANS8BDG6 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | ANS8BDG6 | 6 | rad | 0.0362164974 | 0.601382 | 0.0360195935 | 0.0544263 |
| A | HBQ8 | 1 | in | 3.03065281e-05 | 1.02176 | 3.01533746e-05 | 0.511249 |
| A | HBQ8 | 2 | in | 0.107539706 | 0.518311 | 0.107539706 | 0.518311 |
| A | HBQ8 | 3 | in | 0.436072707 | 0.919395 | 0.428902864 | 0.739906 |
| A | HBQ8 | 4 | in | 0.00303439144 | 11.0149 | 0.00303439144 | 11.0149 |
| A | HBQ8 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | HBQ8 | 6 | rad | 0.0362164974 | 0.601382 | 0.0360195935 | 0.0544263 |
| A | MACQ8D | 1 | in | 3.03065281e-05 | 1.02176 | 3.01533746e-05 | 0.511249 |
| A | MACQ8D | 2 | in | 0.107539706 | 0.518311 | 0.107539706 | 0.518311 |
| A | MACQ8D | 3 | in | 0.433254719 | 0.267234 | 0.426322073 | 1.33717 |
| A | MACQ8D | 4 | in | 0.0030258426 | 11.2656 | 0.0030258426 | 11.2656 |
| A | MACQ8D | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| A | MACQ8D | 6 | rad | 0.0362164751 | 0.60132 | 0.0360176414 | 0.049004 |
| E | NASTRAN | 1 | in | 3.06351394e-05 | 2.11713 | 3.03175711e-05 | 1.05857 |
| E | NASTRAN | 2 | in | 0.107494935 | 0.559727 | 0.107494935 | 0.559727 |
| E | NASTRAN | 3 | in | 0.440287083 | 1.89472 | 0.427932769 | 0.964414 |
| E | NASTRAN | 4 | in | 0.00334176095 | 2.00115 | 0.00334176095 | 2.00115 |
| E | NASTRAN | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| E | NASTRAN | 6 | rad | 0.0365954377 | 1.65399 | 0.036247056 | 0.686267 |
| E | SIMOT6 | 1 | in | 3.06341535e-05 | 2.11385 | 3.03170764e-05 | 1.05692 |
| E | SIMOT6 | 2 | in | 0.107494853 | 0.559803 | 0.107494853 | 0.559803 |
| E | SIMOT6 | 3 | in | 0.440287083 | 1.89472 | 0.427932769 | 0.964414 |
| E | SIMOT6 | 4 | in | 0.00334176095 | 2.00115 | 0.00334176095 | 2.00115 |
| E | SIMOT6 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| E | SIMOT6 | 6 | rad | 0.0365954377 | 1.65399 | 0.036247056 | 0.686267 |
| E | MITC6 | 1 | in | 3.06341535e-05 | 2.11385 | 3.03170764e-05 | 1.05692 |
| E | MITC6 | 2 | in | 0.107494853 | 0.559803 | 0.107494853 | 0.559803 |
| E | MITC6 | 3 | in | 0.442402065 | 2.38419 | 0.431028724 | 0.247923 |
| E | MITC6 | 4 | in | 0.00334879966 | 1.79473 | 0.00334879966 | 1.79473 |
| E | MITC6 | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| E | MITC6 | 6 | rad | 0.0365524553 | 1.5346 | 0.0362315588 | 0.643219 |
| E | MH6T | 1 | in | 3.06341535e-05 | 2.11385 | 3.03170764e-05 | 1.05692 |
| E | MH6T | 2 | in | 0.107494853 | 0.559803 | 0.107494853 | 0.559803 |
| E | MH6T | 3 | in | 0.443781018 | 2.70331 | 0.429846406 | 0.521545 |
| E | MH6T | 4 | in | 0.0032863114 | 3.62723 | 0.0032863114 | 3.62723 |
| E | MH6T | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| E | MH6T | 6 | rad | 0.0367502645 | 2.08407 | 0.0363349169 | 0.930325 |
| E | REZAIEE | 1 | in | 3.06341535e-05 | 2.11385 | 3.03170764e-05 | 1.05692 |
| E | REZAIEE | 2 | in | 0.107494853 | 0.559803 | 0.107494853 | 0.559803 |
| E | REZAIEE | 3 | in | 0.442402065 | 2.38419 | 0.431028724 | 0.247923 |
| E | REZAIEE | 4 | in | 0.00334879966 | 1.79473 | 0.00334879966 | 1.79473 |
| E | REZAIEE | 5 | in | 0.000899999985 | 1.71777e-06 | 0.000899999985 | 1.71777e-06 |
| E | REZAIEE | 6 | rad | 0.0365524553 | 1.5346 | 0.0362315588 | 0.643219 |


## Solver status and limitations

Both variants have 11 readable OP2 files and six subcases each: 264 signed tip observations, 132 aggregated measurements, and108 Nastran/MYSTRAN quantity pairs in total. Every GRID displacement is finite; explicit SPC residuals are zero; hashes match the prepared DAT; physical signatures match across solvers within each boundary variant. No failed output is converted to zero.

Nastran uses BAILOUT=−1 under the earlier user instruction. All four Nastran runs retain warning4698 (matrix/factor conditioning or pivot ratios); values are warning-qualified, not an unconditional PASS. AUTOSPC is enabled explicitly on both solvers, so source-assigned restraints and numerical automatic-DOF treatment must not be confused. Successful termination does not establish physical accuracy. MYSTRAN pilot executions terminated normally.

All numerical results come from native OP2; F06 is used only for diagnostics. The two extra Nastran boundary variants were executed directly using the existing run_nastran.bat in their separate directories, successfully without a license failure. Future Nastran work should first try that configured batch runner before asking the user to run it.

Earlier regular nx1…64 curves used a different mesh family and are not the current 2-002 comparison. They remain historical research evidence, not an active convergence requirement. No B/C/D results are relabelled as A or E.

## Reference and reproducibility

[CSI Example 2-002](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-002.pdf), supplied pages and the two user-supplied thick S2K files. S2K engineering nodes, connectivity, material, force/moment vectors and case scales were audited. No S2K program/license metadata is redistributed. Original files are preserved; internal Python and raw solver outputs remain local.

## Input archives

[Fixed-mesh A/Q8 and E/T6 DAT](DAT/CSI_2_002_A_Q8_E_T6_FIXED_MESH_DAT_20261020.zip): 22 DAT, source-corner and midpoint-SPC126 variants separated. No Python or solver output.

[Benchmark index](README.md) · [All comparisons](comparison_all.md)
