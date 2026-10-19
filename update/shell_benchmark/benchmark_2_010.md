# CSI2-010B — cylinder under internal pressure, fixed source mesh

20 October2026. **One source24×48 parent mesh, no refinement series.** Five Q8 and four T6 MYSTRAN formulations are compared with Nastran using native OP2, local joint-coordinate transforms and explicit pressure integration checks. Displacement agreement and all-point stress recovery are assessed separately.

## Source and physical model

User-supplied Example2-010b-thick.s2k. Radius60 in, height200 in, thickness1 in; E=29,000 kip/in², G=11,153.8461538462 kip/in², nu0.3. Source units are **kip–in**, not lb–in; pressure is1 kip/in². Source Shell-Thick; PSHELL MID1=MID2=MID3=1 with transverse shear factor5/6. No end-cap pressure or artificial axial force is added.

All1200 source vertex coordinates/IDs and1152 source quads are retained. Q8 enriches to1152 CQUAD8 /3552 nodes. T6 splits each quad into two CTRIA6, yielding2304 elements /4704 nodes. Shared midsides are arithmetic source-edge midpoints, preserving the48-sided polygon facets. They are not projected to a true cylindrical surface. Original source axial-coordinate spacing is preserved, including its stored rounding.

## Local restraints, including midsides

At every node, local axis1 is radial outward, axis2 is circumferential tangent, axis3 is global Z. Node-specific CORD2R is GRID CD, while geometry remains in basic CP=0. Source joint-local angles were checked against the coordinates.

Local **U2 is restrained at every original and new quadratic node**. At the base, local U3 is additionally restrained:96 base vertices/midsides per family. Local radial U1 remains free and rotations are not clamped. This extends the stated source symmetry/kinematic conditions to the quadratic edge trace. AUTOSPC is enabled explicitly; no BAILOUT or extra rotational restraint is used.

## Internal pressure and native load audit

Source bottom-face positive pressure is mapped to outward facet-normal PLOAD4, magnitude1 kip/in², explicit basic CID0 direction. Each source face is planar; both T6 children share its outward normal. No nodal-force rescaling is used.

For a flat Q8 facet under constant pressure, consistent corner weights are−pA/12 and midsides+pA/3. For T6, corner weights are zero and midsides+pA/3. Contributions are assembled across neighboring elements. Native OP2 OLOAD, rotated from GRID CD to basic, matches this expected nodal distribution: maximum difference below1.87e−6 kip for Q8 and9.33e−7 kip for T6. This verifies the pressure RHS for these models; it does not validate every curved PLOAD4 case or every recovery table.

## Independent targets and fixed observations

With p=1 kip/in², R=60 in, E=29000 kip/in², h=1 in and H=200 in:

- Ur=pR²/(Eh)=+0.124137931… in, printed CSI target **+0.12414 in**.
- UZ(top)=−pR nu H/(Eh)=−0.124137931… in, printed **−0.12414 in**.
- Circumferential hoop stress sigma11=pR/h=**60 kip/in² (ksi)**.

Primary radial observation: original GRID577 at(60,0,100) in. Primary top vertical observation: original GRID1153 at(60,0,200) in. Native OUG translations are in GRID CD; local U1 gives Ur and local U3 gives UZ. Basic vectors are also checked through the CD transform. No node-ID-only assumption, moving maximum, or nearest constrained corner is used.

All1200 original source vertices are additionally exported per model to check the 'any joint' radial statement, and all48 original top vertices are retained for vertical ranges. Quadratic chord-midpoint radii differ from60, so they are not silently pooled with the source-radius vertex statistic.

## Stress coordinate/recovery audit

Native OES local sigmaXX/sigmaYY/tauXY is transformed into the actual flat facet's circumferential tangent. Q8 uses its center derivative basis; T6 uses the first-corner-edge basis, with the element normal and second in-plane axis completing the frame. No local component is blindly equated with global X stress. Both fibers are retained.

The E1 center membrane indicator is the mean of its two native fiber hoop stresses. It is reported as a specific diagnostic, **not a substitute for the CSI 'anywhere' stress condition**. All native center/corner fiber values are separately inspected; their full min/max and worst independent error are reported. No GPSTRESS nodal averaging, extrapolated smoothing or Python analytical stress is substituted for native OES.

## Geometry qualification

The source48-sided polygon has facet-mid radius a=60cos(3.75°)=59.8715354 in, a0.2141077% radial-location difference from the ideal cylinder. For the polygon membrane response, Ur≈p a R/(Eh), UZ(top)≈−p a nu H/(Eh), and hoop stress≈p a/h. The target remains the ideal-cylinder independent value; it is not adjusted to force a match. The printed displacement target introduces a small additional rounding difference, so a0.2158% printed-target error is consistent with the approximately0.2141% geometric effect.

CSI SAP2000 ModelB thick reports Ur=0.12387 in, UZ(top)=−0.12387 in and hoop stress59.87 ksi. These are published solver evidence, not the independent denominators and not native Nastran/MYSTRAN outputs.

## Single-mesh error bars

The fourth panel deliberately shows the worst native fiber/corner hoop error. A small center membrane error does not establish uniform stress throughout the model.

![Q8 fixed-mesh errors](figures/2_010B_Q8_fixed_mesh_errors.png)

![T6 fixed-mesh errors](figures/2_010B_T6_fixed_mesh_errors.png)

## Fixed-point displacement results

| Family | Solver/formulation | Ur (in) | Ur error (%) | TopUZ (in) | UZ error (%) |
|---|---|---|---|---|---|
| Q8 | NASTRAN | 0.123872139 | 0.215774 | -0.123872139 | 0.215774 |
| Q8 | ANS8BDG6 | 0.123872139 | 0.215774 | -0.123872139 | 0.215774 |
| Q8 | HBQ8 | 0.123872139 | 0.215774 | -0.123872139 | 0.215774 |
| Q8 | MACQ8D | 0.123872139 | 0.215774 | -0.123872139 | 0.215774 |
| Q8 | MITC8 | 0.123872139 | 0.215774 | -0.123872139 | 0.215774 |
| Q8 | SIMOQ8 | 0.122419544 | 1.3859 | -0.123820923 | 0.25703 |
| T6 | NASTRAN | 0.123872124 | 0.215786 | -0.123873614 | 0.214585 |
| T6 | MH6T | 0.123872101 | 0.215804 | -0.123859234 | 0.226169 |
| T6 | MITC6 | 0.123857528 | 0.227543 | -0.123873912 | 0.214345 |
| T6 | REZAIEE | 0.123857528 | 0.227543 | -0.123873912 | 0.214345 |
| T6 | SIMOT6 | 0.123872139 | 0.215774 | -0.123875 | 0.213469 |


## Native hoop-stress recovery ranges

| Family | Solver/formulation | E1 center membrane (ksi) | Center error (%) | All center/corner fiber min (ksi) | Max (ksi) | Worst fiber error (%) |
|---|---|---|---|---|---|---|
| Q8 | NASTRAN | 59.8715363 | 0.214106 | 59.8715363 | 59.8715363 | 0.214106 |
| Q8 | ANS8BDG6 | 59.8715363 | 0.214106 | 59.8715363 | 59.8715363 | 0.214106 |
| Q8 | HBQ8 | 59.8715363 | 0.214106 | 59.8715363 | 59.8715363 | 0.214106 |
| Q8 | MACQ8D | 59.8715363 | 0.214106 | 59.8715363 | 59.8715363 | 0.214106 |
| Q8 | MITC8 | 59.8715363 | 0.214106 | 59.8715363 | 59.8715363 | 0.214106 |
| Q8 | SIMOQ8 | 60.3397312 | 0.566219 | 58.5656548 | 60.6623039 | 2.39058 |
| T6 | NASTRAN | 59.861269 | 0.231218 | 46.7889824 | 73.0469513 | 22.0184 |
| T6 | MH6T | 59.8220234 | 0.296628 | 38.1436844 | 80.6615753 | 36.4272 |
| T6 | MITC6 | 59.8266716 | 0.288881 | 35.4798012 | 84.1514206 | 40.867 |
| T6 | REZAIEE | 59.8266716 | 0.288881 | 35.4798012 | 84.1514206 | 40.867 |
| T6 | SIMOT6 | 59.8678684 | 0.220219 | 47.2344666 | 72.6783981 | 21.2759 |


## Original source-vertex displacement ranges

| Family | Formulation | Ur minimum (in) | Ur maximum (in) | TopUZ minimum (in) | TopUZ maximum (in) |
|---|---|---|---|---|---|
| Q8 | NASTRAN | 0.123872139 | 0.123872139 | -0.123872139 | -0.123872139 |
| Q8 | ANS8BDG6 | 0.123872139 | 0.123872139 | -0.123872139 | -0.123872139 |
| Q8 | HBQ8 | 0.123872139 | 0.123872139 | -0.123872139 | -0.123872139 |
| Q8 | MACQ8D | 0.123872139 | 0.123872139 | -0.123872139 | -0.123872139 |
| Q8 | MITC8 | 0.123872139 | 0.123872139 | -0.123872139 | -0.123872139 |
| Q8 | SIMOQ8 | 0.122173443 | 0.122424707 | -0.123820923 | -0.123820923 |
| T6 | NASTRAN | 0.123804525 | 0.12396317 | -0.123873614 | -0.123873614 |
| T6 | MH6T | 0.122931987 | 0.124036275 | -0.123859234 | -0.123859234 |
| T6 | MITC6 | 0.123771578 | 0.124042757 | -0.123873912 | -0.123873912 |
| T6 | REZAIEE | 0.123771578 | 0.124042757 | -0.123873912 | -0.123873912 |
| T6 | SIMOT6 | 0.123854041 | 0.124034338 | -0.123875 | -0.123875 |


## Findings and limitations

Nastran Q8 radial and top displacement errors are about0.2158% against the printed ideal reference; its native hoop stress is essentially uniform59.87154 ksi (0.2141%). Q8 MITC8/ANS8BDG6/HBQ8/MACQ8D show the same polygon response to output precision. SIMOQ8 radial error is about1.386%, top error0.257%, with native hoop nonuniformity up to2.39% error; that difference remains visible.

T6 Nastran and MYSTRAN fixed-point displacements are close to the ideal reference at about0.21–0.23%, but native fiber/corner stress recovery is not uniform. Worst hoop errors are approximately22.02% for Nastran,21.28% SIMOT6,40.87% MITC6/REZAIEE, and36.43% MH6T. These values are not hidden by averaging the two fibers or choosing a favorable center. The 'sigma11 anywhere' condition is therefore not established for these T6 outputs. A separate recovery/edge-behavior study would be needed to identify the cause; no formulation or recovery code was changed to improve the table.

All11 runs terminate with readable, finite native OP2 and no numbered fatal/warning in the diagnostic scan. There are13200 original-vertex displacement observations; all native center/corner stress rows were inspected, while the local numeric release retains both centroid fibers and complete range/extremum metadata. All assigned local U2/baseU3 residuals are zero; input hashes and matched physical signatures are unchanged. OLOAD is verified node by node. Successful execution, displacement accuracy and pointwise stress uniformity are separate conclusions.

Both Nastran cases were run directly with the configured run_nastran.bat, without a license failure or BAILOUT. Original DAT/S2K and raw OP2 stay local. This is a one-mesh source-model test, not a convergence study or a global element ranking.

## Reference and reproducibility

[CSI Example2-010](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-010.pdf), supplied verification pages and user-supplied ModelB thick S2K. The source is kip–in and retains the documented pressure/local restraints. Source program/license metadata and internal Python are not published.

[Executed fixed-source Q8/T6 DAT](DAT/CSI_2_010B_Q8_T6_FIXED_SOURCE_DAT_20261020.zip):11 inputs, no solver outputs or Python. [Complete benchmark DAT archive](DAT/quadratic_shell_benchmark_DAT.zip).

[Benchmark index](README.md) · [All comparisons](comparison_all.md)
