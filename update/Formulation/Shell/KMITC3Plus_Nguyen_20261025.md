# CTRIAR KMITC3+: Nguyen2024 research reconstruction

Status: **EXPERIMENTAL / HOLD**, explicit selection `PARAM,TRIARTYP,KMITC3+`. This replaces the retired MITC3+HB option by owner request. Old MITC3+HB aliases cause an explicit input error; they are not silently mapped to a different formulation. DKMT18 remains the CTRIAR default.

This is a literal native port of the supplied complete Python research class `MITC3Plus_Nguyen2024_CompositeCurvature`, not the historical Allman-only experiment and not a claim of paper-author code equivalence. The source paper is Son H. Nguyen, Computers & Structures295 (2024),107272, DOI10.1016/j.compstruc.2024.107272. The supplied paper's Eqs24-32 were inspected; equivalence of every subdomain integration detail remains a research question.

## Discrete formulation

1. Right-handed flat facet with18 external nodal degrees of freedom and two internal bubble rotations; nodal director smoothing is not consumed.
2. Modified Allman membrane: `Bm=Bm_CST-(2/3)*H*Gm`, integrated at three triangle Gauss points.
3. Independent true-drill rank stabilizer: `K_Add=G*h*A*c^T*c`, with `c*u=mean(rz)-0.5*(v,x-u,y)`. No additional Lee weak-spin penalty is superposed.
4. Compatible rotational bubble curvature and six-point Lee tied transverse shear.
5. P1 curvature projection: `f=[1,r,s]`, `F=integral(f*f^T)dA`, `Q=integral(f*Bb)dA`, and `Bb_projected=f^T*inv(F)*Q`. With the supplied three-point projection rule this is exactly linear interpolation through the three Gauss samples; the implementation uses equivalent interpolation weights rather than explicitly inverting F.
6. Projected curvature stiffness integrated at the centroid of each of three barycenter-to-edge subtriangles, exactly as supplied Python reconstruction.
7. Condense the two bubbles with the same Schur complement for stiffness and recovery; pressure, consistent Allman mean-thermal loads, gradient loads and local stress operators are assembled within this branch.
8. Initial-stress matrix uses the Allman membrane resultants at the three Gauss points and all three translational slope blocks. No global eigensolver changes.

## Verified

Native/Python K agrees to <5e-16 on skew, rotated, cyclic and reversed triangles; six rigid modes. Stress parity relative error<7e-7 using printed native displacement/stress tables. Mixed CTRIA3 Lee + CTRIAR KMITC3+ gives correct matrix and element IDs in F06 and OP2. OP2 currently uses the existing CTRIA3 table layout for both triangular card families; input element IDs distinguish them.

All71 native benchmark cases complete in the qualified run, including2-002,2-004,2-005,2-007,2-010,2-017. Direct Python/native static parity on exact source models is within1.7e-7 relative printed-displacement error. On coarse2-017, elastic K agrees to1.38e-16, geometric Kg to4.49e-14, and preload/first critical factor agree.

## HOLD limitations

- Original128-triangle hemisphere node1 result/reference=0.25566649, versus Lee0.99192798. This severe coarse-mesh regression is retained openly. Four-child benchmark meshes have different results and are not interchangeable with the original DAT.
- 2-002 torsion n24 ratio1.12966745 versus Lee0.98445718.
- Higher buckling modes are not certified. Coarse2-017 native six reported factors are131.5733,685.5907,686.1437,686.2421,713.3178,1062.787. Direct independent generalized eigenanalysis of matched K/Kg gives131.5733113,1106.8927132,2739.7882026,3756.7243028,3756.7243028,3958.0125874. First agrees; modes2-6 do not, despite small native reported residuals. Reduced KLL/KLLD print audits also match the Python restriction. Cause remains open; no shared solver modification is made.
- Native2-013 restrained thermal stress extrema differ from Lee/control and require a separate Allman/drilling-boundary recovery audit. No blanket thermal-stress certification.
- Laminates, arbitrary nodal thickness, exact curved-shell edge stress, and full paper quadrature equivalence are not certified.

This branch is provided for explicit research use; successful solver termination does not remove these HOLD conditions. See [audit](MITC3_Lee_KMITC3_AUDIT_20261025.md).
