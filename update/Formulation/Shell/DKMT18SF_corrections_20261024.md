# DKMT18 pressure correction and opt-in DKMT18SF research branch

24 October 2026. Applies to CTRIAR/DKMT18. The SF branch is experimental and is not enabled by default.

## Changes

Remove duplicate face-pressure assembly that doubled the nodal pressure force vector. Nodal point loads and the baseline elastic stiffness are unchanged by this correction.

Add an opt-in shear reconstruction using adjacent edge-tangent systems to recover nodal shear vectors, then barycentric interpolation at the existing three stiffness points. The SF expression is phi/(1+phi)*Au; the original Adelta inverse scaling is accounted for algebraically. Membrane, bending, fictitious stabilization, directors, mass, pressure, thermal and geometric-stiffness expressions are retained. Recovery calls use the centroid. No shared solver change is included.

```text
PARAM,TRIARTYP,DKMT18
DEBUG,241,1
```

DEBUG 241 enables SF. With it zero or absent, the original shear reconstruction remains. This is an element-local research switch, not an official new PARAM selector.

## Validation

The pressure correction matches Python nodal forces and leaves stiffness unchanged; finest 2-005 pressure errors are 0.5423% clamped and 0.0715% simple. The original 2-002 LC4 torsion accuracy failure is retained as a baseline; drilling energy is zero there.

For SF, exact-DAT native/Python KGG relative differences are 3.58e-14 at nx2 and 2.58e-16 at nx24. Flag-off nx2 retains baseline agreement (3.57e-14). All six load-case displacement vectors agree within 1.82e-7 relative, limited by F06 precision. On the original nx24/ny1 mesh, SF UY=0.002475087 in versus 0.00341 in reference; it does not cure every coarse mesh. Two-axis native refinement reaches UY=0.0033967025 in at 192x16 (12,288 triangles), error 0.38996%, matching Python within 1.27e-8 relative.

The opt-in 71-case run completed 70 configurations. The six-mode 2-016c shear case still fails residual mode6=4.249713e-6 against tolerance 1e-6; the prior pressure-corrected baseline had a similar failure. Failed output is not accepted as a valid result. No tolerance/backend changes were made to bypass it.

Accuracy is not uniformly improved. On 2-010 12x24, the retained native_fiber_principal_max_min error increases from 25.3538% to 67.6590%. Curved and mode metrics also change. Through-thickness temperature-gradient loading remains unresolved. These limits preclude default promotion.

Builds replace only the relevant object in the prior configured archive and are not full clean rebuilds. Local evidence is retained under investigations/dkmt18_joint and mitc3 dkmt18 upgrade/review_20261024/native_dkmt18sf in the development workspace. Binaries, run outputs and original benchmark publication values are not changed by this source publication.
