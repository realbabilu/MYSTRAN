# CBEAM Validation 04: 2D Portal Frame

Reference worked example:

- [Portal Frame Worked Example](https://www.valuedes.co.uk/portal-frame-worked-example.html)

Source problem statement:

- portal frame span `4 m`
- height `3 m`
- section `80 x 80 x 5 SHS`
- steel `S355`
- bases fixed
- horizontal load `10 kN` at the top corner

Input deck:

- [cbeam_portal_frame_2d.dat](E:/mystran17/mystran/Binaries/cbeam_portal_frame_2d.dat)

Output file:

- [cbeam_portal_frame_2d.F06](E:/mystran17/mystran/Binaries/cbeam_portal_frame_2d.F06)

## Reference values from the source page

The worked example reports:

- base column bending moment: about `8.86 kN m` hand, `8.87 kN m` FEA
- top corner column bending moment: about `6.14 kN m` hand, `6.13 kN m` FEA
- vertical base reaction: `3.07 kN` at each base
- horizontal base reaction: `5.00 kN` at each base

The source page also reports a principal stress comparison, but the current MYSTRAN `CBEAM` section-stress recovery is not yet a trustworthy benchmark for SHS stress comparison, so this validation focuses on:

- frame reactions
- end moments
- overall 2D frame response

## MYSTRAN setup notes

- model is constrained to 2D frame behavior in the `XY` plane
- top nodes keep only `UX, UY, RZ` active
- bases are fully fixed
- three `CBEAM` elements are used:
  - left column
  - top beam
  - right column

## Comparison

From [cbeam_portal_frame_2d.F06](E:/mystran17/mystran/Binaries/cbeam_portal_frame_2d.F06):

Base reactions:

- left base node `1001`
  - horizontal reaction `Fx = -4.998102 kN`
  - vertical reaction `Fy = -3.066363 kN`
  - base moment `Mz = 8.863288 kN m`
- right base node `1004`
  - horizontal reaction `Fx = -5.001898 kN`
  - vertical reaction `Fy = 3.066363 kN`
  - base moment `Mz = 8.871261 kN m`

Element engineering-force moments:

- left column base end moment: `8.863288 kN m`
- left column top corner moment: `6.131017 kN m`
- beam left end moment: `6.131017 kN m`
- beam right end moment: `6.134434 kN m`
- right column top corner moment: `6.134434 kN m`
- right column base end moment: `8.871261 kN m`

Reference comparison against the worked example:

- base horizontal reactions
  - reference: `5.00 kN` at each base
  - MYSTRAN: `4.998102 kN` and `5.001898 kN`
- base vertical reactions
  - reference: `3.07 kN` at each base
  - MYSTRAN: `3.066363 kN` at each base
- base moments
  - reference: `8.86` to `8.87 kN m`
  - MYSTRAN: `8.863288` and `8.871261 kN m`
- top corner moments
  - reference: `6.13` to `6.14 kN m`
  - MYSTRAN: `6.131017` and `6.134434 kN m`

Small numerical differences are at rounding level only.

## Conclusion

This check is intended to validate whether the current `CBEAM` implementation behaves correctly for a small 2D frame, not just isolated single-member tests.

For this portal-frame benchmark, the current `CBEAM` matches the published worked-example reaction forces and end moments very closely.

That is a strong validation that:

- the element axial-bending coupling in a frame is behaving correctly
- member-to-member joint transfer is working correctly
- the 2D frame behavior from the 3D `CBEAM` element is working as intended when planar constraints are applied
