# CTRIA3 MITC3+: Lee interior tying and consistent spin

Owner-authorized default update, 25 October 2026. Syntax: `PARAM,TRIA3TYP,MITC3+`. Historical diagnostic operators remain accessible only through `DEBUG,242,2`; no automatic fallback occurs. F06 identifies the selected branch.

For a right-handed facet triad, use nodal degrees of freedom `[u,v,w,rx,ry,rz]` and two internal rotational bubbles. The CST membrane is retained. With nodal gradient matrix G, the compatible curvature contains theta_y,x, -theta_x,y, and theta_y,y-theta_x,x. The physical bubble gradient uses Jinv, matching the supplied Python reference on skew triangles.

1. Bubble function: `b=27*r*s*(1-r-s)`.
2. Interior shear tying: six Lee sample points, with project tying distance `d=1e-4`.
3. Uncondensed energy: `K20=integral(Bm^T*A*Bm + Bb^T*D*Bb + Bs^T*T*Bs)dA + Kspin`, with three symmetric Gauss points.
4. Continuum spin: `omega_z=0.5*(v,x-u,y)`.
5. Project drilling stabilization: `Kspin=0.01*G*h*integral((omega_z-rz)^2)dA` in the energy convention used by the supplied Python class. The coefficient is project-specific, not a universal paper parameter.
6. Bubble condensation: `K18=Kcc-Kcb*inv(Kbb)*Kbc`; recovery uses the same bubble back-substitution. The old DKMT18 recovery substitution is removed.

Constant face pressure and the condensed gradient load are retained. Membrane thermal assembly precedes the bending-gradient contribution so its initialization cannot erase that load. Existing thermal-recovery DT(4) correction applies to triangular shells.

## Buckling configuration

For the tested six-mode 2-016 shear benchmark explicitly set:

```nastran
PARAM,TRIA3TYP,MITC3+
PARAM,DARPACK,4
EIGRL,1,,,6
```

F06 reports the configuration requirement and current DARPACK. No default, backend, eigen residual tolerance or requested mode count is changed. DARPACK2 mode6 residual=4.8927350789982573e-6 (known fail). DARPACK4 mode6 residual=1.3077452361049327e-10 (six pass); DARPACK4/6 modal MAC=1.0 to reported precision.

## Validation and limits

Eight combined Lee/Nguyen skew/rotated/cyclic/reversed matrix tests agree with supplied Python to <5e-16; six rigid modes retained. Lee static parity verified for all six 2-002 n24 load cases, supplied pressure-cylinder 2-010, buckling preload 2-017 and original 128-triangle hemisphere, including GRID CD transformations. Native71-case run completes with explicit DARPACK4 on2-016.

Limits remain: constrained shear patch response about1.13408, CST membrane under-response on2-002 LC2/LC5, and documented trade-offs on roof/coarse meshes. Completion is not universal accuracy acceptance. See [audit](MITC3_Lee_KMITC3_AUDIT_20261025.md).
