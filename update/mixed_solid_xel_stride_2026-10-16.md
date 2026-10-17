# Mixed-solid XEL leading dimension correction

Observed in the user's 5-002 eight-solid model after correcting its blank PSOLID integration order: CHEXA8 EAS9 helper reported invalid geometry in recovery. Full XEL is allocated for MELGP=20; passing it to an explicit x(8,3) dummy reinterprets the coordinate stride incorrectly. Homogeneous tests had ELGP=MELGP and did not expose this.

Production HEX8_EAS9_STATE calls now pass XEL(1:8,1:3), PENTA6_ANSEAS_STATE passes XEL(1:6,1:3), and Bedrosian recovery/Kg passes XEL(1:13,1:3). Fortran creates a correctly packed input temporary as needed. Shape functions, material properties, enrichment and global DOFs are unchanged.

The production-linked HEXA8, CPENTA6 ANSEAS and CPYRAM13 probes now allocate XEL with 20 rows while populating their own element node counts, exercising the actual stride failure. Their direct helper calls also use explicit slices. The existing reference comparisons and full solver regressions are retained. 5-002_psolid_fixed.dat finishes all three subcases with no *ERROR; mixed 5-001 eight-patch interior displacement comparison passes all eight topologies at 1e-7 tolerance.

Separate input fix: 5-002_psolid_fixed.dat assigns each topology a distinct PSOLID with IN=2 or IN=3 as required. No original user DAT was overwritten. Source helper correction and input correction are independent.
