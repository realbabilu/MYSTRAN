# NEWSOLID Snippets

These snippets are the log markers requested for the NEWSOLID integration. Full copied files are in `files/`.

## PARAM Switch

```fortran
! --- newsolid_add begin ---
ELSE IF (CHRPARM == 'NEWSOLID') THEN
   SOLIDTYP = 'NEWSOLID'
ELSE IF (CHRPARM == 'EAS') THEN
   SOLIDTYP = 'NEWSOLID'
! --- newsolid_add end ---
```
## CHEXA8 / CPENTA6 Condensed EAS Routing

```fortran
! --- newsolid_add begin ---
USE_EAS9_NEWSOLID = ((SOLIDTYP == 'NEWSOLID') .AND. linear_element_guard)
IF (USE_EAS9_NEWSOLID) THEN
   ! Build standard strain modes, enhanced strain modes, and condense:
   ! KE = Kuu - Kua * inverse(Kaa) * Kau
END IF
! --- newsolid_add end ---
```

## CTETRA4 Smooth Assembly Routing

```fortran
! --- newsolid_add begin ---
IF (SOLIDTYP == 'NEWSOLID') THEN
   ! Assemble CTETRA4 smooth alpha 0.9 nodal-patch stiffness through
   ! CTETRA4S_SMOOTH_ASSEMBLY instead of changing legacy TETRA.f90.
END IF
! --- newsolid_add end ---
```

## CPYRA NEWSOLID Routing

```fortran
! --- newsolid_add begin ---
IF ((TYPE == 'PYRA5   ') .AND. (SOLIDTYP == 'NEWSOLID')) THEN
   CALL PYRA5_EAS54_STIFFNESS(...)
ELSE IF (TYPE == 'PYRA14  ') THEN
   CALL PYRA14_LIU_COMPOSITE_STIFFNESS(...)
END IF
! --- newsolid_add end ---
```
