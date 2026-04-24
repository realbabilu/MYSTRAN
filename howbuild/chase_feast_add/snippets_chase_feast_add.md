# CHASE/FEAST Change Snippets

Below are core change snippets with explicit markers.

## `mystran/CMakeLists.txt`

```cmake
! !--- CHASE and FEAST --- begin!
# FEAST/ChASE wrapper sources use preprocessor guards for optional native hooks.
set_source_files_properties(
  "${CMAKE_SOURCE_DIR}/LK4/EIG_LANCZOS_FEAST.f90"
  "${CMAKE_SOURCE_DIR}/LK4/EIG_LANCZOS_CHASE.f90"
  PROPERTIES COMPILE_OPTIONS "-cpp"
)
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/Modules/PARAMS.f90`

```fortran
! !--- CHASE and FEAST --- begin!
CHARACTER(8*BYTE)               :: LANCMETH         = 'ARPACK  '
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/LK1/L1A-BD/BD_PARAM.f90`

```fortran
! !--- CHASE and FEAST --- begin!
ELSE IF (JCARD(2)(1:8) == 'LANCMETH ') THEN
   IF (JCARD(3)(1:8) == 'ARPACK  ') THEN
      LANCMETH = 'ARPACK  '
   ELSE IF (JCARD(3)(1:8) == 'FEAST   ') THEN
      LANCMETH = 'FEAST   '
   ELSE IF (JCARD(3)(1:8) == 'CHASE   ') THEN
      LANCMETH = 'CHASE   '
   ENDIF
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/LK4/LINK4.f90`

```fortran
! !--- CHASE and FEAST --- begin!
IF (((LANCMETH(1:6) == 'FEAST ') .OR. (LANCMETH(1:6) == 'CHASE ')) .AND. (EIG_METH(1:7) /= 'LANCZOS')) THEN
   WRITE(ERR,4910) LANCMETH, EIG_METH
ELSE
   IF (LANCMETH(1:6) == 'FEAST ') THEN
      CALL EIG_LANCZOS_FEAST
   ELSE IF (LANCMETH(1:6) == 'CHASE ') THEN
      CALL EIG_LANCZOS_CHASE
   ELSE
      CALL EIG_LANCZOS_ARPACK
   ENDIF
ENDIF
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/LK4/EIG_LANCZOS_CHASE.f90`

```fortran
! !--- CHASE and FEAST --- begin!
SUBROUTINE EIG_LANCZOS_CHASE
...
CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE NATIVE (GENERALIZED SYMMETRIC SPARSE)')
...
CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - CHASE FALLBACK (ARPACK LANCZOS)')
CALL EIG_LANCZOS_ARPACK
END SUBROUTINE EIG_LANCZOS_CHASE
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/LK4/EIG_LANCZOS_FEAST.f90`

```fortran
! !--- CHASE and FEAST --- begin!
SUBROUTINE EIG_LANCZOS_FEAST
...
UPLO = 'F'
...
MREG_FLOOR = MAX(EPS1, ABS(AVG_MDIAG_POS)*1.0D-10)
...
CALL LINK_MESSAGE('SOLVE FOR EIGENVALS/VECTORS - FEAST FALLBACK (ARPACK LANCZOS)')
CALL EIG_LANCZOS_ARPACK
END SUBROUTINE EIG_LANCZOS_FEAST
! !--- CHASE and FEAST --- end!
```

## `mystran/Source/USE_IFs/*.f90` interfaces

```fortran
! !--- CHASE and FEAST --- begin!
USE EIG_LANCZOS_FEAST_USE_IFs
USE EIG_LANCZOS_CHASE_USE_IFs
! !--- CHASE and FEAST --- end!
```
