# RCM BANDED Add-on (rcm_add)

## Quick run modes
- Banded only:
  - `PARAM,SOLLIB,BANDED`
  - `PARAM,GRIDSEQ,INPUT`
  - `PARAM,BANDEDOPT,N`
- RCM only (no BANDIT):
  - `PARAM,SOLLIB,BANDED`
  - `PARAM,GRIDSEQ,RCM`
  - `PARAM,BANDEDOPT,N`
- Add-on auto gate:
  - `PARAM,SOLLIB,BANDED`
  - `PARAM,GRIDSEQ,INPUT`
  - `PARAM,BANDEDOPT,Y`

## Supporting files in this snapshot
- `mystran/Source/Modules/PARAMS.f90`
- `mystran/Source/LK1/L1A-BD/BD_PARAM.f90`
- `mystran/dev_docs/v18_banded_reorder_strategy.md`
- `mystran/examples/*.dat`

## Snippet locations (begin/end markers)
### 1) `mystran/Source/Modules/PARAMS.f90`
```fortran
! !--- RCM BANDED ADD-ON --- begin!
      CHARACTER(  8*BYTE)      :: GRIDSEQ        = 'INPUT   '! Method for sequencing grids:
!                                                              BANDIT for bandit auto grid swquencing
!                                                              GRID for grid numerical order
!                                                              INPUT for grid input order
!                                                              RCM for reserved/add-on RCM path (currently mapped to INPUT flow)
      CHARACTER(  1*BYTE)      :: SEQQUIT        =    'N'    !*'Y', 'N' indicator to stop processing if G.P. auto sequencing failed
      CHARACTER(  1*BYTE)      :: SEQPRT         =    'N'    !*'Y', 'N' indicator to print SEQGP card images from bandit
! !--- RCM BANDED ADD-ON --- end!
```

### 2) `mystran/Source/LK1/L1A-BD/BD_PARAM.f90` (BANDEDOPT parse)
```fortran
! !--- RCM BANDED ADD-ON --- begin!
      ELSE IF (JCARD(2)(1:8) == 'BANDEDOP') THEN
         PARNAM = 'BANDEDOPT'
         CALL YES_NO_CHECK(CARD, JCARD, CHRPARM, PARNAM, BANDEDOPT)
! !--- RCM BANDED ADD-ON --- end!
```

### 3) `mystran/Source/LK1/L1A-BD/BD_PARAM.f90` (GRIDSEQ accepts RCM)
```fortran
! !--- RCM BANDED ADD-ON --- begin!
      ELSE IF (JCARD(2)(1:8) == 'GRIDSEQ ') THEN
         ...
            ELSE IF (CHRPARM(1:3) == 'RCM') THEN
               GRIDSEQ = 'RCM     '
         ...
! !--- RCM BANDED ADD-ON --- end!
```

### 4) `mystran/Source/USE_IFs/SEQ_PROC_USE_IFs.f90` (new interface import)
```fortran
      USE OURTIM_Interface
      USE GET_ARRAY_ROW_NUM_Interface
      USE GET_ELGP_Interface
      USE OUTA_HERE_Interface
```

### 5) `mystran/Source/LK1/L1B/SEQ_PROC.f90` (in-core RCM hook)
```fortran
      LOGICAL                         :: DO_RCM
      LOGICAL                         :: RCM_OK
...
      IF ((GRIDSEQ(1:3) == 'RCM') .OR. ((BANDEDOPT == 'Y') .AND. (SOLLIB == 'BANDED  ') .AND. (GRIDSEQ(1:6) /= 'BANDIT'))) THEN
         DO_RCM = .TRUE.
      ENDIF
...
      IF (DO_RCM) THEN
         CALL RCM_SEQ_PROC ( R_GSEQ, RCM_OK )
         IF (RCM_OK) THEN
            NSEQ = 0
            WRITE(ERR,102) NGRID, NSEQ_SAVE, GRIDSEQ, BANDEDOPT, SOLLIB
         ELSE
            WRITE(ERR,103) GRIDSEQ, BANDEDOPT, SOLLIB
         ENDIF
      ENDIF
```
