! --- codex_mod\bug_rbe2_dy_dz ---
! RBE2 DELTA_0 bug: dy and dz were swapped in the first row/column pair.
! Pre-patch path read RGRID(*,3) where this benchmark needed RGRID(*,2),
! and the paired term read RGRID(*,2) where this benchmark needed RGRID(*,3).

! pre-patch code
! DELTA_0(1,2) =  (RGRID(GRID_ID_ROW_NUM_D,3) - RGRID(GRID_ID_ROW_NUM_I,3))
! DELTA_0(1,3) = -(RGRID(GRID_ID_ROW_NUM_D,2) - RGRID(GRID_ID_ROW_NUM_I,2))

! post-patch code
      DELTA_0(1,2) =  (RGRID(GRID_ID_ROW_NUM_D,2) - RGRID(GRID_ID_ROW_NUM_I,2))
      DELTA_0(1,3) = -(RGRID(GRID_ID_ROW_NUM_D,3) - RGRID(GRID_ID_ROW_NUM_I,3))
      DELTA_0(2,1) = -DELTA_0(1,2)
      DELTA_0(2,3) =  (RGRID(GRID_ID_ROW_NUM_D,1) - RGRID(GRID_ID_ROW_NUM_I,1))
      DELTA_0(3,1) = -DELTA_0(1,3)
      DELTA_0(3,2) = -DELTA_0(2,3)
