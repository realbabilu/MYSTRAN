! #################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.
      SUBROUTINE MITC8_LEGACY ( OPT, INT_ELEM_ID )

! Calculates, or calls subr's to calculate, quadrilateral element matrices:

!  1) ME        = element mass matrix                  , if OPT(1) = 'Y'
!  2) PTE       = element thermal load vectors         , if OPT(2) = 'Y'
!  3) SEi, STEi = element stress data recovery matrices, if OPT(3) = 'Y'
!  4) KE        = element linea stiffness matrix       , if OPT(4) = 'Y'
!  5) PPE       = element pressure load matrix         , if OPT(5) = 'Y'
!  6) KED       = element differen stiff matrix calc   , if OPT(6) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_ORDER_GAUSS, MAX_STRESS_POINTS, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, DT, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, ELGP, ES, KE, EM, ET, BE1, BE2, BE3,       &
                                         PHI_SQ, FCONV, EPROP, SHELL_STR_ANGLE, ME, MASS_PER_UNIT_AREA, PPE, PRESS, PTE, TREF,   &
                                         XEB, XEL
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, FOUR
      USE PARAMS, ONLY                :  TSTM_DEF, COUPMASS
      USE MITC_STUF, ONLY             :  GP_RS

      USE MITC_INITIALIZE_Interface
      USE ORDER_GAUSS_Interface
      USE OUTA_HERE_Interface
      USE PLANE_COORD_TRANS_21_Interface
      USE MATL_TRANSFORM_MATRIX_Interface
      USE MATMULT_FFF_Interface
      USE MATMULT_FFF_T_Interface
      USE MITC_SHAPE_FUNCTIONS_Interface
      USE MATMULX_FFF_Interface
      USE MATMULX_FFF_T_Interface
      USE MITC_DETJ_Interface
      USE MITC8_B_Interface
      USE MITC8_CARTESIAN_LOCAL_BASIS_Interface
      USE MITC8_ELEMENT_CS_BASIS_Interface
      USE MITC_ELASTICITY_Interface

      USE QUADRATIC_SURFACE_PRESSURE_Interface
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MITC8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG), PARAMETER        :: IORD_IJ = 3       ! Integration order for stiffness matrix
      INTEGER(LONG), PARAMETER        :: IORD_K = 2        ! Integration order for stiffness matrix in thickness direction
      INTEGER(LONG), PARAMETER        :: IORD_STRESS_Q8 = 2! Gauss integration order for stress/strain recovery matrices
      INTEGER(LONG)                   :: I,J,K,L,M,JSUB    ! DO loop indices
      INTEGER(LONG)                   :: STR_PT_NUM        ! Stress recovery point number

      REAL(DOUBLE), PARAMETER         :: BETA_DRILL = 1.0D-4
      REAL(DOUBLE)                    :: HH_IJ(MAX_ORDER_GAUSS) ! Gauss weights for integration in in-layer directions
      REAL(DOUBLE)                    :: SS_IJ(MAX_ORDER_GAUSS) ! Gauss abscissa's for integration in in-layer directions
      REAL(DOUBLE)                    :: HH_K(MAX_ORDER_GAUSS)  ! Gauss weights for integration in thickness direction
      REAL(DOUBLE)                    :: SS_K(MAX_ORDER_GAUSS)  ! Gauss abscissa's for integration in thickness direction
      REAL(DOUBLE)                    :: R, S, T                ! Isoparametric coordinates of a point
      REAL(DOUBLE)                    :: BI(6,6*ELGP)      ! Strain-displ matrix for this element for one Gauss point
      REAL(DOUBLE)                    :: BI1(6,6*ELGP)     ! Strain-displ matrix for this element for one Gauss point bottom
      REAL(DOUBLE)                    :: BI2(6,6*ELGP)     ! Strain-displ matrix for this element for one Gauss point top
      REAL(DOUBLE)                    :: DUM1(6,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: DUM2(6*ELGP,6*ELGP)    ! Intermediate matrix
      REAL(DOUBLE)                    :: INTFAC            ! An integration factor (constant multiplier for the Gauss integration)
      REAL(DOUBLE)                    :: DETJ              ! Jacobian determinant
      REAL(DOUBLE)                    :: E(6,6)            ! Elasticity matrix in the material coordinate system.
      REAL(DOUBLE)                    :: EE(6,6)           ! Elasticity matrix in the cartesian local coordinate system.
      REAL(DOUBLE)                    :: BDRILL(1,6*ELGP)  ! Drilling strain-displacement matrix
      REAL(DOUBLE)                    :: LOCAL_BASIS(3,3)  ! Cartesian local basis
      REAL(DOUBLE)                    :: ELEMENT_BASIS(3,3)! Element coordinate system basis
      REAL(DOUBLE)                    :: M_1DOF(ELGP,ELGP) ! Consistent translational mass matrix with 1 DOF per node.
      REAL(DOUBLE)                    :: PSH(ELGP)         ! Shape functions
      REAL(DOUBLE)                    :: DPSHG(2,ELGP)     ! Shape function derivatives
      REAL(DOUBLE)                    :: DENSITY           ! Mass density
      REAL(DOUBLE)                    :: MASS_ELEM         ! Total translational element mass
      REAL(DOUBLE)                    :: MASS_NODE         ! Lumped translational mass per node
      REAL(DOUBLE)                    :: UNIT_PPE(6*ELGP)  ! Pressure load vector for unit pressure
      REAL(DOUBLE)                    :: UNIT_PTE(6*ELGP)  ! Thermal load vector for unit temperature change
      REAL(DOUBLE)                    :: DXDR(3), DXDS(3), SURF_VEC(3)
      REAL(DOUBLE)                    :: XL(3)
      REAL(DOUBLE)                    :: ZL(3)
      REAL(DOUBLE)                    :: XE(3)
      REAL(DOUBLE)                    :: CROSS_XLE(3)
      REAL(DOUBLE)                    :: GDRILL
      REAL(DOUBLE)                    :: CTE(6), THERMAL_STRAIN(6), TBAR, MATL_AXES_ROTATE, E3(6,6), T66(6,6), DUM66(6,6)
      REAL(DOUBLE)                    :: CLB(3,3), TRANSFORM(3,3)

! **********************************************************************************************************************************

! COORDINATE SYSTEMS
! ==================
!
! Cartesian local
!  e^_1, e^_2, e^_3 in Bathe but they are oriented differently there.
!  e^_3 is the midsurface normal which is the same as the director vector for MITC8 because it doesn't support SNORM.
!  Used for strain in the strain-displacement matrix and the material elasticity matrix is transformed to this to integrate KE.
!  Defined the same way as the zero THETA material coordinate system in Siemens SimCenter. This definition is used because it
!  has uniform orientation on distorted flat elements and is only non-uniform as needed to accomodate out-of-plane curvature.
!  The uniformity allows stress to be interpolated and extrapolated to different locations conveniently.
!  Orthogonal
!
! Element local (Nastran definition)
!  x_l, y_l, z_l
!  Used for element stress, strain, and force outputs
!  Defined the same way as MSC (element coordinate system) and SimCenter (local coordinate system).
!  Defined by x_l being the bisection of the R, S isoparametric basis vectors rotated about the normal by -45 degrees.
!  Orthogonal
!
! XEL element (internal use)
!  Used for the grid point DOFs of the strain-displacement and the element stiffness matrices.
!  Used for extrapolating stress and strain from Gauss points to corners.
!  Grid point coordinates stored in XEL are in a coordinate system which is flat, with the normal being the cross product
!  of vectors from grid points 1-3 and 2-4. The x axis is an arbitrary direction in this plane. The flat coordinate system
!  is used for grid point coordinates for extrapolating stress because the polynomial curve fit code to extrapolate stress/strain
!  from Gauss points to corners is only 2D.
!  Orthogonal
!
! Material
!  Used for material elasticity read from the input file.
!  Currently, this is the same as the cartesian local coordinate system. To allow non-isotropic materials, it
!  should find the angle between the two systems's x axes at each integration point and rotate the material
!  elasticity matrix about that when building the stiffness matrix KE.
!  Orthogonal
!
! Isoparametric (natural)
!  R, S, T in code. r_1, r_2, r_3 in Bathe.
!  Each coordinate has range [-1,1].
!  T is parallel to the director vector.
!  Not orthogonal
!
! Covariant
!  g_r, g_s, g_t. g_1, g_2, g_3 in Bathe.
!  Parallel to the isoparametric coordinates but scaled by the element size. Eg. |g_t| = half thickness.
!  Not orthogonal
!
! Contravariant
!  g^r, g^s, g^t. g^1, g^2, g^3 in Bathe.
!  Contravariant to the covariant.
!  Not orthogonal
!

! **********************************************************************************************************************************

! Initialize
      PHI_SQ  = ONE                                        ! Not used for this element
      CALL MITC_INITIALIZE ()


      IF (PCOMP_PROPS == 'Y') THEN
        WRITE(ERR,*) ' *ERROR: Code not written for composite material with QUAD8'
        WRITE(F06,*) ' *ERROR: Code not written for composite material with QUAD8'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )
      ENDIF



! **********************************************************************************************************************************
! Generate the mass matrix for this element.

      IF (OPT(1) == 'Y') THEN
         M_1DOF = ZERO
         MASS_ELEM = ZERO
         DENSITY = MASS_PER_UNIT_AREA/EPROP(1)

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         CALL ORDER_GAUSS ( IORD_K , SS_K , HH_K  )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               DO K=1,IORD_K
                  R = SS_IJ(I)
                  S = SS_IJ(J)
                  T = SS_K(K)
                  CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )
                  DETJ = MITC_DETJ ( R, S, T )
                  INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)
                  MASS_ELEM = MASS_ELEM + DENSITY*INTFAC
                  DO L=1,ELGP
                     DO M=1,ELGP
                        M_1DOF(L,M) = M_1DOF(L,M) + PSH(L)*PSH(M)*DENSITY*INTFAC
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         ME = ZERO
         IF ((SOL_NAME(1:5) == 'MODES') .AND. (COUPMASS > 0)) THEN
            DO L=1,ELGP
               DO M=1,ELGP
                  DO K=1,3
                     ME(6*(L-1)+K,6*(M-1)+K) = M_1DOF(L,M)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            MASS_NODE = MASS_ELEM/REAL(ELGP,DOUBLE)
            DO L=1,ELGP
               DO K=1,3
                  ME(6*(L-1)+K,6*(L-1)+K) = MASS_NODE
               ENDDO
            ENDDO
         ENDIF
      ENDIF



! **********************************************************************************************************************************
! Calculate element thermal loads.

      IF (OPT(2) == 'Y') THEN
         E = MITC_ELASTICITY()
         UNIT_PTE(:) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         CALL ORDER_GAUSS ( IORD_K , SS_K , HH_K  )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               DO K=1,IORD_K
                  R = SS_IJ(I)
                  S = SS_IJ(J)
                  T = SS_K(K)
                  CLB = MITC8_CARTESIAN_LOCAL_BASIS ( R, S )
                  MATL_AXES_ROTATE = -ATAN2(CLB(2,1), CLB(1,1))
                  CALL PLANE_COORD_TRANS_21 ( MATL_AXES_ROTATE, TRANSFORM, SUBR_NAME )
                  CALL MATL_TRANSFORM_MATRIX ( TRANSFORM, T66 )
                  T66 = TRANSPOSE(T66)
                  CALL MATMULT_FFF   ( E  , T66   , 6, 6, 6, DUM66 )
                  CALL MATMULT_FFF_T ( T66 , DUM66 , 6, 6, 6, E3    )
                  CTE(:) = ALPVEC(:,1)
                  CTE(4:6) = CTE(4:6) / TWO
                  CTE = MATMUL(TRANSPOSE(T66), CTE)
                  CTE(4:6) = CTE(4:6) * TWO
                  DETJ = MITC_DETJ ( R, S, T )
                  INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)
                  CALL MITC8_B ( R, S, T, .TRUE., .TRUE., BI )
                  CALL MATMULX_FFF_T ( BI, E3, 6, 6*ELGP, 6, DUM1 )
                  THERMAL_STRAIN = MATMUL(E3, CTE)
                  UNIT_PTE = UNIT_PTE + MATMUL( TRANSPOSE(BI), THERMAL_STRAIN ) * INTFAC
               ENDDO
            ENDDO
         ENDDO

         DO JSUB=1,SIZE(PTE,2)
            TBAR = ZERO
            DO J=1,ELGP
               TBAR = TBAR + DT(J,JSUB)
            ENDDO
            TBAR = TBAR / REAL(ELGP,DOUBLE) - TREF(1)
            PTE(1:6*ELGP,JSUB) = UNIT_PTE(1:6*ELGP) * TBAR
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! BE1 matrix (3 x 48) for membrane strain/stress/force data recovery.
! BE2 matrix (3 x 48) for bending strain/stress/force data recovery.
! BE3 matrix (2 x 48) for transverse shear strain/stress/force data recovery.
! All calculated at Gauss points and not center.
! The displacements are in basic coordinates and the strains are in element coordinates.


! There's a possible bug where the numbering of the element grid points (eg. 1-2-3-4 vs 2-3-4-1) affects the stress/strain/elforce,
! when the element is distorted. This includes von Mises which should be invariant to rotation.
!
! - It doesn't occur if we skip the extrapolation from Gauss points to corners so Gauss point strains are probably OK.
! - It makes no differences if the cartesian local coordiante system is changed to be G1G2 or the bisected diagonals projected
!   onto the surface everywhere or like Siemens material coordinates with an intermediate reference plane.
! - It works OK if the cartesian local coordinate system is defined using the same vector for each element, eg. (1,0,0), instead
!   of G1G2 with either the Siemens or direct projection. However, this won't generalize to elements in any orientation and might
!   just be hiding the problem.
! - The problem is probably the extrapolation from Gauss points to grid points. Maybe it should be done in covariant coordinates
!   the way strains are interpolated by MITC. Somehow.

      IF (OPT(3) == 'Y') THEN

         STR_PT_NUM = 1

         CALL ORDER_GAUSS ( IORD_STRESS_Q8, SS_IJ, HH_IJ )

         DO I=1,IORD_STRESS_Q8
            DO J=1,IORD_STRESS_Q8

               STR_PT_NUM = STR_PT_NUM + 1

               R = SS_IJ(I)
               S = SS_IJ(J)

               CALL MITC8_B( R, S, -ONE, .TRUE., .TRUE., BI1)
               CALL MITC8_B( R, S, +ONE, .TRUE., .TRUE., BI2)

                                                  ! Membrane strain is the average of the strains at the two t points.
               BE1(1,:,STR_PT_NUM) = (BI2(1,:) + BI1(1,:)) / TWO           ! xx
               BE1(2,:,STR_PT_NUM) = (BI2(2,:) + BI1(2,:)) / TWO           ! yy
               BE1(3,:,STR_PT_NUM) = (BI2(4,:) + BI1(4,:)) / TWO           ! xy

                                                  ! Curvature is (strain_top - strain_bottom) / thickness
                                                  ! To allow grid point thicknesses, this should be the thickness
                                                  ! interpolated at the Gauss point.
               BE2(1,:,STR_PT_NUM) = (BI2(1,:) - BI1(1,:)) / EPROP(1)      ! xx
               BE2(2,:,STR_PT_NUM) = (BI2(2,:) - BI1(2,:)) / EPROP(1)      ! yy
               BE2(3,:,STR_PT_NUM) = (BI2(4,:) - BI1(4,:)) / EPROP(1)      ! xy

                                                  ! Transverse shear strain. Note reversed order of rows.
               BE3(1,:,STR_PT_NUM) = (BI2(6,:) + BI1(6,:)) / TWO           ! zx
               BE3(2,:,STR_PT_NUM) = (BI2(5,:) + BI1(5,:)) / TWO           ! yz

            ENDDO
         ENDDO

                                                           ! Find angle of the element coordinate system's x axis from
                                                           ! the cartesian local coordinate system's x axis at each
                                                           ! corner.
                                                           ! This will be used to transform stress and strain to the
                                                           ! element coordinate system after extrapolating to corners.
         DO STR_PT_NUM=2,5

            R = GP_RS(1, STR_PT_NUM - 1)
            S = GP_RS(2, STR_PT_NUM - 1)

            LOCAL_BASIS = MITC8_CARTESIAN_LOCAL_BASIS( R, S )
            XL = LOCAL_BASIS(:,1)                          ! X axis of cartesian local basis
            ZL = LOCAL_BASIS(:,3)                          ! Normal
            ELEMENT_BASIS = MITC8_ELEMENT_CS_BASIS( R, S )
            XE = ELEMENT_BASIS(:,1)                        ! X axis of element coordinate system

            CALL CROSS( XL, XE, CROSS_XLE )
            SHELL_STR_ANGLE( STR_PT_NUM ) = ATAN2(DOT_PRODUCT( ZL, CROSS_XLE ), DOT_PRODUCT( XL, XE ))

         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate element stiffness matrix KE.

      IF(OPT(4) == 'Y') THEN

! Based on
! MITC4 paper "A continuum mechanics based four-node shell element for general nonlinear analysis"
!   by Dvorkin and Bathe
! MITC8 paper "A FORMULATION OF GENERAL SHELL ELEMENTS-THE USE OF MIXED INTERPOLATION OF TENSORIAL COMPONENTS"
!   by Dvorkin and Bathe, 1986


         ! K = int( [B]^T [EE] [B] dV )
         !    dV = |det(J)|dr ds dt
         ! K = int( [B]^T [EE] [B] |det(J)| dr ds dt )
         !
         ! [EE] is material elasticity matrix in cartesian local coordinates
         ! stress = [EE] * strain
         ! strain = [B] * displacement
         ! K is in the basic coordinate system

         E = MITC_ELASTICITY()
         GDRILL = E(4,4)

         KE(1:6*ELGP,1:6*ELGP) = ZERO

         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         CALL ORDER_GAUSS ( IORD_K, SS_K, HH_K )

         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               DO K=1,IORD_K
                  R = SS_IJ(I)
                  S = SS_IJ(J)
                  T = SS_K(K)
                  CALL MITC8_B( R, S, T, .TRUE., .TRUE., BI)

                  ! For non-isotropic materials, this should be rotated from the material coordinate system to the cartesian local
                  ! coordinate system here. The rotation angle may be different at each Gauss point.
                  EE(:,:) = E(:,:)

                  CALL MATMULX_FFF ( EE, BI, 6, 6, 6*ELGP, DUM1 )
                  CALL MATMULX_FFF_T ( BI, DUM1, 6, 6*ELGP, 6*ELGP, DUM2 )
                  DETJ = MITC_DETJ ( R, S, T )
                  INTFAC = DETJ*HH_IJ(I)*HH_IJ(J)*HH_K(K)
                  KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + DUM2(:,:)*INTFAC
                  CALL MITC8_DRILL_B ( R, S, BDRILL )
                  CALL MATMULT_FFF_T ( BDRILL, BDRILL, 1, 6*ELGP, 6*ELGP, DUM2 )
                  KE(1:6*ELGP,1:6*ELGP) = KE(1:6*ELGP,1:6*ELGP) + BETA_DRILL*GDRILL*DUM2(:,:)*INTFAC
               ENDDO
            ENDDO
         ENDDO



      ENDIF


! **********************************************************************************************************************************
! Determine element pressure loads

      IF (OPT(5) == 'Y') THEN

         UNIT_PPE = ZERO
         CALL ORDER_GAUSS ( IORD_IJ, SS_IJ, HH_IJ )
         DO I=1,IORD_IJ
            DO J=1,IORD_IJ
               R = SS_IJ(I)
               S = SS_IJ(J)
               CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH, DPSHG )
               DXDR = ZERO
               DXDS = ZERO
               DO L=1,ELGP
                  DXDR(:) = DXDR(:) + DPSHG(1,L)*XEB(L,:)
                  DXDS(:) = DXDS(:) + DPSHG(2,L)*XEB(L,:)
               ENDDO
               CALL CROSS ( DXDR, DXDS, SURF_VEC )
               INTFAC = HH_IJ(I)*HH_IJ(J)
               DO L=1,ELGP
                  UNIT_PPE(6*(L-1)+1) = UNIT_PPE(6*(L-1)+1) + PSH(L)*SURF_VEC(1)*INTFAC
                  UNIT_PPE(6*(L-1)+2) = UNIT_PPE(6*(L-1)+2) + PSH(L)*SURF_VEC(2)*INTFAC
                  UNIT_PPE(6*(L-1)+3) = UNIT_PPE(6*(L-1)+3) + PSH(L)*SURF_VEC(3)*INTFAC
               ENDDO
            ENDDO
         ENDDO
         DO J=1,SIZE(PPE,2)
            PPE(1:6*ELGP,J) = PPE(1:6*ELGP,J) + UNIT_PPE(1:6*ELGP)*PRESS(3,J)
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Calculate linear differential stiffness matrix

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN

        WRITE(ERR,*) ' *ERROR: Code not written for QUAD8 differential stiffness matrix'
        WRITE(F06,*) ' *ERROR: Code not written for QUAD8 differential stiffness matrix'
        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
        FATAL_ERR = FATAL_ERR + 1
        CALL OUTA_HERE ( 'Y' )

      ENDIF






      RETURN

! **********************************************************************************************************************************

      CONTAINS

      SUBROUTINE MITC8_DRILL_B ( R, S, BDOUT )

      REAL(DOUBLE), INTENT(IN)        :: R, S
      REAL(DOUBLE), INTENT(OUT)       :: BDOUT(1,6*ELGP)

      INTEGER(LONG)                   :: II
      REAL(DOUBLE)                    :: PSH_D(ELGP)
      REAL(DOUBLE)                    :: DPSHG_D(2,ELGP)
      REAL(DOUBLE)                    :: CLB_D(3,3)
      REAL(DOUBLE)                    :: XI_LOC(ELGP), ETA_LOC(ELGP)
      REAL(DOUBLE)                    :: J11, J12, J21, J22, DET2
      REAL(DOUBLE)                    :: DNDX, DNDY

      BDOUT = ZERO
      CALL MITC_SHAPE_FUNCTIONS ( R, S, PSH_D, DPSHG_D )
      CLB_D = MITC8_CARTESIAN_LOCAL_BASIS ( R, S )

      DO II=1,ELGP
         XI_LOC(II)  = DOT_PRODUCT( XEL(II,:), CLB_D(:,1) )
         ETA_LOC(II) = DOT_PRODUCT( XEL(II,:), CLB_D(:,2) )
      ENDDO

      J11 = DOT_PRODUCT( DPSHG_D(1,:), XI_LOC  )
      J12 = DOT_PRODUCT( DPSHG_D(1,:), ETA_LOC )
      J21 = DOT_PRODUCT( DPSHG_D(2,:), XI_LOC  )
      J22 = DOT_PRODUCT( DPSHG_D(2,:), ETA_LOC )
      DET2 = J11*J22 - J12*J21
      IF (DABS(DET2) < 1.0D-14) RETURN

      DO II=1,ELGP
         DNDX = ( J22*DPSHG_D(1,II) - J12*DPSHG_D(2,II) ) / DET2
         DNDY = (-J21*DPSHG_D(1,II) + J11*DPSHG_D(2,II) ) / DET2
         BDOUT(1,6*(II-1)+1) = -0.5D0 * DNDY
         BDOUT(1,6*(II-1)+2) =  0.5D0 * DNDX
         BDOUT(1,6*(II-1)+6) =  PSH_D(II)
      ENDDO

      END SUBROUTINE MITC8_DRILL_B

! **********************************************************************************************************************************

      END SUBROUTINE MITC8_LEGACY

! #################################################################################################################################
! CQUAD8 MITC8 v3 shell for PARAM,QUAD8TYP,MITC8.

      SUBROUTINE MITC8 ( OPT, INT_ELEM_ID )

! Ported from:
!   MITC8_ShellElement_v3.py / MITC8_ShellElement_v3 (static stiffness and recovery).
! Adaptation: modified field metric, standard geometry area; not a paper reproduction.
! Geometric stiffness uses MITC8_LEGACY;
! these auxiliary paths are not validated against the v3 static reference.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MAX_STRESS_POINTS, SOL_NAME, NSNORM
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  ALPVEC, BGRID, GRID_ID, GRID_SNORM, SNORM, Q8_POINT_BASIS, DT, EID, ELGP, KE, KED, ME, BE1, BE2, BE3, EPROP, MASS_PER_UNIT_AREA, PPE,   &
                                         PRESS, PTE, SHELL_A, SHELL_D, SHELL_T, TREF, UEL, XEB, NUM_EMG_FATAL_ERRS,           &
                                         PCOMP_PROPS
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO
      USE PARAMS, ONLY                :  COUPMASS
      USE ELMDIS_Interface
      USE OUTA_HERE_Interface

      USE QUADRATIC_SURFACE_MASS_Interface
      USE QUADRATIC_SURFACE_PRESSURE_Interface
      USE MODEL_STUF, ONLY : UEB, KED
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MITC8'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)
      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID

      INTEGER(LONG), PARAMETER        :: NNODE = 8
      INTEGER(LONG), PARAMETER        :: NDOF  = 48
      REAL(DOUBLE), PARAMETER         :: KT_DRILL = 1.0D-4
      INTEGER(LONG)                   :: I, J, K, L, IA, JSUB, GP
      REAL(DOUBLE)                    :: XYZ(8,3), XY8(8,2), E1F(3), E2F(3), E3F(3)
      REAL(DOUBLE)                    :: GP3(3), W3(3), R, S, WT, DETJ, THICK, TBAR, GVAL
      REAL(DOUBLE)                    :: BM(3,NDOF), BB(3,NDOF), BS(2,NDOF), BD(1,NDOF)
      REAL(DOUBLE)                    :: M1(8,8), N8(8), DNDX(8), DNDY(8), MASS_ELEM, MASS_NODE
      REAL(DOUBLE)                    :: UNIT_PPE(NDOF), UNIT_PTE(NDOF), DXDR(3), DXDS(3), SURF_VEC(3)
      REAL(DOUBLE)                    :: CTE(3), THERMAL_RESULTANT(3), CDRILL
      REAL(DOUBLE)                    :: SIG0(2,2), KG8(8,8), STRAIN0(3), N0V(3)
      REAL(DOUBLE) :: NORMALS(8,3),NORMAL_SIGN,FIELD_COEFF(8),RR(9),SS(9),E1OUT(3),E2OUT(3),E3OUT(3),DN8(2,8)
      LOGICAL :: USE_ANS
      CHARACTER(1*BYTE) :: LEGACY_OPT(6)
      CHARACTER(8),PARAMETER :: ANSSHEAR="TENSOR6"
      INTEGER(LONG)                   :: KI, KJ

      REAL(DOUBLE) :: UNIT_PTG(48), TEMP_GRAD_SIGN, TEMP_NORMAL(3), TEMP_G1(3), TEMP_G2(3)

      IF (ELGP /= NNODE) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9001) SUBR_NAME, EID, ELGP
         WRITE(F06,9001) SUBR_NAME, EID, ELGP
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (PCOMP_PROPS == 'Y') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,MITC8'
         WRITE(F06,*) ' *ERROR: Code not written for composite material with PARAM,QUAD8TYP,MITC8'
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      CALL LOAD_BASIC_COORDS_Q8(XYZ)
      CALL FIXED_FRAME_Q8(XYZ, E1F, E2F, E3F)
      CALL BUILD_LOCAL_XY_Q8(XYZ, E1F, E2F, XY8)

      CALL MI_FIELD_COEFF_Q8(XY8)
      CALL MI_SET_NORMAL_SIGN_Q8(XYZ)
      CALL MI_CALC_NODAL_NORMALS_Q8(XYZ,NORMALS)
      USE_ANS=.FALSE.

      GP3 = (/-DSQRT(3.0D0/5.0D0), ZERO, DSQRT(3.0D0/5.0D0)/)
      W3  = (/5.0D0/9.0D0, 8.0D0/9.0D0, 5.0D0/9.0D0/)
      THICK = EPROP(1)
      GVAL = SHELL_T(1,1)
      CDRILL = KT_DRILL*GVAL

! All active options, including geometric stiffness, bypass the legacy branch.
      LEGACY_OPT=OPT
      LEGACY_OPT(1)='N'
      LEGACY_OPT(2)='N'
      LEGACY_OPT(3)='N'
      LEGACY_OPT(4)='N'
      LEGACY_OPT(5)='N'
      LEGACY_OPT(6)='N'
      IF (ANY(LEGACY_OPT == 'Y')) CALL MITC8_LEGACY(LEGACY_OPT,INT_ELEM_ID)

      IF (OPT(1) == 'Y') THEN
         CALL MI_SHAPE_Q8(ZERO,ZERO,N8,DN8)
         N8(1:4)=N8(1:4)+0.25D0
         N8(5:8)=N8(5:8)-0.5D0
         CALL QUADRATIC_SURFACE_MASS(8,XYZ,N8)
      ENDIF

      IF (OPT(2) == 'Y') THEN
! TEMPP1 gradient is along connectivity normal; directors use the normalized frame.
         TEMP_G1=XYZ(2,:)-XYZ(1,:)
         TEMP_G2=XYZ(3,:)-XYZ(1,:)
         TEMP_NORMAL(1)=TEMP_G1(2)*TEMP_G2(3)-TEMP_G1(3)*TEMP_G2(2)
         TEMP_NORMAL(2)=TEMP_G1(3)*TEMP_G2(1)-TEMP_G1(1)*TEMP_G2(3)
         TEMP_NORMAL(3)=TEMP_G1(1)*TEMP_G2(2)-TEMP_G1(2)*TEMP_G2(1)
         TEMP_GRAD_SIGN=SIGN(ONE,DOT_PRODUCT(TEMP_NORMAL,NORMALS(1,:)))
         CTE=(/ALPVEC(1,1),ALPVEC(2,1),ALPVEC(4,1)/)
         UNIT_PTE=ZERO
         UNIT_PTG=ZERO
         DO I=1,3
            DO J=1,3
               R=GP3(I)
               S=GP3(J)
               WT=W3(I)*W3(J)
               CALL MI_BM_Q8_AT(XYZ,R,S,BM,DETJ)
               CALL MI_BB_Q8_AT(XYZ,NORMALS,R,S,BB,DETJ)
               CALL MI_GEOM_SHAPE_Q8(R,S,N8,DN8)
               DXDR=MATMUL(DN8(1,:),XYZ)
               DXDS=MATMUL(DN8(2,:),XYZ)
               CALL MI_CROSS3(DXDR,DXDS,SURF_VEC)
               WT=WT*MI_VNORM(SURF_VEC)
               UNIT_PTE=UNIT_PTE+MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,CTE))*WT
! eps(z)=Bm*u-z*Bb*u: bending thermal force has the physical minus sign.
               UNIT_PTG=UNIT_PTG-MATMUL(TRANSPOSE(BB),MATMUL(SHELL_D,CTE))*WT
            ENDDO
         ENDDO
         DO JSUB=1,SIZE(PTE,2)
            TBAR=SUM(DT(1:ELGP,JSUB))/REAL(ELGP,DOUBLE)-TREF(1)
            PTE(1:48,JSUB)=UNIT_PTE*TBAR+UNIT_PTG*TEMP_GRAD_SIGN*DT(ELGP+1,JSUB)
         ENDDO
      ENDIF

      IF (OPT(3) == 'Y') THEN
         RR=(/ZERO,-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO,-ONE/)
         SS=(/ZERO,-ONE,-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO/)
         DO GP=1,9
            R=RR(GP)
            S=SS(GP)
            CALL MI_BM_Q8_AT(XYZ,R,S,BM,DETJ)
            CALL MI_BB_Q8_AT(XYZ,NORMALS,R,S,BB,DETJ)
            CALL MI_BS_Q8_AT(XYZ,NORMALS,R,S,BS,DETJ)
            BE1(1:3,1:48,GP)=BM
! Physical fiber strain is membrane-z*Bb*u; Bb is +Hessian(w).
            BE2(1:3,1:48,GP)=BB
            BE3(1:2,1:48,GP)=BS
            CALL MI_LOCAL_BASIS_AT_Q8(XYZ,R,S,E1OUT,E2OUT,E3OUT,DETJ)
            Q8_POINT_BASIS(1,:,GP)=E1OUT
            Q8_POINT_BASIS(2,:,GP)=E2OUT
            Q8_POINT_BASIS(3,:,GP)=E3OUT
         ENDDO
      ENDIF

      IF (OPT(4) == 'Y') THEN
         KE=ZERO
! Python beta_drill=1e-4; remove PSHELL shear correction from G*h.
         CDRILL=KT_DRILL*SHELL_T(1,1)/(5.0D0/6.0D0)
         DO I=1,3
            DO J=1,3
               R=GP3(I)
               S=GP3(J)
               CALL MI_BM_Q8_AT(XYZ,R,S,BM,DETJ)
               CALL MI_BB_Q8_AT(XYZ,NORMALS,R,S,BB,DETJ)
               CALL MI_BS_Q8_AT(XYZ,NORMALS,R,S,BS,DETJ)
               CALL MI_BDRILL_Q8_AT(XYZ,NORMALS,R,S,BD,DETJ)
! Area is always standard Q8 geometry, independently of modified field metric.
               CALL MI_GEOM_SHAPE_Q8(R,S,N8,DN8)
               DXDR=MATMUL(DN8(1,:),XYZ)
               DXDS=MATMUL(DN8(2,:),XYZ)
               CALL MI_CROSS3(DXDR,DXDS,SURF_VEC)
               WT=W3(I)*W3(J)*MI_VNORM(SURF_VEC)
               KE=KE+WT*(MATMUL(TRANSPOSE(BM),MATMUL(SHELL_A,BM)) &
                       +MATMUL(TRANSPOSE(BB),MATMUL(SHELL_D,BB)) &
                       +MATMUL(TRANSPOSE(BS),MATMUL(SHELL_T,BS)) &
                       +CDRILL*MATMUL(TRANSPOSE(BD),BD))
            ENDDO
         ENDDO
      ENDIF

      IF (OPT(5) == 'Y') THEN
         CALL MI_SHAPE_Q8(ZERO,ZERO,N8,DN8)
         N8(1:4)=N8(1:4)+0.25D0
         N8(5:8)=N8(5:8)-0.5D0
         CALL QUADRATIC_SURFACE_PRESSURE(INT_ELEM_ID,8,XYZ,N8)
      ENDIF

      IF ((OPT(6) == 'Y') .AND. (LOAD_ISTEP > 1)) THEN
         CALL NATIVE_MEMBRANE_KG
      ENDIF

      RETURN

 9001 FORMAT(' *ERROR: ',A,' expects ELGP=8 for element ',I8,' but got ',I8)

      CONTAINS

! Flat-shell membrane initial-stress stiffness, tension-positive resultants.
! Four-point Gauss/Duffy; standard geometry area, active translation field.
! Mechanical linear reference state only; no director or follower tangent.
      SUBROUTINE NATIVE_MEMBRANE_KG
      REAL(DOUBLE) :: GX(4),GW(4),RG,SG,WG,NVAL(8),DG(2,8),DF(2,8)
      REAL(DOUBLE) :: TG(2,3),TF(2,3),CV(3),AREA,AJ,MT(2,2),INV(2,2),DETMT
      REAL(DOUBLE) :: E1(3),E2(3),E3(3),GRAD(2,8),BMG(3,48),NV(3),SIG(2,2),BLOCK(8,8)
      REAL(DOUBLE) :: SCALE_GEOM,NORMAL(3)
      INTEGER(LONG) :: IG,JG,IN,JN,ID
      GX=(/-0.8611363115940526D0,-0.3399810435848563D0, &
            0.3399810435848563D0,0.8611363115940526D0/)
      GW=(/0.3478548451374538D0,0.6521451548625461D0, &
           0.6521451548625461D0,0.3478548451374538D0/)
      TG(1,:)=XYZ(2,:)-XYZ(1,:)
      TG(2,:)=XYZ(3,:)-XYZ(1,:)
      NORMAL=(/TG(1,2)*TG(2,3)-TG(1,3)*TG(2,2), &
                TG(1,3)*TG(2,1)-TG(1,1)*TG(2,3), &
                TG(1,1)*TG(2,2)-TG(1,2)*TG(2,1)/)
      AREA=SQRT(SUM(NORMAL*NORMAL))
      SCALE_GEOM=MAXVAL(ABS(XYZ-SPREAD(XYZ(1,:),1,8)))
      IF (AREA <= 1.0D-14) THEN
         CALL KG_GEOMETRY_ERROR
      ENDIF
      NORMAL=NORMAL/AREA
      DO IN=1,8
         IF (ABS(DOT_PRODUCT(XYZ(IN,:)-XYZ(1,:),NORMAL)) > 1.0D-9*MAX(SCALE_GEOM,ONE)) THEN
            CALL KG_GEOMETRY_ERROR
         ENDIF
      ENDDO
      CALL ELMDIS
      BLOCK=ZERO
      DO IG=1,4
         DO JG=1,4
               RG=GX(IG)
               SG=GX(JG)
               WG=GW(IG)*GW(JG)
               CALL MI_GEOM_SHAPE_Q8(RG,SG,NVAL,DG)
               TG=MATMUL(DG,XYZ)
               CV=(/TG(1,2)*TG(2,3)-TG(1,3)*TG(2,2), &
                     TG(1,3)*TG(2,1)-TG(1,1)*TG(2,3), &
                     TG(1,1)*TG(2,2)-TG(1,2)*TG(2,1)/)
               AREA=SQRT(SUM(CV*CV))
               CALL MI_SHAPE_Q8(RG,SG,NVAL,DF)
               TF=MATMUL(DF,XYZ)
               MT=MATMUL(TF,TRANSPOSE(TF))
               DETMT=MT(1,1)*MT(2,2)-MT(1,2)*MT(2,1)
               IF (AREA <= 1.0D-14 .OR. DETMT <= 1.0D-30) CALL KG_GEOMETRY_ERROR
               INV(1,1)=MT(2,2)/DETMT
               INV(2,2)=MT(1,1)/DETMT
               INV(1,2)=-MT(1,2)/DETMT
               INV(2,1)=INV(1,2)
               CALL MI_LOCAL_BASIS_AT_Q8(XYZ,RG,SG,E1,E2,E3,AJ)
               GRAD(1,:)=MATMUL(MATMUL(INV,MATMUL(TF,E1)),DF)
               GRAD(2,:)=MATMUL(MATMUL(INV,MATMUL(TF,E2)),DF)
               CALL MI_BM_Q8_AT(XYZ,RG,SG,BMG,AJ)
               NV=MATMUL(SHELL_A,MATMUL(BMG,UEB(1:48)))
               SIG(1,:)=(/NV(1),NV(3)/)
               SIG(2,:)=(/NV(3),NV(2)/)
               BLOCK=BLOCK+MATMUL(TRANSPOSE(GRAD),MATMUL(SIG,GRAD))*AREA*WG
         ENDDO
      ENDDO
      BLOCK=(BLOCK+TRANSPOSE(BLOCK))/TWO
      KED(1:48,1:48)=ZERO
      DO IN=1,8
         DO JN=1,8
            DO ID=1,3
               KED(6*(IN-1)+ID,6*(JN-1)+ID)=BLOCK(IN,JN)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE NATIVE_MEMBRANE_KG

      SUBROUTINE KG_GEOMETRY_ERROR
      WRITE(ERR,*) ' *ERROR: native Q8/T6 membrane buckling requires nondegenerate flat geometry. EID=',EID
      WRITE(F06,*) ' *ERROR: native Q8/T6 membrane buckling requires nondegenerate flat geometry. EID=',EID
      FATAL_ERR=FATAL_ERR+1
      NUM_EMG_FATAL_ERRS=NUM_EMG_FATAL_ERRS+1
      CALL OUTA_HERE('Y')
      END SUBROUTINE KG_GEOMETRY_ERROR



      SUBROUTINE LOAD_BASIC_COORDS_Q8 ( XYZOUT )
      REAL(DOUBLE), INTENT(OUT) :: XYZOUT(8,3)
      INTEGER(LONG) :: II, JJ
      DO II=1,8
         DO JJ=1,3
            XYZOUT(II,JJ) = XEB(II,JJ)
         ENDDO
      ENDDO
      END SUBROUTINE LOAD_BASIC_COORDS_Q8

      SUBROUTINE SHAPE_Q8_STD_DERIVS ( XI, ETA, NVAL, DN )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), DN(2,8)
      NVAL(1) = 0.25D0*(ONE-XI)*(ONE-ETA)*(-XI-ETA-ONE)
      NVAL(2) = 0.25D0*(ONE+XI)*(ONE-ETA)*( XI-ETA-ONE)
      NVAL(3) = 0.25D0*(ONE+XI)*(ONE+ETA)*( XI+ETA-ONE)
      NVAL(4) = 0.25D0*(ONE-XI)*(ONE+ETA)*(-XI+ETA-ONE)
      NVAL(5) = 0.5D0*(ONE-XI*XI)*(ONE-ETA)
      NVAL(6) = 0.5D0*(ONE+XI)*(ONE-ETA*ETA)
      NVAL(7) = 0.5D0*(ONE-XI*XI)*(ONE+ETA)
      NVAL(8) = 0.5D0*(ONE-XI)*(ONE-ETA*ETA)

      DN(1,1) = 0.25D0*(ONE-ETA)*(2.0D0*XI+ETA)
      DN(1,2) = 0.25D0*(ONE-ETA)*(2.0D0*XI-ETA)
      DN(1,3) = 0.25D0*(ONE+ETA)*(2.0D0*XI+ETA)
      DN(1,4) = 0.25D0*(ONE+ETA)*(2.0D0*XI-ETA)
      DN(1,5) = -XI*(ONE-ETA)
      DN(1,6) = 0.5D0*(ONE-ETA*ETA)
      DN(1,7) = -XI*(ONE+ETA)
      DN(1,8) = -0.5D0*(ONE-ETA*ETA)

      DN(2,1) = 0.25D0*(ONE-XI)*(2.0D0*ETA+XI)
      DN(2,2) = 0.25D0*(ONE+XI)*(2.0D0*ETA-XI)
      DN(2,3) = 0.25D0*(ONE+XI)*(2.0D0*ETA+XI)
      DN(2,4) = 0.25D0*(ONE-XI)*(2.0D0*ETA-XI)
      DN(2,5) = -0.5D0*(ONE-XI*XI)
      DN(2,6) = -(ONE+XI)*ETA
      DN(2,7) = 0.5D0*(ONE-XI*XI)
      DN(2,8) = -(ONE-XI)*ETA
      END SUBROUTINE SHAPE_Q8_STD_DERIVS

      SUBROUTINE SHAPE_Q8_STD ( XI, ETA, NVAL, G1, G2, XYZN, DETJ )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA, XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), G1(3), G2(3), DETJ
      REAL(DOUBLE) :: DN(2,8), GV(3)
      CALL SHAPE_Q8_STD_DERIVS(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL CROSS3(G1, G2, GV)
      DETJ = VNORM(GV)
      END SUBROUTINE SHAPE_Q8_STD

      SUBROUTINE FIXED_FRAME_Q8 ( XYZN, E1OUT, E2OUT, E3OUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: E1OUT(3), E2OUT(3), E3OUT(3)
      REAL(DOUBLE) :: C(4,3), TMP(3), NM
      C = XYZN(1:4,:)
      E1OUT = 0.5D0*((C(2,:)-C(1,:)) + (C(3,:)-C(4,:)))
      NM = VNORM(E1OUT)
      IF (NM > 1.0D-15) THEN
         E1OUT = E1OUT/NM
      ELSE
         E1OUT = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL CROSS3(C(3,:)-C(1,:), C(4,:)-C(2,:), E3OUT)
      NM = VNORM(E3OUT)
      IF (NM > 1.0D-15) THEN
         E3OUT = E3OUT/NM
      ELSE
         E3OUT = (/ZERO, ZERO, ONE/)
      ENDIF
      CALL CROSS3(E3OUT, E1OUT, E2OUT)
      NM = VNORM(E2OUT)
      IF (NM > 1.0D-15) E2OUT = E2OUT/NM
      CALL CROSS3(E2OUT, E3OUT, TMP)
      NM = VNORM(TMP)
      IF (NM > 1.0D-15) E1OUT = TMP/NM
      END SUBROUTINE FIXED_FRAME_Q8

      SUBROUTINE BUILD_LOCAL_XY_Q8 ( XYZN, E1IN, E2IN, XYOUT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), E1IN(3), E2IN(3)
      REAL(DOUBLE), INTENT(OUT) :: XYOUT(8,2)
      REAL(DOUBLE) :: ORG(3)
      INTEGER(LONG) :: II
      ORG = 0.25D0*(XYZN(1,:) + XYZN(2,:) + XYZN(3,:) + XYZN(4,:))
      DO II=1,8
         XYOUT(II,1) = DOT_PRODUCT(XYZN(II,:) - ORG, E1IN)
         XYOUT(II,2) = DOT_PRODUCT(XYZN(II,:) - ORG, E2IN)
      ENDDO
      END SUBROUTINE BUILD_LOCAL_XY_Q8

      FUNCTION L3 ( X, XN, K1 ) RESULT(VAL)
      REAL(DOUBLE), INTENT(IN) :: X, XN(3)
      INTEGER(LONG), INTENT(IN):: K1
      REAL(DOUBLE) :: VAL
      INTEGER(LONG) :: J1
      REAL(DOUBLE) :: NUM, DEN
      NUM = ONE
      DEN = ONE
      DO J1=1,3
         IF (J1 /= K1) THEN
            NUM = NUM*(X - XN(J1))
            DEN = DEN*(XN(K1) - XN(J1))
         ENDIF
      ENDDO
      IF (DABS(DEN) > 1.0D-14) THEN
         VAL = NUM/DEN
      ELSE
         VAL = ZERO
      ENDIF
      END FUNCTION L3

      SUBROUTINE CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN)  :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      END SUBROUTINE CROSS3

      FUNCTION VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(MAX(ZERO, DOT_PRODUCT(V,V)))
      END FUNCTION VNORM

      SUBROUTINE MI_GEOM_SHAPE_Q8 ( XI, ETA, NVAL, DN )
      REAL(DOUBLE), INTENT(IN)  :: XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: NVAL(8), DN(2,8)
      NVAL(1) = 0.25D0*(ONE-XI)*(ONE-ETA)*(-XI-ETA-ONE)
      NVAL(2) = 0.25D0*(ONE+XI)*(ONE-ETA)*( XI-ETA-ONE)
      NVAL(3) = 0.25D0*(ONE+XI)*(ONE+ETA)*( XI+ETA-ONE)
      NVAL(4) = 0.25D0*(ONE-XI)*(ONE+ETA)*(-XI+ETA-ONE)
      NVAL(5) = 0.5D0*(ONE-XI*XI)*(ONE-ETA)
      NVAL(6) = 0.5D0*(ONE+XI)*(ONE-ETA*ETA)
      NVAL(7) = 0.5D0*(ONE-XI*XI)*(ONE+ETA)
      NVAL(8) = 0.5D0*(ONE-XI)*(ONE-ETA*ETA)

      DN(1,1) = 0.25D0*(ONE-ETA)*(TWO*XI+ETA)
      DN(1,2) = 0.25D0*(ONE-ETA)*(TWO*XI-ETA)
      DN(1,3) = 0.25D0*(ONE+ETA)*(TWO*XI+ETA)
      DN(1,4) = 0.25D0*(ONE+ETA)*(TWO*XI-ETA)
      DN(1,5) = -XI*(ONE-ETA)
      DN(1,6) = 0.5D0*(ONE-ETA*ETA)
      DN(1,7) = -XI*(ONE+ETA)
      DN(1,8) = -0.5D0*(ONE-ETA*ETA)

      DN(2,1) = 0.25D0*(ONE-XI)*(TWO*ETA+XI)
      DN(2,2) = 0.25D0*(ONE+XI)*(TWO*ETA-XI)
      DN(2,3) = 0.25D0*(ONE+XI)*(TWO*ETA+XI)
      DN(2,4) = 0.25D0*(ONE-XI)*(TWO*ETA-XI)
      DN(2,5) = -0.5D0*(ONE-XI*XI)
      DN(2,6) = -(ONE+XI)*ETA
      DN(2,7) = 0.5D0*(ONE-XI*XI)
      DN(2,8) = -(ONE-XI)*ETA
      END SUBROUTINE MI_GEOM_SHAPE_Q8

      SUBROUTINE MI_CALC_NODAL_NORMALS_Q8 ( XYZN, NORMS )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3)
      REAL(DOUBLE), INTENT(OUT) :: NORMS(8,3)
      REAL(DOUBLE) :: RS(8,2), NVAL(8), DN(2,8), G1(3), G2(3), N(3), NM, SN(3)
      INTEGER(LONG) :: II, BIDX, ISN
      RS(1,:) = (/-ONE, -ONE/)
      RS(2,:) = (/ ONE, -ONE/)
      RS(3,:) = (/ ONE,  ONE/)
      RS(4,:) = (/-ONE,  ONE/)
      RS(5,:) = (/ZERO, -ONE/)
      RS(6,:) = (/ ONE, ZERO/)
      RS(7,:) = (/ZERO,  ONE/)
      RS(8,:) = (/-ONE, ZERO/)
      DO II=1,8
         CALL MI_GEOM_SHAPE_Q8(RS(II,1), RS(II,2), NVAL, DN)
         G1 = MATMUL(DN(1,:), XYZN)
         G2 = MATMUL(DN(2,:), XYZN)
         CALL MI_CROSS3(G1, G2, N)
         NM = MI_VNORM(N)
         IF (NM <= 1.0D-12) THEN
            CALL MI_GEOM_SHAPE_Q8(ZERO, ZERO, NVAL, DN)
            G1 = MATMUL(DN(1,:), XYZN)
            G2 = MATMUL(DN(2,:), XYZN)
            CALL MI_CROSS3(G1, G2, N)
            NM = MI_VNORM(N)
         ENDIF
         IF (NM > 1.0D-15) THEN
            NORMS(II,:) = NORMAL_SIGN*N/NM
         ELSE
            NORMS(II,:) = (/ZERO, ZERO, ONE/)
         ENDIF
         IF (NSNORM > 0 .AND. ALLOCATED(GRID_SNORM) .AND. ALLOCATED(SNORM)) THEN
            DO ISN=1,NSNORM
            IF (SNORM(ISN,1) /= GRID_ID(BGRID(II))) CYCLE
            BIDX = 0
            IF (II <= SIZE(BGRID)) BIDX = BGRID(II)
            IF ((BIDX > 0) .AND. (BIDX <= SIZE(GRID_SNORM,1))) THEN
               SN = GRID_SNORM(BIDX,:)
               NM = MI_VNORM(SN)
               IF (NM > 1.0D-15) THEN
                  SN = SN/NM
                  IF (DOT_PRODUCT(SN, NORMS(II,:)) < ZERO) SN = -SN
                  NORMS(II,:) = SN
               ENDIF
            ENDIF
            ENDDO
         ENDIF
      ENDDO
      END SUBROUTINE MI_CALC_NODAL_NORMALS_Q8

      SUBROUTINE MI_LOCAL_BASIS_AT_Q8 ( XYZN, XI, ETA, E1, E2, E3, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: E1(3), E2(3), E3(3), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), G3(3), NM
      CALL MI_SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL MI_CROSS3(G1, G2, G3)
      JAC = MI_VNORM(G3)
      IF (JAC > 1.0D-15) THEN
         E3 = NORMAL_SIGN*G3/JAC
      ELSE
         E3 = (/ZERO, ZERO, ONE/)
      ENDIF
      NM = MI_VNORM(G1)
      IF (NM > 1.0D-15) THEN
         E1 = G1/NM
      ELSE
         E1 = (/ONE, ZERO, ZERO/)
      ENDIF
      CALL MI_CROSS3(E3, E1, E2)
      NM = MI_VNORM(E2)
      IF (NM > 1.0D-15) THEN
         E2 = E2/NM
      ELSE
         E2 = (/ZERO, ONE, ZERO/)
      ENDIF
      END SUBROUTINE MI_LOCAL_BASIS_AT_Q8

      SUBROUTINE MI_COV_MAP_Q8 ( XYZN, XI, ETA, G1, G2, C, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: G1(3), G2(3), C(4), JAC
      REAL(DOUBLE) :: E1(3), E2(3), E3(3), NVAL(8), DN(2,8)
      REAL(DOUBLE) :: A11, A22, A12, DET, AI11, AI22, AI12, GC1(3), GC2(3)
      CALL MI_SHAPE_Q8(XI, ETA, NVAL, DN)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      CALL MI_LOCAL_BASIS_AT_Q8(XYZN, XI, ETA, E1, E2, E3, JAC)
      A11 = DOT_PRODUCT(G1,G1)
      A22 = DOT_PRODUCT(G2,G2)
      A12 = DOT_PRODUCT(G1,G2)
      DET = A11*A22 - A12*A12
      IF (DABS(DET) <= 1.0D-30) THEN
         C = ZERO
         RETURN
      ENDIF
      AI11 =  A22/DET
      AI22 =  A11/DET
      AI12 = -A12/DET
      GC1 = AI11*G1 + AI12*G2
      GC2 = AI12*G1 + AI22*G2
      ! Python _covariant_maps order:
      !   C = (c1_1, c1_2, c2_1, c2_2)
      ! where c1_* maps the first physical in-plane basis direction and
      ! c2_* maps the second physical in-plane basis direction.
      C(1) = DOT_PRODUCT(E1, GC1)
      C(2) = DOT_PRODUCT(E2, GC1)
      C(3) = DOT_PRODUCT(E1, GC2)
      C(4) = DOT_PRODUCT(E2, GC2)
      END SUBROUTINE MI_COV_MAP_Q8

      SUBROUTINE MI_TENSOR_PHYS_Q8 ( V11, V22, V12, C, VOUT )
      REAL(DOUBLE), INTENT(IN)  :: V11(3), V22(3), V12(3), C(4)
      REAL(DOUBLE), INTENT(OUT) :: VOUT(3,3)
      VOUT(1,:) = C(1)*C(1)*V11 + C(3)*C(3)*V22 + TWO*C(1)*C(3)*V12
      VOUT(2,:) = C(2)*C(2)*V11 + C(4)*C(4)*V22 + TWO*C(2)*C(4)*V12
      VOUT(3,:) = TWO*(C(1)*C(2)*V11 + C(3)*C(4)*V22 + (C(1)*C(4)+C(2)*C(3))*V12)
      END SUBROUTINE MI_TENSOR_PHYS_Q8


! Matched-Gauss ILS: transport sampled physical strain tensors through basic
! coordinates before projecting into the evaluation-point physical basis.
      SUBROUTINE MI_BM_Q8_AT(XYZN,R,S,BMOUT,JAC)
      REAL(DOUBLE),INTENT(IN) :: XYZN(8,3),R,S
      REAL(DOUBLE),INTENT(OUT) :: BMOUT(3,48),JAC
      REAL(DOUBLE) :: A,N(8),DN(2,8),RR(8),SS(8),E1(3),E2(3),E3(3)
      REAL(DOUBLE) :: P1(3),P2(3),P3(3),BJ(3,48),T(3,3),GT(3,3,48),DUMMY
      INTEGER(LONG) :: II,JJ,KK,DD
      A=ONE/SQRT(3.0D0)
      RR=(/-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO,-ONE/)
      SS=(/-ONE,-ONE,ONE,ONE,-ONE,ZERO,ONE,ZERO/)
      CALL MI_GEOM_SHAPE_Q8(R/A,S/A,N,DN)
      GT=ZERO
      DO II=1,8
         CALL MI_BM_DIRECT_Q8_AT(XYZN,A*RR(II),A*SS(II),BJ,DUMMY)
         CALL MI_LOCAL_BASIS_AT_Q8(XYZN,A*RR(II),A*SS(II),P1,P2,P3,DUMMY)
         DO DD=1,48
            DO JJ=1,3
               DO KK=1,3
                  GT(JJ,KK,DD)=GT(JJ,KK,DD)+N(II)*(P1(JJ)*P1(KK)*BJ(1,DD) &
                       +P2(JJ)*P2(KK)*BJ(2,DD) &
                       +0.5D0*(P1(JJ)*P2(KK)+P2(JJ)*P1(KK))*BJ(3,DD))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL MI_LOCAL_BASIS_AT_Q8(XYZN,R,S,E1,E2,E3,JAC)
      DO DD=1,48
         T=GT(:,:,DD)
         BMOUT(1,DD)=DOT_PRODUCT(E1,MATMUL(T,E1))
         BMOUT(2,DD)=DOT_PRODUCT(E2,MATMUL(T,E2))
         BMOUT(3,DD)=TWO*DOT_PRODUCT(E1,MATMUL(T,E2))
      ENDDO
      END SUBROUTINE MI_BM_Q8_AT

      SUBROUTINE MI_BM_DIRECT_Q8_AT ( XYZN, XI, ETA, BMOUT, JAC )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), XI, ETA
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), C(4), VP(3,3)
      INTEGER(LONG) :: II, COL
      CALL MI_SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL MI_COV_MAP_Q8(XYZN, XI, ETA, G1, G2, C, JAC)
      BMOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         CALL MI_TENSOR_PHYS_Q8(DN(1,II)*G1, DN(2,II)*G2, 0.5D0*(DN(1,II)*G2 + DN(2,II)*G1), C, VP)
         BMOUT(1:3,COL+1:COL+3) = VP
      ENDDO
      END SUBROUTINE MI_BM_DIRECT_Q8_AT

      SUBROUTINE MI_BB_Q8_AT ( XYZN, NORMS, XI, ETA, BBOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BBOUT(3,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), C(4), VP(3,3)
      REAL(DOUBLE) :: T1(3), T2(3), T0(3), CG1(3), CG2(3)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL MI_SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL MI_COV_MAP_Q8(XYZN, XI, ETA, G1, G2, C, JAC)
      T1 = MATMUL(DN(1,:), NORMS_LOC)
      T2 = MATMUL(DN(2,:), NORMS_LOC)
      BBOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         T0 = NORMS_LOC(II,:)
         CALL MI_TENSOR_PHYS_Q8(DN(1,II)*T1, DN(2,II)*T2, 0.5D0*(DN(1,II)*T2 + DN(2,II)*T1), C, VP)
         BBOUT(1:3,COL+1:COL+3) = -VP
         CALL MI_CROSS3(G1, T0, CG1)
         CALL MI_CROSS3(G2, T0, CG2)
         CALL MI_TENSOR_PHYS_Q8(DN(1,II)*CG1, DN(2,II)*CG2, 0.5D0*(DN(1,II)*CG2 + DN(2,II)*CG1), C, VP)
         BBOUT(1:3,COL+4:COL+6) = VP
      ENDDO
      END SUBROUTINE MI_BB_Q8_AT

      SUBROUTINE MI_BS_Q8_AT(XYZN,NORMS,XI,ETA,BSOUT,JAC)
      REAL(DOUBLE),INTENT(IN) :: XYZN(8,3),NORMS(8,3),XI,ETA
      REAL(DOUBLE),INTENT(OUT) :: BSOUT(2,48),JAC
      REAL(DOUBLE) :: G1(3),G2(3),C(4)
      CALL MI_COV_MAP_Q8(XYZN,XI,ETA,G1,G2,C,JAC)
      CALL MI_ANS_BS_Q8(XYZN,NORMS,XI,ETA,C,BSOUT)
      END SUBROUTINE MI_BS_Q8_AT

      SUBROUTINE MI_BDRILL_Q8_AT ( XYZN, NORMS, XI, ETA, BDOUT, JAC, NORMS_EXT )
      REAL(DOUBLE), INTENT(IN)  :: XYZN(8,3), NORMS(8,3), XI, ETA
      REAL(DOUBLE), INTENT(IN), OPTIONAL :: NORMS_EXT(8,3)
      REAL(DOUBLE), INTENT(OUT) :: BDOUT(1,48), JAC
      REAL(DOUBLE) :: NVAL(8), DN(2,8), G1(3), G2(3), E1(3), E2(3), E3(3)
      REAL(DOUBLE) :: A(2,2), AINV(2,2), DLOC(2,8), DX(8), DY(8)
      REAL(DOUBLE) :: NORMS_LOC(8,3)
      INTEGER(LONG) :: II, COL
      NORMS_LOC = NORMS
      IF (PRESENT(NORMS_EXT)) NORMS_LOC = NORMS_EXT
      CALL MI_SHAPE_Q8(XI, ETA, NVAL, DN)
      CALL MI_LOCAL_BASIS_AT_Q8(XYZN, XI, ETA, E1, E2, E3, JAC)
      G1 = MATMUL(DN(1,:), XYZN)
      G2 = MATMUL(DN(2,:), XYZN)
      A(1,1) = DOT_PRODUCT(G1,G1)
      A(1,2) = DOT_PRODUCT(G1,G2)
      A(2,1) = DOT_PRODUCT(G1,G2)
      A(2,2) = DOT_PRODUCT(G2,G2)
      CALL MI_INV2(A, AINV)
      DLOC = MATMUL(AINV, DN)
      DX=DLOC(1,:)*DOT_PRODUCT(G1,E1)+DLOC(2,:)*DOT_PRODUCT(G2,E1)
      DY=DLOC(1,:)*DOT_PRODUCT(G1,E2)+DLOC(2,:)*DOT_PRODUCT(G2,E2)
      BDOUT = ZERO
      DO II=1,8
         COL = (II-1)*6
         BDOUT(1,COL+1:COL+3) = 0.5D0*(DX(II)*E2 - DY(II)*E1)
         BDOUT(1,COL+4:COL+6) = BDOUT(1,COL+4:COL+6) - NVAL(II)*E3
      ENDDO
      END SUBROUTINE MI_BDRILL_Q8_AT

      SUBROUTINE MI_SET_NORMAL_SIGN_Q8(XYZN)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3)
      REAL(DOUBLE) :: N(8),DN(2,8),NC(3),G1(3),G2(3)
      CALL MI_GEOM_SHAPE_Q8(ZERO,ZERO,N,DN)
      G1=MATMUL(DN(1,:),XYZN)
      G2=MATMUL(DN(2,:),XYZN)
      CALL MI_CROSS3(G1,G2,NC)
      NORMAL_SIGN=ONE
      IF (NC(3) < -1.0D-6*MI_VNORM(NC)) NORMAL_SIGN=-ONE
      END SUBROUTINE

      SUBROUTINE MI_ANS_ROWS_Q8(XYZN,NORMS,R,S,EM,ES)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S
      REAL(DOUBLE), INTENT(OUT) :: EM(3,48),ES(2,48)
      REAL(DOUBLE) :: N(8),DN(2,8),G1(3),G2(3),T0(3),CV(3)
      INTEGER(LONG) :: II,COL
      CALL MI_SHAPE_Q8(R,S,N,DN)
      G1=MATMUL(DN(1,:),XYZN)
      G2=MATMUL(DN(2,:),XYZN)
      T0=MATMUL(N,NORMS)
      EM=ZERO
      ES=ZERO
      DO II=1,8
         COL=6*(II-1)
         EM(1,COL+1:COL+3)=DN(1,II)*G1
         EM(2,COL+1:COL+3)=DN(2,II)*G2
         EM(3,COL+1:COL+3)=0.5D0*(DN(1,II)*G2+DN(2,II)*G1)
         ES(1,COL+1:COL+3)=DN(1,II)*T0
         ES(2,COL+1:COL+3)=DN(2,II)*T0
         CALL MI_CROSS3(NORMS(II,:),G1,CV)
         ES(1,COL+4:COL+6)=N(II)*CV
         CALL MI_CROSS3(NORMS(II,:),G2,CV)
         ES(2,COL+4:COL+6)=N(II)*CV
      ENDDO
      END SUBROUTINE

      SUBROUTINE MI_ANS_INTERP_Q8(XYZN,NORMS,R,S,COMP,ISMEM,ROW)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S
      INTEGER(LONG), INTENT(IN) :: COMP
      LOGICAL, INTENT(IN) :: ISMEM
      REAL(DOUBLE), INTENT(OUT) :: ROW(48)
      REAL(DOUBLE) :: A,B,H,PA(2),PB(3),LR(2),LS(2),QR(3),QS(3),EM(3,48),ES(2,48),P,Q,W
      INTEGER(LONG) :: II,JJ,NJ
      A=ONE/DSQRT(3.0D0)
      B=DSQRT(0.6D0)
      H=A
      IF (ISMEM .AND. COMP /= 3) H=ONE
      PA=(/-H,H/)
      PB=(/-B,ZERO,B/)
      LR=(/0.5D0*(ONE-R/H),0.5D0*(ONE+R/H)/)
      LS=(/0.5D0*(ONE-S/H),0.5D0*(ONE+S/H)/)
      QR=(/R*(R-B)/(TWO*B*B),ONE-R*R/(B*B),R*(R+B)/(TWO*B*B)/)
      QS=(/S*(S-B)/(TWO*B*B),ONE-S*S/(B*B),S*(S+B)/(TWO*B*B)/)
      IF (.NOT.ISMEM .AND. ANSSHEAR == 'BDG4') THEN
         PB(1:2)=(/-ONE,ONE/)
         QR(1:2)=(/0.5D0*(ONE-R),0.5D0*(ONE+R)/)
         QS(1:2)=(/0.5D0*(ONE-S),0.5D0*(ONE+S)/)
      ENDIF
      ROW=ZERO
      NJ=3
      IF (COMP == 3 .OR. (.NOT.ISMEM .AND. ANSSHEAR == 'BDG4')) NJ=2
      DO II=1,2
         DO JJ=1,NJ
            IF (COMP == 1) THEN
               P=PA(II); Q=PB(JJ); W=LR(II)*QS(JJ)
            ELSE IF (COMP == 2) THEN
               P=PB(JJ); Q=PA(II); W=LS(II)*QR(JJ)
            ELSE
               P=PA(II); Q=PA(JJ); W=LR(II)*LS(JJ)
            ENDIF
            CALL MI_ANS_ROWS_Q8(XYZN,NORMS,P,Q,EM,ES)
            IF (ISMEM) THEN
               ROW=ROW+W*EM(COMP,:)
            ELSE
               ROW=ROW+W*ES(COMP,:)
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE

      SUBROUTINE MI_ANS_BM_Q8(XYZN,R,S,C,BMOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),R,S,C(4)
      REAL(DOUBLE), INTENT(OUT) :: BMOUT(3,48)
      REAL(DOUBLE) :: ROWS(3,48),NORMS(8,3)
      INTEGER(LONG) :: II
      NORMS=ZERO
      DO II=1,3
         CALL MI_ANS_INTERP_Q8(XYZN,NORMS,R,S,II,.TRUE.,ROWS(II,:))
      ENDDO
      BMOUT(1,:)=C(1)**2*ROWS(1,:)+C(3)**2*ROWS(2,:)+TWO*C(1)*C(3)*ROWS(3,:)
      BMOUT(2,:)=C(2)**2*ROWS(1,:)+C(4)**2*ROWS(2,:)+TWO*C(2)*C(4)*ROWS(3,:)
      BMOUT(3,:)=TWO*(C(1)*C(2)*ROWS(1,:)+C(3)*C(4)*ROWS(2,:)+(C(1)*C(4)+C(3)*C(2))*ROWS(3,:))
      END SUBROUTINE

      SUBROUTINE MI_ANS_BS_Q8(XYZN,NORMS,R,S,C,BSOUT)
      REAL(DOUBLE), INTENT(IN) :: XYZN(8,3),NORMS(8,3),R,S,C(4)
      REAL(DOUBLE), INTENT(OUT) :: BSOUT(2,48)
      REAL(DOUBLE) :: ROWS(2,48)
      INTEGER(LONG) :: II
      DO II=1,2
         CALL MI_ANS_INTERP_Q8(XYZN,NORMS,R,S,II,.FALSE.,ROWS(II,:))
      ENDDO
      BSOUT(1,:)=C(1)*ROWS(1,:)+C(3)*ROWS(2,:)
      BSOUT(2,:)=C(2)*ROWS(1,:)+C(4)*ROWS(2,:)
      END SUBROUTINE

      SUBROUTINE MI_CROSS3 ( A, B, C )
      REAL(DOUBLE), INTENT(IN)  :: A(3), B(3)
      REAL(DOUBLE), INTENT(OUT) :: C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      END SUBROUTINE MI_CROSS3

      FUNCTION MI_VNORM ( V ) RESULT(NM)
      REAL(DOUBLE), INTENT(IN) :: V(3)
      REAL(DOUBLE) :: NM
      NM = DSQRT(MAX(ZERO, DOT_PRODUCT(V,V)))
      END FUNCTION MI_VNORM

      SUBROUTINE MI_INV2 ( A, AINV )
      REAL(DOUBLE), INTENT(IN)  :: A(2,2)
      REAL(DOUBLE), INTENT(OUT) :: AINV(2,2)
      REAL(DOUBLE) :: DET
      DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      IF (DABS(DET) < 1.0D-20) THEN
         AINV = ZERO
         AINV(1,1) = ONE
         AINV(2,2) = ONE
      ELSE
         AINV(1,1) =  A(2,2)/DET
         AINV(1,2) = -A(1,2)/DET
         AINV(2,1) = -A(2,1)/DET
         AINV(2,2) =  A(1,1)/DET
      ENDIF
      END SUBROUTINE MI_INV2


      SUBROUTINE MI_FIELD_COEFF_Q8(XY)
      REAL(DOUBLE),INTENT(IN) :: XY(8,2)
      REAL(DOUBLE) :: D(4),V(2),W(2),DEN
      INTEGER(LONG) :: II,JJ,KK,MM
      DO II=1,4
         JJ=MOD(II,4)+1
         MM=MOD(II+2,4)+1
         V=XY(JJ,:)-XY(II,:)
         W=XY(MM,:)-XY(II,:)
         D(II)=V(1)*W(2)-V(2)*W(1)
      ENDDO
      DO II=1,4
         KK=MOD(II+1,4)+1
         MM=MOD(II+2,4)+1
         DEN=D(II)+D(KK)
         FIELD_COEFF(II)=-0.25D0+(D(II)-D(KK))/(8.0D0*DEN)
         FIELD_COEFF(II+4)=0.5D0+(D(MM)-D(II))/(4.0D0*DEN)
      ENDDO
      END SUBROUTINE

      SUBROUTINE MI_SHAPE_Q8(R,S,N,DN)
      REAL(DOUBLE),INTENT(IN) :: R,S
      REAL(DOUBLE),INTENT(OUT) :: N(8),DN(2,8)
      REAL(DOUBLE) :: LR(3),LS(3),DR(3),DS(3),N9
      INTEGER(LONG) :: II,IR(8),IS(8)
      LR=(/0.5D0*R*(R-ONE),ONE-R*R,0.5D0*R*(R+ONE)/)
      LS=(/0.5D0*S*(S-ONE),ONE-S*S,0.5D0*S*(S+ONE)/)
      DR=(/R-0.5D0,-TWO*R,R+0.5D0/)
      DS=(/S-0.5D0,-TWO*S,S+0.5D0/)
      IR=(/1,3,3,1,2,3,2,1/)
      IS=(/1,1,3,3,1,2,3,2/)
      N9=(ONE-R*R)*(ONE-S*S)
      DO II=1,8
         N(II)=LR(IR(II))*LS(IS(II))+N9*FIELD_COEFF(II)
         DN(1,II)=DR(IR(II))*LS(IS(II))-TWO*R*(ONE-S*S)*FIELD_COEFF(II)
         DN(2,II)=LR(IR(II))*DS(IS(II))-TWO*S*(ONE-R*R)*FIELD_COEFF(II)
      ENDDO
      END SUBROUTINE

      END SUBROUTINE MITC8
