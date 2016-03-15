! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  Johns Hopkins University
!
!  VICAR3Dp, a parallelized version of VICAR3D.
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Rajat Mittal
!  Fady Najjar
!
!  Other contributing programmers:
!     Haibo Dong
!     Haoxiang Luo
!     Meliha Bozkurttas
!     Rupesh Babu K. A.
!     Xudong Zheng 
!     Reza Ghias 
!     S. A. Mohsen Karimian
!
!  Filename: BOUNDARY_SET.PAR.F90
!  Latest Modification: Jan 20, 2010 (ver. P1.5.5)
!  Made by X.Zheng
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     set_boundary()
!     set_outer_iup()
!     set_internal_boundary()
!     set_outflow_area()
!     set_iblank_slow(nBody_Begin,nBody_End)
!     find_isolatedIblank(fix)
!     find_iblankHoles(fix)
!     find_iblankHoles_outerboundary(fix,iBody,iblanktmp)
!     identify_freshcells()
!     set_iblank_general_body()
!     get_scanRange(scanRange)
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE set_boundary(bsCommTime)

! ------------------------------------------------------
!  This subroutine sets the inner and outer boundaries.
! ------------------------------------------------------

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE
    
!   Parameters
!   ----------
    REAL(KIND=CGREAL),INTENT (OUT) :: bsCommTime



    bsCommTime=0.d0   ! For Now!!!!

!   set outer boundary
!   ------------------
    IF (monitorON) WRITE(STDOUT,'(5X,A)') 'Outer Boundary'
    CALL set_outer_iup()  !Only periodic is left...

!   set internal boundary
!   ---------------------
    IF ( boundary_formulation /= NO_INTERNAL_BOUNDARY ) THEN

      IF (monitorON) WRITE(STDOUT,'(5X,A)') 'Internal Boundary'
      CALL set_internal_boundary()!..............................COMPLETE(SAMK)

!!      IF ( boundary_formulation == SSM_METHOD .and. flow_type == POTENTIAL_FLOW) THEN
!!        CALL set_internal_potential_flag()
!!      END IF ! boundary_formulation

    ENDIF ! boundary_formulation 

!   calculate outflow area
!   ----------------------
    CALL set_outflow_area()!.....................................COMPLETE(SAMK)

END SUBROUTINE set_boundary
!---------------------------------------------------------------------



SUBROUTINE set_outer_iup()

! ----------------------------------------------------------------------
!  This subroutine initializes the outer boundaries. For all boundaries
!  xup and xum are equal to one.
! ----------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE



    ium = 0.0_CGREAL
    jum = 0.0_CGREAL
    kum = 0.0_CGREAL
    iup = 0.0_CGREAL
    jup = 0.0_CGREAL
    kup = 0.0_CGREAL

    iumm = 0.0_CGREAL    !
    jumm = 0.0_CGREAL    !
    kumm = 0.0_CGREAL    !
    iupp = 0.0_CGREAL    !
    jupp = 0.0_CGREAL    !
    kupp = 0.0_CGREAL    !   Added by Rupesh (used for 2nd upwinding)

!   West and East Boundaries
!   ------------------------
    IF (myCoords(1)==0) THEN
      ium( 1,:,:) = 1.0_CGREAL
      iumm(1,:,:) = 1.0_CGREAL  !
      iumm(2,:,:) = 1.0_CGREAL  !
    END IF
    IF (myCoords(1)==Np(1)-1) THEN
      iup( nxc  ,:,:) = 1.0_CGREAL
      iupp(nxc  ,:,:) = 1.0_CGREAL  !   Added by Rupesh (used for 2nd upwinding)   
      iupp(nxc-1,:,:) = 1.0_CGREAL  !   
    END IF

!   South and North Boundaries
!   --------------------------
    IF (myCoords(2)==0) THEN
      jum( :,1,:) = 1.0_CGREAL
      jumm(:,1,:) = 1.0_CGREAL  !
      jumm(:,2,:) = 1.0_CGREAL  !
    END IF
    IF (myCoords(2)==Np(2)-1) THEN
      jup( :,nyc  ,:) = 1.0_CGREAL
      jupp(:,nyc  ,:) = 1.0_CGREAL  !   Added by Rupesh (used for 2nd upwinding)
      jupp(:,nyc-1,:) = 1.0_CGREAL  !
    END IF

!   Bottom and Top Boundaries
!   -------------------------
    IF (nDim == DIM_3D) THEN
      kum(:,:,1  ) = 1.0_CGREAL
      kup(:,:,nzc) = 1.0_CGREAL

      kumm(:,:,1    ) = 1.0_CGREAL  !
      kumm(:,:,2    ) = 1.0_CGREAL  !
      kupp(:,:,nzc  ) = 1.0_CGREAL  !   Added by Rupesh (used for 2nd upwinding)
      kupp(:,:,nzc-1) = 1.0_CGREAL  !   
    ENDIF

! new changes for the periodic condition
! --------------------------------------
!!    IF (bcx1 == BC_TYPE_PERIODIC .OR.   &
!!        bcy1 == BC_TYPE_PERIODIC .OR.   &
!!        bcz1 == BC_TYPE_PERIODIC  ) THEN 
!!!!     it_solver_type == IT_SOLVER_TYPE_MG) THEN

!!      CALL remove_up_um
!!    END IF

    IF (bcx1 == BC_TYPE_PERIODIC .OR.   &                   !Ehsan added for periodic BC
        bcy1 == BC_TYPE_PERIODIC .OR.   &                  
        bcz1 == BC_TYPE_PERIODIC  ) THEN                    
      CALL remove_up_um                                     
    END IF                                                  


END SUBROUTINE set_outer_iup
!---------------------------------------------------------------------



SUBROUTINE set_internal_boundary()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER :: i, j, k ,iBody

    num_fresh  = 0
    num_dead   = 0

!   Initialize values
!   -----------------
    DO k=0,nzc
    DO j=myJLL,myJUL
    DO i=myILL,myIUL
      bodyNumOld(i,j,k)    = bodyNum(i,j,k)*iblank(i,j,k)
      bodyNum(i,j,k) = 0
      fresh_cell(i,j,k) = 0
      ghostCellMark(i,j,k) = 0
    ENDDO
    ENDDO
    ENDDO
   


!  set iblank at current time step 
!  -------------------------------
    SELECT CASE(body_type)

    CASE(CANONICAL)
      IF(nBody_solid /=0) THEN

        IF (monitorON) THEN
          WRITE(STDOUT,'(7X,A)') 'Solid Bodies'
          WRITE(STDOUT,'(7X,A)') '============'
          WRITE(STDOUT,*)
        END IF

        iblank(:,:,:)  = iblank_solid(:,:,:)

        SELECT CASE(iblankFast)

        CASE(IBLANK_SLOW)
          IF (ImtheBOSS) WRITE(STDOUT,*) 'IBLANK_SLOW is no more supported in the parallel version.'
          CALL flow_stop()
          STOP

        CASE(IBLANK_FAST)

         CALL set_iblank_canonical_body_fast(1, nBody_Solid)
        END SELECT ! iblankFast

        CALL find_IsolatedIblank(1)
!        CALL find_iblankHoles(1)

        iblank_solid(:,:,:) = iblank(:,:,:)

        CALL identify_ghostcells_solid()
        CALL set_internal_iup_solid()

        iblank(:,:,:) = 0

      END IF

      IF(nBody_membrane /=0) THEN

        IF (monitorON) THEN
          WRITE(STDOUT,*)
          WRITE(STDOUT,'(7X,A)') 'Membranes'
          WRITE(STDOUT,'(7X,A)') '========='
          WRITE(STDOUT,*)
        END IF


        DO iBody=nBody_solid+1, nBody

          iblank(:,:,:) = iblank_memb(:,:,:,iBody-nBody_solid)

          SELECT CASE(iblankFast)

          CASE(IBLANK_SLOW)
            IF (ImtheBOSS) WRITE(STDOUT,*) 'IBLANK_SLOW is no more supported in the parallel version.'
            CALL flow_stop()
            STOP

          CASE(IBLANK_FAST)
            CALL set_iblank_canonical_body_fast(iBody,iBody)
          END SELECT ! iblankFast

          CALL find_IsolatedIblank(1)
!          CALL find_iblankHoles(1)

          iblank_memb(:,:,:,iBody-nBody_solid) = iblank(:,:,:)

          CALL identify_ghostcells_membrane(iBody)
          CALL set_internal_iup_membrane(iBody)
          CALL identify_gc_membrane_final(iBody)

          iblank(:,:,:) = 0
  
          IF (monitorON) WRITE(STDOUT,*)
        ENDDO
      END IF

      CALL identify_freshcells()

      iblank(:,:,:) = iblank_solid(:,:,:)

      IF ( boundary_formulation == GCM_METHOD ) THEN
        CALL GCM_set_internal_boundary()
      END IF ! boundary_formulation

    CASE(GENERAL)
      CALL set_iblank_general_body()

    END SELECT ! body_type

    IF ( boundary_formulation == SSM_METHOD ) CALL SSM_set_internal_area()

END SUBROUTINE set_internal_boundary
!---------------------------------------------------------------------



SUBROUTINE set_outflow_area()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
    USE scalar
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     

    IMPLICIT NONE

    INTEGER :: i, j, k, m



!   compute areas of external boundaries
!   ------------------------------------
    area_left    = 0.0_CGREAL
    area_right   = 0.0_CGREAL
    area_bot     = 0.0_CGREAL
    area_top     = 0.0_CGREAL
    area_back    = 0.0_CGREAL
    area_front   = 0.0_CGREAL
    outflow_area = 0.0_CGREAL
    outflow_area_FP = 0.0_CGREAL
    prim_left    = 0.0_CGREAL

    IF (myCoords(1)==0) THEN
      DO k=1,nzc
      DO j=1,nyc
        area_left  = area_left  + dy(L2GJ(j))*dz(k)*REAL(1-iblank(1  ,j,k),KIND=CGREAL)

        prim_left  = prim_left  + REAL(1-iblank(1,j,k),KIND=CGREAL) *       &
                                  SQRT(dy(L2GJ(j))*dy(L2GJ(j))*(jum(1,j,k)+jup(1,j,k))  &
                                      +dz(     k )*dz(     k )*(kum(1,j,k)+kup(1,j,k)))

      ENDDO
      ENDDO
    END IF
    
    IF (myCoords(1)==Np(1)-1) THEN
      DO k=1,nzc
      DO j=1,nyc
        area_right = area_right + dy(L2GJ(j))*dz(k)*REAL(1-iblank(nxc,j,k),KIND=CGREAL)
      ENDDO
      ENDDO
    END IF

    IF (myCoords(2)==0) THEN
      DO k=1,nzc
      DO i=1,nxc
        area_bot   = area_bot   + dx(L2GI(i))*dz(k)*REAL(1-iblank(i,1  ,k),KIND=CGREAL)
      ENDDO
      ENDDO
    END IF

    IF (myCoords(2)==Np(2)-1) THEN
      DO k=1,nzc
      DO i=1,nxc
        area_top   = area_top   + dx(L2GI(i))*dz(k)*REAL(1-iblank(i,nyc,k),KIND=CGREAL)
      ENDDO
      ENDDO
    END IF

    DO j=1,nyc
    DO i=1,nxc
      area_back  = area_back  + dx(L2GI(i))*dy(L2GJ(j))*REAL(1-iblank(i,j,1  ),KIND=CGREAL)
      area_front = area_front + dx(L2GI(i))*dy(L2GJ(j))*REAL(1-iblank(i,j,nzc),KIND=CGREAL)
    ENDDO
    ENDDO

#   ifdef MPI
      CALL par_getSumReal(area_left )
      CALL par_getSumReal(area_right)
      CALL par_getSumReal(area_bot  )
      CALL par_getSumReal(area_top  )
      CALL par_getSumReal(area_back )
      CALL par_getSumReal(area_front)
      CALL par_getSumReal(prim_left )
#   endif

!!   calculate inflow area without using IF statements
!!   -------------------------------------------------
!    inflow_area  = area_left

!   calculate outflow area without using IF statements
!   --------------------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    IF (bcx1 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_left
    IF (bcx2 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_right
    IF (bcy1 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_bot
    IF (bcy2 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_top
    IF (bcz1 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_back
    IF (bcz2 == BC_TYPE_ZERO_GRADIENT) outflow_area = outflow_area + area_front
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
    !IF (bcx1 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_left
    !IF (bcx2 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_right
    !IF (bcy1 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_bot
    !IF (bcy2 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_top
    !IF (bcz1 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_back
    !IF (bcz2 == BC_TYPE_ZERO_GRADIENT .OR. VDVactive == 1) outflow_area_FP = outflow_area_FP + area_front
    outflow_area_FP = area_left+area_right+area_bot+area_top+area_back+area_front
END SUBROUTINE set_outflow_area
!---------------------------------------------------------------------



SUBROUTINE find_isolatedIblank(fix)

! -------------------------------------------------------------------------
!  This subroutine checks for isolated iblanked cells and eliminates them.
!  such isolated cells can be formed due to roundoff error.
! -------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: fix

    INTEGER :: iMin,iMax,jMin,jMax,kMin,kMax
    INTEGER :: i,j,k,scanRange,sumIblank

    REAL(KIND=CGREAL) :: myTime



!    scanRange = 1

!   ----------------------------------------------------------------------------------
!    Searching for isolated blank cells

!    SAMK: Using myImin, myImax, myJmin and myJmax will for a safety net in searching 
!          for isolated iblanks. This method requires a communication phase for the
!          outermost layer of the domain.
!   ----------------------------------------------------------------------------------
    DO k=1,nzc
    DO j=myJmin,myJmax
    DO i=myImin,myImax

      IF (iBlank(i,j,k) == 1) THEN
        scanRange = 1
        iMin = MAX(myILL,i-scanRange)
        iMax = MIN(myIUL,i+scanRange)
        jMin = MAX(myJLL,j-scanRange)
        jMax = MIN(myJUL,j+scanRange)
        kMin = MAX(1    ,k-scanRange)
        kMax = MIN(nzc  ,k+scanRange)

        sumIblank = 0    

        sumIblank = SUM(iBlank(iMin:iMax, &
                               jMin:jMax, &
                               kMin:kMax))

        IF (nDim == DIM_2D) sumIblank = sumIblank/2
        
        IF (sumIblank <= 1) THEN
          WRITE(ifuParLog,'(7X,A,3I5)') 'Found Isolated blank Cell at ',L2GI(i),L2GJ(j),k
          IF (fix == 1) THEN
            WRITE(ifuParLog,'(7X,A,3I5)') 'This cell is unblanked'
            iBlank(i,j,k)  = 0
            bodyNum(i,j,k) = 0
          ENDIF
        ENDIF 
      ENDIF

    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

#   ifdef MPI
!     Communicating the outmost layer of the subdomain for iBlank.
!     ------------------------------------------------------------
      CALL par_comm_outermost_int(iblank ,myTime)
      CALL par_comm_outermost_int(bodyNum,myTime)
#   endif

   END SUBROUTINE find_isolatedIblank
!------------------------------------------------------------------------------



SUBROUTINE find_iblankHoles(fix)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: fix

    INTEGER :: i, j, k
    
    LOGICAL :: withinBoundary

    REAL(KIND=CGREAL) :: myTime



!   finding holes
!   -------------
    IF ( boundary_formulation == SSM_METHOD ) THEN

      DO k=1,nzc
      DO j=myJmin,myJmax
      DO i=myImin,myImax
      
        IF (i>=1 .AND. i<=nxc .AND. &
            j>=1 .AND. j<=nyc) THEN
          withinBoundary=.TRUE.
        ELSE
          withinBoundary=.FALSE.
        END IF

        IF (iBlank(i,j,k) == 0) THEN

          IF ( (iblank(i+1,j,k) == 1 .OR. L2GI(i) == nxc_GLBL) .AND. (iblank(i-1,j,k) == 1 .OR. L2GI(i) == 1) ) THEN
           WRITE(ifuParLog,'(7X,A,3I5)') 'flow hole-1 at',L2GI(i),L2GJ(j),k
           IF (fix == 1) THEN
             WRITE(ifuParLog,'(7X,A)') 'closing hole'
             iblank(i,j,k) = 1
             IF (fresh_cell(i,j,k) == 1 .AND. withinBoundary) THEN
               fresh_cell(i,j,k) = 0
             ENDIF
           ENDIF
          ENDIF

          IF ( (iblank(i,j+1,k) == 1 .OR. L2GJ(j) == nyc_GLBL) .AND. (iblank(i,j-1,k) == 1 .OR. L2GJ(j) == 1) ) THEN
           WRITE(ifuParLog,'(7X,A,3I5)') 'flow hole-2 at',L2GI(i),L2GJ(j),k
           IF (fix == 1) THEN
             WRITE(ifuParLog,'(7X,A)') 'closing hole'
             iblank(i,j,k) = 1
             IF (fresh_cell(i,j,k) == 1 .AND. withinBoundary) THEN
               fresh_cell(i,j,k) = 0
             ENDIF
           ENDIF
          ENDIF

          IF ( (iblank(i,j,k+1) == 1 .OR. k == nzc) .AND. (iblank(i,j,k-1) == 1 .OR. k == 1) ) THEN
           WRITE(ifuParLog,'(7X,A,3I5)') 'flow hole-3 at',L2GI(i),L2GJ(j),k
           IF (fix == 1) THEN
             WRITE(ifuParLog,'(7X,A)') 'closing hole'
             iblank(i,j,k) = 1
             IF (fresh_cell(i,j,k) == 1 .AND. withinBoundary) THEN
               fresh_cell(i,j,k) = 0
             ENDIF
           ENDIF
          ENDIF

        ENDIF

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

#     ifdef MPI
!       Communicating the outmost layer of the subdomain for iBlank.
!       ------------------------------------------------------------
        CALL par_comm_outermost_int(iblank    ,myTime)
        CALL par_comm_outermost_int(fresh_cell,myTime)
#     endif

    ENDIF

END SUBROUTINE find_iblankHoles
!---------------------------------------------------------------------

SUBROUTINE find_iblankHoles_outerboundary(fix,iBody,iblanktmp)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: fix,iBody
    INTEGER , DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT(INOUT) :: iblanktmp

    INTEGER              :: i,j,k

!   finding holes
!   -------------
    IF ( boundary_formulation == GCM_METHOD ) RETURN

!   First check for holes on outer boundaries (e.g. i=1 and nx_GLBL -1)
!   Then look for holes in the domain interior

    DO k=1,nz-1
    DO j=myJmin, myJmax
    DO i=myImin, myImax
     IF(L2GI(i) < 1 .OR. L2GI(i) > nx_GLBL-1 .OR. L2GJ(j) < 1 .OR. L2GJ(j) > ny_GLBL-1) CYCLE
      IF (iblanktmp(i,j,k) == 0 ) THEN

!       i-boundary..............
        IF (    ( L2GI(i) == nx_GLBL-1 .AND. iblanktmp(i-1,j,k) == 1 ) &
           .OR. ( L2GI(i) == 1         .AND. iblanktmp(i+1,j,k) == 1 ) ) THEN
!         PRINT*,'flow hole-1 at',i,j,k,'at processor ', myrank
         IF (fix == 1) THEN
!           PRINT*,'closing outer boundary hole'
           iblanktmp(i,j,k) = 1
           IF( L2GI(i) == nx_GLBL-1 ) bodynum(i,j,k) = bodynum(i-1,j,k)
           IF( L2GI(i) == 1 )          bodynum(i,j,k) = bodynum(i+1,j,k)
         ENDIF
        ENDIF ! L2GI

!       j-boundary..............
        IF (    ( L2GJ(j) == ny_GLBL-1 .AND. iblanktmp(i,j-1,k) == 1 ) &
           .OR. ( L2GJ(j) == 1         .AND. iblanktmp(i,j+1,k) == 1 ) ) THEN
!          PRINT*,'flow hole-2 at',i,j,k,'processor', myrank
          IF (fix == 1) THEN
!            PRINT*,'closing outer boundary hole'
            iblanktmp(i,j,k) = 1
            IF( L2GJ(j) == ny_GLBL-1 ) bodynum(i,j,k) = bodynum(i,j-1,k)
            IF( L2GJ(j) == 1 )         bodynum(i,j,k) = bodynum(i,j+1,k)
          ENDIF
        ENDIF ! L2GJ

!       k-boundary..............
        IF (    ( k == nz-1 .AND. iblanktmp(i,j,k-1) == 1 ) &
           .OR. ( k == 1    .AND. iblanktmp(i,j,k+1) == 1 ) ) THEN
!         PRINT*,'flow hole-3 at',i,j,k, 'at processor ', myrank
         IF (fix == 1) THEN
!           PRINT*,'closing outer boundary hole'
           iblanktmp(i,j,k) = 1
           IF(k == nz-1) bodynum(i,j,k) = bodynum(i,j,k-1)
           IF(k == 1)    bodynum(i,j,k) = bodynum(i,j,k+1)
         ENDIF
        ENDIF ! k

      ENDIF ! iblanktmp

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!   Fix holes for membranes
!    Used only when Membrane body is only present.

    IF (unstruc_surface_type(iBody) == MEMBRANE) THEN

      DO k=1,nz-1
      DO j=myJmin, myJmax
      DO i=myImin, myImax
       IF (iblanktmp(i,j,k) == 1 ) THEN

!       i-boundary..............
        IF (    ( L2GI(i) == nx_GLBL-1 .AND. iblanktmp(i-1,j,k) == 0 ) &
           .OR. ( L2GI(i) == 1         .AND. iblanktmp(i+1,j,k) == 0 ) ) THEN
!         PRINT*,'flow hole-1 at',i,j,k,'at processor ', myrank
         IF (fix == 1) THEN
!           PRINT*,'closing outer boundary hole'
           iblanktmp(i,j,k) = 1
           IF( L2GI(i) == nx_GLBL-1) bodynum(i,j,k) = bodynum(i-1,j,k)
           IF( L2GI(i) == 1        ) bodynum(i,j,k) = bodynum(i+1,j,k)
         ENDIF
        ENDIF ! L2GI

!       j-boundary..............
        IF (    ( L2GJ(j) == ny_GLBL-1 .AND. iblanktmp(i,j-1,k) == 0 ) &
           .OR. ( L2GJ(j) == 1         .AND. iblanktmp(i,j+1,k) == 0 ) ) THEN
!          PRINT*,'flow hole-2 at',i,j,k,'at processor ', myrank
          IF (fix == 1) THEN
!            PRINT*,'closing outer boundary hole'
            iblanktmp(i,j,k) = 1
            IF(L2GJ(j) == ny_GLBL-1) bodynum(i,j,k) = bodynum(i,j-1,k)
            IF(L2GJ(j) == 1        ) bodynum(i,j,k) = bodynum(i,j+1,k)
          ENDIF
        ENDIF ! L2GJ

!       k-boundary..............
        IF (    ( k == nz-1 .AND. iblanktmp(i,j,k-1) == 0 ) &
           .OR. ( k == 1    .AND. iblanktmp(i,j,k+1) == 0 ) ) THEN
!         PRINT*,'flow hole-3 at',i,j,k,'at processor ', myrank
         IF (fix == 1) THEN
!           PRINT*,'closing outer boundary hole'
           iblanktmp(i,j,k) = 1
           IF(k == nz-1) bodynum(i,j,k) = bodynum(i,j,k-1)
           IF(k == 1   ) bodynum(i,j,k) = bodynum(i,j,k+1)
         ENDIF
        ENDIF ! k

       ENDIF ! iblanktmp

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k
    ENDIF ! Membrane

END SUBROUTINE find_iblankHoles_outerboundary
!---------------------------------------------------------------------


SUBROUTINE identify_freshcells()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER :: i, j, k
    
    LOGICAL :: withinBoundary

    num_fresh = 0

    WRITE(ifuParLog,*)
    WRITE(ifuParLog,'(A,I10)') 'List of fresh cells at time step ', ntime
    WRITE(ifuParLog,'(A    )') '==========================================='
    WRITE(ifuParLog,'(3X,A )') '   i       j       k    Body Number'
    WRITE(ifuParLog,'(3X,A )') '------- ------- ------- -----------'
    
    DO k=1,nzc
    DO j=myJLL,myJUL
    DO i=myILL,myIUL
      
        IF (i>=1 .AND. i<=nxc .AND. &
            j>=1 .AND. j<=nyc) THEN
          withinBoundary=.TRUE.
        ELSE
          withinBoundary=.FALSE.
        END IF

        IF (fresh_cell(i,j,k) == 1) THEN
           IF (ghostCellMark(i,j,k)/=1 .AND. ghostCellSolid(i,j,k)==0 ) THEN   
              fresh_cell(i,j,k) = 0
           ELSE IF (withinBoundary) THEN
            WRITE(ifuParLog,'(4X,I5,3X,I5,3X,I5,5X,I5)') L2GI(i),L2GJ(j),k, bodyNum(i,j,k)
            num_fresh = num_fresh + 1
           ENDIF
        ENDIF
    ENDDO
    ENDDO
    ENDDO

    WRITE(ifuParLog,'(A    )') '==========================================='
    WRITE(ifuParLog,*)
    
    num_fresh = SUM(fresh_cell(1:nxc,1:nyc,1:nzc))

#   ifdef MPI
      CALL par_getSumInteger(num_fresh)
#   endif

    IF (monitorON) THEN 
      WRITE(STDOUT,'(7X,A,I8)') 'Number of Actual      Fresh Cells = ',num_fresh
      WRITE(STDOUT,*)
    ENDIF ! ntime

END SUBROUTINE identify_freshcells
!---------------------------------------------------------------------



SUBROUTINE set_iblank_general_body()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k



    DO k=1,nzc
      DO j=1,L2GJ(myJLL)-1
        DO i=1,nxc_GLBL
          READ(ifuIblankIn,*)
        ENDDO
      ENDDO

      DO j=myJLL,myJUL
        DO i=1,L2GI(myILL)-1
          READ(ifuIblankIn,*)
        ENDDO

        DO i=myILL,myIUL
          READ(ifuIblankIn,*) iblank(i,j,k)
        ENDDO

        DO i=L2GI(myIUL)+1,nxc_GLBL
          READ(ifuIblankIn,*)
        ENDDO
      ENDDO

      DO j=L2GJ(myJUL)+1,nyc_GLBL
        DO i=1,nxc_GLBL
          READ(ifuIblankIn,*)
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE set_iblank_general_body
!---------------------------------------------------------------------
