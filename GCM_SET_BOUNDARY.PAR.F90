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
!     Fady Najjar
!     Rajat Mittal
!
!  Other contributing programmers:
!     Haibo Dong
!     Haoxiang Luo
!     Meliha Bozkurttas
!     Reza Ghias 
!     Rupesh Babu K. A.
!     S. A. Mohsen Karimian
!     Xudong Zheng 
!
!  Filename: GCM_SET_BOUNDARY.PAR.F90
!  Latest Modification: Dec, 29 2010 (PAT 2.1.0)
!  by Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     identify_cell_bodyNum()
!     identify_ghostcells_solid()
!     identify_ghostcells_membrane(iBody)
!     GCM_set_internal_boundary()
!     GCM_Calc_VanMatrixDN()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1

# define G2LI(i)      i-(myIs-1)
# define G2LJ(j)      j-(myJs-1)



SUBROUTINE identify_cell_bodynum()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!   loop variables
!   --------------
    INTEGER :: iBody
    INTEGER :: i, j, k, m
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax

!   local variables
!   ---------------
    INTEGER :: iC, jC, kC



    
    DO iBody=1,nBody
      DO m = 1,nPtsBodyMarker(iBody)
        iC = 0
        jC = 0
        kC = 0

        DO i = 0, nxc_GLBL
          IF ((xc(i) <= xBodyMarker(iBody,m) .AND. xc(i+1) > xBodyMarker(iBody,m))) THEN
            iC = i
          ENDIF ! xc
        ENDDO ! i

        DO j = 0, nyc_GLBL
          IF ((yc(j) <= yBodyMarker(iBody,m) .AND. yc(j+1) > yBodyMarker(iBody,m))) THEN
            jC = j
          ENDIF ! yc
        ENDDO ! j

        IF (body_dim(iBody) == BODY_DIM2) THEN
          kC = 1 
        ELSE
          DO k = 0, nzc
            IF ((zc(k) <= zBodyMarker(iBody,m) .AND. zc(k+1) > zBodyMarker(iBody,m))) THEN
              kC = k
            ENDIF ! zc
          ENDDO ! k
        ENDIF

!       Projecting on the subdomain (SAMK)
!       ----------------------------------
        iMin=G2LI(iC  )
        iMax=G2LI(iC+1)
        jMin=G2LJ(jC  )
        jMax=G2LJ(jC+1)

!       Cycle if element is outside the subdomain (SAMK)
!       ------------------------------------------------
        IF ((iMin>myIUL) .OR. (iMax<myILL)) CYCLE
        IF ((jMin>myJUL) .OR. (jMax<myJLL)) CYCLE

        iMin=MAX(iMin,myILL)
        iMax=MIN(iMax,myIUL)
        jMin=MAX(jMin,myJLL)
        jMax=MIN(jMax,myJUL)
        kMin=MAX(kC  ,1    )
        kMax=MIN(kC+1,nzc  )

        bodyNum(iMin:iMAX,jMin:jMax,kMin:kMax) = iBody
      ENDDO
    ENDDO
 
END SUBROUTINE identify_cell_bodynum
!---------------------------------------------------------------------



SUBROUTINE identify_ghostcells_solid()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!   loop variables
!   --------------
    INTEGER :: i, j, k

!   local variables
!   ---------------
    INTEGER :: iM, iP, jM, jP, kM, kP
    
    REAL(KIND=CGREAL) :: myTime


      iblank(:,:,:) = iblank_solid(:,:,:)
      ghostCellMark(:,:,:)  = 0
      ghostCellSolid(:,:,:) = 0

!     Mark boundary cells (ie. cells inside body which have at least one neighbor in fluid)
!     -------------------------------------------------------------------------------------
      DO k = 1     , nzc
      DO j = myJmin, myJmax
      DO i = myImin, myImax

        iM = MAX(i-1,myILL)
        iP = MIN(i+1,myIUL)
        jM = MAX(j-1,myJLL)
        jP = MIN(j+1,myJUL)
        kM = MAX(k-1,1)
        kP = MIN(k+1,nzc)

        IF ( iblank(i,j,k)/=iblank(iM,j,k) .OR. &
             iblank(i,j,k)/=iblank(iP,j,k) .OR. &
             iblank(i,j,k)/=iblank(i,jM,k) .OR. &
             iblank(i,j,k)/=iblank(i,jP,k) .OR. &
             iblank(i,j,k)/=iblank(i,j,kM) .OR. &
             iblank(i,j,k)/=iblank(i,j,kP)      ) THEN

           ghostCellSolid(i,j,k) = 1
           IF (iblank(i,j,k) == 1) ghostCellMark(i,j,k) = 1
        ENDIF

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

#     ifdef MPI
!       Communicating the outmost layer of the subdomain for ghost cells.
!       -----------------------------------------------------------------
        CALL par_comm_outermost_int(ghostCellSolid,myTime)
        CALL par_comm_outermost_int(ghostCellMark ,myTime)
#     endif

END SUBROUTINE identify_ghostcells_solid 
!---------------------------------------------------------------------



SUBROUTINE identify_ghostcells_membrane(iBody)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

    INTEGER, INTENT(IN):: iBody

!   loop variables
!   --------------
    INTEGER :: i, j, k, m

!   local variables
!   ---------------
    INTEGER :: iG, jG, kG, nbdr
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iM, iP, jM, jP, kM, kP
    INTEGER :: iC, jC, kC
    INTEGER :: cElemG

    REAL(KIND=CGREAL) :: xBIT,yBIT,zBIT
    REAL(KIND=CGREAL) :: xGC,yGC,zGC
    REAL(KIND=CGREAL) :: myTime



!   Find all cells that contain a IB node
!   -------------------------------------
    ghostCellMemb = 0
    boundCell     = 0

    DO m = 1, nPtsBodyMarker(iBody)

      iC = -1
      jC = -1
      kC = -1

      DO i = 1,nxc_GLBL
        IF ( ( xc(i) <= xBodyMarker(iBody,m) .AND. xc(i+1) > xBodyMarker(iBody,m) ) ) THEN
          iC = i
        ENDIF ! xc
      ENDDO ! i

      DO j = 1,nyc_GLBL
        IF ( ( yc(j) <= yBodyMarker(iBody,m) .AND. yc(j+1) > yBodyMarker(iBody,m) ) ) THEN
          jC = j
        ENDIF ! yc
      ENDDO ! j

      IF (body_dim(iBody) == BODY_DIM2) THEN
        kC = 1 
      ELSE 
        DO k = 1,nzc
          IF ( ( zc(k) <= zBodyMarker(iBody,m) .AND. zc(k+1) > zBodyMarker(iBody,m) ) ) THEN
            kC = k
          ENDIF ! zc
        ENDDO ! k
      ENDIF

      IF ( iC == -1 .OR. jC == -1 .OR. kC == -1 ) THEN
        PRINT*,'IC,JC,KC',IC,JC,KC
        PRINT*,'iBody,m',iBody,m
        PRINT*,'xBodyMarker,yBodyMarker,zBodyMarker',xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
        CALL flow_stop
        STOP
      ENDIF
          
!     Projecting on the subdomain (SAMK)
!     ----------------------------------
      iMin=G2LI(iC  )
      iMax=G2LI(iC+1)
      jMin=G2LJ(jC  )
      jMax=G2LJ(jC+1)

!     Cycle if element is outside the subdomain (SAMK)
!     ------------------------------------------------
      IF ((iMin>myIUL) .OR. (iMax<myILL)) CYCLE
      IF ((jMin>myJUL) .OR. (jMax<myJLL)) CYCLE

      iMin=MAX(iMin,myILL)
      iMax=MIN(iMax,myIUL)
      jMin=MAX(jMin,myJLL)
      jMax=MIN(jMax,myJUL)
      kMin=MAX(kC  ,1    )
      kMax=MIN(kC+1,nzc  )
      
      boundCell(iMin:iMax,jMin:jMax,kMin:kMax) = 1

    ENDDO

!   Mark boundary cells (ie. cells inside body which have at least one neighbor in fluid)
!   -------------------------------------------------------------------------------------
    DO k = 1, nzc
    DO j = myJmin, myJmax
    DO i = myImin, myImax

      iM = MAX(i-1,myILL)
      iP = MIN(i+1,myIUL)
      jM = MAX(j-1,myJLL)
      jP = MIN(j+1,myJUL)
      kM = MAX(k-1,1  )
      kP = MIN(k+1,nzc)

      IF (boundCell(i,j,k) > 0) THEN

        IF ( ( iblank(i,j,k)+iblank(iM,j,k) == 1 .AND. boundCell(iM,j,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(iP,j,k) == 1 .AND. boundCell(iP,j,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,jM,k) == 1 .AND. boundCell(i,jM,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,jP,k) == 1 .AND. boundCell(i,jP,k) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,j,kM) == 1 .AND. boundCell(i,j,kM) == 1 )  .OR.  &
             ( iblank(i,j,k)+iblank(i,j,kP) == 1 .AND. boundCell(i,j,kP) == 1 )      )  THEN
          ghostCellMemb(i,j,k) = 1
          bodyNum(i,j,k)= iBody
        ENDIF

      ENDIF

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

     DO m = 1,nPtsBodyMarker(iBody)

       IF (edge_node(iBody,m) == 1) THEN

         iC = -1
         jC = -1
         kC = -1

         DO i = 1,nxc_GLBL
           IF ( ( xc(i) <= xBodyMarker(iBody,m)  .AND. xc(i+1) > xBodyMarker(iBody,m) ) ) THEN
             iC = i
           ENDIF ! xc
         ENDDO ! i
         DO j = 1,nyc_GLBL
           IF ( ( yc(j) <= yBodyMarker(iBody,m)  .AND. yc(j+1) > yBodyMarker(iBody,m) ) ) THEN
             jC = j
           ENDIF ! yc
         ENDDO ! j
         IF (body_dim(iBody) == BODY_DIM2) THEN
            kC = 1
         ELSE
           DO k = 1,nzc
             IF ( ( zc(k) <= zBodyMarker(iBody,m)  .AND. zc(k+1) > zBodyMarker(iBody,m) ) ) THEN
               kC = k
             ENDIF ! zc
           ENDDO ! k
         ENDIF

         IF ( iC == -1 .OR. jC == -1 .OR. kC == -1 ) THEN
           write(*,*) xc(0), yc(0), zc(0)
           write(*,*)'IC,JC,KC',IC,JC,KC
           write(*,*)
           write(*,*)'iBody,m',iBody,m
           write(*,*)
           write(*,*)'xBMark,yBMark,zBMark',xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
!           CALL abort_vicar3d(30)
           CALL flow_stop
           STOP
         ENDIF

      iMin=G2LI(iC-1)
      iMax=G2LI(iC+2)
      jMin=G2LJ(jC-1)
      jMax=G2LJ(jC+2)

      IF ((iMin>myIUL) .OR. (iMax<myILL)) CYCLE
      IF ((jMin>myJUL) .OR. (jMax<myJLL)) CYCLE

      iMin=MAX(iMin,myILL)
      iMax=MIN(iMax,myIUL)
      jMin=MAX(jMin,myJLL)
      jMax=MIN(jMax,myJUL)
      kMin=MAX(kC-1  ,1    )
      kMax=MIN(kC+2,nzc  )

      boundCell(iMin:iMax,jMin:jMax,kMin:kMax) = 2

       ENDIF

     ENDDO

#   ifdef MPI
!     Communicating the outmost layer of the subdomain for ghost cells.
!     -----------------------------------------------------------------
      CALL par_comm_outermost_int(ghostCellMemb,myTime)
      CALL par_comm_outermost_int(bodyNum,myTime)
#   endif
   
!   Determine which ghost cells are real and calculate body intercept and closest element for these ghost cells.
!   ------------------------------------------------------------------------------------------------------------
!    IF (boundary_formulation == GCM_METHOD) THEN
      nbdr  = 0

      DO k = 1    , nzc
      DO j = myJLL, myJUL
      DO i = myILL, myIUL

        IF ( ghostCellMemb(i,j,k) == 1 .AND. boundCell(i,j,k) == 2 )   THEN
          iG    = i
          jG    = j
          kG    = k
          xGC   = xc(L2GI(iG))
          yGC   = yc(L2GJ(jG))
          zGC   = zc(kG)

          nbdr = nbdr + 1
     
          CALL GCM_Calc_BodyIntercept_Unstruc( iG, jG, kG, xGC, yGC, zGC, xBIT, yBIT, zBIT, cElemG )

          IF ( cElemG == -1 )  THEN    ! did not find an intercept
            ghostCellMemb(i,j,k) = 0
            bodyNum(i,j,k) = 0
          ENDIF
        ENDIF ! ghostCellMemb

      ENDDO ! i 
      ENDDO ! j
      ENDDO ! k

      IF (monitorON) WRITE(STDOUT,'(7X,A,I5)') 'GCM_set_internal_boundary: nGhost_MEMB = ',nbdr

!    ENDIF ! boundary_formulation

    DO k=1,nzc
    DO j=myJLL,myJUL
    DO i=myILL,myIUL
      ghostCellMark(i,j,k) = ghostCellMemb(i,j,k) + ghostCellMark(i,j,k)
      IF ( ghostCellMark(i,j,k) > 1 ) THEN
        PRINT*,'ghostCellMark is > 1;   aborting.'
        PRINT*, i,j,k,myRank
        CALL flow_stop
        STOP
      ENDIF 
      
      IF ( ghostCellMark(i,j,k) /= 1 .AND. bodyNum(i,j,k)==iBody) THEN
        write(*,*) 'Error in bodyNum at', i,j,k,iBody
        stop
      END IF
    ENDDO
    ENDDO
    ENDDO

#ifdef MPI
!     Communicating the outmost layer of the subdomain for ghost cells.
!     -----------------------------------------------------------------
      CALL par_comm_outermost_int(bodyNum,myTime)
      CALL par_comm_outermost_int(ghostCellMemb,myTime)
      CALL par_comm_outermost_int(ghostCellMark ,myTime)
#endif

END SUBROUTINE identify_ghostcells_membrane
!---------------------------------------------------------------------



SUBROUTINE GCM_set_internal_boundary()

! ---------------------------------------------------------------------------------
!  Given a set of marker points, this subroutine computes
!   1) normal intercepts and associated geometrical info. from ghost nodes to body
!   2) Image point location
!   3) Weights in stencil for computing values at the image points
! ---------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays
	USE implicit_coupling_parameters !   Added by Rajneesh  
    
    IMPLICIT NONE

!   loop variables
!   --------------
    INTEGER :: i, j, k, iBody

!   local variables
!   ---------------
    INTEGER :: ibdr
    INTEGER :: iG, jG, kG, iCIndx, jCIndx, kCIndx
    INTEGER :: iGG, jGG, kGG
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iRange, jRange, kRange
    INTEGER :: node1, node2, node3

    REAL(KIND=CGREAL) :: xGC, yGC, zGC, xBI, yBI, zBI, xIP, yIP, zIP, xBIN, yBIN, zBIN
    REAL(KIND=CGREAL) :: dsIntercept

    REAL(KIND=CGREAL) :: coeffGCMDirc(iRowMax), coeffGCMNeum(iRowMax)



    ibdr = 0
    DO k = 1, nzc
    DO j = 1, nyc
    DO i = 1, nxc
      IF ( ghostCellMark(i,j,k) == 1 )   THEN
        ibdr = ibdr + 1
      ENDIF ! ghostCellMark    
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k
      
    nGhost = ibdr         ! total number of ghost points

    WRITE(ifuParLog,'(A,I8)') 'Number of Ghost Cells = ',nGhost 

!   Deallocate arrays pertinent to Ghost Cells 
!   ------------------------------------------

!    IF(( nread /= 0.and.ntime > ntime_start).OR.(nread == 0 .and. ntime > ntime_start+1)) THEN
	
!   Added by Rajneesh
    IF(( nread /= 0.and.ntime > ntime_start).OR.(nread == 0 .and. ntime > ntime_start+1)&    
	.OR.(implicit_coupling_flag == 1 .and. ntime>ntime_start .and. kimplicit.gt.1) 	)THEN	 	
      CALL GCM_DeallocateGhostCellArrays !........................COMPLETE(SAMK)
    ENDIF ! ntime
   

!   Allocate Arrays pertinent to Ghost Cells
!   ----------------------------------------
    CALL GCM_AllocateGhostCellArrays()!..........................COMPLETE(SAMK)

!   Set appropriate values for iGhost and jGhost by doing search
!   ------------------------------------------------------------
    ibdr = 0
    DO k = 1, nzc
    DO j = 1, nyc
    DO i = 1, nxc
      IF ( ghostCellMark(i,j,k) == 1 )   THEN
        ibdr = ibdr + 1
        iGhost(ibdr) = i
        jGhost(ibdr) = j
        kGhost(ibdr) = k
      ENDIF ! ghostCellMark
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

!   Find marker point closest to each boundary node 
!   and compute normal at closest marker point
!   -----------------------------------------------
    DO ibdr = 1, nGhost
      iCellIndex(ibdr) = -100
      jCellIndex(ibdr) = -100
      kCellIndex(ibdr) = -100

      iG = iGhost(ibdr)
      jG = jGhost(ibdr)
      kG = kGhost(ibdr)
      
      iGG = L2GI(iG)
      jGG = L2GJ(jG)
      kGG = kG
        
      xGC = xc(iGG)
      yGC = yc(jGG)
      zGC = zc(kGG)

      iBody = bodyNum(iG,jG,kG)

      CALL GCM_Calc_BodyIntercept_Unstruc( iG, jG, kG, xGC, yGC, zGC, xBI, yBI, zBI, closestElementGC(ibdr) )

!     Extract coordinates of Body Intercept
!     -------------------------------------
      xBodyIntercept(ibdr) = xBI
      yBodyIntercept(ibdr) = yBI
      zBodyIntercept(ibdr) = zBI

!     Get length of normal
!     --------------------
      dsIntercept = SQRT( (xGC-xBI)**2 + (yGC-yBI)**2 + (zGC-zBI)**2 )

!     check intercept length against cell size
!     if intercept is longer than cell diagonal then potential sign of  oblique intecept.
!     -----------------------------------------------------------------------------------
      IF ( dsIntercept > SQRT(dxc(iGG)**2 + dyc(jGG)**2 + dzc(kGG)**2) ) THEN
        PRINT*,' Normal intercept might not be correct'
        PRINT*,ibdr
        PRINT*,iG,jG,kG
        PRINT*,dsIntercept,SQRT(dxc(iGG)**2 + dyc(jGG)**2 + dzc(kGG)**2)
        PRINT*,dsIntercept/SQRT(dxc(iGG)**2 + dyc(jGG)**2 + dzc(kGG)**2)
        PRINT*,'Check fort.198 for more info.'

        node1   = triElemNeig(iBody,1,closestElementGC(ibdr))
        node2   = triElemNeig(iBody,2,closestElementGC(ibdr))
        node3   = triElemNeig(iBody,3,closestElementGC(ibdr))

        IF (dsIntercept > 2.0_CGREAL*SQRT(dxc(iGG)**2 + dyc(jGG)**2 + dzc(kGG)**2)) THEN
          PRINT*,'Intercept is too long'
          CALL write_dump()
          CALL flow_Stop
          STOP
        ENDIF

      ENDIF

         
!     Now compute location of probe-tip (Image Point)
!     Equation of 3D line  (parametric form)
!     SAMK :
!     Xip = Xgc + probelength . (Xbi-Xgc)
!     Capital letters indicate VECTORS
!     -----------------------------------------------
      xImagePoint(ibdr) = xGC + (xBI-xGC)*probeLengthNormalized
      yImagePoint(ibdr) = yGC + (yBI-yGC)*probeLengthNormalized
      zImagePoint(ibdr) = zGC + (zBI-zGC)*probeLengthNormalized

      probeLength(ibdr) = dsIntercept*probeLengthNormalized

      xIP = xImagePoint(ibdr)
      yIP = yImagePoint(ibdr)
      zIP = zImagePoint(ibdr)

!     Search for the lower left cell with active coefficient for the ghost cell of interest
!     -------------------------------------------------------------------------------------

!     Base range on probeLength
!     -------------------------
      iRange  = NINT(probeLength(ibdr)/dx(iGG)) +1
      jRange  = NINT(probeLength(ibdr)/dy(jGG)) +1
      kRange  = NINT(probeLength(ibdr)/dz(kGG)) +1
      
      iMin = iGhost(ibdr)-iRange
      iMax = iGhost(ibdr)+iRange
   
      iMin = MAX(iMin, 1-Ngl)
      IF (myCoords(1)==0) THEN
        iMin = MAX(iMin,0)    ! note image point is allowed to be between xc(0) and x(1)
      ELSE IF (iMin<1-Ngl) THEN
        WRITE(STDOUT,*) 'Probe length exceeds the parallel ghost layer thickness: iMin = ',iMin
        CALL flow_stop
        STOP
      END IF
!debug temp
      iMax = MIN(iMax, nxc+Ngl)
      IF (myCoords(1)==Np(1)-1) THEN
        iMax = MIN(iMax,nxc+1)   ! note image point is allowed to be between x(nx) and xc(nxc+1)
      ELSE IF (iMax>nxc+Ngl) THEN
        WRITE(STDOUT,*) 'Probe length exceeds the parallel ghost layer thickness: iMax = ',iMax, 'nxc+Ngl = ', nxc+Ngl
        WRITE(STDOUT,*)
        CALL flow_stop
        STOP
      END IF

      DO i = iMin,iMax

        IF ( ( xc(L2GI(i)) <= xIP .AND. xc(L2GI(i+1)) > xIP ) ) THEN
          iCellIndex(ibdr) = i
!          EXIT
        ENDIF ! xc

      ENDDO ! i

      jMin = jGhost(ibdr)-jRange
      jMax = jGhost(ibdr)+jRange

      jMin = MAX(jMin, 1-Ngl)
      IF (myCoords(2)==0) THEN
        jMin = MAX(jMin,0)
      ELSE IF (jMin<1-Ngl) THEN
        WRITE(STDOUT,*) 'Probe length exceeds the parallel ghost layer thickness: jMin = ',jMin
        CALL flow_stop
        STOP
      END IF

      Jmax = MIn(Jmax, nyc+Ngl)
      IF (myCoords(2)==Np(2)-1) THEN
        jMax = MIN(jMax,nyc+1)
      ELSE IF (jMax>nyc+Ngl) THEN
        WRITE(STDOUT,*) 'Probe length exceeds the parallel ghost layer thickness: jMax = ',jMax, 'nyc+Ngl = ', nyc+Ngl
        CALL flow_stop
        STOP
      END IF

      DO j = jMin,jMax

        IF ( ( yc(L2GJ(j)) <= yIP .AND. yc(L2GJ(j+1)) > yIP ) ) THEN
          jCellIndex(ibdr) = j
          EXIT
        ENDIF ! xc

      ENDDO ! j

      kMin = kGhost(ibdr)-kRange
      kMax = kGhost(ibdr)+kRange

      kMin = MAX(kMin,0)
      kMax = MIN(kMax,nzc+1)

      DO k = kMin,kMax

        IF ( ( zc(k) <= zIP .AND. zc(k+1) > zIP ) ) THEN
          kCellIndex(ibdr) = k
          EXIT
        ENDIF ! xc

      ENDDO ! k

      IF ( iCellIndex(ibdr) ==-100 .OR. jCellIndex(ibdr) ==-100 .OR. kCellIndex(ibdr) ==-100 ) THEN
        PRINT*,'Failed to Find four nodes surrounding an image point'
        PRINT*,ibdr
        PRINT*,iG, jG, kG
        PRINT*,xgc,ygc,zgc
        PRINT*,xBI,yBI,zBI
        PRINT*,xIP,yIP,zIP
        PRINT*,'Aborting Run'
        CALL flow_stop
        STOP
      ENDIF

!     Perform bilinear interpolation
!     ------------------------------
      iCIndx = iCellIndex(ibdr)
      jCIndx = jCellIndex(ibdr)
      kCIndx = kCellIndex(ibdr)
 
      xBIN = triElemNormx(iBody,closestElementGC(ibdr))
      yBIN = triElemNormy(iBody,closestElementGC(ibdr))
      zBIN = triElemNormz(iBody,closestElementGC(ibdr))
	  
	  xBodyInterceptNorm(ibdr) = xBIN
	  yBodyInterceptNorm(ibdr) = yBIN
	  zBodyInterceptNorm(ibdr) = zBIN
	  
      CALL GCM_Calc_vanMatrixDN( iG, jG, kG, iCIndx, jCIndx, kCIndx,             &
                                 xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                 coeffGCMDirc, coeffGCMNeum      )
 
      coeffGCMD(1:iRowMax,ibdr) = coeffGCMDirc(1:iRowMax)
      coeffGCMN(1:iRowMax,ibdr) = coeffGCMNeum(1:iRowMax)

    ENDDO ! ibdr

END SUBROUTINE GCM_set_internal_boundary
!---------------------------------------------------------------------



SUBROUTINE GCM_Calc_VanMatrixDN( iG, jG, kG, iCIndex,jCIndex, kCIndex,      &
                                 xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                 coeffGCMDirc, coeffGCMNeum)
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE gcm_arrays

    IMPLICIT NONE

!   parameters variables
!   --------------------
    INTEGER, INTENT(IN)            :: iG, jG, kG, iCIndex,jCIndex, kCIndex
    REAL(KIND=CGREAL), INTENT(IN ) :: xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN
    REAL(KIND=CGREAL), INTENT(OUT) :: coeffGCMDirc(iRowMax), coeffGCMNeum(iRowMax)

!   loop variables
!   --------------
    INTEGER :: i, j, k, iRow
    INTEGER :: info

!   local variables
!   ---------------
    REAL(KIND=CGREAL) :: xC1, xC2, xC3, xN1, xN2, xN3

!*****************************************************************************************
  
!   |-------|-------|---/---|-------|--         N : Nth ghost point
!   |   ii  |  iii  |  *    |       |           * : markers
!   |   0...|...O   | / .   |   .   |           O : other nodes used in bilinear interpolation
!   |   .   |   .   |*      |       |           + : probe tip (Image Point) 
!   |---.---|--+.---/-------|-------|--
!   |   .   |   .  *|       |       |
!   |   0...|. .O / |   N   |   .   |
!   |   i   |  iv*  |       |       |
!   |-------| --/ --|-------|-------|--

! interpolant      U = a X X X  + b X X  + c X X  + d X X
!                         1 2 3      1 2      1 3      2 3
!
!                    + e X  + f X + g X  + h
!                         1      2     3
!
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!     iv  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  iv ]
!
!
!   Van Matrix For Dirichlet conditions at Intersection Point (N)
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!      N  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  N  ]
!
!   Van Matrix For Neumann conditions at Intersection point (N)
!    B1 = n_x, B2 = n_y (components of normal vectors)
!    F_m = value of normal derivative 
!
!         [  X X           X     X   1  ]  [   ]     [     ] 
!      i  [   1 2           1     2     ]  [ a ]     [ U   ]
!         [                             ]  [   ]     [  i  ]
!         [  X X           X     X   1  ]  [   ]     [     ]
!      ii [   1 2           1     2     ]  [ b ]     [ U   ]
!         [                             ]  [   ]  =  [  ii ]
!     iii [  X X           X     X   1  ]  [   ]     [     ]
!         [   1 2           1     2     ]  [ c ]     [ U   ]
!         [                             ]  [   ]     [  iii]
!      N  [  B X  + B X    B     B   0  ]  [   ]     [     ]
!         [   1 2    2  1   1     2     ]  [ d ]     [ F   ]
!         [                             ]  [   ]     [  N  ]
!


    DO iRow= 1, iRowMax
      i  = iCIndex + incI(iRow)
      j  = jCIndex + incJ(iRow)
      k  = kCIndex + incK(iRow)

!     Construct Vandermonde Matrices
!     ------------------------------

!     Check if Ghost node is part of cell formation
!     ---------------------------------------------
      IF ( i/=iG .OR. j /= jG  .OR. k/= kG) THEN

!       Normal interpolation
!       --------------------
        xC1 = xc(L2GI(i))
        xC2 = yc(L2GJ(j))
        xC3 = zc(k)

!       Dirichlet conditions for velocity field
!       ---------------------------------------
        vanMatrixD(iRow,1) = xC1*xC2*xC3
        vanMatrixD(iRow,2) = xC1*xC2
        vanMatrixD(iRow,3) = xC1*xC3
        vanMatrixD(iRow,4) = xC2*xC3
        vanMatrixD(iRow,5) = xC1
        vanMatrixD(iRow,6) = xC2
        vanMatrixD(iRow,7) = xC3
        vanMatrixD(iRow,8) = 1.0_CGREAL

!       Neumann conditions for pressure field
!       -------------------------------------
        vanMatrixN(iRow,1) = xC1*xC2*xC3
        vanMatrixN(iRow,2) = xC1*xC2
        vanMatrixN(iRow,3) = xC1*xC3
        vanMatrixN(iRow,4) = xC2*xC3
        vanMatrixN(iRow,5) = xC1
        vanMatrixN(iRow,6) = xC2
        vanMatrixN(iRow,7) = xC3
        vanMatrixN(iRow,8) = 1.0_CGREAL

      ELSE

!       Correct For Ghost node part of cell formation, switch to Body Intercept point
!       -----------------------------------------------------------------------------
        xC1 = xBI
        xC2 = yBI
        xC3 = zBI
        xN1 = xBIN
        xN2 = yBIN
        xN3 = zBIN

        vanMatrixD(iRow,1) = xC1*xC2*xC3
        vanMatrixD(iRow,2) = xC1*xC2
        vanMatrixD(iRow,3) = xC1*xC3
        vanMatrixD(iRow,4) = xC2*xC3
        vanMatrixD(iRow,5) = xC1
        vanMatrixD(iRow,6) = xC2
        vanMatrixD(iRow,7) = xC3
        vanMatrixD(iRow,8) = 1.0_CGREAL

        vanMatrixN(iRow,1) = xN1*xC2*XC3 + xN2*xC1*XC3 + xN3*XC1*XC2
        vanMatrixN(iRow,2) = xN1*xC2 + xN2*xC1
        vanMatrixN(iRow,3) = xN1*xC3 + xN3*xC1
        vanMatrixN(iRow,4) = xN2*xC3 + xN3*xC2
        vanMatrixN(iRow,5) = xN1
        vanMatrixN(iRow,6) = xN2
        vanMatrixN(iRow,7) = xN3
        vanMatrixN(iRow,8) = 0.0_CGREAL
      ENDIF ! i
    ENDDO ! iRow

!   Compute inverse of Vandermonde Matrices
!   ---------------------------------------
    CALL DGETRF(8, 8, vanMatrixD,8,iPvt, info) 
    CALL DGETRI(8, vanMatrixD,8,iPvt,work, 8, info) 

    CALL DGETRF(8, 8, vanMatrixN,8,iPvt, info)
    CALL DGETRI(8, vanMatrixN,8,iPvt,work, 8, info)

!   Load Coeff-Matrices
!   -------------------
    DO iRow = 1, iRowMax
      coeffGCMDirc(iRow) = vanMatrixD(1,iRow)*xIP*yIP*zIP  &
                         + vanMatrixD(2,iRow)*xIP*yIP      &
                         + vanMatrixD(3,iRow)*xIP*zIP      &
                         + vanMatrixD(4,iRow)*yIP*zIP      &
                         + vanMatrixD(5,iRow)*xIP          &
                         + vanMatrixD(6,iRow)*yIP          &
                         + vanMatrixD(7,iRow)*zIP          &
                         + vanMatrixD(8,iRow)

      coeffGCMNeum(iRow) = vanMatrixN(1,iRow)*xIP*yIP*zIP  &
                         + vanMatrixN(2,iRow)*xIP*yIP      &
                         + vanMatrixN(3,iRow)*xIP*zIP      &
                         + vanMatrixN(4,iRow)*yIP*zIP      &
                         + vanMatrixN(5,iRow)*xIP          &
                         + vanMatrixN(6,iRow)*yIP          &
                         + vanMatrixN(7,iRow)*zIP          &
                         + vanMatrixN(8,iRow)
    ENDDO ! iRow 

END SUBROUTINE GCM_Calc_VanMatrixDN
!---------------------------------------------------------------------

SUBROUTINE identify_gc_membrane_final(iBody)

! Subroutine removes final straggler ghost cells based on iup-ium

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN):: iBody

!... loop variables

    INTEGER :: i,j,k

!... local variables
    REAL(KIND=CGREAL):: myTime

    SELECT CASE( body_dim(iBody) )

      CASE(BODY_DIM2)
        DO j = myJmin, myJmax
        DO i = myImin, myImax
          k = 1
          IF (ghostCellMemb(i,j,k) == 1) THEN

            IF (   ium(i+1,j,k) == 1 .OR. iup(i-1,j,k)== 1  .OR. &
                   jum(i,j+1,k) == 1 .OR. jup(i,j-1,k)== 1  ) THEN
              DO k = 1,nz-1
                ghostCellMemb(i,j,k) = 1
              ENDDO
            ELSE
              DO k = 1,nzc
                ghostCellMemb(i,j,k) = 0
                ghostCellMark(i,j,k) = 0
                bodyNum(i,j,k)       = 0
              ENDDO
            ENDIF
          ENDIF

        ENDDO ! i
        ENDDO ! j

      CASE(BODY_DIM3)
        DO k = 1, nzc
        DO j = myJmin, myJmax
        DO i = myImin, myImax

          IF (ghostCellMemb(i,j,k) == 1) THEN

            IF (   ium(i+1,j,k) == 1 .OR. iup(i-1,j,k)== 1  .OR. &
                   jum(i,j+1,k) == 1 .OR. jup(i,j-1,k)== 1  .OR. &
                   kum(i,j,k+1) == 1 .OR. kup(i,j,k-1)== 1  ) THEN
              ghostCellMemb(i,j,k) = 1
            ELSE
              !PRINT*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
              !PRINT*,'Removing GC_FINAL AT ',i,j,k
              ghostCellMemb(i,j,k) = 0
              ghostCellMark(i,j,k) = 0
              bodyNum(i,j,k)       = 0
            ENDIF
          ENDIF

        ENDDO ! i
        ENDDO ! j
        ENDDO ! k

      END SELECT ! BODY_DIM

#ifdef MPI
!     Communicating the outmost layer of the subdomain for ghost cells.
!     -----------------------------------------------------------------
      CALL par_comm_outermost_int(bodyNum,myTime)
      CALL par_comm_outermost_int(ghostCellMemb,myTime)
      CALL par_comm_outermost_int(ghostCellMark,myTime)
#endif

END SUBROUTINE identify_gc_membrane_final
!---------------------------------------------------------------------

!==================================================================================
   SUBROUTINE set_internal_iup_GCM_memb()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
   

    IMPLICIT NONE

    INTEGER           :: i,j,k

    DO k=1,nzc
    DO j=myJmin,myJmax
    DO i=myImin,myImax
      IF (ghostCellMark(i,j,k) == 1 ) THEN
        IF (iblank(i+1,j,k) /= iblank(i,j,k) .AND. ghostCellMark(i+1,j,k) == 1) iup(i,j,k)=1.0_CGREAL
        IF (iblank(i-1,j,k) /= iblank(i,j,k) .AND. ghostCellMark(i-1,j,k) == 1) ium(i,j,k)=1.0_CGREAL
        IF (iblank(i,j+1,k) /= iblank(i,j,k) .AND. ghostCellMark(i,j+1,k) == 1) jup(i,j,k)=1.0_CGREAL
        IF (iblank(i,j-1,k) /= iblank(i,j,k) .AND. ghostCellMark(i,j-1,k) == 1) jum(i,j,k)=1.0_CGREAL
        IF (iblank(i,j,k+1) /= iblank(i,j,k) .AND. ghostCellMark(i,j,k+1) == 1) kup(i,j,k)=1.0_CGREAL
        IF (iblank(i,j,k-1) /= iblank(i,j,k) .AND. ghostCellMark(i,j,k-1) == 1) kum(i,j,k)=1.0_CGREAL       
      ENDIF
    ENDDO
    ENDDO
    ENDDO

   END SUBROUTINE set_internal_iup_GCM_memb

!------------------------------------------------------------

