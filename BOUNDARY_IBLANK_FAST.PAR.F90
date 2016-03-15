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
!  Filename: BOUNDARY_IBLANK_FAST.PAR.F90
!  Latest Modification: Jan 20, 2010 (ver. P1.5.5)
!  Made by X. Zheng
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     set_iblank_canonical_body_fast()
!     search_vertex_dotNorm(iBody,iCell,jCell,kCell,closestElement,dotNorm)
!     check_BIInsideTriangle()
!     calc_BIOutsideTriangle()
!     calc_crossProduct(r,s,cross_product)
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1

# define G2LI(i)      i-(myIs-1)
# define G2LJ(j)      j-(myJs-1)



SUBROUTINE set_iblank_canonical_body_fast(nBody_Begin,nBody_End)

! -----------------------------------------------------------------------------------------
!  In this subroutine iBlank is determined for each LOCAL cell including the ghost layers.
!  Please note that some loops are performed globally and some locally. (SAMK)
! -----------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

# ifdef MPI
    use mpi
# endif

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nbody_Begin, nbody_End

    INTEGER, PARAMETER :: UNDECIDED = 100000, DECIDED = -UNDECIDED

    !INTEGER, DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1) :: iblankUndecided, iblankTemp

    INTEGER :: i, j, k, iBody, iElem
    INTEGER :: iCellMinG, iCellMaxG, jCellMinG, jCellMaxG
    INTEGER :: iCellMin,  iCellMax,  jCellMin,  jCellMax, kCellMin, kCellMax
    INTEGER :: sumBlank
    INTEGER :: iclock1,iclock2,iclock3,iclock4,iclock5,iclock6,clock_rate
    INTEGER :: iErr
    INTEGER :: CoronaSize
    INTEGER :: mVert1,mVert2,mVert3
    INTEGER :: kMin,kMax
    INTEGER :: cElement
    INTEGER :: nxct
    
    LOGICAL :: DecidedYet, InExtension

    REAL(KIND=CGREAL) :: xBoundMax,xBoundMin
    REAL(KIND=CGREAL) :: yBoundMax,yBoundMin
    REAL(KIND=CGREAL) :: zBoundMax,zBoundMin
    REAL(KIND=CGREAL) :: dotNorm
    REAL(KIND=CGREAL) :: myTime

    REAL(KIND=CGREAL), DIMENSION(3) :: xVert,yVert,zVert

!   Allocate local memory
    ALLOCATE(iblankUndecided(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for iblankUndecided'
      CALL flow_stop
      STOP
    ENDIF ! iErr

    ALLOCATE(iblankTemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for iblankTemp'
      CALL flow_stop
      STOP
    ENDIF ! iErr

!   Initialize values
!   -----------------
    InExtension=.FALSE.
    nxct=0

!   Checking the outflow extension zone (SAMK)
!   ------------------------------------------
    DO i = 1, nxc_GLBL
      IF (xc(i)>=Xext .AND. nxct==0) nxct=G2LI(i)
    ENDDO
    IF (nxct<=MyILL) THEN
      InExtension=.TRUE.
      IF (ntime == ntime_start) THEN
        WRITE(ifuParLog,'(A)') 'This sub-domain completely falls into the outflow extension zone!'
        WRITE(ifuParLog,*)
      END IF
    END IF
    IF (monitorON) THEN
      WRITE(STDOUT,'(7X,A,I5,A,I5)') 'Body number(s) ', nbody_Begin,' to ',nbody_End
      WRITE(STDOUT,'(7X,A)') '-----------------------------'
    END IF

!   Initialize iblank to undecided value over all bodies in domain
!   --------------------------------------------------------------
    SELECT CASE (nDim)
    CASE (DIM_2D)
      iblankUndecided(myILL:myIUL,myJLL:myJUL,1) = UNDECIDED

    CASE (DIM_3D)
      iblankUndecided(myILL:myIUL,myJLL:myJUL,1:nzc) = UNDECIDED
    END SELECT ! nDim

!   Loop over all bodies in domain and body surface elements
!   --------------------------------------------------------
    CALL system_clock(iclock1)

    DO iBody=nbody_Begin,nbody_End

!     Define size of corona of cells around body for which iblank will be determined from first principles
!     ----------------------------------------------------------------------------------------------------
      IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
        coronaSize = 1  ! SAMK: 3 layers of corona may not be necessary for membranes
      ELSE 
        coronaSize = 1
      ENDIF

      DO iElem=1,totNumTriElem(iBody)

!       Extract cell indices and vertices of triangular element
!       -------------------------------------------------------
        mVert1    = triElemNeig(iBody,1,iElem)
        mVert2    = triElemNeig(iBody,2,iElem)
        mVert3    = triElemNeig(iBody,3,iElem)

        xVert(1)  = xBodyMarker(iBody,mVert1)
        xVert(2)  = xBodyMarker(iBody,mVert2)
        xVert(3)  = xBodyMarker(iBody,mVert3)

        yVert(1)  = yBodyMarker(iBody,mVert1)
        yVert(2)  = yBodyMarker(iBody,mVert2)
        yVert(3)  = yBodyMarker(iBody,mVert3)

        zVert(1)  = zBodyMarker(iBody,mVert1)
        zVert(2)  = zBodyMarker(iBody,mVert2)
        zVert(3)  = zBodyMarker(iBody,mVert3)

!       cycle if element is outside computational domain
!       ------------------------------------------------
        IF ((MAXVAL(xVert(1:3)) <= x(1      )) .OR. &
            (MINVAL(xVert(1:3)) >= x(nx_GLBL)) .OR. &
            (MAXVAL(yVert(1:3)) <= y(1      )) .OR. &
            (MINVAL(yVert(1:3)) >= y(ny_GLBL)) .OR. &
            (MAXVAL(zVert(1:3)) <= z(1      )) .OR. &
            (MINVAL(zVert(1:3)) >= z(nz     ))) CYCLE

!       -------------------------------------------------------------------------------
!        Find all cells within the bounding box of the vertices for each element
!
!        SAMK: Note that even if the element is not in the subdomain, since coronaSize
!              can be 3, the searching layer might fall into the subdomain. Therefore,
!              the search must be performed globally and then the results must be
!              limited to the subdomains.
!       -------------------------------------------------------------------------------
        xBoundMin = MINVAL(xVert(1:3))
        xBoundMax = MAXVAL(xVert(1:3))

        yBoundMin = MINVAL(yVert(1:3))
        yBoundMax = MAXVAL(yVert(1:3))

        zBoundMin = MINVAL(zVert(1:3))
        zBoundMax = MAXVAL(zVert(1:3))

        iCellMinG = 0
        iCellMaxG = 0
        jCellMinG = 0
        jCellMaxG = 0
        kCellMin  = 0
        kCellMax  = 0

!       i-direction 
!       -----------
        DO i = 0, nxc_GLBL
          IF ( (xc(i)-xBoundMin)*(xc(i+1)-xBoundMin) < 0.0_CGREAL ) &
            iCellMinG = (i  )-(coronaSize-1)
          IF ( (xc(i)-xBoundMax)*(xc(i+1)-xBoundMax) < 0.0_CGREAL ) &
            iCellMaxG = (i+1)+(coronaSize-1)
        ENDDO
        IF (iCellMinG <= 0 )                             iCellMinG = 1
        IF (iCellMaxG == 0 .OR. iCellMaxG >= nxc_GLBL+1) iCellMaxG = nxc_GLBL

!       Projecting on the subdomain (SAMK)
!       ----------------------------------
        iCellMin=G2LI(iCellMinG)
        iCellMax=G2LI(iCellMaxG)

!       Cycle if element is outside the subdomain (SAMK)
!       ------------------------------------------------
        IF (iCellMax<myILL .OR. iCellMin>myIUL) CYCLE

!       Limit the search to within the subdomain
!       ----------------------------------------
        iCellMin=MAX(iCellMin,myILL)
        iCellMax=MIN(iCellMax,myIUL)

!       j-direction 
!       -----------
        DO j = 0, nyc_GLBL
          IF ( (yc(j)-yBoundMin)*(yc(j+1)-yBoundMin) < 0.0_CGREAL ) &
            jCellMinG = (j  )-(coronaSize-1)
          IF ( (yc(j)-yBoundMax)*(yc(j+1)-yBoundMax) < 0.0_CGREAL ) &
            jCellMaxG = (j+1)+(coronaSize-1)
        ENDDO
        IF (jCellMinG <= 0 )                             jCellMinG = 1
        IF (jCellMaxG <= 0 .OR. jCellMaxG >= nyc_GLBL+1) jCellMaxG = nyc_GLBL

!       Projecting on the subdomain
!       ---------------------------
        jCellMin=G2LJ(jCellMinG)
        jCellMax=G2LJ(jCellMaxG)

!       Cycle if element is outside the subdomain (SAMK)
!       ------------------------------------------------
        IF (jCellMax<myJLL .OR. jCellMin>myJUL) CYCLE
          
!       Limit the search to within the subdomain
!       ----------------------------------------
        jCellMin=MAX(jCellMin,myJLL)
        jCellMax=MIN(jCellMax,myJUL)

!       k-direction 
!       ------------
        SELECT CASE (nDim)
        CASE (DIM_2D)
          kCellMin = 1
          kCellMax = nzc

        CASE (DIM_3D)
          kMin = 0
          kMax = nzc+1
          DO k = kMin, kMax-1
            IF ( (zc(k)-zBoundMin)*(zc(k+1)-zBoundMin) < 0.0_CGREAL ) &
              kCellMin = (k+1)-coronaSize
            IF ( (zc(k)-zBoundMax)*(zc(k+1)-zBoundMax) < 0.0_CGREAL ) &
              kCellMax = k+coronaSize
          ENDDO
          IF (kCellMin <= 0 ) kCellMin = 1
          IF (kCellMax == 0 .OR. kCellMax >= nzc+1) kCellMax = nzc

        END SELECT ! nDim

!       Cycle if element is outside the subdomain (SAMK)
!       ------------------------------------------------
        IF (kCellMax<1 .OR. kCellMin>nzc) CYCLE
          
!       Set iblankUndecided to NEGATIVE undecided value for the ring of cells extracted
!       -------------------------------------------------------------------------------
        iblankUndecided(iCellMin:iCellMax, &
                        jCellMin:jCellMax, &
                        kCellMin:kCellMax) = DECIDED

      END DO ! iElem
    END DO ! iBody
    CALL system_clock(iclock3)

!   Loop over cells whose iblank is DECIDED for unstructured surfaces invoking dot normal algorithm 
!   -----------------------------------------------------------------------------------------------
    CALL system_clock(iclock4)
    SELECT CASE(nDim)
      CASE (DIM_2D)
        kMin=1
        kMax=1

      CASE (DIM_3D)
        kMin=1
        kMax=nzc
    END SELECT ! nDim

!   Initialize IblankTemp

!    iblankTemp = 0

    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
     iblankTemp(i,j,k) = 0
    END DO
    END DO
    END DO
  
    DO k=kMin,kMax
    DO j=myJLL,myJUL

!     ---------------------------------------------------------------------------------------------------
!      SAMK: Make sure that, at each row, at least one cell is decided. Then, a ray-tracing search
!            only in X direction will be enough. This will handle the subdomain unawareness problem too!
!     ---------------------------------------------------------------------------------------------------

!     Check if the subdomain is completely in the extension zone (SAMK)
!     -----------------------------------------------------------------
      IF (InExtension .AND. Extended_Outflow==1 .AND. unstruc_surface_type(nbody_Begin)==SOLID_BODY) THEN
        i=nxct-1
        iblankUndecided(MyILL,j,k)=DECIDED
        
!       ------------------------------------------------------------------
!        Search for vertex closest to cells in band
!        Extract the elements that share that vertex
!        Drop the normal and find BI point
!        Check if normal intercept lies inside closest triangular element
!        Authors: Fady and Rajat Oct 18, 2005
!       ------------------------------------------------------------------
        DO iBody = nbody_Begin, nbody_End
            
          CALL search_vertex_dotNorm(iBody,i,j,k,cElement,dotNorm)

          IF (dotNorm >= 0.0_CGREAL) THEN
            iblankTemp(MyILL,j,k) = 1
            IF (unstruc_surface_type(iBody) == SOLID_BODY) bodyNum(MyILL,j,k) = iBody 
            EXIT
          ELSE IF ((iblank(MyILL,j,k) == 1 .AND. unstruc_surface_type(iBody) == SOLID_BODY)) THEN

!           SAMK: Not a provisional fresh cell
!           ----------------------------------
            bodyNum(MyILL,j,k) = bodyNumOld(MyILL,j,k) 
          ENDIF
        ENDDO
      ELSE
        DecidedYet=.FALSE.

        DO i=myILL,myIUL
          IF ( iblankUndecided(i,j,k) == DECIDED ) THEN
            DecidedYet=.TRUE.

!           ------------------------------------------------------------------
!            Search for vertex closest to cells in band
!            Extract the elements that share that vertex
!            Drop the normal and find BI point
!            Check if normal intercept lies inside closest triangular element
!            Authors: Fady and Rajat Oct 18, 2005
!           ------------------------------------------------------------------
            DO iBody = nbody_Begin, nbody_End

              CALL search_vertex_dotNorm(iBody,i,j,k,cElement,dotNorm)
    
              IF (dotNorm >= 0.0_CGREAL) THEN
                iblankTemp(i,j,k) = 1
                IF (unstruc_surface_type(iBody) == SOLID_BODY) bodyNum(i,j,k) = iBody 
                EXIT
              ELSE IF ((iblank(i,j,k) == 1 .AND. unstruc_surface_type(iBody) == SOLID_BODY)) THEN

!               SAMK: Not a provisional fresh cell
!               ----------------------------------
                bodyNum(i,j,k) = bodyNumOld(i,j,k) 
              ENDIF
            ENDDO
          END IF ! iblankUndecided
        ENDDO ! i
       
        IF ((unstruc_surface_type(nBody_Begin) == SOLID_BODY) .AND. (.NOT. DecidedYet)) THEN
          i=myILL
          iblankUndecided(i,j,k)=DECIDED

!         ------------------------------------------------------------------
!          Search for vertex closest to cells in band
!          Extract the elements that share that vertex
!          Drop the normal and find BI point
!          Check if normal intercept lies inside closest triangular element
!          Authors: Fady and Rajat Oct 18, 2005
!         ------------------------------------------------------------------
          DO iBody = nbody_Begin, nbody_End

            CALL search_vertex_dotNorm(iBody,i,j,k,cElement,dotNorm)

            IF (dotNorm >= 0.0_CGREAL) THEN
              iblankTemp(i,j,k) = 1
              bodyNum(i,j,k) = iBody
              EXIT
            ELSE IF ((iblank(i,j,k) == 1)) THEN
!             SAMK: Not a provisional fresh cell
!             ----------------------------------
              bodyNum(i,j,k) = bodyNumOld(i,j,k) 
            ENDIF
          ENDDO
        END IF ! iblankUndecided
      END IF ! InExtension
        
    ENDDO ! j
    ENDDO ! k
    CALL system_clock(iclock5)

!DEBUG
!    write(ifuParLog,*)'VARIABLES="X","Y","Z","IBLANK", "decided"'
!    write(ifuParLog,*)'ZONE F=POINT, I=',myIUL-myILL+1,', J=',myJUL-myJLL+1,' ,K=',nzc
!    do k=1,nzc
!    do j=myJLL,myJUL
!    do i=myILL,myIUL
!       write(ifuParLog,'(3(3X,1PE12.5),3(3X,I10))')xc(L2GI(i)),yc(L2GJ(j)),zc(k),iblankUndecided(i,j,k),iblanktemp(i,j,k),bodyNum(i,j,k)
!    enddo
!    enddo
!    enddo
!return
!call flow_stop
!stop
!DEBUG

    IF (unstruc_surface_type(nBody_Begin) == SOLID_BODY) THEN

!     --------------------------------------------------------------------
!      Set undecided iblank values outside the ring of cells for body
!      by searching horizontal, similar to a ray tracing routine. 
!      Set iblank value at grid cell by searching for first DECIDED VALUE
!      Move along i-direction, j-direction, then k-direction
!     --------------------------------------------------------------------
      CALL system_clock(iclock6)

      DO k=kMin,kMax
      DO j=myJLL,myJUL
      DO i=myILL,myIUL
        IF ( iblankUndecided(i,j,k) == DECIDED ) THEN
          iblankTemp(myILL,j,k) = iblankTemp(i,j,k)
          IF ( iblankTemp(i,j,k) == 1) THEN
            bodyNum(myILL,j,k) = bodyNum(i,j,k)
          END IF
          EXIT
        END IF
      END DO ! i
      END DO ! j
      END DO ! k

!     At this point the left face of each sub-domain is DECIDED.
!     No need for ray tracing in other direcions. (SAMK)
!     ----------------------------------------------------------

      DO k=kMin,kMax
      DO j=myJLL,myJUL
        IF (Extended_Outflow==1) THEN

!         Prior to extension zone (SAMK)
!         ------------------------------
          DO i=myILL+1,MIN(nxct-1,MyIUL)
            IF ( iblankUndecided(i,j,k) == UNDECIDED ) THEN
              iblankTemp(i,j,k) = iblankTemp(i-1,j,k)
              IF ( iblankTemp(i,j,k) == 1) THEN
                bodyNum(i,j,k) = bodyNum(i-1,j,k)
              END IF
            END IF
          END DO ! i

!         In the extension zone (SAMK)
!         ----------------------------
          DO i=MAX(nxct,MyILL+1),myIUL
            iblankTemp(i,j,k) = iblankTemp(i-1,j,k)
            IF ( iblankTemp(i,j,k) == 1) THEN
              bodyNum(i,j,k) = bodyNum(i-1,j,k)
            END IF
          END DO ! i

        ELSE

!         No extension zone
!         -----------------
          DO i=myILL+1,myIUL
            IF ( iblankUndecided(i,j,k) == UNDECIDED ) THEN
              iblankTemp(i,j,k) = iblankTemp(i-1,j,k)
              IF ( iblankTemp(i,j,k) == 1) THEN
                bodyNum(i,j,k) = bodyNum(i-1,j,k)
              END IF
            END IF
          END DO ! i
  
        END IF !Extended_Outflow
      END DO ! j
      END DO ! k

    END IF  !unstruc_surface_type

    CALL find_iblankHoles_outerboundary(1,nbody_begin,iblankTemp)

!   Extend iblank for 2D simulations 
!   --------------------------------
    IF (nDim == DIM_2D) THEN
      DO k=kMin+1,nzc
      DO j=myJLL,myJUL
      DO i=myILL,myIUL
        iblankTemp(i,j,k)      = iblankTemp(i,j,kMin)
        iblankUndecided(i,j,k) = iblankUndecided(i,j,kMin)
        bodyNum(i,j,k) = bodyNum(i,j,kMin)
      END DO !i
      END DO !j
      END DO !k
    END IF ! nDim

!   Extend iblank to the outflow extension Zone
!   -------------------------------------------
!!    DO WHILE (Extended_GLBL==0)
!!      IF (InExtension) THEN
!!      ELSE  ! InExtension
!!        IF (Extended==0) THEN
!!          DO k=1,nzc
!!          DO j=myJLL,MyJUL
!!          DO i=nxc_xt,MyIUL
!!        iblankTemp(i,j,k)      = iblankTemp(i,j,kMin)
!!        iblankUndecided(i,j,k) = iblankUndecided(i,j,kMin)
!!        IF ( bodyNum(i,j,kMin) /= 0 ) bodyNum(i,j,k) = bodyNum(i,j,kMin)
!!          END DO
!!          END DO
!!          END DO
!!        END IF
!!      END IF
!!  END DO

    CALL par_comm_outermost_int(bodyNum,myTime)

    IF (unstruc_surface_type(nBody_Begin) == MEMBRANE) THEN
      DO k=1,nzc
      DO j=myJLL,myJUL
      DO i=myILL,myIUL
        IF ( iblankUndecided(i,j,k) == DECIDED ) THEN
          IF ( iblankTemp(i,j,k) == 1 ) THEN
            IF ( iblank(i,j,k) == 0 .AND. ntime > 1 ) THEN
              fresh_cell(i,j,k) = 1
              num_fresh         = num_fresh+1
              WRITE(ifuFreshCellOut,*) ntime,i,j,k,'   --- provisional fresh cell'
            ENDIF
          ELSE
            IF ( iblank(i,j,k) == 1 .AND. ntime > 1 ) THEN
              fresh_cell(i,j,k) = 1
              num_fresh         = num_fresh+1
              WRITE(ifuFreshCellOut,*) ntime,i,j,k,'   --- provisional fresh cell'
            ENDIF
          ENDIF
        ENDIF
        iblank(i,j,k) = iblankTemp(i,j,k)
      END DO ! i
      END DO ! j
      END DO   ! k
    ELSE

      DO k=1,nzc
      DO j=myJLL,myJUL
      DO i=myILL,myIUL
        IF  ( iblankTemp(i,j,k) == 1 ) THEN
          IF ( iblank(i,j,k) == 0  .AND. ntime > 1 ) THEN
             !dddWRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- dead cell'
             num_dead = num_dead+1
          ENDIF
        ELSE
          IF ( iblank(i,j,k) == 1 .AND. ntime > 1 ) THEN
            fresh_cell(i,j,k)=  1
            num_fresh        = num_fresh+1
            !dddWRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- provisional  fresh cell'
          ELSE IF (bodyNum(i,j,k) /= 0) THEN
            WRITE(*,*) bodyNum(i, j, k), iblankTemp(i, j, k), iblank(i, j, k),myrank
            bodyNum(i,j,k)=0
            WRITE(*,'(A,4I7)') 'Error: a bodyNum is assigned to a non-Fress fluid cell',i,j,k, myrank
            CALL flow_stop
            STOP
          ENDIF
        ENDIF
        iblank(i,j,k) = iblankTemp(i,j,k)
      END DO ! i
      END DO ! j
      END DO   ! k
    ENDIF

    sumBlank=SUM(IBLANK(1:nxc,1:nyc,1:nzc))
#   ifdef MPI
      CALL par_getSumInteger(num_dead)
      CALL par_getSumInteger(num_fresh)
      CALL par_getSumInteger(sumBlank)
#   endif

    IF (monitorON) THEN
      WRITE(STDOUT,'(7X,A,I8)') 'Number of Dead              Cells = ',num_dead
      WRITE(STDOUT,'(7X,A,I8)') 'Number of Provisional Fresh Cells = ',num_fresh
      WRITE(STDOUT,'(7X,A,I8)') 'Number of BLANK             Cells = ',sumBlank
    ENDIF ! ntime

    CALL system_clock(iclock2,clock_rate)
    IF ( monitorON )  THEN
      WRITE(STDOUT,'(7X,A,1PE15.7)') 'CPU Time for Fast Initial Iblank = ',&
        REAL(iclock2-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(STDOUT,'(9X,A,1PE15.7)') 'CPU Time for Initial Setup     = ',&
        REAL(iclock3-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(STDOUT,'(9X,A,1PE15.7)') 'CPU Time for dotNorm           = ',&
        REAL(iclock5-iclock4,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
      WRITE(STDOUT,'(9X,A,1PE15.7)') 'CPU Time for Filling Iblank    = ',&
         REAL(iclock2-iclock6,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
    END IF ! ntime

!DEBUG
 !   write(111,*)'VARIABLES="X","Y","Z","IBLANK", "BodyNum"'
 !   write(111,*)'ZONE F=POINT, I=',myIUL-myILL+1,', J=',myJUL-myJLL+1,' ,K=',1!nzc
 !   do k=1,1!nzc
 !   do j=myJLL,myJUL
 !   do i=myILL,myIUL
 !      write(111,'(3(3X,1PE12.5),3(3X,I10))')xc(L2GI(i)),yc(L2GJ(j)),zc(k),iblank(i,j,k),bodyNum(i,j,k)
 !   enddo
 !   enddo
 !   enddo
 ! stop
!    write(*,*) ifuParLog
!DEBUG


!  Deallocate local memory

   DEALLOCATE(iblankUndecided, STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'set_iblank_canonical_body_fast: Memory Deallocation Error for iblankUndecided'
      CALL flow_stop
      STOP
    ENDIF ! iErr

   DEALLOCATE(iblankTemp, STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'set_iblank_canonical_body_fast: Memory Deallocation Error for iblankUndecided'
      CALL flow_stop
      STOP
    ENDIF ! iErr

   END SUBROUTINE set_iblank_canonical_body_fast
!------------------------------------------------------------------------------



SUBROUTINE search_vertex_dotNorm(iBody,iCell,jCell,kCell,closestElement,dotNorm)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!   parameters
!   ----------
    INTEGER, INTENT(IN)             :: iBody,iCell,jCell,kCell
    INTEGER, INTENT(OUT)            :: closestElement
    REAL(KIND=CGREAL) , INTENT(OUT) :: dotNorm

!   loop variables
!   --------------
    INTEGER :: m,n,nc
    
!   local variables
!   ---------------
    INTEGER, PARAMETER :: NSIZE = 1000
    INTEGER, PARAMETER :: MSIZE = 20

    INTEGER,           DIMENSION(:),ALLOCATABLE :: NeighElemInd
    REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE :: distMarkerS

    INTEGER :: iErr,nBodySelect,numNeighElement
    INTEGER :: elemInd,node1,node2,node3,nMarker
    INTEGER :: nCheck
    INTEGER :: shortestProbe

    INTEGER,DIMENSION(1)     :: iDummy(1)
    INTEGER,DIMENSION(MSIZE) :: cElement,closestVert

    REAL(KIND=CGREAL) :: xCell,yCell,zCell
    REAL(KIND=CGREAL) :: dsIntercept,xM,yM,zM
    REAL(KIND=CGREAL) :: areaDiffMin
    REAL(KIND=CGREAL) :: distBIElem, distBIElemMin
    REAL(KIND=CGREAL) :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL) :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL) :: area123,area124,area234,area314
    REAL(KIND=CGREAL) :: semiPerimeter123,semiPerimeter124,semiPerimeter234,semiPerimeter314
    REAL(KIND=CGREAL) :: epsiArea,areaDiff
    REAL(KIND=CGREAL) :: xBI, yBI, zBI
    REAL(KIND=CGREAL) :: xBITemp, yBITemp, zBITemp

    REAL(KIND=CGREAL),DIMENSION(3)     :: xVert, yVert, zVert
    REAL(KIND=CGREAL),DIMENSION(MSIZE) :: xBIT,yBIT,zBIT,dist



    xCell = xc(L2GI(iCell))
    yCell = yc(L2GJ(jCell))
    zCell = zc(     kCell)

    dotNorm = -1.0_CGREAL

    nMarker = nPtsBodyMarker(iBody)

!   nCheck:  Number of closesest nodes to check
!   and also CPU time for finding body intercept.
!   ---------------------------------------------
    nCheck = 1

    IF (nCheck > MSIZE) THEN
      PRINT*,'nCheck in GCM_calc_bodyIntercept_Unstruc is limited to', MSIZE
      PRINT*,'Increase array size'
      CALL flow_stop
      STOP
    ENDIF

!   Allocate local array 
!   --------------------
    ALLOCATE(distMarkerS(nMarker),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Allocation Error for distMarker'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    ALLOCATE(NeighElemInd(NSIZE),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Allocation Error for NeighElemInd'
      CALL flow_stop
      STOP
    ENDIF ! ierr 

!   Get closestMarker for cell 
!   --------------------------
    DO m = 1, nMarker
      xM = xBodyMarker(iBody,m)
      yM = yBodyMarker(iBody,m)
      zM = zBodyMarker(iBody,m)
      distMarkerS(m) = (xM-xCell)**2 + (yM-yCell)**2 + (zM-zCell)**2
    ENDDO ! m

    DO nc = 1,NCheck
      iDummy                       = MINLOC(distMarkerS(1:nMarker))
      closestVert(nc)              = iDummy(1)
      distMarkerS(closestVert(nc)) = 1.0E20_CGREAL
    ENDDO


!   Find elements that share closest node/marker
!   --------------------------------------------
    DO nc = 1, NCheck

    numNeighElement = 0
    DO m=1,totNumTriElem(iBody)
      IF ( triElemNeig(iBody,1,m) == closestVert(nc) .OR. &
           triElemNeig(iBody,2,m) == closestVert(nc) .OR. &
           triElemNeig(iBody,3,m) == closestVert(nc)      ) THEN
        numNeighElement               = numNeighElement + 1
        NeighElemInd(numNeighElement) = m
      ENDIF
    ENDDO ! m

!   Trap error if array NeighElemenInd overflows 
!   --------------------------------------------
    IF ( numNeighElement > NSIZE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Overflow Error for NeighElemInd'
      WRITE(STDOUT,*) ' Allocated size = ',NSIZE
      WRITE(STDOUT,*) ' Current size   = ',numNeighElement
      WRITE(STDOUT,*) ' Aborting Run'
      CALL flow_stop
      STOP
    ENDIF ! NeighElemInd

!   Check which element contains normal intercept
!   ---------------------------------------------
    distBIElemMin = 1.0E+16_CGREAL
    areaDiffMin   = 1.0E+16_CGREAL
    epsiArea      = 1.0E-4_CGREAL

    closestElement = 0

    DO n = 1,numNeighElement
      elemInd = NeighElemInd(n)

      node1   = triElemNeig(iBody,1,elemInd)
      node2   = triElemNeig(iBody,2,elemInd)
      node3   = triElemNeig(iBody,3,elemInd)

!     Check if BI inside the triangle of the surface element through area differences
!     -------------------------------------------------------------------------------    
      CALL check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xCell,yCell,zCell,xBITemp,yBITemp,zBITemp,area123,areaDiff)

!     --------------------------------------------------------------------
!      Select closest Elem and BI coordinates:
!      If BI falls inside the element use that
!      Else Base the selection on the minimum distance 
!      between BI and either the norm to closest side or vertices of side
!     --------------------------------------------------------------------
      IF ( DABS(areaDiff) < epsiArea*area123) THEN
        xBI = xBITemp
        yBI = yBITemp
        zBI = zBITemp
        closestElement = elemInd
        dist(nc) = 0.0_CGREAL
        GOTO 999
      ELSE
        CALL calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,xBITemp,yBITemp,zBITemp,distBIElem, area123, areaDiff)

        IF (distBIElem <= distBIElemMin) THEN
          distBIElemMin = distBIElem
          closestElement = elemInd
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp
        ENDIF ! distBIElem

        dist(nc) = distBIElemMin

      ENDIF ! areaDiff
    ENDDO ! n

999 CONTINUE

    xBIT(nc) = xBI
    yBIT(nc) = yBI
    zBIT(nc) = zBI
    cElement(nc) = closestElement
    ENDDO  ! nc

    iDummy         = MINLOC(dist(1:nCheck))
    shortestProbe  = iDummy(1)
    closestElement = cElement(shortestProbe)

!   Perform the dot product with element that has shortest distance or area
!   -----------------------------------------------------------------------
    dotNorm = (xCell - triElemCentx(iBody,closestElement)) &
                      *triElemNormx(iBody,closestElement)  &
             +(yCell - triElemCenty(iBody,closestElement)) &
                      *triElemNormy(iBody,closestElement)  &
             +(zCell - triElemCentz(iBody,closestElement)) &
                      *triElemNormz(iBody,closestElement)

!   Deallocate local array 
!   ----------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)      
    DEALLOCATE(NeighElemInd,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Deallocation Error for NeighElemInd,NSIZE=',NSIZE      
      STOP   
      CALL flow_stop
    ENDIF ! ierr 
    
    DEALLOCATE(distMarkerS, STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Deallocation Error for distMarkerS,nMarker=',nMarker      
      STOP   
      CALL flow_stop
    ENDIF ! ierr  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
END SUBROUTINE search_vertex_dotNorm
!---------------------------------------------------------------------


SUBROUTINE check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xCell,yCell,zCell,&
                                  xBITemp,yBITemp,zBITemp,area123,areaDiff)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN) :: elemInd,iBody,node1,node2,node3
    REAL(KIND=CGREAL), INTENT(IN)  :: xCell,yCell,zCell
    REAL(KIND=CGREAL), INTENT(OUT) :: xBITemp,yBITemp,zBITemp,area123,areaDiff

!... loop variables

    INTEGER :: iside

!... local variables

    REAL(KIND=CGREAL)        :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL)        :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL)        :: area124,area234,area314
    REAL(KIND=CGREAL)        :: semiPerimeter123,semiPerimeter124
    REAL(KIND=CGREAL)        :: semiPerimeter234,semiPerimeter314
	
    REAL(KIND=CGREAL)        :: denom, num1, num2
    REAL(KIND=CGREAL)        :: VX, VY, VZ, PX, PY, PZ, RR, Rx, Ry, Rz
    REAL(KIND=CGREAL)        :: UX, UY, UZ, WX, WY, WZ, SS, TT, PIx, PIy, PIz, tol	


! ******************************************************************************
! equation of plane (note our normals are unit normals)
     
!  n  x + n  y + n  z + planeConst = 0
!   x         y      z
! ******************************************************************************

    planeConst =- triElemNormx(iBody,elemInd)*xBodyMarker(iBody,node1) &
                - triElemNormy(iBody,elemInd)*yBodyMarker(iBody,node1) &
                - triElemNormz(iBody,elemInd)*zBodyMarker(iBody,node1)
       
! ******************************************************************************
! Compute coordinates of normal intercept
!      
! Consider point Po = (xo,yo,zo)  not on plane
!   and  point   P1 = (x1,y1,z1)  on the plane
!                                                               ^
! equation of line through Po normal to plane is  P(s) = Po + s n
!      
! normal distance from Po to Plane is given by
!            ^                               ^         ^ ^
!            n. ( P(s) - P1 ) = 0   => so = -n.(Po-P1)/n.n                     ^ ^
!                                         = -(n xo + n yo + n zo + planeConst)/n.n
!                                              x      y      z
!
!                                                 ^
!   subsequently normal intersection point = Po + so n  
!                   
! ******************************************************************************

     distanceToPlane = -(  triElemNormx(iBody,elemInd)*xCell  &
                         + triElemNormy(iBody,elemInd)*yCell  &
                         + triElemNormz(iBody,elemInd)*zCell  &
                         + planeConst )

     xBITemp = xCell + triElemNormx(iBody,elemInd)*distanceToPlane
     yBITemp = yCell + triElemNormy(iBody,elemInd)*distanceToPlane
     zBITemp = zCell + triElemNormz(iBody,elemInd)*distanceToPlane
	 
!*******************************************************************************
! Check whether BI on the triangular element or not
!
! New Method: source ( http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm )
!
!*******************************************************************************

	 
	tol = 1.e-10
	
	RX = xBITemp
	RY = yBITemp
	RZ = zBITemp
	
	UX = xBodyMarker(iBody,node2) - xBodyMarker(iBody,node1)
	UY = yBodyMarker(iBody,node2) - yBodyMarker(iBody,node1)
	Uz = zBodyMarker(iBody,node2) - zBodyMarker(iBody,node1)
	
	VX = xBodyMarker(iBody,node3) - xBodyMarker(iBody,node1)
	VY = yBodyMarker(iBody,node3) - yBodyMarker(iBody,node1)
	Vz = zBodyMarker(iBody,node3) - zBodyMarker(iBody,node1)	
	
	WX = RX - xBodyMarker(iBody,node1)
	WY = RY - yBodyMarker(iBody,node1)
	Wz = RZ - zBodyMarker(iBody,node1)		
	
	denom = (UX*VX+UY*VY+UZ*VZ)**2 - (UX*UX+UY*UY+UZ*UZ)*(VX*VX+VY*VY+VZ*VZ)
	
	num1 = (UX*VX+UY*VY+UZ*VZ)*(WX*VX+WY*VY+WZ*VZ) - (VX*VX+VY*VY+VZ*VZ)*(WX*UX+WY*UY+WZ*UZ)
	
	num2 = (UX*VX+UY*VY+UZ*VZ)*(WX*UX+WY*UY+WZ*UZ) - (UX*UX+UY*UY+UZ*UZ)*(WX*VX+WY*VY+WZ*VZ)
	
	SS = num1/denom
	TT = num2/denom
	
	IF( (SS>=(0.0-tol)) .and. (TT>=(0.0-tol)) .and. (SS+TT<=(1.0+tol)) ) THEN  !! Inside Triangle
	
    area123  = 1.0
	areadiff = 0.0
	
	ELSE   
	
	area123  = 1.0
	areadiff = 1.0
	
    ENDIF

END SUBROUTINE check_BIInsideTriangle 
!---------------------------------------------------------------------



SUBROUTINE check_BIInsideTriangle_old(iBody,elemInd,node1,node2,node3,xCell,yCell,zCell,&
                                  xBITemp,yBITemp,zBITemp,area123,areaDiff)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!   parameters
!   ----------
    INTEGER, INTENT(IN) :: elemInd,iBody,node1,node2,node3

    REAL(KIND=CGREAL), INTENT(IN)  :: xCell,yCell,zCell
    REAL(KIND=CGREAL), INTENT(OUT) :: xBITemp,yBITemp,zBITemp,area123,areaDiff

!   loop variables
!   --------------
    INTEGER :: iside

!   local variables
!   ---------------
    REAL(KIND=CGREAL) :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL) :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL) :: area124,area234,area314
    REAL(KIND=CGREAL) :: semiPerimeter123,semiPerimeter124
    REAL(KIND=CGREAL) :: semiPerimeter234,semiPerimeter314




! ******************************************************************************
! equation of plane (note our normals are unit normals)
     
!  n  x + n  y + n  z + planeConst = 0
!   x         y      z
! ******************************************************************************

    planeConst =- triElemNormx(iBody,elemInd)*xBodyMarker(iBody,node1) &
                - triElemNormy(iBody,elemInd)*yBodyMarker(iBody,node1) &
                - triElemNormz(iBody,elemInd)*zBodyMarker(iBody,node1)
       
! ******************************************************************************
! Compute coordinates of normal intercept
!      
! Consider point Po = (xo,yo,zo)  not on plane
!   and  point   P1 = (x1,y1,z1)  on the plane
!                                                               ^
! equation of line through Po normal to plane is  P(s) = Po + s n
!      
! normal distance from Po to Plane is given by
!            ^                               ^         ^ ^
!            n. ( P(s) - P1 ) = 0   => so = -n.(Po-P1)/n.n                     ^ ^
!                                         = -(n xo + n yo + n zo + planeConst)/n.n
!                                              x      y      z
!
!                                                 ^
!   subsequently normal intersection point = Po + so n  
!                   
! ******************************************************************************

     distanceToPlane = -(  triElemNormx(iBody,elemInd)*xCell  &
                         + triElemNormy(iBody,elemInd)*yCell  &
                         + triElemNormz(iBody,elemInd)*zCell  &
                         + planeConst )

     xBITemp = xCell + triElemNormx(iBody,elemInd)*distanceToPlane
     yBITemp = yCell + triElemNormy(iBody,elemInd)*distanceToPlane
     zBITemp = zCell + triElemNormz(iBody,elemInd)*distanceToPlane

! ******************************************************************************
! Check to see if normal intercept lies inside the closest trianglular element
!               3 
!               *  .
!              /  \   .
!             /    \    .
!            /      \    * 4=BI
!           /        \  .
!         1*__________*2
!         
! Basic Idea :  IF [ AREA(124) + AREA(234) + AREA(314) ] > AREA(123) THEN  POINT(4) is
! outside triangle (123)
!
! using Heron formula for area of triangle
! AREA(123) = SQRT[ S * ( S - S12) * (S - S23) * (S - S31) ]
! S = 0.5*(S12 + S23 + S31) 
! ******************************************************************************

     side12 =  SQRT( (xBodyMarker(iBody,node2)-xBodyMarker(iBody,node1))**2  &
                    +(yBodyMarker(iBody,node2)-yBodyMarker(iBody,node1))**2  &
                    +(zBodyMarker(iBody,node2)-zBodyMarker(iBody,node1))**2  )
     side23 =  SQRT( (xBodyMarker(iBody,node3)-xBodyMarker(iBody,node2))**2  &
                    +(yBodyMarker(iBody,node3)-yBodyMarker(iBody,node2))**2  &
                    +(zBodyMarker(iBody,node3)-zBodyMarker(iBody,node2))**2  )
     side31 =  SQRT( (xBodyMarker(iBody,node1)-xBodyMarker(iBody,node3))**2  &
                    +(yBodyMarker(iBody,node1)-yBodyMarker(iBody,node3))**2  &
                    +(zBodyMarker(iBody,node1)-zBodyMarker(iBody,node3))**2  )
     side14 =  SQRT( (xBITemp-xBodyMarker(iBody,node1))**2  &
                    +(yBITemp-yBodyMarker(iBody,node1))**2  &
                    +(zBITemp-zBodyMarker(iBody,node1))**2  )
     side24 =  SQRT( (xBITemp-xBodyMarker(iBody,node2))**2  &
                    +(yBITemp-yBodyMarker(iBody,node2))**2  &
                    +(zBITemp-zBodyMarker(iBody,node2))**2  )
     side34 =  SQRT( (xBITemp-xBodyMarker(iBody,node3))**2  &
                    +(yBITemp-yBodyMarker(iBody,node3))**2  &
                    +(zBITemp-zBodyMarker(iBody,node3))**2  )

     semiPerimeter123 = 0.5_CGREAL*(side12 + side23 + side31)
     semiPerimeter124 = 0.5_CGREAL*(side12 + side24 + side14)
     semiPerimeter234 = 0.5_CGREAL*(side23 + side24 + side34)
     semiPerimeter314 = 0.5_CGREAL*(side31 + side34 + side14)

     area123       = SQRT( DABS(semiPerimeter123*(semiPerimeter123-side12) &
                                           *(semiPerimeter123-side23) &
                                           *(semiPerimeter123-side31))   )
   
     area124       = SQRT( DABS(semiPerimeter124*(semiPerimeter124-side12) &
                                           *(semiPerimeter124-side24) &
                                           *(semiPerimeter124-side14))   )
    
     area234       = SQRT( DABS(semiPerimeter234*(semiPerimeter234-side23) &
                                           *(semiPerimeter234-side24) &
                                           *(semiPerimeter234-side34))   )

     area314       = SQRT( DABS(semiPerimeter314*(semiPerimeter314-side31) &
                                           *(semiPerimeter314-side34) &
                                           *(semiPerimeter314-side14))   )

     areaDiff  = area124 + area234 + area314 - area123

END SUBROUTINE check_BIInsideTriangle_old 
!---------------------------------------------------------------------



SUBROUTINE calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,  &
                                  xBITemp,yBITemp,zBITemp,distBIElem, area1, aread)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN)            :: elemInd,iBody,node1,node2,node3
    REAL(KIND=CGREAL), INTENT(IN)  :: xBITemp,yBITemp,zBITemp
    REAL(KIND=CGREAL), INTENT(OUT) :: distBIElem
    REAL(KIND=CGREAL) :: area1, aread

!... loop variables

    INTEGER :: iside, i

!... local variables

    INTEGER :: isideSelect,mside
    INTEGER :: nodeSelect1,nodeSelect2

    REAL(KIND=CGREAL) :: aCrossbVectMagn,dotVal
    REAL(KIND=CGREAL) :: distIntBI,distIntCG,distIntMin
    REAL(KIND=CGREAL) :: distNorm,distNode1BINorm,distNode2BINorm
    REAL(KIND=CGREAL) :: distVert1BI,distVert2BI
    REAL(KIND=CGREAL) :: magnitude12,magnitudeBICG,magnitude,projectedLength
    REAL(KIND=CGREAL) :: node12x,node12y,node12z
    REAL(KIND=CGREAL) :: vec01x,vec01y,vec01z
    REAL(KIND=CGREAL) :: vec12x,vec12y,vec12z
    REAL(KIND=CGREAL) :: xBINorm,yBINorm,zBINorm
    REAL(KIND=CGREAL) :: xCG,yCG,zCG
    REAL(KIND=CGREAL), DIMENSION(3) :: aVect,bVect,cVect,aCrossbVect,cCrossbVect
    REAL(KIND=CGREAL), DIMENSION(3) :: xInt,yInt,zInt
    REAL(KIND=CGREAL), DIMENSION(3) :: xVert,yVert,zVert
    REAL(KIND=CGREAL), DIMENSION(3) :: vect1,vect2,vect3,vect4
    REAL(KIND=CGREAL), DIMENSION(3,3) :: vectInt


          
    xVert(1)  = xBodyMarker(iBody,node1)
    yVert(1)  = yBodyMarker(iBody,node1)
    zVert(1)  = zBodyMarker(iBody,node1)

    xVert(2)  = xBodyMarker(iBody,node2)
    yVert(2)  = yBodyMarker(iBody,node2)
    zVert(2)  = zBodyMarker(iBody,node2)

    xVert(3)  = xBodyMarker(iBody,node3)
    yVert(3)  = yBodyMarker(iBody,node3)
    zVert(3)  = zBodyMarker(iBody,node3)

    xCG       = triElemCentx(iBody,elemInd)
    yCG       = triElemCenty(iBody,elemInd)
    zCG       = triElemCentz(iBody,elemInd)

! ============================================================================
!   Construct Intersection points between 
!     line linking BI and Centroid of Surface Element and the triangle sides
!     L1: BI-CG, L2:Sides of Vertices
!
!   use formula for intersection point between 2 co-planar lines from Mathworld
!   http://mathworld.wolfram.com/Line-LineIntersect.html
!
!             x4
!             *
!             |
!             |
!    x1       | Int      x2
!    *--------*----------*
!             |
!             |  
!             | 
!             * x3
!
! ============================================================================
   
    vect1(1:3) = (/xBITemp,yBITemp,zBITemp/)
    vect2(1:3) = (/xCG,yCG,zCG/)

    aVect(1:3) = vect2(1:3)-vect1(1:3)

    DO iside = 1,3
      mside = iside +1
      IF(iside == 3) mside = 1

      vect3(1:3) =(/xVert(iside),yVert(iside),zVert(iside)/)
      vect4(1:3) =(/xVert(mside),yVert(mside),zVert(mside)/)

      bVect(1:3)  = vect4(1:3) -vect3(1:3)
      cVect(1:3)  = vect3(1:3) -vect1(1:3)

      call calc_crossProduct(aVect,bVect,aCrossbVect)
      call calc_crossProduct(cVect,bVect,cCrossbVect)

      aCrossbVectMagn = aCrossbVect(1)**2 + aCrossbVect(2)**2 +aCrossbVect(3)**2
      dotVal = DOT_PRODUCT(cCrossbVect,aCrossbVect)

      vectInt(1:3,iside) = vect1(1:3) + aVect(1:3) *dotVal/aCrossbVectMagn
    END DO ! iside 
      
! ============================================================================
!   Choose closest intersection point lying between BI and CG
!     Normalsize value with L1
! ============================================================================

    magnitudeBICG = SQRT( (vect1(1)-vect2(1))**2 &
                        + (vect1(2)-vect2(2))**2 &
                        + (vect1(3)-vect2(3))**2 )

    distIntMin = 1.0E+16_CGREAL
    isideSelect = -1000

    DO iside = 1,3
      distIntBI = SQRT( (vect1(1)-vectInt(1,iside))**2 &
                      + (vect1(2)-vectInt(2,iside))**2 &
                      + (vect1(3)-vectInt(3,iside))**2 )/magnitudeBICG

      distIntCG = SQRT( (vect2(1)-vectInt(1,iside))**2 &
                      + (vect2(2)-vectInt(2,iside))**2 &
                      + (vect2(3)-vectInt(3,iside))**2 )/magnitudeBICG

      IF ( distIntBI <= 1.0_CGREAL .AND. distIntCG <= 1.0_CGREAL) THEN
        distIntMin  = DMIN1(distIntBI,distIntCG)
        isideSelect = iside
      END IF ! distIntBI
    END DO ! iside 

! ============================================================================
!   Trap error for isideSelect 
! ============================================================================
   
     IF ( isideSelect < 0 ) THEN
      WRITE(STDOUT,*) &
       'calc_BIdistMin: Incorrect selection of iside (Should be either 1, 2 or 3'
      WRITE(STDOUT,*) &   
       '                default value selected = ',isideSelect
      WRITE(STDOUT,*) myrank
      OPEN(22222, file = 'outBItriangle.err')
      WRITE(22222,*)'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(22222,*)'VARIABLES= "X","Y","Z"'
      WRITE(22222,*) 'ZONE T="unstruc" N= 6 E= 2 F=FEPOINT  ET=TRIANGLE'
      DO i = 1, 3
       WRITE(22222,*) xVert(i), yVert(i), zVert(i)
      ENDDO !i
       WRITE(22222,*) xBITemp,yBITemp,zBITemp
       WRITE(22222,*) xCG,yCG,zCG
       WRITE(22222,*) xCG,yCG,zCG
       WRITE(22222,*) 1, 2, 3
       WRITE(22222,*) 4, 5, 6       	
      CLOSE(22222)
      WRITE(STDOUT,*) area1, aread
      STOP
     END IF ! isideSelect 

! ============================================================================
!   Select appropriate vertices from isideSelect 
! ============================================================================
   
     SELECT CASE(isideSelect)
       CASE(1)
         nodeSelect1 = 1
         nodeSelect2 = 2
       CASE(2)
         nodeSelect1 = 2
         nodeSelect2 = 3
       CASE(3)
         nodeSelect1 = 3
         nodeSelect2 = 1
     END SELECT ! isideSelect

! ============================================================================
!   Drop normals from BI to selected side 
!    and find coordinates of intersection point
!
!   unit vector from node 1 to node 2
!   use formula for distance between point and line from Mathworld
!   http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
! 
!    x1                  x2
!    *-------------------*
!             |                              |(x2-x1) x (x1-x0)|
!             | d                       d =  --------------------
!             |                                   |x2-x1|
!             * x0
!   x0: BI
! ============================================================================

     vec12x = xVert(nodeSelect2) - xVert(nodeSelect1)
     vec12y = yVert(nodeSelect2) - yVert(nodeSelect1)
     vec12z = zVert(nodeSelect2) - zVert(nodeSelect1)

     magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)

     vec01x = xVert(nodeSelect1) - xBITemp
     vec01y = yVert(nodeSelect1) - yBITemp
     vec01z = zVert(nodeSelect1) - zBITemp
     
     distNorm = SQRT(  (vec12y*vec01z - vec12z*vec01y)**2  &
                     + (vec12z*vec01x - vec01z*vec12x)**2  &
                     + (vec12x*vec01y - vec01x*vec12y)**2  )/magnitude12

! ============================================================================
!    Project vector BI-node1 onto node12 to find body intercept point
! ============================================================================

     node12x = xVert(nodeSelect2) - xVert(nodeSelect1)
     node12y = yVert(nodeSelect2) - yVert(nodeSelect1)
     node12z = zVert(nodeSelect2) - zVert(nodeSelect1)

     magnitude = SQRT(node12x**2 + node12y**2 + node12z**2)

     node12x = node12x/magnitude
     node12y = node12y/magnitude
     node12z = node12z/magnitude

     projectedLength = (xBITemp - xVert(nodeSelect1))*node12x  &
                      +(yBITemp - yVert(nodeSelect1))*node12y  &
                      +(zBITemp - zVert(nodeSelect1))*node12z

     xBINorm = xVert(nodeSelect1) + projectedLength*node12x
     yBINorm = yVert(nodeSelect1) + projectedLength*node12y
     zBINorm = zVert(nodeSelect1) + projectedLength*node12z

! ============================================================================
!    Determine distance between BINorm and vertices of selected side.
!     If normal point lies inside the side, select that distance.
!     If it lies outside, find the minimum distance with vertices
!     Use normalized length
! ============================================================================

      distNode1BINorm = SQRT( (xVert(nodeSelect1)-xBINorm)**2 &
                            + (yVert(nodeSelect1)-yBINorm)**2 &
                            + (zVert(nodeSelect1)-zBINorm)**2 )/magnitude

      distNode2BINorm = SQRT( (xVert(nodeSelect2)-xBINorm)**2 &
                            + (yVert(nodeSelect2)-yBINorm)**2 &
                            + (zVert(nodeSelect2)-zBINorm)**2 )/magnitude

      IF ( distNode1BINorm <= 1.0_CGREAL .AND. distNode2BINorm <= 1.0_CGREAL) THEN
        distBIElem  =  distNorm
      ELSE
        distVert1BI = SQRT( (xVert(nodeSelect1)-xBITemp)**2 &
                          + (yVert(nodeSelect1)-yBITemp)**2 &
                          + (zVert(nodeSelect1)-zBITemp)**2 )

        distVert2BI = SQRT( (xVert(nodeSelect2)-xBITemp)**2 &
                          + (yVert(nodeSelect2)-yBITemp)**2 &
                          + (zVert(nodeSelect2)-zBITemp)**2 )

        IF (distVert1BI <= distVert2BI) THEN
           distBIElem  = distVert1BI
        ELSE
           distBIElem  = distVert2BI
        ENDIF

      END IF ! distNode1BINorm

END SUBROUTINE calc_BIOutsideTriangle 
!---------------------------------------------------------------------



SUBROUTINE calc_crossProduct(r,s,cross_product)
   
    USE global_parameters
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN)  :: r,s
    REAL(KIND=CGREAL), DIMENSION(3), INTENT(OUT) :: cross_product 

    INTEGER :: component,i,j



    DO component = 1,3
      i = MODULO(component,3) + 1
      j = MODULO(i,3) + 1
      cross_product(component) = r(i)*s(j) - s(i)*r(j)
    END DO ! component 

END SUBROUTINE calc_crossProduct
!---------------------------------------------------------------------
