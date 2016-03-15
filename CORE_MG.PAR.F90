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
!  Fady Najjar
!  Rajat Mittal
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
!  Filename: CORE_MG.PAR.F90
!  Latest Modification: Oct, 20 2010 (ver. PAT 1.2.0)
!  by JHSeo
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!
!  Purpose: Routines to solve the PPE with Multigrid approach 
!           using V, W and F-cycles
!
!  Copyright: (c) 2007 by the George Washington University
!
!  Originally written by H. Dong
!  Reformulated by F. Najjar and R. Mittal
!
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     MG_initial_domain()
!     MG_initial_grid()
!     MG_initial() 
!     MG_Allocate_Memory()
!     MG_Solver(var,rr,compTime,commTime)
!     MG_solv_x()
!     MG_solv_y()
!     MG_solv_z()
!     MG_Nullify_Arrays(dir)
!     MG_Prepare_Iblank(rr, NC_dim)
!     MG_FINE_TO_COARSE_X(var,rr,levBegin,nlev)
!     MG_COARSE_TO_FINE_X(var,rr,levEnd,nlev)
!     MG_FINE_TO_COARSE_Y(var,rr,levBegin,nlev)
!     MG_COARSE_TO_FINE_Y(var,rr,levEnd,nlev)
!     MG_FINE_TO_COARSE_Z(var,rr,levBegin,nlev)
!     MG_COARSE_TO_FINE_Z(var,rr,levEnd,nlev)
!     MG_Prepare_BC(var)
!     MG_Prepare(var,IL,JL,KL)
!     MG_itsolv(var,r,nLevX,nLevY,nLevZ,IL,JL,KL)
!     MG_residual(var,rrr,nLevX,nLevY,nLevZ,IL,JL,KL,resCheckFlag)
!     MG_Restrict(var,r,nlev,NC_dim)
!     MG_Prolong(phi,correc,nlev,NC_dim)
!------------------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1
# define L2GI_MG(l,i) MGX(l)%myIs+i-1
# define L2GJ_MG(l,j) MGY(l)%myJs+j-1



SUBROUTINE MG_initial_domain()

! -------------------------------------------------------------
!  Purpose: Compute MG levels, initialize MG grid-based arrays
!           Calculate grid numbers at each level.
! -------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER, PARAMETER :: MG_LEVEL_MAX=10

    INTEGER:: ii, jj , kk, ilevel



!   ---------------------------------------------------
!    Calculate levels along each direction separately.
!    NOTE: No node retiring for now
!   ---------------------------------------------------
    IF (mgLevels_X==0) THEN
      ii = nxc_GLBL
      DO ilevel = 2, MG_LEVEL_MAX
        ii = (ii+1)/2 
        IF (ii>max(nCoarseGrids,Np(1))) CYCLE
        EXIT 
      END DO ! ilevel
      mgLevels_X = ilevel
    END IF ! mgLevels_X

    IF (mgLevels_Y==0) THEN
      jj = nyc_GLBL
      DO ilevel = 2, MG_LEVEL_MAX
        jj = (jj+1)/2
        IF (jj>max(nCoarseGrids,Np(2))) CYCLE 
        EXIT  
      END DO ! ilevel 
      mgLevels_Y = ilevel
    END IF ! mgLevels_Y

    IF (mgLevels_Z==0 .and. ndim==DIM_3D) THEN
      kk = nzc
      DO ilevel = 2, MG_LEVEL_MAX
        kk = (kk+1)/2
        IF (kk>nCoarseGrids) CYCLE 
        EXIT  
      END DO ! ilevel
      mgLevels_Z = ilevel
    ELSE IF (ndim==DIM_2D) THEN
        mgLevels_Z = 1
    END IF ! mgLevels_Z
    
    if (monitorON) then
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(7X,A,1X,I5)') 'mgLevels_X = ',mgLevels_X 
      WRITE(STDOUT,'(7X,A,1X,I5)') 'mgLevels_Y = ',mgLevels_Y 
      IF (ndim==DIM_3D) WRITE(STDOUT,'(7X,A,1X,I5)') 'mgLevels_Z = ',mgLevels_Z 
      WRITE(STDOUT,*)
    end if ! monitorON

!   Allocate MG  grid-based arrays (phase I)
!   ----------------------------------------
    if (monitorON) WRITE(STDOUT,'(7X,A)') 'Allocating MG Memory, Phase I: MG Structures'
    if (monitorON) WRITE(STDOUT,*)
    CALL MG_Allocate_Memory(1)
!   ----------------------------------------------------------
!    Calculate the grid numbers and mesh sizes at each level.
!   ----------------------------------------------------------

!   X direction
!   -----------
    ilevel=1
    MGX(ilevel)%nxc_GLBL = nxc_GLBL
    MGX(ilevel)%nx_GLBL  = nx_GLBL

    DO ilevel = 2, mgLevels_X
      MGX(ilevel)%nxc_GLBL = (MGX(ilevel-1)%nxc_GLBL+1)/2
      MGX(ilevel)%nx_GLBL  =  MGX(ilevel  )%nxc_GLBL+1
    END DO ! ilevel

!   Y direction
!   -----------
    ilevel=1
    MGY(ilevel)%nyc_GLBL = nyc_GLBL
    MGY(ilevel)%ny_GLBL  = ny_GLBL

    DO ilevel = 2, mgLevels_Y
      MGY(ilevel)%nyc_GLBL = (MGY(ilevel-1)%nyc_GLBL+1)/2 
      MGY(ilevel)%ny_GLBL  =  MGY(ilevel  )%nyc_GLBL+1
    END DO ! n

!   Z direction
!   -----------
    ilevel=1
    MGZ(ilevel)%nzc = nzc
    MGZ(ilevel)%nz  = nz

    DO ilevel = 2, mgLevels_Z
      MGZ(ilevel)%nzc = (MGZ(ilevel-1)%nzc+1)/2 
      MGZ(ilevel)%nz  =  MGZ(ilevel  )%nzc+1
    END DO ! ilevel

!   Print output to screen
!   ----------------------
    if (monitorON) then
      WRITE(STDOUT,'(7X,A)') 'Global MG domain information:'
      WRITE(STDOUT,'(7X,A)') '============================='
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(7X,A)') 'nxc_MG:'
      WRITE(STDOUT,'(7X,A)') 'ilevel   nxc_MG   nx_MG'
      WRITE(STDOUT,'(7X,A)') '------   ------   -----'
      DO ilevel = 1, mgLevels_X
        WRITE(STDOUT,'(7X,I6,3X,I6,3X,I5)') ilevel, MGX(ilevel)%nxc_GLBL, MGX(ilevel)%nx_GLBL
      END DO ! ilevel
      WRITE(STDOUT,*)
    
      WRITE(STDOUT,'(7X,A)') 'nyc_MG:'
      WRITE(STDOUT,'(7X,A)') 'ilevel   nyc_MG   ny_MG'
      WRITE(STDOUT,'(7X,A)') '------   ------   -----'
      DO ilevel = 1, mgLevels_Y
        WRITE(STDOUT,'(7X,I6,3X,I6,3X,I5)') ilevel, MGY(ilevel)%nyc_GLBL, MGY(ilevel)%ny_GLBL
      END DO ! ilevel
      WRITE(STDOUT,*)

      IF (ndim == DIM_3D) THEN
        WRITE(STDOUT,'(7X,A)') 'nzc_mg '
        WRITE(STDOUT,'(7X,A)') 'ilevel   nzc_MG   nz_MG'
        WRITE(STDOUT,'(7X,A)') '------   ------   -----'
        DO ilevel = 1, mgLevels_Z
          WRITE(STDOUT,'(7X,I6,3X,I6,3X,I5)') ilevel, MGZ(ilevel)%nzc, MGZ(ilevel)%nz
        END DO ! ilevel
        WRITE(STDOUT,*)
      END IF  ! ndim   
    end if  ! monitorON

!   ---------------------------------------------------------------------------
!    Domain decomposition at each level (Building up from the coarsest level).
!   ---------------------------------------------------------------------------

!   X direction
!   -----------
    WRITE(ifuParLog,'(A)') "Details of MG domain decomposition, X direction"
    WRITE(ifuParLog,'(A)') "==============================================="
    WRITE(ifuParLog,'(A)') "ilevel    nxc     Is      Ie  "
    WRITE(ifuParLog,'(A)') "------  ------  ------  ------"

    ilevel=mgLevels_X

#   ifdef MPI
      call par_decomp(MGX(ilevel)%nxc_GLBL,ICOORD,MGX(ilevel)%nxc,MGX(ilevel)%nx,MGX(ilevel)%myIs,MGX(ilevel)%myIe)
#   else
      MGX(ilevel)%nxc=MGX(ilevel)%nxc_GLBL
      MGX(ilevel)%nx =MGX(ilevel)%nx_GLBL

      MGX(ilevel)%myIs=1
      MGX(ilevel)%myIe=MGX(ilevel)%nxc_GLBL
#   endif
    WRITE(ifuParLog,'(4(I6,2X))') ilevel, MGX(ilevel)%nxc, MGX(ilevel)%myIs, MGX(ilevel)%myIe

    DO ilevel=mgLevels_X-1,1,-1
      MGX(ilevel)%myIs=2*MGX(ilevel+1)%myIs-1
      MGX(ilevel)%myIe=MIN(2*MGX(ilevel+1)%myIe,MGX(ilevel)%nxc_GLBL)
      
      MGX(ilevel)%nxc=MGX(ilevel)%myIe-MGX(ilevel)%myIs+1
      MGX(ilevel)%nx =MGX(ilevel)%nxc+1

      WRITE(ifuParLog,'(4(I6,2X))') ilevel, MGX(ilevel)%nxc, MGX(ilevel)%myIs, MGX(ilevel)%myIe 
    END DO
    WRITE(ifuParLog,*)

!   Y direction
!   -----------
    WRITE(ifuParLog,'(A)') "Details of MG domain decomposition, Y direction"
    WRITE(ifuParLog,'(A)') "==============================================="
    WRITE(ifuParLog,'(A)') "ilevel    nyc     Js      Je  "
    WRITE(ifuParLog,'(A)') "------  ------  ------  ------"

    ilevel=mgLevels_Y

#   ifdef MPI
      call par_decomp(MGY(ilevel)%nyc_GLBL,JCOORD,MGY(ilevel)%nyc,MGY(ilevel)%ny,MGY(ilevel)%myJs,MGY(ilevel)%myJe)
#   else
      MGY(ilevel)%nyc=MGY(ilevel)%nyc_GLBL
      MGY(ilevel)%ny =MGY(ilevel)%ny_GLBL

      MGY(ilevel)%myJs=1
      MGY(ilevel)%myJe=MGY(ilevel)%nyc_GLBL
#   endif
    WRITE(ifuParLog,'(4(I6,2X))') ilevel, MGY(ilevel)%nyc, MGY(ilevel)%myJs, MGY(ilevel)%myJe

    DO ilevel=mgLevels_Y-1,1,-1
      MGY(ilevel)%myJs=2*MGY(ilevel+1)%myJs-1
      MGY(ilevel)%myJe=MIN(2*MGY(ilevel+1)%myJe,MGY(ilevel)%nyc_GLBL)
      
      MGY(ilevel)%nyc=MGY(ilevel)%myJe-MGY(ilevel)%myJs+1
      MGY(ilevel)%ny =MGY(ilevel)%nyc+1

      WRITE(ifuParLog,'(4(I6,2X))') ilevel, MGY(ilevel)%nyc, MGY(ilevel)%myJs, MGY(ilevel)%myJe 
    END DO
    WRITE(ifuParLog,*)

!   Allocate MG  grid-based arrays (phase II)
!   ----------------------------------------
    if (monitorON) WRITE(STDOUT,'(5X,A)') 'Allocating MG Memory, Phase II: MG Variables'
    if (monitorON) WRITE(STDOUT,*)
    CALL MG_Allocate_Memory(2)

END SUBROUTINE MG_initial_domain
!---------------------------------------------------------------------



SUBROUTINE MG_initial_grid()

! ----------------------------------------------
!  Purpose: Calculate mesh sizes at each level.
! ----------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

    INTEGER :: i, j, k, ii, jj, kk, ilevel



!   -------------------------------------------------------------------
!    Calculate GLOBAL domain variables at each level.
!    Cannot emphasize enough on word "GLOBAL"!
!    Compute xc, yc, zc, dxcinv, dycinv, dzcinv etc for MG structures.
!   -------------------------------------------------------------------

!   Level 1
!   -------
    MGX(1)%x(0:nx_GLBL+1)       =      x(0:nx_GLBL +1)
    MGX(1)%xc(0:nxc_GLBL+1)     =     xc(0:nxc_GLBL+1)
    MGX(1)%dxinv(0:nxc_GLBL+1)  =  dxinv(0:nxc_GLBL+1)
    MGX(1)%dxcinv(0:nxc_GLBL+1) = dxcinv(0:nxc_GLBL+1)

    MGY(1)%y(0:ny_GLBL+1)       =      y(0:ny_GLBL +1) 
    MGY(1)%yc(0:nyc_GLBL+1)     =     yc(0:nyc_GLBL+1) 
    MGY(1)%dyinv(0:nyc_GLBL+1)  =  dyinv(0:nyc_GLBL+1)
    MGY(1)%dycinv(0:nyc_GLBL+1) = dycinv(0:nyc_GLBL+1)

    MGZ(1)%z(0:nz+1)            =      z(0:nz +1)
    MGZ(1)%zc(0:nzc+1)          =     zc(0:nzc+1) 
    MGZ(1)%dzinv(0:nzc+1)       =  dzinv(0:nzc+1)
    MGZ(1)%dzcinv(0:nzc+1)      = dzcinv(0:nzc+1)

!   Remaining Levels
!   ----------------

!   a. computing the node coordinates
!   ---------------------------------

!   x direction
!   -----------
    DO ilevel = 2, mgLevels_X

      MGX(ilevel)%x(0)=x(0)
      DO i=1,MGX(ilevel)%nx_GLBL-1
        ii = 2*i-1
        MGX(ilevel)%x(i) = MGX(ilevel-1)%x(ii)
      ENDDO ! i
      MGX(ilevel)%x(MGX(ilevel)%nx_GLBL  )=x(nx_GLBL  )
      MGX(ilevel)%x(MGX(ilevel)%nx_GLBL+1)=x(nx_GLBL+1)

      DO i=1,MGX(ilevel)%nxc_GLBL
        MGX(ilevel)%xc(i)    = 0.5_CGREAL*(MGX(ilevel)%x(i+1)+MGX(ilevel)%x(i))
      ENDDO ! i
    !This implementation has been commented out to match serial code 14.4.3 
     ! MGX(ilevel)%xc(0) = 2.0_CGREAL * MGX(ilevel)%x(1) - MGX(ilevel)%xc(2)
     ! MGX(ilevel)%xc(MGX(ilevel)%nx_GLBL) = 2.0_CGREAL * MGX(ilevel)%x(MGX(ilevel)%nx_GLBL) - MGX(ilevel)%xc(MGX(ilevel)%nxc_GLBL)
       MGX(ilevel)%xc(0) = MGX(ilevel-1)%xc(0) 
       MGX(ilevel)%xc(MGX(ilevel)%nx_GLBL) = MGX(ilevel-1)%xc(MGX(ilevel-1)%nx_GLBL)
 
      DO i=0,MGX(ilevel)%nxc_GLBL+1
        MGX(ilevel)%dxinv(i) = 1.0_CGREAL/(MGX(ilevel)%x(i+1)-MGX(ilevel)%x(i))
      ENDDO ! i


      DO i= 1,MGX(ilevel)%nxc_GLBL+1
        MGX(ilevel)%dxcinv(i) = 1.0_CGREAL/(MGX(ilevel)%xc(i)-MGX(ilevel)%xc(i-1))
      ENDDO ! i

    END DO ! ilevel

!   y direction
!   -----------
    DO ilevel = 2, mgLevels_Y

      MGY(ilevel)%y(0)=y(0)
      DO j=1,MGY(ilevel)%ny_GLBL-1
        jj = 2*j-1
        MGY(ilevel)%y(j) = MGY(ilevel-1)%y(jj)
      ENDDO ! j
      MGY(ilevel)%y(MGY(ilevel)%ny_GLBL  )=y(ny_GLBL)
      MGY(ilevel)%y(MGY(ilevel)%ny_GLBL+1)=y(ny_GLBL+1)

      DO j=1,MGY(ilevel)%nyc_GLBL
        MGY(ilevel)%yc(j)    = 0.5_CGREAL*(MGY(ilevel)%y(j+1)+MGY(ilevel)%y(j))
      ENDDO ! j
      !MGY(ilevel)%yc(0) = 2.0_CGREAL * MGY(ilevel)%y(1) - MGY(ilevel)%yc(2)
      !MGY(ilevel)%yc(MGY(ilevel)%ny_GLBL) = 2.0_CGREAL * MGY(ilevel)%y(MGY(ilevel)%ny_GLBL) - MGY(ilevel)%yc(MGY(ilevel)%nyc_GLBL)
       MGY(ilevel)%yc(0) = MGY(ilevel-1)%yc(0)    
       MGY(ilevel)%yc(MGY(ilevel)%ny_GLBL) = MGY(ilevel-1)%yc(MGY(ilevel-1)%ny_GLBL)  

      DO j=0,MGY(ilevel)%nyc_GLBL+1
        MGY(ilevel)%dyinv(j) = 1.0_CGREAL/(MGY(ilevel)%y(j+1)-MGY(ilevel)%y(j))
      ENDDO ! j


      DO j= 1,MGY(ilevel)%nyc_GLBL+1
        MGY(ilevel)%dycinv(j) = 1.0_CGREAL/(MGY(ilevel)%yc(j)-MGY(ilevel)%yc(j-1))
      ENDDO ! j

    END DO ! ilevel

!   z direction
!   -----------
    DO ilevel = 2, mgLevels_Z

      MGZ(ilevel)%z(0)=z(0)
      DO k=1,MGZ(ilevel)%nz-1
        kk = 2*k-1
        MGZ(ilevel)%z(k) = MGZ(ilevel-1)%z(kk)
      ENDDO ! k
      MGZ(ilevel)%z(MGZ(ilevel)%nz  )=z(nz)
      MGZ(ilevel)%z(MGZ(ilevel)%nz+1)=z(nz+1)

      DO k=1,MGZ(ilevel)%nzc
        MGZ(ilevel)%zc(k)    = 0.5_CGREAL*(MGZ(ilevel)%z(k+1)+MGZ(ilevel)%z(k))
      ENDDO
      !MGZ(ilevel)%zc(0) = 2.0_CGREAL * MGZ(ilevel)%z(1) - MGZ(ilevel)%zc(2)
      !MGZ(ilevel)%zc(MGZ(ilevel)%nzc+1) = 2.0_CGREAL * MGZ(ilevel)%z(MGZ(ilevel)%nz) - MGZ(ilevel)%zc(MGZ(ilevel)%nzc) 
      MGZ(ilevel)%zc(0)=MGZ(ilevel-1)%zc(0)
      MGZ(ilevel)%zc(MGZ(ilevel)%nzc+1) = MGZ(ilevel-1)%zc(MGZ(ilevel-1)%nzc+1)

      DO k=0,MGZ(ilevel)%nzc+1
        MGZ(ilevel)%dzinv(k) = 1.0_CGREAL/(MGZ(ilevel)%z(k+1)-MGZ(ilevel)%z(k))
      ENDDO


      DO k= 1,MGZ(ilevel)%nzc+1
        MGZ(ilevel)%dzcinv(k) = 1.0_CGREAL/(MGZ(ilevel)%zc(k)-MGZ(ilevel)%zc(k-1))
      ENDDO

    END DO ! ilevel

!   Allocate MG  grid-based arrays (phase II)
!   ----------------------------------------
    if (monitorON) WRITE(STDOUT,'(5X,A)') 'Allocating MG Memory, Phase III: MG Arrays'
    CALL MG_Allocate_Memory(3)

#   ifdef MPI
!     Generate MG communication pattern
!     -------------------------------
      if (monitorON) WRITE(STDOUT,'(5X,A)') 'Generating MG communication pattern'
      CALL par_init_MG_domain
#   endif
    if (monitorON) WRITE(STDOUT,*)

END SUBROUTINE MG_initial_grid
!---------------------------------------------------------------------



SUBROUTINE MG_Allocate_Memory(phase1)

! ---------------------------------------------------------------------------
!  This subroutine allocates MG arrays relevant to the coarser level meshes.
! ---------------------------------------------------------------------------
   
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: phase1

!   Local variables
!   ---------------
    INTEGER :: ilevel, iErr



    SELECT CASE(phase1)
    
    CASE(1)

!     ------------------------
!      Allocate MG structures
!     ------------------------
      ALLOCATE(MGX(mgLevels_X), STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for MGX'
        STOP
      ENDIF ! iErr

      ALLOCATE(MGY(mgLevels_Y), STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for MGY'
        STOP
      ENDIF ! iErr

      ALLOCATE(MGZ(mgLevels_Z), STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for MGZ'
        STOP
      ENDIF ! iErr

    CASE(2)

!     -----------------------
!      Allocate MG variables
!     -----------------------
      DO ilevel=1,mgLevels_X
        ALLOCATE(MGX(ilevel)%x(0:MGX(ilevel)%nx_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for x_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGX(ilevel)%xc(0:MGX(ilevel)%nxc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for xc_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGX(ilevel)%dxinv(0:MGX(ilevel)%nxc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dxinv_MG'
          STOP
        ENDIF ! iErr 

        ALLOCATE(MGX(ilevel)%dxcinv(0:MGX(ilevel)%nxc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dxcinv_MG',ilevel,MGX(ilevel)%nxc_GLBL
          STOP
        ENDIF ! iErr  
      END DO

!     --------------------------------------------------------------
      DO ilevel=1,mgLevels_Y
        ALLOCATE(MGY(ilevel)%y(0:MGY(ilevel)%ny_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for y_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGY(ilevel)%yc(0:MGY(ilevel)%nyc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for yc_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGY(ilevel)%dyinv(0:MGY(ilevel)%nyc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dyinv_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGY(ilevel)%dycinv(0:MGY(ilevel)%nyc_GLBL+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dycinv_MG'
          STOP
        ENDIF ! iErr 
      END DO

!     --------------------------------------------------------------
      DO ilevel=1,mgLevels_Z
        ALLOCATE(MGZ(ilevel)%z(0:MGZ(ilevel)%nz+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for z_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGZ(ilevel)%zc(0:MGZ(ilevel)%nzc+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for zc_MG'
          STOP
        ENDIF ! iErr

        ALLOCATE(MGZ(ilevel)%dzinv(0:MGZ(ilevel)%nzc+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dzinv_MG'
          STOP
        ENDIF ! iErr 

        ALLOCATE(MGZ(ilevel)%dzcinv(0:MGZ(ilevel)%nzc+1),STAT=iErr)
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Allocate_Memory: Memory Allocation Error for dzcinv_MG'
          STOP
        ENDIF ! iErr 
      END DO

!     --------------------------------------------------------------
      ALLOCATE(ium_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for ium_MG'
        STOP
      ENDIF ! iErr 

      ALLOCATE(iup_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for iup_MG'
        STOP
      ENDIF ! iErr

      ALLOCATE(jum_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for jum_MG'
        STOP
      ENDIF ! iErr

      ALLOCATE(jup_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for jup_MG'
        STOP
      ENDIF ! iErr

      ALLOCATE(kum_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for kum_MG'
        STOP
      ENDIF ! iErr

      ALLOCATE(kup_MG(MGX(1)%nxc,MGY(1)%nyc,MGZ(1)%nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for kup_MG'
        STOP
      ENDIF ! iErr

      ALLOCATE(iblank_MG(1-Ngl:MGX(1)%nxc+Ngl,1-Ngl:MGY(1)%nyc+Ngl,0:MGZ(1)%nzc+1),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'MG_Allocate_Memory: Memory Allocation Error for iblank_MG'
        STOP
      ENDIF ! iErr
    
    CASE(3)

!     --------------------
!      Allocate MG arrays
!     --------------------
      ilevel=1

      ALLOCATE(MGX(ilevel)%rhs(MGX(ilevel)%nxc,nyc,nzc),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGX(ilevel)%res(MGX(ilevel)%nxc,nyc,nzc),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGX(ilevel)%phi(1-Ngl:MGX(ilevel)%nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGX(ilevel)%iblank(1-Ngl:MGX(ilevel)%nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
        STOP
      ENDIF ! ierr

      DO ilevel=2,mgLevels_X
        ALLOCATE(MGX(ilevel)%rhs(MGX(ilevel)%nxc,nyc,nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGX(ilevel)%res(MGX(ilevel)%nxc,nyc,nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGX(ilevel)%phi(0:MGX(ilevel)%nxc+1,0:nyc+1,0:nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGX(ilevel)%iblank(0:MGX(ilevel)%nxc+1,0:nyc+1,0:nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr
      END DO

!     --------------------------------------------------------------
      ilevel=1

      ALLOCATE(MGY(ilevel)%rhs(nxc,MGY(ilevel)%nyc,nzc),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGY(ilevel)%res(nxc,MGY(ilevel)%nyc,nzc),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGY(ilevel)%phi(1-Ngl:nxc+Ngl,1-Ngl:MGY(ilevel)%nyc+Ngl,0:nzc+1),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
        STOP
      ENDIF ! ierr

      ALLOCATE(MGY(ilevel)%iblank(1-Ngl:nxc+Ngl,1-Ngl:MGY(ilevel)%nyc+Ngl,0:nzc+1),STAT=ierr )
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
        STOP
      ENDIF ! ierr

      DO ilevel=2,mgLevels_Y
        ALLOCATE(MGY(ilevel)%rhs(nxc,MGY(ilevel)%nyc,nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
          'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGY(ilevel)%res(nxc,MGY(ilevel)%nyc,nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGY(ilevel)%phi(0:nxc+1,0:MGY(ilevel)%nyc+1,0:nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGY(ilevel)%iblank(0:nxc+1,0:MGY(ilevel)%nyc+1,0:nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr
      END DO

!     --------------------------------------------------------------
      IF (ndim==DIM_3D) THEN
        ilevel=1
        ALLOCATE(MGZ(ilevel)%rhs(nxc,nyc,MGZ(ilevel)%nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGZ(ilevel)%res(nxc,nyc,MGZ(ilevel)%nzc),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGZ(ilevel)%phi(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:MGZ(ilevel)%nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr

        ALLOCATE(MGZ(ilevel)%iblank(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:MGZ(ilevel)%nzc+1),STAT=ierr )
        IF ( iErr /= ERR_NONE ) THEN
          WRITE(STDOUT,*) &
            'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
          STOP
        ENDIF ! ierr

        DO ilevel=2,mgLevels_Z
          ALLOCATE(MGZ(ilevel)%rhs(nxc,nyc,MGZ(ilevel)%nzc),STAT=ierr )
          IF ( iErr /= ERR_NONE ) THEN
            WRITE(STDOUT,*) &
              'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
            STOP
          ENDIF ! ierr

          ALLOCATE(MGZ(ilevel)%res(nxc,nyc,MGZ(ilevel)%nzc),STAT=ierr )
          IF ( iErr /= ERR_NONE ) THEN
            WRITE(STDOUT,*) &
              'MG_Memory_Allocation: Memory Allocation Error for rhs_MG'
            STOP
          ENDIF ! ierr

          ALLOCATE(MGZ(ilevel)%phi(0:nxc+1,0:nyc+1,0:MGZ(ilevel)%nzc+1),STAT=ierr )
          IF ( iErr /= ERR_NONE ) THEN
            WRITE(STDOUT,*) &
              'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
            STOP
          ENDIF ! ierr

          ALLOCATE(MGZ(ilevel)%iblank(0:nxc+1,0:nyc+1,0:MGZ(ilevel)%nzc+1),STAT=ierr )
          IF ( iErr /= ERR_NONE ) THEN
            WRITE(STDOUT,*) &
              'MG_Memory_Allocation: Memory Allocation Error for phi_MG'
            STOP
          ENDIF ! ierr
        END DO
      END IF

    END SELECT

END SUBROUTINE MG_Allocate_Memory
!---------------------------------------------------------------------



SUBROUTINE MG_Solver(var,rr,compTime,commTime)

! --------------------------------------------------------------
!  Purpose: Main driver for MG code
!           calling semi-coarsening MG method for each direction
!
!  Input: 
!         var:  pressure values
!         rr:   values at the right-hand side
!  
!  Output: 
!         var: storing the approximation of the solution.
! --------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
#   ifdef MPI
      USE mpi
#   endif

    IMPLICIT NONE

!   Parameters
!   ---------- 
    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(kind=CGREAL),INTENT (OUT)   :: compTime  , commTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime
    REAL(kind=CGREAL) :: totalTime  

    INTEGER :: clock1, clock2, clock_rate
    INTEGER :: ierr, i,j,k
!DEBUG
!    CHARACTER*20              :: fname1



!    WRITE(fname1,"('q.',I3.3,'.',I7.7)") myRank,ntime

!    OPEN(UNIT=70,FILE=fname1,FORM='UNFORMATTED')
!DEBUG

 
#   ifdef MPI
      startTime=MPI_WTIME()
#   else
      CALL system_clock(clock1)
#   endif

    commTime=0.d0

!   MG in X direction
!   -----------------
    CALL MG_solv_x(var,rr,commTime)!......................................................ONLY_PER(SAMK)
 
!   MG in Y direction
!   -----------------
   CALL MG_solv_y(var,rr,commTime)!......................................................ONLY_PER(SAMK)

!   MG in Z direction
!   -----------------
    IF ( nDim==DIM_3D ) THEN
      CALL MG_solv_z(var,rr,commTime)!....................................................ONLY_PER(SAMK)
    END IF ! nDim

#   ifdef MPI
      endTime=MPI_WTIME()
      totalTime=endTime-startTime
#   else
      CALL system_clock(clock2, clock_rate)
      totalTime=REAL(clock2-clock1)/REAL(clock_rate)
#   endif
    compTime=totalTime-commTime

END SUBROUTINE MG_Solver
!------------------------------------------------------------------------------



SUBROUTINE MG_Solv_x(var,rr,commTime)

! ---------------------------------------------------------------
!  Purpose: apply semi-coarsening MG method for x direction
!
!  Input:
!         var: input pressure values
!         rr:  values at the right-hand side
!
!  Output: 
!         var: storing the approximation of the solution.
!
!
!  V cycle:                                 F cycle  W cycle:
!  O           O     O                         *    O
!   \         /       \                       /    /
!    O       O         O             O       *    O
!     \     /           \           / \     /    / 
!      O   O             O    O    O   O   *    O
!       \ /               \  / \  /     \ / \  /
!        O                 O    O        O   O
! ---------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
 
    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    INTEGER, PARAMETER :: V_CYCLE=1, W_CYCLE=2, F_CYCLE=3



!   Allocate memory, prepare iblank
!   -------------------------------
    CALL MG_Nullify_Arrays(ICOORD)!..................................................COMPLETE(SAMK)
    CALL MG_Prepare_Iblank(rr, ICOORD)!..............................................COMPLETE(SAMK)

!   Coarsening is done in each direction separately.
!   Then LSOR method is invoked for that level.
!
!   mgCycle = 1-V-Cycle
!   mgCycle = 2-W cycle
!   mgCycle = 3-F cycle
!   -------------------
    SELECT CASE(mgCycleX)
    
    CASE(V_CYCLE) 
      CALL MG_FINE_TO_COARSE_X(var, rr, 1, mgLevels_X, commTime)!.........................ONLY-PER(SAMK)
!DEBUG
!  return
!DEBUG
      CALL MG_COARSE_TO_FINE_X(var, rr, 1, mgLevels_X, commTime)!.........................ONLY-PER(SAMK)

    CASE(W_CYCLE) 
      CALL MG_FINE_TO_COARSE_X(var, rr, 1           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK) 
      CALL MG_COARSE_TO_FINE_X(var, rr, 2           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_X(var, rr, 2           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_X(var, rr, 1           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      
    CASE(F_CYCLE) 
      CALL MG_FINE_TO_COARSE_X(var, rr, 1           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_X(var, rr, mgLevels_X-1, mgLevels_X, commTime)!..............ONLY-PER(SAMK) 
      CALL MG_COARSE_TO_FINE_X(var, rr, 2           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_X(var, rr, 2           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_X(var, rr, 1           , mgLevels_X, commTime)!..............ONLY-PER(SAMK)

    END SELECT ! mgCycleX

END SUBROUTINE MG_Solv_x
!---------------------------------------------------------------------



SUBROUTINE MG_Solv_y(var,rr,commTime)

! ---------------------------------------------------------------
!  Purpose: apply semi-coarsening MG method for y direction
!
!  Input:
!         var: input pressure values
!         rr:  values at the right-hand side
!
!  Output: 
!         var: storing the approximation of the solution.
!
!
!  V cycle:                                 F cycle  W cycle:
!  O           O     O                         *    O
!   \         /       \                       /    /
!    O       O         O             O       *    O
!     \     /           \           / \     /    / 
!      O   O             O    O    O   O   *    O
!       \ /               \  / \  /     \ / \  /
!        O                 O    O        O   O
! ---------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
 
    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    INTEGER, PARAMETER :: V_CYCLE=1, W_CYCLE=2, F_CYCLE=3



!   Allocate memory, prepare iblank
!   -------------------------------
    CALL MG_Nullify_Arrays(JCOORD)
    CALL MG_Prepare_Iblank(rr, JCOORD)

!   Coarsening is done in each direction separately.
!   Then LSOR method is invoked for that level.
!
!   mgCycle = 1-V-Cycle
!   mgCycle = 2-W cycle
!   mgCycle = 3-F cycle
!   -------------------
    SELECT CASE(mgCycleY)
    
    CASE(V_CYCLE)
      CALL MG_FINE_TO_COARSE_Y(var, rr, 1, mgLevels_Y, commTime)!.........................ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, 1, mgLevels_Y, commTime)!.........................ONLY-PER(SAMK)
      
    CASE(W_CYCLE)    
      CALL MG_FINE_TO_COARSE_Y(var, rr, 1           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, 2           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_Y(var, rr, 2           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, 1           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)

    CASE(F_CYCLE)    
      CALL MG_FINE_TO_COARSE_Y(var, rr, 1           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_Y(var, rr, mgLevels_Y-1, mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, 2           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_FINE_TO_COARSE_Y(var, rr, 2           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)
      CALL MG_COARSE_TO_FINE_Y(var, rr, 1           , mgLevels_Y, commTime)!..............ONLY-PER(SAMK)

    END SELECT ! mgCycleY

END SUBROUTINE MG_Solv_y
!---------------------------------------------------------------------



SUBROUTINE MG_Solv_z(var,rr,commTime)

! ---------------------------------------------------------------
!  Purpose: apply semi-coarsening MG method for z direction
!
!  Input:
!         var: input pressure values
!         rr:  values at the right-hand side
!
!  Output: 
!         var: storing the approximation of the solution.
!
!
!  V cycle:                                 F cycle  W cycle:
!  O           O     O                         *    O
!   \         /       \                       /    /
!    O       O         O             O       *    O
!     \     /           \           / \     /    / 
!      O   O             O    O    O   O   *    O
!       \ /               \  / \  /     \ / \  /
!        O                 O    O        O   O
! ---------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
 
    IMPLICIT NONE

!   Parameters

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables

    INTEGER, PARAMETER :: V_CYCLE=1, W_CYCLE=2, F_CYCLE=3



!   Allocate memory, prepare iblank
!   -------------------------------
    CALL MG_Nullify_Arrays(KCOORD)
    CALL MG_Prepare_Iblank(rr, KCOORD)

!   Coarsening is done in each direction separately.
!   Then LSOR method is invoked for that level.
!
!   mgCycle = 1-V-Cycle
!   mgCycle = 2-W cycle
!   mgCycle = 3-F cycle
!   -------------------

    SELECT CASE(mgCycleZ)

    CASE(V_CYCLE)
      CALL MG_FINE_TO_COARSE_Z(var, rr, 1, mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, 1, mgLevels_Z, commTime)
    
    CASE(W_CYCLE) 
      CALL MG_FINE_TO_COARSE_Z(var, rr, 1           , mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_FINE_TO_COARSE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, 2           , mgLevels_Z, commTime)
      CALL MG_FINE_TO_COARSE_Z(var, rr, 2           , mgLevels_Z, commTime) 
      CALL MG_COARSE_TO_FINE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_FINE_TO_COARSE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, 1           , mgLevels_Z, commTime)

    CASE(F_CYCLE) 
      CALL MG_FINE_TO_COARSE_Z(var, rr, 1           , mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_FINE_TO_COARSE_Z(var, rr, mgLevels_Z-1, mgLevels_Z, commTime)
      CALL MG_COARSE_TO_FINE_Z(var, rr, 2           , mgLevels_Z, commTime)
      CALL MG_FINE_TO_COARSE_Z(var, rr, 2           , mgLevels_Z, commTime) 
      CALL MG_COARSE_TO_FINE_Z(var, rr, 1           , mgLevels_Z, commTime)

    END SELECT ! mgCycleX

END SUBROUTINE MG_Solv_z
!---------------------------------------------------------------------



SUBROUTINE MG_Nullify_Arrays(dir)

! -----------------------------------------------------------------------------
!  Purpose: Allocate memory specifically for  MG_solv_x, MG_solv_y, MG_solv_z 
! -----------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER,INTENT (IN)  :: dir

!   Local variables
!   ---------------
    INTEGER :: ilevel



!   Nullify MG arrays
!   -----------------
    SELECT CASE(dir)
    CASE(ICOORD)
      DO ilevel=1,mgLevels_X
        MGX(ilevel)%rhs=0.0_CGREAL
        MGX(ilevel)%phi=0.0_CGREAL

        MGX(ilevel)%iblank=0
      END DO

    CASE(JCOORD)
      DO ilevel=1,mgLevels_Y
        MGY(ilevel)%rhs=0.0_CGREAL
        MGY(ilevel)%phi=0.0_CGREAL

        MGY(ilevel)%iblank=0
      END DO

    CASE(KCOORD)
      DO ilevel=1,mgLevels_Z
        MGZ(ilevel)%rhs=0.0_CGREAL
        MGZ(ilevel)%phi=0.0_CGREAL

        MGZ(ilevel)%iblank=0
      END DO
    END SELECT

END SUBROUTINE MG_Nullify_Arrays
!---------------------------------------------------------------------



SUBROUTINE MG_Prepare_Iblank(rr, NC_dim)

!------------------------------------------------------------
!  Purpose: Compute rhs_MG values for each level
!
!  Input: 
!         rr:   values at the right-hand side
!         NC_dim: x, y, z direction
!  
!  Output: 
!         rhs_MG
!------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE MG_parameters
    USE MG_arrays
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN)           :: NC_dim
    REAL(KIND=CGREAL),INTENT (IN) :: rr(nxc,nyc,nzc)

!   Local variables
!   ---------------
    INTEGER :: i, j, k
    INTEGER :: ilevel
    INTEGER :: ii, iim, jj, jjm, kk, kkm
    
    REAL(KIND=CGREAL) :: myTime



    SELECT CASE (NC_dim)

!   x-direction
!   -----------
    CASE(ICOORD)
      MGX(1)%rhs(:,:,:) = rr(:,:,:)
      MGX(1)%iblank(:,:,:) = iblank(:,:,:)

      DO ilevel = 2, mgLevels_X  
        DO k = 1, nzc
        DO j = 0, nyc+1
        DO i = 1, MGX(ilevel)%nxc
          ii  = MIN(2*i,MGX(ilevel-1)%nxc)
          iim = 2*i-1

          MGX(ilevel)%iblank(i,j,k) = MGX(ilevel-1)%iblank(iim,j,k)*MGX(ilevel-1)%iblank(ii,j,k) 
        END DO ! i
        END DO ! j
        END DO ! k
        
      END DO

#     ifdef MPI
        CALL par_comm_iblankX(mgLevels_X,myTime)
#     endif

!   y-direction
!   -----------
    CASE(JCOORD)
      MGY(1)%rhs(:,:,:) = rr(:,:,:) 
      MGY(1)%iblank(:,:,:) = iblank(:,:,:)

      DO ilevel = 2, mgLevels_Y
        DO k = 1, nzc
        DO i = 0, nxc+1
        DO j = 1, MGY(ilevel)%nyc
          jj  = MIN(2*j,MGY(ilevel-1)%nyc)
          jjm = 2*j-1

          MGY(ilevel)%iblank(i,j,k) = MGY(ilevel-1)%iblank(i,jjm,k)*MGY(ilevel-1)%iblank(i,jj,k)
        END DO ! j
        END DO ! i
        END DO ! k
      END DO ! n

#     ifdef MPI
        CALL par_comm_iblankY(mgLevels_Y,myTime)
#     endif

!   z-direction
!   -----------
    CASE(KCOORD)
      MGZ(1)%rhs(:,:,:) = rr(:,:,:) 
      MGZ(1)%iblank(:,:,:) = iblank(:,:,:)

      DO ilevel = 2, mgLevels_Z
        DO j = 0, nyc+1
        DO i = 0, nxc+1
        DO k = 1, MGZ(ilevel)%nzc
          kk  = MIN(2*k,MGZ(ilevel-1)%nzc)
          kkm = 2*k-1

          MGZ(ilevel)%iblank(i,j,k) = MGZ(ilevel-1)%iblank(i,j,kkm)*MGZ(ilevel-1)%iblank(i,j,kk)
        END DO ! k
        END DO ! j
        END DO ! i
      END DO ! n
    END SELECT ! NC_dim 

END SUBROUTINE MG_Prepare_Iblank
!---------------------------------------------------------------------



SUBROUTINE MG_FINE_TO_COARSE_X(var, rr, levBegin, levEnd, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply coarsening on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levBegin: beginning level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the x-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
    USE grid_arrays
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levBegin, levEnd

    REAL(KIND=CGREAL),INTENT (IN)    :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    logical :: resCheckFlag

    INTEGER :: iIter
    INTEGER :: ilevel, i, j, k, totIter
    INTEGER :: MGgl
    
    REAL(KIND=CGREAL) :: myTime

    REAL(KIND=CGREAL),DIMENSION(nxc,nyc,nzc) :: r, r1



!   Initialize variable
!   -------------------
    resCheckFlag = .true.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin, levEnd
      IF (ilevel == 1) THEN

!       ilevel Parameters
!       -----------------
        totIter = iterFinest
        MGgl=Ngl

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare_BC(var)!.....................................................COMPLETE(SAMK)

!       Prepare MG Arrays
!       -----------------
        MGX(1)%phi(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)=var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
        MGX(1)%rhs(1:nxc,1:nyc,1:nzc) = rr(1:nxc,1:nyc,1:nzc)

!       Compute ghost cells and impose compatibility for _first_level with GCM
!       ----------------------------------------------------------------------
        IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)        
          SELECT CASE(POSSION_INDEX) 
          CASE(index_pressure)
            CALL set_outer_ghost_pres(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)    
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                    
          END SELECT
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation 

      ELSE ! ilevel > 1

!       ilevel Parameters
!       -----------------
        IF (ilevel == mgLevels_X) THEN
          totIter = iterCoarsest
        ELSE 
          totIter = 1
        END IF ! ilevel==mgLevels_X
        MGgl=1

!       Prepare MG Arrays
!       -----------------
        MGX(ilevel)%phi(:,:,:) = 0.0_CGREAL

        do k=1,            nzc
        do j=0,            nyc+1
        do i=0,MGX(ilevel)%nxc+1
          iblank_MG(i,j,k) = MGX(ilevel)%iblank(i,j,k)
        end do
        end do
        end do

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGX(ilevel)%phi,MGX(ilevel)%nxc,nyc,nzc)!.........................ONLY_PER(SAMK)

      END IF ! ilevel==1

!DEBUG
!  WRITE(70) nxc,myIs,myIe,nyc,MyJs,myJe,nzc
!  do k=1,nzc
!  do j=1,nyc
!  do i=1,nxc
!    WRITE(70) MGX(ilevel)%phi(i,j,k), MGX(ilevel)%rhs(i,j,k), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
!  end do
!  end do
!  end do
!  close(70)
!  return
!DEBUG

!     Invoke Line SOR
!     ---------------
      DO iIter = 1,totIter
        CALL MG_itsolv(MGX(ilevel)%phi,MGX(ilevel)%rhs,ilevel,1,1,MGX(ilevel)%nxc,nyc,nzc,MGgl,myTime)!.......ONLY_PER(SAMK)
      END DO ! iIter

!DEBUG
!  WRITE(70) nxc,myIs,myIe,nyc,MyJs,myJe,nzc
!  do k=1,nzc
!  do j=1,nyc
!  do i=1,nxc
!    WRITE(70) MGX(ilevel)%phi(i,j,k), MGX(ilevel)%rhs(i,j,k), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
!  end do
!  end do
!  end do
!  close(70)
!  return
!DEBUG

!     For all levels except the coarsest, compute residuals and apply restriction (ilevel ==> ilevel+1).
!     --------------------------------------------------------------------------------------------------
      IF ( ilevel /= levEnd ) THEN
        CALL MG_Residual(MGX(ilevel)%phi,MGX(ilevel)%rhs,MGX(ilevel)%res,ilevel,1,1,  &
                         MGX(ilevel)%nxc,nyc,nzc,MGgl,resCheckFlag)!......................ONLY_PER(SAMK)

        CALL MG_Restrict(ilevel,ICOORD)!.............................................COMPLETE(SAMK)

      END IF ! ilevel /= levEnd

    END DO  ! ilevel

!  open(123456, file = 'debugvar.dat')
!  write(123456, *) 'VARIABLES="X","Y","var"'
!  write(123456,*)'ZONE F=POINT, I=',(nx-1)/2,', J=',ny-1
!   k = 1
!   do j=1,ny-1
!   do i=1,(nx-1)/2
!   write(123456, *)MGX(2)%xc(i),yc(j),MGX(2)%phi(i,j,k)
!   enddo
!   enddo
!  close(123456)
!  CALL mpi_finalize(i)
!  stop


END SUBROUTINE MG_FINE_TO_COARSE_X
!---------------------------------------------------------------------



SUBROUTINE MG_COARSE_TO_FINE_X(var, rr, levEnd, levBegin, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply refinement on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levEnd: end level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the x-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     

   USE grid_arrays 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levEnd, levBegin

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    INTEGER :: ilevel, i, j, k, iIter, totIter, MGgl

    LOGICAL :: resCheckFlag

    REAL(KIND=CGREAL) :: myTime



!   Initialize variable
!   -------------------
    resCheckFlag = .false.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin-1, levEnd, -1 
    
      IF (ilevel==1) THEN
        MGgl=Ngl
      ELSE
        MGgl=1
      END IF

      IF (internal_boundary_present==INTR_BOUND_PRESENT) THEN
        DO k = 1,               nzc
        DO j = 0,               nyc+1
        DO i = 0, MGX(ilevel+1)%nxc+1
          iblank_MG(i,j,k) = MGX(ilevel+1)%iblank(i,j,k)
        END DO ! i
        END DO ! j
        END DO ! k
      END IF ! internal_boundary_present

!     Invoke prolongation (ilevel+1 ==> ilevel)
!     -----------------------------------------
      CALL MG_Prolong(ilevel,ICOORD)!................................................COMPLETE(SAMK)
!   open(123456, file = 'debugvar.dat')
!  write(123456, *) 'VARIABLES="X","Y","var"'
!  write(123456,*)'ZONE F=POINT, I=',(nx-1),', J=',ny-1
!   k = 1
!   do j=1,ny-1
!   do i=1,(nx-1)
!   write(123456, *)MGX(1)%xc(i),yc(j),MGX(1)%phi(i,j,k)
!   enddo
!   enddo
!  close(123456)
!  CALL mpi_finalize(i)
!  stop

#     ifdef MPI
        CALL par_comm_MGX(ilevel,myTime)
        commTime=commTime+myTime
#     endif

!     Load variables
!     --------------
      DO k = 1     ,             nzc
      DO j = 1-MGgl,             nyc+MGgl
      DO i = 1-MGgl, MGX(ilevel)%nxc+MGgl
        iblank_MG(i,j,k) = MGX(ilevel)%iblank(i,j,k)
      END DO ! i
      END DO ! j
      END DO ! k

      if (ilevel/=1) then

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGX(ilevel)%phi,MGX(ilevel)%nxc,nyc,nzc)!.........................ONLY-PER(SAMK)

!       Iteration count
!       ---------------
        totIter = iterInter

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGX(ilevel)%phi,MGX(ilevel)%rhs,ilevel,1,1,MGX(ilevel)%nxc,nyc,nzc,MGgl,myTime)!.....ONLY_PER(SAMK)
        END DO ! iIter
           
!       Compute residuals
!       -----------------
!!        IF (infoconv == 1) THEN
!!          CALL MG_Residual(var,r, n, 1, 1, mgrid_I(n), ny-1, nz-1, resCheckFlag)
!!        END IF ! infoConv
      ELSE  

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
       IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)       
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure) 
            CALL set_outer_ghost_pres(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)  
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                     
          END SELECT          
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation

         CALL MG_Prepare_BC(MGX(ilevel)%phi)!.........................................COMPLETE(SAMK)

!       Iteration count
!       ---------------
        totIter = iterFinest

!          open(123456, file = 'debugvar.dat')
!  write(123456, *) 'VARIABLES="X","Y","var","rsh"'
!  write(123456,*)'ZONE F=POINT, I=',(nx-1),', J=',ny-1
!   k = 1
!   do j=1,ny-1
!   do i=1,(nx-1)
!   write(123456, '(4(1x, e14.7))')MGX(1)%xc(i),yc(j),MGX(1)%phi(i,j,k), MGX(1)%rhs(i, j, k)
!   enddo
!   enddo
!  close(123456)
!  CALL mpi_finalize(i)
!  stop

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGX(ilevel)%phi,MGX(ilevel)%rhs,ilevel,1,1,MGX(ilevel)%nxc,nyc,nzc,MGgl,myTime)!.....ONLY_PER(SAMK)
        END DO ! iIter

      END IF ! ilevel==1

    END DO ! ilevel

    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, MGX(1)%nxc+Ngl
          var(i,j,k) = MGX(1)%phi(i,j,k)
    END DO ! i
    END DO ! j
    END DO ! k
  !open(123456, file = 'debugvar.dat')
  !write(123456, *) 'VARIABLES="X","Y","var"'
  !write(123456,*)'ZONE F=POINT, I=',(nx-1),', J=',ny-1
  ! k = 1
  ! do j=1,ny-1
  ! do i=1,(nx-1)
  ! write(123456, *)MGX(1)%xc(i),yc(j),MGX(1)%phi(i,j,k)
  ! enddo
  ! enddo
  !close(123456)
  !CALL mpi_finalize(i)
  !stop

END SUBROUTINE MG_COARSE_TO_FINE_X
!---------------------------------------------------------------------



SUBROUTINE MG_FINE_TO_COARSE_Y(var, rr, levBegin, levEnd, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply coarsening on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levBegin: beginning level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the y-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levBegin, levEnd

    REAL(KIND=CGREAL),INTENT (IN) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN) :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    logical :: resCheckFlag

    INTEGER :: iIter
    INTEGER :: ilevel, i, j, k, totIter
    INTEGER :: MGgl

    REAL(KIND=CGREAL) :: myTime
    
    REAL(KIND=CGREAL),DIMENSION(nxc,nyc,nzc) :: r, r1



!   Initialize variable
!   -------------------
    resCheckFlag = .true.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin, levEnd
      IF (ilevel == 1) THEN

!       ilevel Parameters
!       -----------------
        totIter = iterFinest
        MGgl=Ngl

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare_BC(var)!.....................................................COMPLETE(SAMK)

!       Prepare MG Arrays
!       -----------------
        MGY(1)%phi(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)=var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
        MGY(1)%rhs(1:nxc,1:nyc,1:nzc) = rr(1:nxc,1:nyc,1:nzc)

!       Compute ghost cells and impose compatibility for _first_level with GCM
!       ----------------------------------------------------------------------
        IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)        
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure) 
            CALL set_outer_ghost_pres(MGY(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGY(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGY(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGY(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGY(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGY(1)%phi)!..............................COMPLETE(SAMK)   
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                    
          END SELECT              
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation 

      ELSE ! ilevel > 1

!       ilevel Parameters
!       -----------------
        IF (ilevel == mgLevels_Y) THEN
          totIter = iterCoarsest
        ELSE 
          totIter = 1
        END IF ! ilevel==mgLevels_Y
        MGgl=1

!       Prepare MG Arrays
!       -----------------
        MGY(ilevel)%phi(:,:,:) = 0.0_CGREAL

        do k=1,            nzc
        do j=0,MGY(ilevel)%nyc+1
        do i=0,            nxc+1
          iblank_MG(i,j,k) = MGY(ilevel)%iblank(i,j,k)
        end do
        end do
        end do

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGY(ilevel)%phi,nxc,MGY(ilevel)%nyc,nzc)!.........................ONLY_PER(SAMK)

      END IF ! ilevel==1

!     Invoke Line SOR
!     ---------------
      DO iIter = 1, totIter
        CALL MG_itsolv(MGY(ilevel)%phi,MGY(ilevel)%rhs,1,ilevel,1,nxc,MGY(ilevel)%nyc,nzc,MGgl,myTime)!.......ONLY_PER(SAMK)
      END DO ! iIter
      
!     For all levels except the coarsest, compute residuals and apply restriction.
!     ----------------------------------------------------------------------------
      IF ( ilevel /= levEnd ) THEN
        CALL MG_Residual(MGY(ilevel)%phi,MGY(ilevel)%rhs,MGY(ilevel)%res,1,ilevel,1,  &
                         nxc,MGY(ilevel)%nyc,nzc,MGgl,resCheckFlag)!......................ONLY_PER(SAMK)

        CALL MG_Restrict(ilevel,JCOORD)!.............................................COMPLETE(SAMK)
      END IF ! ilevel /= levEnd

    END DO  ! ilevel

END SUBROUTINE MG_FINE_TO_COARSE_Y
!---------------------------------------------------------------------



SUBROUTINE MG_COARSE_TO_FINE_Y(var, rr, levEnd, levBegin, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply refinement on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levEnd: end level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the y-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     
    
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levEnd, levBegin

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    INTEGER :: ilevel, i, j, k, iIter, totIter, MGgl

    logical :: resCheckFlag

    REAL(KIND=CGREAL) :: myTime



!   Initialize variable
!   -------------------
    resCheckFlag = .false.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin-1, levEnd, -1 
    
      IF (ilevel==1) THEN
        MGgl=Ngl
      ELSE
        MGgl=1
      END IF

      IF (internal_boundary_present==INTR_BOUND_PRESENT) THEN
        DO k = 1     ,               nzc
        DO j = 0, MGY(ilevel+1)%nyc+1
        DO i = 0,               nxc+1
          iblank_MG(i,j,k) = MGY(ilevel+1)%iblank(i,j,k)
        END DO ! i
        END DO ! j
        END DO ! k
      END IF ! internal_boundary_present

!     Invoke prolongation
!     -------------------
      CALL MG_Prolong(ilevel,JCOORD)!................................................COMPLETE(SAMK)
#     ifdef MPI
        CALL par_comm_MGY(ilevel,myTime)
        commTime=commTime+myTime
#     endif

!     Load variables
!     --------------
      DO k = 1     ,             nzc
      DO j = 1-MGgl, MGY(ilevel)%nyc+MGgl
      DO i = 1-MGgl,             nxc+MGgl
        iblank_MG(i,j,k) = MGY(ilevel)%iblank(i,j,k)
      END DO ! i
      END DO ! j
      END DO ! k

      if (ilevel/=1) then

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGY(ilevel)%phi, nxc, MGY(ilevel)%nyc, nzc)!......................ONLY_PER(SAMK)

!       Iteration count
!       ---------------
        totIter = iterInter

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGY(ilevel)%phi,MGY(ilevel)%rhs,1,ilevel,1,nxc,MGY(ilevel)%nyc,nzc,1,myTime)!........ONLY_PER(SAMK)
        END DO ! iIter
           
!       Compute residuals
!       -----------------
!!        IF (infoconv == 1) THEN
!!          CALL MG_Residual(var,r, 1, n, 1, nxc, nyc_MG(ilevel), nzc, resCheckFlag)
!!        END IF ! infoConv
      ELSE

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)        
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
            CALL set_outer_ghost_pres(MGY(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGY(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGY(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGY(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGY(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGY(1)%phi)!..............................COMPLETE(SAMK)  
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                     
          END SELECT 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation

        CALL MG_Prepare_BC(MGY(ilevel)%phi)!.........................................COMPLETE(SAMK)

!       Iteration count
!       ---------------
        totIter = iterFinest

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGY(ilevel)%phi,MGY(ilevel)%rhs,1,ilevel,1,nxc,MGY(ilevel)%nyc,nzc,MGgl,myTime)!.....ONLY_PER(SAMK)
        END DO ! iIter

      END IF ! ilevel==1
    END DO ! ilevel

    DO k = 0, nzc+1
    DO j = 1-Ngl, MGY(1)%nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
      var(i,j,k) = MGY(1)%phi(i,j,k)
    END DO ! i
    END DO ! j
    END DO ! k

END SUBROUTINE MG_COARSE_TO_FINE_Y
!---------------------------------------------------------------------



SUBROUTINE MG_FINE_TO_COARSE_Z(var, rr, levBegin, levEnd, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply coarsening on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levBegin: beginning level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the z-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levBegin, levEnd

    REAL(KIND=CGREAL),INTENT (IN) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN) :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    logical :: resCheckFlag

    INTEGER :: iIter
    INTEGER :: ilevel, i, j, k, totIter
    INTEGER :: MGgl

    REAL(KIND=CGREAL) :: myTime
    
    REAL(KIND=CGREAL),DIMENSION(nxc,nyc,nzc) :: r, r1



!   Initialize variable
!   -------------------
    resCheckFlag = .true.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin, levEnd
      IF (ilevel == 1) THEN

!       ilevel Parameters
!       -----------------
        totIter = iterFinest
        MGgl=Ngl

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare_BC(var)!.....................................................COMPLETE(SAMK)

!       Prepare MG Arrays
!       -----------------
        MGZ(1)%phi(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)=var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
        MGZ(1)%rhs(1:nxc,1:nyc,1:nzc) = rr(1:nxc,1:nyc,1:nzc)

!       Compute ghost cells and impose compatibility for _first_level with GCM
!       ----------------------------------------------------------------------
        IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)        
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
            CALL set_outer_ghost_pres(MGZ(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGZ(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGZ(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGZ(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGZ(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGZ(1)%phi)!..............................COMPLETE(SAMK) 
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                      
          END SELECT           
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation 

      ELSE ! ilevel > 1

!       ilevel Parameters
!       -----------------
        IF (ilevel == mgLevels_Z) THEN
          totIter = iterCoarsest
        ELSE 
          totIter = 1
        END IF ! ilevel==mgLevels_Z
        MGgl=1

!       Prepare MG Arrays
!       -----------------
        MGZ(ilevel)%phi(:,:,:) = 0.0_CGREAL

        do k=1,MGZ(ilevel)%nzc
        do j=0,            nyc+1
        do i=0,            nxc+1
          iblank_MG(i,j,k) = MGZ(ilevel)%iblank(i,j,k)
        end do
        end do
        end do

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGZ(ilevel)%phi,nxc,nyc,MGZ(ilevel)%nzc)!.........................ONLY_PER(SAMK)

      END IF ! ilevel==1

!     Invoke Line SOR
!     ---------------
      DO iIter = 1, totIter
        CALL MG_itsolv(MGZ(ilevel)%phi,MGZ(ilevel)%rhs,1,1,ilevel,nxc,nyc,MGZ(ilevel)%nzc,MGgl,myTime)!.......ONLY_PER(SAMK)
      END DO ! iIter

!     For all levels except coarsest one, compute residuals and apply restriction.
!     ----------------------------------------------------------------------------
      IF ( ilevel /= levEnd ) THEN
        CALL MG_Residual(MGZ(ilevel)%phi,MGZ(ilevel)%rhs,MGZ(ilevel)%res,1,1,ilevel,  &
                         nxc,nyc,MGZ(ilevel)%nzc,MGgl,resCheckFlag)!......................ONLY_PER(SAMK)

        CALL MG_Restrict(ilevel,KCOORD)!.............................................COMPLETE(SAMK)
 
      END IF ! ilevel /= levEnd
      
    END DO  ! ilevel

END SUBROUTINE MG_FINE_TO_COARSE_Z
!---------------------------------------------------------------------



SUBROUTINE MG_COARSE_TO_FINE_Z(var, rr, levEnd, levBegin, commTime)

! ---------------------------------------------------------------------
!  Purpose: Apply refinement on the variables between levBegin to nLev
!
!  Input: 
!         rr:  right-hand side
!         var: pressure variable
!         levEnd: end level
!         nLev: end level
!  
!  Output: 
!         var: pressure variable
!
!  Notes:
!        1. This is done _only_ for the z-direction
! ---------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     
    

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: levEnd, levBegin

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: rr(nxc,nyc,nzc)
    REAL(KIND=CGREAL),INTENT (INOUT) :: commTime

!   Local variables
!   ---------------
    INTEGER :: ilevel, i, j, k, iIter, totIter, MGgl

    logical :: resCheckFlag

    REAL(KIND=CGREAL) :: myTime



!   Initialize variable
!   -------------------
    resCheckFlag = .false.
    
!   Loop over levels
!   ----------------
    DO ilevel = levBegin-1, levEnd, -1
    
      IF (ilevel==1) THEN
        MGgl=Ngl
      ELSE
        MGgl=1
      END IF

      IF (internal_boundary_present==INTR_BOUND_PRESENT) THEN
        DO k = 1     , MGZ(ilevel+1)%nzc
        DO j = 0,               nyc+1
        DO i = 0,               nxc+1
          iblank_MG(i,j,k) = MGZ(ilevel+1)%iblank(i,j,k)
        END DO ! i
        END DO ! j
        END DO ! k
      END IF ! internal_boundary_present



!     Invoke prolongation
!     -------------------
      CALL MG_Prolong(ilevel,KCOORD)!................................................COMPLETE(SAMK)
#     ifdef MPI
        CALL par_comm_MGZ(ilevel,myTime)
        commTime=commTime+myTime
#     endif

!     Load variables
!     --------------
      DO k = 1     , MGZ(ilevel)%nzc
      DO j = 1-MGgl,             nyc+MGgl
      DO i = 1-MGgl,             nxc+MGgl
        iblank_MG(i,j,k) = MGZ(ilevel)%iblank(i,j,k)
      END DO ! i
      END DO ! j
      END DO ! k

      if (ilevel/=1) then

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
        CALL MG_Prepare(MGZ(ilevel)%phi, nxc, nyc, MGZ(ilevel)%nzc)!......................ONLY_PER(SAMK)

!       Iteration count
!       ---------------
        totIter = iterInter

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGZ(ilevel)%phi,MGZ(ilevel)%rhs,1,1,ilevel,nxc,nyc,MGZ(ilevel)%nzc,1,myTime)!........ONLY_PER(SAMK)
        END DO ! iIter
           
!       Compute residuals
!       -----------------
!!        IF (infoconv == 1) THEN
!!          CALL MG_Residual(var,r, n, 1, 1, mgrid_I(n), ny-1, nz-1, resCheckFlag)
!!        END IF ! infoConv
      ELSE  

!       Prepare ium, iup, jum, jup, kum, kup
!       ------------------------------------
          IF ( boundary_formulation == GCM_METHOD ) THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)          
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
            CALL set_outer_ghost_pres(MGZ(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_p_set_bc_internal(MGZ(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_p_compatibility(MGZ(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_scalar) 
            CALL set_outer_ghost_scalar(MGZ(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_scalar_set_bc_internal(MGZ(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_scalar_compatibility(MGZ(1)%phi)!..............................COMPLETE(SAMK)
          CASE(index_potential) 
            CALL set_outer_ghost_potential(MGX(1)%phi)!..........................................ONLY_PER(SAMK)
            CALL GCM_potential_set_bc_internal(MGX(1)%phi)!....................................COMPLETE(SAMK)
            CALL GCM_enforce_potential_compatibility(MGX(1)%phi)!..............................COMPLETE(SAMK)                       
          END SELECT  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
        ENDIF ! boundary_formulation

        CALL MG_Prepare_BC(MGZ(ilevel)%phi)

!       Iteration count
!       ---------------
        totIter = iterFinest

!       Invoke Line SOR
!       ---------------
        DO iIter = 1, totIter
          CALL MG_itsolv(MGZ(ilevel)%phi,MGZ(ilevel)%rhs,1,1,ilevel,nxc,nyc,MGZ(ilevel)%nzc,MGgl,myTime)!.....ONLY_PER(SAMK)
        END DO ! iIter

      END IF ! ilevel==1
    END DO ! ilevel

    DO k = 0, MGz(1)%nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
      var(i,j,k) = MGz(1)%phi(i,j,k)
    END DO ! i
    END DO ! j
    END DO ! k

END SUBROUTINE MG_COARSE_TO_FINE_Z
!---------------------------------------------------------------------



SUBROUTINE MG_Prepare_BC(var)

! ----------------------------------------------------------------------------------
!  Purpose: Prepare the boundary conditions for _FIRST_ level only
!
!  Output: 
!         ium_MG, iup_MG, jum_MG, jup_MG, kum_MG, kup_MG
!
!  Notes:
!   1. With new unified formulation, there is no need to supply bc for lower levels
! ----------------------------------------------------------------------------------
 
    USE flow_parameters
    USE boundary_arrays
    USE MG_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)



!   Load values for first level
!   ---------------------------
    iblank_MG(:,:,:) = iblank(:,:,:)
 
    !ium_MG(:,:,:) = ium(:,:,:)
    !iup_MG(:,:,:) = iup(:,:,:)
    !jum_MG(:,:,:) = jum(:,:,:)
    !jup_MG(:,:,:) = jup(:,:,:)
    !kum_MG(:,:,:) = kum(:,:,:)
    !kup_MG(:,:,:) = kup(:,:,:)
!lingxiao change
    ium_MG(1:nxc,1:nyc,1:nzc) = ium(1:nxc,1:nyc,1:nzc)
    iup_MG(1:nxc,1:nyc,1:nzc) = iup(1:nxc,1:nyc,1:nzc)
    jum_MG(1:nxc,1:nyc,1:nzc) = jum(1:nxc,1:nyc,1:nzc)
    jup_MG(1:nxc,1:nyc,1:nzc) = jup(1:nxc,1:nyc,1:nzc)
    kum_MG(1:nxc,1:nyc,1:nzc) = kum(1:nxc,1:nyc,1:nzc)
    kup_MG(1:nxc,1:nyc,1:nzc) = kup(1:nxc,1:nyc,1:nzc)
    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
SELECT CASE(POSSION_INDEX) 
CASE(index_pressure)

    IF ((myCoords(1)==0      ) .AND. (pbcx1 == PBC_DIRICHLET)) var(0    ,:,:) = 2.0_CGREAL * pppx1 - var(1, :,:)
    IF ((myCoords(1)==Np(1)-1) .AND. (pbcx2 == PBC_DIRICHLET)) var(nxc+1,:,:) = 2.0_CGREAL * pppx2 - var(nxc,:,:)

    IF ((myCoords(2)==0      ) .AND. (pbcy1 == PBC_DIRICHLET)) var(:,0    ,:) = 2.0_CGREAL * pppy1 - var(:,1, :)
    IF ((myCoords(2)==Np(2)-1) .AND. (pbcy2 == PBC_DIRICHLET)) var(:,nyc+1,:) = 2.0_CGREAL * pppy2 - var(:, nyc,:)

    IF ( pbcz1 == PBC_DIRICHLET )                              var(:,:,0    ) = 2.0_CGREAL * pppz1 - var(:,:,0)
    IF ( pbcz2 == PBC_DIRICHLET )                              var(:,:,nzc+1) = 2.0_CGREAL * pppz2 - var(:,:,nzc)
    

CASE(index_scalar)

    IF ((myCoords(1)==0      ) .AND. (Xbcx1 == XBC_DIRICHLET)) var(0    ,:,:) = 2.0_CGREAL * xxxx1 - var(1, :,:)
    IF ((myCoords(1)==Np(1)-1) .AND. (Xbcx2 == XBC_DIRICHLET)) var(nxc+1,:,:) = 2.0_CGREAL * xxxx2 - var(nxc,:,:)

    IF ((myCoords(2)==0      ) .AND. (Xbcy1 == XBC_DIRICHLET)) var(:,0    ,:) = 2.0_CGREAL * xxxy1 - var(:,1, :)
    IF ((myCoords(2)==Np(2)-1) .AND. (Xbcy2 == XBC_DIRICHLET)) var(:,nyc+1,:) = 2.0_CGREAL * xxxy2 - var(:, nyc,:)

    IF ( Xbcz1 == XBC_DIRICHLET )                              var(:,:,0    ) = 2.0_CGREAL * xxxz1 - var(:,:,0)
    IF ( Xbcz2 == XBC_DIRICHLET )                              var(:,:,nzc+1) = 2.0_CGREAL * xxxz2 - var(:,:,nzc)
    
    
CASE(index_potential)

    IF ((myCoords(1)==0      ) .AND. (Phibcx1 == PhiBC_DIRICHLET)) var(0    ,:,:) = 2.0_CGREAL * phix1 - var(1, :,:)
    IF ((myCoords(1)==Np(1)-1) .AND. (Phibcx2 == PhiBC_DIRICHLET)) var(nxc+1,:,:) = 2.0_CGREAL * phix2 - var(nxc,:,:)

    IF ((myCoords(2)==0      ) .AND. (Phibcy1 == PhiBC_DIRICHLET)) var(:,0    ,:) = 2.0_CGREAL * phiy1 - var(:,1, :)
    IF ((myCoords(2)==Np(2)-1) .AND. (Phibcy2 == PhiBC_DIRICHLET)) var(:,nyc+1,:) = 2.0_CGREAL * phiy2 - var(:, nyc,:)

    IF ( Phibcz1 == PhiBC_DIRICHLET )                              var(:,:,0    ) = 2.0_CGREAL * phiz1 - var(:,:,0)
    IF ( Phibcz2 == PhiBC_DIRICHLET )                              var(:,:,nzc+1) = 2.0_CGREAL * phiz2 - var(:,:,nzc)    

END SELECT   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

END SUBROUTINE MG_Prepare_BC
!---------------------------------------------------------------------



SUBROUTINE MG_Prepare(var,IL,JL,KL) 

! -----------------------------------------------------------------------------------
!  Purpose: Preparing the grids information for smoothing.
!
!  Input: ilevel     -- current level
!         IL, JL, KL -- the number of grid points in x, y, z direction respectively
!                       at current level (used for semi-coarsening MG method)
!
!  Output: all grids information such as ium_mg, iup_mg etc. 
! -----------------------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_arrays

    IMPLICIT NONE
 
!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: IL, JL, KL

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(0:IL+1,0:JL+1,0:KL+1)

!   Local variables
!   ---------------
    INTEGER :: i, j, k
    INTEGER :: ttemp



!   -------------------------------------------------------
!   Compute ium, iup, jum, jup, kum, kup for coarser levels
!   Notes:
!    1. This kernel is _only_ relevant for solid bodies
!    2. Probably broken for membranes
!   -------------------------------------------------------

    IF ( internal_boundary_present == 1 ) THEN 
      DO k = 1, KL  

!       Compute ium and iup
!       -------------------
        DO j = 1, JL 
        DO i = 1, IL
          ttemp = ((iblank_MG(i-1,j,k)-iblank_MG(i,j,k))+1)/2
          ium_MG(i,j,k) = REAL(ttemp, KIND=CGREAL)

          ttemp = ((iblank_MG(i+1,j,k)-iblank_MG(i,j,k))+1)/2
          iup_MG(i,j,k) = REAL(ttemp, KIND=CGREAL) 
        END DO ! i
        END DO ! j

!       Compute jum and jup
!       -------------------
        DO j = 1, JL
        DO i = 1, IL
          ttemp = ((iblank_MG(i,j-1,k)-iblank_MG(i,j,k))+1)/2
          jum_MG(i,j,k) = REAL(ttemp, KIND=CGREAL)

          ttemp = ((iblank_MG(i,j+1,k)-iblank_MG(i,j,k))+1)/2
          jup_MG(i,j,k) = REAL(ttemp, KIND=CGREAL)
        END DO ! i
        END DO ! j

      END DO !  K

!     Compute kum and kup
!     -------------------
      DO j = 1, JL
      DO i = 1, IL
        DO k = 1, KL
          ttemp = ((iblank_MG(i,j,k-1)-iblank_MG(i,j,k))+1)/2
          kum_MG(i,j,k) = REAL(ttemp, KIND=CGREAL)

          ttemp = ((iblank_MG(i,j,k+1)-iblank_MG(i,j,k))+1)/2
          kup_MG(i,j,k) = REAL(ttemp, KIND=CGREAL)
        END DO ! k
      END DO ! i
      END DO ! j

    END IF ! internal_boundary_present

!   -------------------------------------------------------
!   Modify ium, iup, jum, jup, kum, kup
!   for periodic and Dirichlet pressure boundary conditions
!   -------------------------------------------------------

!   x-direction
!   -----------
!!    IF ( bcx1 == BC_TYPE_PERIODIC .AND.  &
!!         bcx2 == BC_TYPE_PERIODIC .AND.  &
!!         IL==nx-1) THEN
!!      ium_MG(1 ,1:JL,1:KL) = 0.0_CGREAL
!!      iup_MG(IL,1:JL,1:KL) = 0.0_CGREAL
!!    ELSE
!!      ium_MG(1 ,1:JL,1:KL) = 1.0_CGREAL
!!      iup_MG(IL,1:JL,1:KL) = 1.0_CGREAL
!!    END IF ! bcx1

    IF (myCoords(1)==0) THEN
      IF ( pbcx1 == PBC_DIRICHLET ) THEN 
        ium_MG(1,1:JL,1:KL) = 0.0_CGREAL
        var   (0,1:JL,1:KL) = 0.0_CGREAL
      ELSE
        ium_MG(1,1:JL,1:KL) = 1.0_CGREAL
      END IF ! pbcx1
    END IF

    IF (myCoords(1)==Np(1)-1) THEN
      IF ( pbcx2 == PBC_DIRICHLET ) THEN
        iup_MG(IL  ,1:JL,1:KL) = 0.0_CGREAL
        var   (IL+1,1:JL,1:KL) = 0.0_CGREAL
      ELSE
        iup_MG(IL  ,1:JL,1:KL) = 1.0_CGREAL
      END IF ! pbcx2
    END IF

!   y-direction
!   -----------
!!    IF ( bcy1 == BC_TYPE_PERIODIC .AND.  &
!!         bcy2 == BC_TYPE_PERIODIC .AND.  &
!!         JL==ny-1) THEN   
!!      jum_MG(1:IL,1 ,1:KL) = 0.0_CGREAL
!!      jup_MG(1:IL,JL,1:KL) = 0.0_CGREAL
!!    ELSE
!!      jum_MG(1:IL,1 ,1:KL) = 1.0_CGREAL
!!      jup_MG(1:IL,JL,1:KL) = 1.0_CGREAL
!!    END IF ! bcy1

    IF (myCoords(2)==0) THEN
      IF ( pbcy1 == PBC_DIRICHLET ) THEN 
        jum_MG(1:IL,1,1:KL) = 0.0_CGREAL
           var(1:IL,0,1:KL) = 0.0_CGREAL
      ELSE
        jum_MG(1:IL,1,1:KL) = 1.0_CGREAL
      END IF ! pbcy1
    END IF

    IF (myCoords(2)==Np(2)-1) THEN
      IF ( pbcy2 == PBC_DIRICHLET ) THEN
        jup_MG(1:IL,JL  ,1:KL) = 0.0_CGREAL
        var   (1:IL,JL+1,1:KL) = 0.0_CGREAL
      ELSE
        jup_MG(1:IL,JL  ,1:KL) = 1.0_CGREAL
      END IF ! pbcy2 
    END IF
    
!   z-direction
!   -----------
!!    IF ( bcz1 == BC_TYPE_PERIODIC .AND.  &
!!         bcz2 == BC_TYPE_PERIODIC .AND.  &
!!         KL==nz-1) THEN   
!!      kum_MG(1:IL,1:JL,1)  = 0.0_CGREAL
!!      kup_MG(1:IL,1:JL,KL) = 0.0_CGREAL
!!    ELSE
!!      kum_MG(1:IL,1:JL,1)  = 1.0_CGREAL
!!      kup_MG(1:IL,1:JL,KL) = 1.0_CGREAL
!!    END IF ! bcz1

    IF ( bcz1 == BC_TYPE_PERIODIC .AND.  &           !Ehsan added for periodic BC
         bcz2 == BC_TYPE_PERIODIC .AND.  &          
        KL==nzc ) THEN                             
      kum_MG(1:IL,1:JL,1)  = 0.0_CGREAL          
      kup_MG(1:IL,1:JL,KL) = 0.0_CGREAL         
    ELSE                                       
      kum_MG(1:IL,1:JL,1)  = 1.0_CGREAL    
      kup_MG(1:IL,1:JL,KL) = 1.0_CGREAL        
    END IF ! bcz1                           


    IF ( pbcz1 == PBC_DIRICHLET ) THEN 
      kum_MG(1:IL,1:JL,1) = 0.0_CGREAL
      var   (1:IL,1:JL,0) = 0.0_CGREAL
!    else                                           !Ehsan removed
!      kum_MG(1:IL,1:JL,1) = 1.0_CGREAL             !Ehsan removed
    END IF ! pbcz1

    IF ( pbcz2 == PBC_DIRICHLET ) THEN
      kup_MG(1:IL,1:JL,KL  ) = 0.0_CGREAL
      var   (1:IL,1:JL,KL+1) = 0.0_CGREAL
!    else                                           !Ehsan removed
!      kup_MG(1:IL,1:JL,KL  ) = 1.0_CGREAL          !Ehsan removed
    END IF ! pbcy2 

END SUBROUTINE MG_Prepare
!---------------------------------------------------------------------



SUBROUTINE MG_itsolv(var,r,nLevX,nLevY,nLevZ,IL,JL,KL,glN,myCommTime)
! ---------------------------------------------------------------------
!  Purpose: Solve system using a Line-SOR with an option for Red-Black
! ---------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_parameters
    USE boundary_arrays
    USE solver_arrays
    USE flow_arrays
    USE MG_arrays
	USE cutcell_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 	

    IMPLICIT NONE

!   parameters
!   ----------
    INTEGER, INTENT(IN)  :: nLevX, nLevY, nLevZ, IL, JL, KL, glN

    REAL(KIND=CGREAL),INTENT (INOUT) :: var(1-glN:IL+glN,1-glN:JL+glN,0:KL+1)
    REAL(KIND=CGREAL),INTENT (IN)    :: r(IL,JL,KL)
    REAL(KIND=CGREAL),INTENT (  OUT) :: myCommTime

!   Local variables
!   ---------------
    INTEGER :: i  , j  , k
    INTEGER :: iG , jG , kG

    INTEGER :: iinit, jinit, kinit, Ncolor   
  
    

    REAL(KIND=CGREAL) :: mgLevelFlag

    REAL(KIND=CGREAL), DIMENSION (IL) :: kill4GCMxm, kill4GCMxp
    REAL(KIND=CGREAL), DIMENSION (JL) :: kill4GCMym, kill4GCMyp
    REAL(KIND=CGREAL), DIMENSION (KL) :: kill4GCMzm, kill4GCMzp
    
    REAL(kind=CGREAL) :: myTime
!DEBUG
!    REAL(KIND=CGREAL) :: tmpvar(1-glN:IL+glN,1-glN:JL+glN,0:KL+1)
!DEBUG



!   -------------------------------------------------------------------------
!   Initialize variables
!   Notes:
!    1. mgLevelFlag activates pGhost and boundPresSource for first level only
!    2. for MG Level higher than 1, PPE system is solved using SSM whether
!       GCM or SSM method is invoked.
!   -------------------------------------------------------------------------

    IF ( nLevX*nLevY*nLevZ == 1 ) THEN
      mgLevelFlag = 1.0_CGREAL
    ELSE
      mgLevelFlag = 0.0_CGREAL
    ENDIF

!   Initializing Neumann/Dirichlet BC flags (SAMK/RM)
!   -------------------------------------------------
    Kill4GCMxm(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMxp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    IF (myCoords(1)==0      ) Kill4GCMxm(1 ) = 1.0_CGREAL
    IF (myCoords(1)==Np(1)-1) Kill4GCMxp(IL) = 1.0_CGREAL

    Kill4GCMym(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMyp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    IF (myCoords(2)==0      ) Kill4GCMym(1 ) = 1.0_CGREAL
    IF (myCoords(2)==Np(2)-1) Kill4GCMyp(JL) = 1.0_CGREAL

    Kill4GCMzm(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMzp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMzm(1 ) = 1.0_CGREAL
    Kill4GCMzp(KL) = 1.0_CGREAL
    
    myCommTime=0.d0
    myTime=0.d0

!   -------------------------------------
!    Impose periodic boundary conditions
!   -------------------------------------
!!    IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN
!!      IF (bcx1 == BC_TYPE_PERIODIC .OR. &
!!          bcy1 == BC_TYPE_PERIODIC .OR. &
!!          bcz1 == BC_TYPE_PERIODIC) THEN
!!       CALL apply_periodic_pres(var)
!!      END IF
!!    ENDIF

    IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN       !Ehsan added for periodic BC
      IF (bcx1 == BC_TYPE_PERIODIC .OR. &                
          bcy1 == BC_TYPE_PERIODIC .OR. &              
          bcz1 == BC_TYPE_PERIODIC) THEN               
          do i=0,IL+1                               
            do j=0,JL+1                             
              var(i,j,0) = var(i,j,KL)             
              var(i,j,KL+1) = var(i,j,1)              
            enddo                                  
          enddo                                 
      END IF                                        
    ENDIF                                   
         
!   -------------------------------
!    Line solve in the x-direction
!   -------------------------------
    DO Ncolor = 1, TNcolorX
      kinit = (Ncolor-1)*(ndim-DIM_2D)+ 1
      jinit = Ncolor 

      DO k = kinit, KL, kStep
        kG=k

        DO j = jinit, JL, jStep
          jG=MGY(nLevY)%L2GJ(j)

          DO i=1, IL
            iG=MGX(nLevX)%L2GI(i)

            amx(i) = MGX(nlevX)%dxcinv(iG  )*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - ium_mg(i,j,k)*Kill4GCMxm(i))*(AREA_W(i,j,k)*mglevelflag + 1.-mglevelflag)
            apx(i) = MGX(nlevX)%dxcinv(iG+1)*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - iup_mg(i,j,k)*Kill4GCMxp(i))*(AREA_E(i,j,k)*mglevelflag + 1.-mglevelflag)
            acx(i) = - ( amx(i) + apx(i) )

            amy(j) = MGY(nlevY)%dycinv(jG  )*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jum_mg(i,j,k)*Kill4GCMym(j))*(AREA_S(i,j,k)*mglevelflag + 1.-mglevelflag)
            apy(j) = MGY(nlevY)%dycinv(jG+1)*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jup_mg(i,j,k)*Kill4GCMyp(j))*(AREA_N(i,j,k)*mglevelflag + 1.-mglevelflag)
            acy(j) = - ( amy(j) + apy(j) )

            amz(k) = MGZ(nlevZ)%dzcinv(kG  )*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kum_mg(i,j,k)*Kill4GCMzm(k))*KillFor2D*(AREA_B(i,j,k)*mglevelflag + 1.-mglevelflag)
            apz(k) = MGZ(nlevZ)%dzcinv(kG+1)*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kup_mg(i,j,k)*Kill4GCMzp(k))*KillFor2D*(AREA_F(i,j,k)*mglevelflag + 1.-mglevelflag)
            acz(k) = - ( amz(k) + apz(k) )
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
            rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k)) &
                              - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k)) &
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k)) &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k)) &
                              + mgLevelFlag *                                    &
                              (- pGhost(i-1,j,k)*amx(i)*ium_MG(i,j,k)*gcmFlag    &
                               - pGhost(i+1,j,k)*apx(i)*iup_MG(i,j,k)*gcmFlag    &
                               - pGhost(i,j-1,k)*amy(j)*jum_MG(i,j,k)*gcmFlag    &
                               - pGhost(i,j+1,k)*apy(j)*jup_MG(i,j,k)*gcmFlag    &
                               - pGhost(i,j,k-1)*amz(k)*kum_MG(i,j,k)*gcmFlag    &
                               - pGhost(i,j,k+1)*apz(k)*kup_MG(i,j,k)*gcmFlag    &
                               + boundPresSource(i-1,j,k)                        &
                               + boundPresSource(i+1,j,k)                        &
                               + boundPresSource(i,j-1,k)                        &
                               + boundPresSource(i,j+1,k)                        &
                               + boundPresSource(i,j,k-1)                        &
                               + boundPresSource(i,j,k+1) )                                
            CASE(index_scalar)
            rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k)) &
                              - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k)) &
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k)) &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k)) &
                              + mgLevelFlag *                                    &
                              (- XGhost(i-1,j,k)*amx(i)*ium_MG(i,j,k)*gcmFlag    &
                               - XGhost(i+1,j,k)*apx(i)*iup_MG(i,j,k)*gcmFlag    &
                               - XGhost(i,j-1,k)*amy(j)*jum_MG(i,j,k)*gcmFlag    &
                               - XGhost(i,j+1,k)*apy(j)*jup_MG(i,j,k)*gcmFlag    &
                               - XGhost(i,j,k-1)*amz(k)*kum_MG(i,j,k)*gcmFlag    &
                               - XGhost(i,j,k+1)*apz(k)*kup_MG(i,j,k)*gcmFlag    &
                               + boundScalarSource(i-1,j,k)                      &
                               + boundScalarSource(i+1,j,k)                      &
                               + boundScalarSource(i,j-1,k)                      &
                               + boundScalarSource(i,j+1,k)                      &
                               + boundScalarSource(i,j,k-1)                      &
                               + boundScalarSource(i,j,k+1) )
            CASE(index_potential)
            rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k)) &
                              - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k)) &
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k)) &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k)) &
                              + mgLevelFlag *                                    &
                              (- PhiGhost(i-1,j,k)*amx(i)*ium_MG(i,j,k)*gcmFlag  &
                               - PhiGhost(i+1,j,k)*apx(i)*iup_MG(i,j,k)*gcmFlag  &
                               - PhiGhost(i,j-1,k)*amy(j)*jum_MG(i,j,k)*gcmFlag  &
                               - PhiGhost(i,j+1,k)*apy(j)*jup_MG(i,j,k)*gcmFlag  &
                               - PhiGhost(i,j,k-1)*amz(k)*kum_MG(i,j,k)*gcmFlag  &
                               - PhiGhost(i,j,k+1)*apz(k)*kup_MG(i,j,k)*gcmFlag  &
                               + boundPhiSource(i-1,j,k)                         &
                               + boundPhiSource(i+1,j,k)                         &
                               + boundPhiSource(i,j-1,k)                         &
                               + boundPhiSource(i,j+1,k)                         &
                               + boundPhiSource(i,j,k-1)                         &
                               + boundPhiSource(i,j,k+1) )                               
            END SELECT                               
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

            amx(i) = amx(i)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL-ium_MG(i,j,k))
            apx(i) = apx(i)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL-iup_MG(i,j,k))

            acx(i) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank_MG(i,j,k),KIND=CGREAL) &
                                           +REAL(  iblank_MG(i,j,k),KIND=CGREAL)                     

            rhs(i) = rhs(i)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
          ENDDO ! i 

!         Impose Periodic BC
!         ------------------
!!          IF ( bcx1 == BC_TYPE_PERIODIC .AND. &
!!               bcx2 == BC_TYPE_PERIODIC .AND. &
!!               nLevX == 1      ) THEN
!!            rhs(1)  = rhs(1) - var(IL,j,k)*amx(1)
!!            rhs(IL) = rhs(IL) - var(1,j,k)*apx(IL)
!!          END IF ! bcx1

!         Impose Subdomain/Dirichlet BC (SAMK)
!         ------------------------------------
          rhs(1)  = rhs(1) -var(0   ,j,k)*amx(1 )
          rhs(IL) = rhs(IL)-var(IL+1,j,k)*apx(IL)


!         Solve TDMA in X Direction
!         -------------------------
          CALL tdma(amx(1:IL),acx(1:IL),apx(1:IL),rhs(1:IL),dummy(1:IL),1,IL)

!          DO i=1,IL
!DEBUG
!!            dummy(i)=(rhs(i)-var(i-1,j,k)*amx(i)-var(i+1,j,k)*apx(i))/acx(i)
!!            tmpvar(i,j,k) = var(i,j,k)*(1.0_CGREAL-0.5d0) + 0.5d0*dummy(i) 
!!            var(i,j,k) = var(i,j,k)*(1.0_CGREAL-0.5d0) + 0.5d0*dummy(i) 
!DEBUG
!            var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega) + omega*dummy(i) 
!          ENDDO ! i

!DEBUG
          IF (myCoords(1)==0) THEN
            var(1,j,k) = var(1,j,k)*(1.0_CGREAL-omega) + omega*dummy(1)
          ELSE
            var(1,j,k) = var(1,j,k)*(1.0_CGREAL-0.75 ) + 0.75 *dummy(1)
          END IF
          DO i=2,IL-1
            var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega) + omega*dummy(i) 
          ENDDO ! i
          IF (myCoords(1)==Np(1)-1) THEN
            var(il,j,k) = var(il,j,k)*(1.0_CGREAL-omega) + omega*dummy(il)
          ELSE
            var(il,j,k) = var(il,j,k)*(1.0_CGREAL-0.75 ) + 0.75 *dummy(il) 
          END IF
!DEBUG
        ENDDO ! j 
      ENDDO ! k

!DEBUG
! var=tmpvar
!DEBUG
#     ifdef MPI      
        CALL par_comm_var(var,IL,JL,KL,glN,myTime)
        myCommTime=myCommTime+myTime
#     endif

    ENDDO ! NcolorX
!DEBUG
! return
!DEBUG

!   -------------------------------------
!    Impose periodic boundary conditions
!   -------------------------------------
!!    IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN
!!      IF (bcx1 == BC_TYPE_PERIODIC .OR. &
!!          bcy1 == BC_TYPE_PERIODIC .OR. &
!!          bcz1 == BC_TYPE_PERIODIC) THEN
!!       CALL apply_periodic_pres(var)
!!      END IF                 
!!    ENDIF 

    IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN              !Ehsan added for Periodic BC
      IF (bcx1 == BC_TYPE_PERIODIC .OR. &                      
          bcy1 == BC_TYPE_PERIODIC .OR. &                     
          bcz1 == BC_TYPE_PERIODIC) THEN                     
         do i=0,IL+1 !1-glN,IL+glN
         do j=0,JL+1 !1-glN,JL+glN
          var(i,j,0) = var(i,j,KL)
          var(i,j,KL+1) = var(i,j,1)
         enddo
         enddo
      END IF                                                 
    ENDIF                                                      

!   -------------------------------
!    Line solve in the y-direction
!   -------------------------------
    DO Ncolor = 1, TNcolorY  
      kinit = (Ncolor-1)*(ndim-DIM_2D)+ 1
      iinit = Ncolor  

      DO k= kinit, KL,kstep  
        kG=k

        DO i= iinit, IL, istep
          iG=MGX(nLevX)%L2GI(i)

          DO j = 1, JL
            jG=MGY(nLevY)%L2GJ(j)

            amx(i) =   MGX(nlevX)%dxcinv(iG  )*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - ium_mg(i,j,k)*Kill4GCMxm(i))*(AREA_W(i,j,k)*mglevelflag + 1.-mglevelflag)
            apx(i) =   MGX(nlevX)%dxcinv(iG+1)*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - iup_mg(i,j,k)*Kill4GCMxp(i))*(AREA_E(i,j,k)*mglevelflag + 1.-mglevelflag)
            acx(i) = - ( amx(i) + apx(i) )

            amy(j) =   MGY(nlevY)%dycinv(jG  )*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jum_mg(i,j,k)*Kill4GCMym(j))*(AREA_S(i,j,k)*mglevelflag + 1.-mglevelflag)
            apy(j) =   MGY(nlevY)%dycinv(jG+1)*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jup_mg(i,j,k)*Kill4GCMyp(j))*(AREA_N(i,j,k)*mglevelflag + 1.-mglevelflag)
            acy(j) = - ( amy(j) + apy(j) )

            amz(k) =   MGZ(nlevZ)%dzcinv(kG  )*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kum_mg(i,j,k)*Kill4GCMzm(k))*KillFor2D*(AREA_B(i,j,k)*mglevelflag + 1.-mglevelflag)
            apz(k) =   MGZ(nlevZ)%dzcinv(kG+1)*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kup_mg(i,j,k)*Kill4GCMzp(k))*KillFor2D*(AREA_F(i,j,k)*mglevelflag + 1.-mglevelflag)
            acz(k) = - ( amz(k) + apz(k) )

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
            rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) )  &  
                              - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) )  & 
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k) )  &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k) )  &
                              + mgLevelFlag *                                      &
                             ( - pGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag      &
                               - pGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag      &
                               - pGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag      &
                               - pGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag      &
                               - pGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag      &
                               - pGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag      &
                               + boundPresSource(i-1,j,k)                          &
                               + boundPresSource(i+1,j,k)                          &
                               + boundPresSource(i,j-1,k)                          &
                               + boundPresSource(i,j+1,k)                          &
                               + boundPresSource(i,j,k-1)                          &
                               + boundPresSource(i,j,k+1) )
                               
          CASE(index_scalar)
            rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) )  &  
                              - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) )  & 
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k) )  &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k) )  &
                              + mgLevelFlag *                                      &
                             ( - XGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag      &
                               - XGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag      &
                               - XGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag      &
                               - XGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag      &
                               - XGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag      &
                               - XGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag      &
                               + boundScalarSource(i-1,j,k)                        &
                               + boundScalarSource(i+1,j,k)                        &
                               + boundScalarSource(i,j-1,k)                        &
                               + boundScalarSource(i,j+1,k)                        &
                               + boundScalarSource(i,j,k-1)                        &
                               + boundScalarSource(i,j,k+1) )    
          CASE(index_potential)
            rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) )  &  
                              - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) )  & 
                              - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum_mg(i,j,k) )  &
                              - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup_mg(i,j,k) )  &
                              + mgLevelFlag *                                      &
                             ( - PhiGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag    &
                               - PhiGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag    &
                               - PhiGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag    &
                               - PhiGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag    &
                               - PhiGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag    &
                               - PhiGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag    &
                               + boundPhiSource(i-1,j,k)                           &
                               + boundPhiSource(i+1,j,k)                           &
                               + boundPhiSource(i,j-1,k)                           &
                               + boundPhiSource(i,j+1,k)                           &
                               + boundPhiSource(i,j,k-1)                           &
                               + boundPhiSource(i,j,k+1) )                                       
          END SELECT
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)                               

            amy(j) = amy(j)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL - jum_mg(i,j,k))
            apy(j) = apy(j)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL - jup_mg(i,j,k))

            acy(j) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank_MG(i,j,k),KIND=CGREAL) &
                                           +REAL(  iblank_MG(i,j,k),KIND=CGREAL)       

            rhs(j) = rhs(j)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
          ENDDO ! j

!         Impose Periodic BC
!         ------------------
!!          IF ( bcy1 == BC_TYPE_PERIODIC .AND. &
!!               bcy2 == BC_TYPE_PERIODIC .AND. &
!!               nLevY == 1       ) THEN
!!            rhs(1)  = rhs(1)  - var(i,JL,k) *amy(1)
!!            rhs(JL) = rhs(JL) - var(i,1,k)  *apy(JL)
!!          END IF ! bcy1

!         Impose Subdomain/Dirichlet BC (SAMK)
!         ------------------------------------
          rhs(1) =rhs(1) -var(i,0   ,k)*amy(1 )
          rhs(JL)=rhs(JL)-var(i,JL+1,k)*apy(JL)

!         Solve TDMA in Y Direction
!         -------------------------
          CALL tdma(amy,acy,apy,rhs(1:JL),dummy(1:JL),1,JL)
           
          IF (myCoords(2)==0) THEN
            var(i,1,k) = var(i,1,k)*(1.0_CGREAL-omega) + omega*dummy(1)          
          ELSE                 
            var(i,1,k) = var(i,1,k)*(1.0_CGREAL-0.75 ) + 0.75 *dummy(1)          
          END IF   

          DO j=2,JL-1
           var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega) + omega*dummy(j)
          ENDDO ! j

          IF (myCoords(2)==Np(2)-1) THEN
            var(i,JL,k) = var(i,JL,k)*(1.0_CGREAL-omega) + omega*dummy(JL)
          ELSE
            var(i,JL,k) = var(i,JL,k)*(1.0_CGREAL-0.75 ) + 0.75 *dummy(JL)
          END IF


        ENDDO ! i
      ENDDO ! k

#     ifdef MPI      
        CALL par_comm_var(var,IL,JL,KL,glN,myTime)
        myCommTime=myCommTime+myTime
#     endif

    ENDDO ! Ncolor

!   -------------------------------
!    Line solve in the z-direction
!   -------------------------------
    IF (ndim == DIM_3D) THEN

!     -------------------------------------
!      Impose periodic boundary conditions
!     -------------------------------------
!!      IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN
!!        IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
!!            bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
!!            bcz1 .EQ. BC_TYPE_PERIODIC) THEN
!!         CALL apply_periodic_pres(var)
!!        END IF                 
!!      ENDIF   

      IF (nLevX==1 .OR. nLevY==1 .OR. nLevZ==1) THEN             !Ehsan added for Periodic BC
        IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &               
            bcy1 .EQ. BC_TYPE_PERIODIC .OR. &                 
            bcz1 .EQ. BC_TYPE_PERIODIC) THEN                  
         do i=0,Il+1 !1-glN,IL+glN
         do j=0,JL+1  !1-glN,JL+glN
         var(i,j,0) = var(i,j,KL)
         var(i,j,KL+1) = var(i,j,1)
         enddo
         enddo
        END IF                                                   !Ehsan added for Periodic BC
      ENDIF                                       

      DO Ncolor = 1, TNcolorZ 
        jinit = Ncolor
        iinit = Ncolor  

        DO j = jinit, JL, jstep 
          jG=MGY(nLevY)%L2GJ(j)

          DO i = iinit, IL, istep 
            iG=MGX(nLevX)%L2GI(i)

            DO k = 1, KL
              kG=k

              amx(i) =   MGX(nlevX)%dxcinv(iG  )*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - ium_mg(i,j,k)*Kill4GCMxm(i))*(AREA_W(i,j,k)*mglevelflag + 1.-mglevelflag)
              apx(i) =   MGX(nlevX)%dxcinv(iG+1)*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - iup_mg(i,j,k)*Kill4GCMxp(i))*(AREA_E(i,j,k)*mglevelflag + 1.-mglevelflag)
              acx(i) = - ( amx(i) + apx(i) )

              amy(j) =   MGY(nlevY)%dycinv(jG  )*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jum_mg(i,j,k)*Kill4GCMym(j))*(AREA_S(i,j,k)*mglevelflag + 1.-mglevelflag)
              apy(j) =   MGY(nlevY)%dycinv(jG+1)*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jup_mg(i,j,k)*Kill4GCMyp(j))*(AREA_N(i,j,k)*mglevelflag + 1.-mglevelflag)
              acy(j) = - ( amy(j) + apy(j) )

              amz(k) =   MGZ(nlevZ)%dzcinv(kG  )*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kum_mg(i,j,k)*Kill4GCMzm(k))*(AREA_B(i,j,k)*mglevelflag + 1.-mglevelflag)
              apz(k) =   MGZ(nlevZ)%dzcinv(kG+1)*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kup_mg(i,j,k)*Kill4GCMzp(k))*(AREA_F(i,j,k)*mglevelflag + 1.-mglevelflag)
              acz(k) = - ( amz(k) + apz(k) )

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
              rhs(k) = r(i,j,k)  - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k) ) &
                                 - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k) ) &
                                 - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) ) &
                                 - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) ) &
                                 + mgLevelFlag *                                     &
                               ( - pGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag      &
                                 - pGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag      &
                                 - pGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag      &
                                 - pGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag      &
                                 - pGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag      &
                                 - pGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag      &
                                 + boundPresSource(i-1,j,k)                          &
                                 + boundPresSource(i+1,j,k)                          &
                                 + boundPresSource(i,j-1,k)                          &
                                 + boundPresSource(i,j+1,k)                          &
                                 + boundPresSource(i,j,k-1)                          &
                                 + boundPresSource(i,j,k+1) )
          CASE(index_scalar)
              rhs(k) = r(i,j,k)  - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k) ) &
                                 - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k) ) &
                                 - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) ) &
                                 - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) ) &
                                 + mgLevelFlag *                                     &
                               ( - XGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag      &
                                 - XGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag      &
                                 - XGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag      &
                                 - XGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag      &
                                 - XGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag      &
                                 - XGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag      &
                                 + boundScalarSource(i-1,j,k)                        &
                                 + boundScalarSource(i+1,j,k)                        &
                                 + boundScalarSource(i,j-1,k)                        &
                                 + boundScalarSource(i,j+1,k)                        &
                                 + boundScalarSource(i,j,k-1)                        &
                                 + boundScalarSource(i,j,k+1) )   
          CASE(index_potential)
              rhs(k) = r(i,j,k)  - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum_mg(i,j,k) ) &
                                 - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup_mg(i,j,k) ) &
                                 - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium_mg(i,j,k) ) &
                                 - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup_mg(i,j,k) ) &
                                 + mgLevelFlag *                                     &
                               ( - PhiGhost(i-1,j,k)*amx(i)*ium_mg(i,j,k)*gcmFlag    &
                                 - PhiGhost(i+1,j,k)*apx(i)*iup_mg(i,j,k)*gcmFlag    &
                                 - PhiGhost(i,j-1,k)*amy(j)*jum_mg(i,j,k)*gcmFlag    &
                                 - PhiGhost(i,j+1,k)*apy(j)*jup_mg(i,j,k)*gcmFlag    &
                                 - PhiGhost(i,j,k-1)*amz(k)*kum_mg(i,j,k)*gcmFlag    &
                                 - PhiGhost(i,j,k+1)*apz(k)*kup_mg(i,j,k)*gcmFlag    &
                                 + boundPhiSource(i-1,j,k)                           &
                                 + boundPhiSource(i+1,j,k)                           &
                                 + boundPhiSource(i,j-1,k)                           &
                                 + boundPhiSource(i,j+1,k)                           &
                                 + boundPhiSource(i,j,k-1)                           &
                                 + boundPhiSource(i,j,k+1) )                                                                 
          END SELECT
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)           
              amz(k) = amz(k)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL - kum_mg(i,j,k))
              apz(k) = apz(k)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)*(1.0_CGREAL - kup_mg(i,j,k))
              acz(k) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank_MG(i,j,k),KIND=CGREAL) &
                                             +REAL(  iblank_MG(i,j,k),KIND=CGREAL)       
              rhs(k) = rhs(k)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
            ENDDO ! k

!           Impose Periodic BC
!           ------------------
!!            IF ( bcz1 == BC_TYPE_PERIODIC .AND. &
!!                 bcz2 == BC_TYPE_PERIODIC .AND. &
!!                 nLevZ == 1       ) THEN
!!               rhs(1)  = rhs(1)  - var(i,j,KL) *amz(1)
!!               rhs(KL) = rhs(KL) - var(i,j,1)  *apz(KL)
!!            END IF !  bcz1

            IF ( bcz1 == BC_TYPE_PERIODIC .AND. &              !Ehsan added for Periodic BC
                 bcz2 == BC_TYPE_PERIODIC .AND. &           
                 nLevZ == 1       ) THEN                      
               rhs(1)  = rhs(1)  - var(i,j,KL) *amz(1)        
               rhs(KL) = rhs(KL) - var(i,j,1)  *apz(KL)   
            ELSE
!           Impose Subdomain/Dirichlet BC (SAMK)
!           ------------------------------------
               rhs(1) =rhs(1) -var(i,j,0)   *amz(1 )                           
               rhs(KL)=rhs(KL)-var(i,j,KL+1)*apz(KL)                          
            END IF !  bcz1                                    

!           Impose Subdomain/Dirichlet BC (SAMK)
!           ------------------------------------
!            rhs(1) =rhs(1) -var(i,j,0)   *amz(1 )                           !Ehsan removed
!            rhs(KL)=rhs(KL)-var(i,j,KL+1)*apz(KL)                           !Ehsan removed


!         Solve TDMA in Z Direction
!         -------------------------
            CALL tdma(amz,acz,apz,rhs(1:KL),dummy(1:KL),1,KL)

            DO k=1,KL
              var(i,j,k) =  var(i,j,k)*(1.0_CGREAL-omega) + omega*dummy(k) 
            ENDDO ! k
          ENDDO ! i 
        ENDDO ! j

#       ifdef MPI      
          CALL par_comm_var(var,IL,JL,KL,glN,myTime)
          myCommTime=myCommTime+myTime
#       endif

      END DO ! Ncolor

    ENDIF ! ndim

END SUBROUTINE MG_itsolv
!---------------------------------------------------------------------



SUBROUTINE MG_residual(var,rhs,res,nLevX,nLevY,nLevZ,IL,JL,KL,glN,resCheckFlag) 

! ---------------------------------------------------------------------------------
!  Purpose: Calculate the residual at each level. 
!           Prepare for the RESTRICTION step in MG method
!
!  Input: all grids information such as ium_mg, iup_mg, iblank_MG
!         nLevX, nLevY, nLevZ: current level in x, y, z direction respectively
!         IL, JL, KL:  the number of grid points in x, y, z direction respectively
!                      at current level (used for semi-coarsening MG method)
!         var: initial guesses
!         rhs: values at the right-hand side
!         resCheckFlag: flag that activates residual checking
!  
!  Output: res -- storing the residual at current level
! ---------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE flow_arrays
    USE MG_arrays
	USE cutcell_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE scalar  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 	

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: nLevX, nLevY, nLevZ, IL, JL, KL, glN
    
    logical, intent(in) :: resCheckFlag

    REAL(KIND=CGREAL),INTENT (INOUT)  :: var(1-glN:IL+glN,1-glN:JL+glN,0:KL+1)     !Ehsan changed IN to INOUT for periodic BC
    REAL(KIND=CGREAL),INTENT (IN)  :: rhs(IL,JL,KL)
    REAL(KIND=CGREAL),INTENT (OUT) :: res(IL,JL,KL)

!   Local variables
!   ---------------
    INTEGER :: i  , j   ,k
    INTEGER :: iG , jG  ,kG

    REAL(KIND=CGREAL) :: bmx,bpx,bcx
    REAL(KIND=CGREAL) :: bmy,bpy,bcy
    REAL(KIND=CGREAL) :: bmz,bpz,bcz
    REAL(KIND=CGREAL) :: bc
    REAL(KIND=CGREAL) :: mgLevelFlag


    REAL(KIND=CGREAL), DIMENSION (IL) :: kill4GCMxm, kill4GCMxp
    REAL(KIND=CGREAL), DIMENSION (JL) :: kill4GCMym, kill4GCMyp
    REAL(KIND=CGREAL), DIMENSION (KL) :: kill4GCMzm, kill4GCMzp

    REAL(KIND=CGREAL) :: rrr, resMax
    REAL(KIND=CGREAL) :: resL(IL,JL,KL)



!   -------------------------------------------------------------------------
!   Initialize variables
!   Notes:
!    1. mgLevelFlag activates pGhost and boundPresSource for first level only
!    2. for MG Level higher than 1, PPE system is solved using SSM whether
!       GCM or SSM method is invoked.
!   -------------------------------------------------------------------------

    resMax = 0.0_CGREAL
    
    IF ( nLevX*nLevY*nLevZ == 1 ) THEN
      mgLevelFlag = 1.0_CGREAL
    ELSE
      mgLevelFlag = 0.0_CGREAL
    END IF

!   Initializing Neumann/Dirichlet BC flags (SAMK/RM)
!   -------------------------------------------------
    Kill4GCMxm(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMxp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    IF (myCoords(1)==0      ) Kill4GCMxm(1 ) = 1.0_CGREAL
    IF (myCoords(1)==Np(1)-1) Kill4GCMxp(IL) = 1.0_CGREAL

    Kill4GCMym(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMyp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    IF (myCoords(2)==0      ) Kill4GCMym(1 ) = 1.0_CGREAL
    IF (myCoords(2)==Np(2)-1) Kill4GCMyp(JL) = 1.0_CGREAL

    Kill4GCMzm(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMzp(: ) = 1.0_CGREAL - gcmFlag*mgLevelFlag*(1-iSSMP*iCC)
    Kill4GCMzm(1 ) = 1.0_CGREAL
    Kill4GCMzp(KL) = 1.0_CGREAL

!   -----------------------------------------
!   Impose periodic boundary conditions
!   Notes: 
!     1. Probably broken-FMN 21 February 2007
!   -----------------------------------------
!!    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. & 
!!        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
!!        bcz1 .EQ. BC_TYPE_PERIODIC) THEN 
!!       CALL enforce_p_periodic(var) 
!!    END IF ! bcx1

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &            !Ehsan added for Periodic BC
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &            
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN         
         do i=0,IL+1 !1-glN,IL+glN              
           do j=0,JL+1 !1-glN,JL+glN           
             var(i,j,0) = var(i,j,KL)        
             var(i,j,KL+1) = var(i,j,1)     
           enddo                         
         enddo                       
    END IF ! bcx1                


!   Compute residual on current MG level
!   ------------------------------------
    DO k = 1, KL
      kG=k
            
      DO j = 1, JL
        jG=MGY(nLevY)%L2GJ(j)

        DO i = 1, IL
          iG=MGX(nLevX)%L2GI(i)

          bmx = MGX(nlevX)%dxcinv(iG  )*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - ium_mg(i,j,k)*Kill4GCMxm(i))*(AREA_W(i,j,k)*mglevelflag + 1.-mglevelflag)
          bpx = MGX(nlevX)%dxcinv(iG+1)*MGX(nlevX)%dxinv(iG)*(1.0_CGREAL - iup_mg(i,j,k)*Kill4GCMxp(i))*(AREA_E(i,j,k)*mglevelflag + 1.-mglevelflag)
          bcx = - ( bmx + bpx )

          bmy = MGY(nlevY)%dycinv(jG  )*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jum_mg(i,j,k)*Kill4GCMym(j))*(AREA_S(i,j,k)*mglevelflag + 1.-mglevelflag)
          bpy = MGY(nlevY)%dycinv(jG+1)*MGY(nlevY)%dyinv(jG)*(1.0_CGREAL - jup_mg(i,j,k)*Kill4GCMyp(j))*(AREA_N(i,j,k)*mglevelflag + 1.-mglevelflag)
          bcy = - ( bmy + bpy )

          bmz = MGZ(nlevZ)%dzcinv(kG  )*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kum_mg(i,j,k)*Kill4GCMzm(k))*KillFor2D*(AREA_B(i,j,k)*mglevelflag + 1.-mglevelflag)
          bpz = MGZ(nlevZ)%dzcinv(kG+1)*MGZ(nlevZ)%dzinv(kG)*(1.0_CGREAL - kup_mg(i,j,k)*Kill4GCMzp(k))*KillFor2D*(AREA_F(i,j,k)*mglevelflag + 1.-mglevelflag)
          bcz = - ( bmz + bpz )

!         Modify the nature of the equation for non-fluid cells
!         -----------------------------------------------------
          bmx = bmx*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
          bpx = bpx*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)

          bmy = bmy*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
          bpy = bpy*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)

          bmz = bmz*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
          bpz = bpz*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)

          bc  = (bcx+bcy+bcz)*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)  &
                             +REAL(  iblank_MG(i,j,k),KIND=CGREAL)                  


!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
          SELECT CASE(POSSION_INDEX)
          CASE(index_pressure)
          rrr = rhs(i,j,k) - var(i,j,k)*bc                                   &
                           - var(i-1,j,k)*bmx*(1.0_CGREAL - ium_mg(i,j,k) )  &
                           - var(i+1,j,k)*bpx*(1.0_CGREAL - iup_mg(i,j,k) )  &
                           - var(i,j-1,k)*bmy*(1.0_CGREAL - jum_mg(i,j,k) )  &
                           - var(i,j+1,k)*bpy*(1.0_CGREAL - jup_mg(i,j,k) )  &
                           - var(i,j,k-1)*bmz*(1.0_CGREAL - kum_mg(i,j,k) )  &
                           - var(i,j,k+1)*bpz*(1.0_CGREAL - kup_mg(i,j,k) )  & 
                           + mgLevelFlag *                                   &
                           ( - pGhost(i-1,j,k)*bmx*ium_mg(i,j,k)*gcmFlag     &
                             - pGhost(i+1,j,k)*bpx*iup_mg(i,j,k)*gcmFlag     &
                             - pGhost(i,j-1,k)*bmy*jum_mg(i,j,k)*gcmFlag     &
                             - pGhost(i,j+1,k)*bpy*jup_mg(i,j,k)*gcmFlag     &
                             - pGhost(i,j,k-1)*bmz*kum_mg(i,j,k)*gcmFlag     &
                             - pGhost(i,j,k+1)*bpz*kup_mg(i,j,k)*gcmFlag     &
                             + boundPresSource(i-1,j,k)                      &
                             + boundPresSource(i+1,j,k)                      &
                             + boundPresSource(i,j-1,k)                      &
                             + boundPresSource(i,j+1,k)                      &
                             + boundPresSource(i,j,k-1)                      &
                             + boundPresSource(i,j,k+1) )    
                             
          CASE(index_scalar)
          rrr = rhs(i,j,k) - var(i,j,k)*bc                                   &
                           - var(i-1,j,k)*bmx*(1.0_CGREAL - ium_mg(i,j,k) )  &
                           - var(i+1,j,k)*bpx*(1.0_CGREAL - iup_mg(i,j,k) )  &
                           - var(i,j-1,k)*bmy*(1.0_CGREAL - jum_mg(i,j,k) )  &
                           - var(i,j+1,k)*bpy*(1.0_CGREAL - jup_mg(i,j,k) )  &
                           - var(i,j,k-1)*bmz*(1.0_CGREAL - kum_mg(i,j,k) )  &
                           - var(i,j,k+1)*bpz*(1.0_CGREAL - kup_mg(i,j,k) )  & 
                           + mgLevelFlag *                                   &
                           ( - XGhost(i-1,j,k)*bmx*ium_mg(i,j,k)*gcmFlag     &
                             - XGhost(i+1,j,k)*bpx*iup_mg(i,j,k)*gcmFlag     &
                             - XGhost(i,j-1,k)*bmy*jum_mg(i,j,k)*gcmFlag     &
                             - XGhost(i,j+1,k)*bpy*jup_mg(i,j,k)*gcmFlag     &
                             - XGhost(i,j,k-1)*bmz*kum_mg(i,j,k)*gcmFlag     &
                             - XGhost(i,j,k+1)*bpz*kup_mg(i,j,k)*gcmFlag     &
                             + boundScalarSource(i-1,j,k)                    &
                             + boundScalarSource(i+1,j,k)                    &
                             + boundScalarSource(i,j-1,k)                    &
                             + boundScalarSource(i,j+1,k)                    &
                             + boundScalarSource(i,j,k-1)                    &
                             + boundScalarSource(i,j,k+1) )      
          CASE(index_potential)
          rrr = rhs(i,j,k) - var(i,j,k)*bc                                   &
                           - var(i-1,j,k)*bmx*(1.0_CGREAL - ium_mg(i,j,k) )  &
                           - var(i+1,j,k)*bpx*(1.0_CGREAL - iup_mg(i,j,k) )  &
                           - var(i,j-1,k)*bmy*(1.0_CGREAL - jum_mg(i,j,k) )  &
                           - var(i,j+1,k)*bpy*(1.0_CGREAL - jup_mg(i,j,k) )  &
                           - var(i,j,k-1)*bmz*(1.0_CGREAL - kum_mg(i,j,k) )  &
                           - var(i,j,k+1)*bpz*(1.0_CGREAL - kup_mg(i,j,k) )  & 
                           + mgLevelFlag *                                   &
                           ( - PhiGhost(i-1,j,k)*bmx*ium_mg(i,j,k)*gcmFlag   &
                             - PhiGhost(i+1,j,k)*bpx*iup_mg(i,j,k)*gcmFlag   &
                             - PhiGhost(i,j-1,k)*bmy*jum_mg(i,j,k)*gcmFlag   &
                             - PhiGhost(i,j+1,k)*bpy*jup_mg(i,j,k)*gcmFlag   &
                             - PhiGhost(i,j,k-1)*bmz*kum_mg(i,j,k)*gcmFlag   &
                             - PhiGhost(i,j,k+1)*bpz*kup_mg(i,j,k)*gcmFlag   &
                             + boundPhiSource(i-1,j,k)                       &
                             + boundPhiSource(i+1,j,k)                       &
                             + boundPhiSource(i,j-1,k)                       &
                             + boundPhiSource(i,j+1,k)                       &
                             + boundPhiSource(i,j,k-1)                       &
                             + boundPhiSource(i,j,k+1) )                                         
          
          END SELECT   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)          
                                               

          resL(i,j,k) = rrr*REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
       
        ENDDO ! i
      ENDDO ! j
    ENDDO ! k

    IF ( resCheckFlag ) res(1:IL,1:JL,1:KL) = resL(1:IL,1:JL,1:KL)

!   ---------------------------------------------
!   Output convergence history
!   Note: 
!    1. This should be turned off for normal runs
!      (set infoconv == 0 in input.dat)
!   ---------------------------------------------
!!    IF (infoconv == 1) THEN
!!      resMax = MAXVAL(ABS(resL))
!!      WRITE(STDOUT,"('MG: residual check : ',1X,3I6,2X,2E19.11)") nLevX, nLevY, nLevZ, resmax
!!    END IF ! infoconv

END SUBROUTINE MG_Residual
!-------------------------------------------------------------------------------



SUBROUTINE MG_Restrict(nlev,NC_dim)

! ------------------------------------------------------------------------------------
!  Purpose: Restriction step in multigrid method.
!
! Input: grids information such as iblank_MG, ghostcellMark_mg.
!        nlev -- current MG level
!        NC_dim  --  could be 1, 2, 3 corresponding to x, y, z direction respectively
!        var  -- input residuals of finer level
!  
! Output:  r   -- output residuals of coaser level
! ------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE MG_arrays
 
    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: nlev, NC_dim

!   Local variables
!   ---------------
    INTEGER  :: i, ii , iim , iip
    INTEGER  ::    iiG, iimG, iipG

    INTEGER  :: j, jj , jjm , jjp
    INTEGER  ::    jjG, jjmG, jjpG

    INTEGER  :: k, kk , kkm , kkp
    INTEGER  ::    kkG, kkmG, kkpG

    REAL(KIND=CGREAL) :: area1, area2, area
    REAL(KIND=CGREAL) :: sum
    real(kind=CGREAL) :: KillForOdd
    
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc) ::  diblank_MG


     
!   Start Restriction
!   -----------------
    SELECT CASE( NC_dim )
    
!   Restrict in x-direction
!   -----------------------
    CASE (ICOORD)
      DO k=1,nzc
      DO j=1,nyc
      DO i=1,MGX(nlev)%nxc
        diblank_MG(i,j,k) = REAL(1-iblank_MG(i,j,k),KIND=CGREAL) 
      END DO ! i
      END DO ! j
      END DO ! k

      DO k=1,nzc
        DO j=1,nyc
          DO i=1,MGX(nlev+1)%nxc-1
            ii  = 2*i
            iim = 2*i-1
            iip = 2*i+1

            iiG  = L2GI_MG(nlev,ii )
            iimG = L2GI_MG(nlev,iim)
            iipG = L2GI_MG(nlev,iip)

            area1 = MGX(nlev)%x(iiG) -MGX(nlev)%x(iimG)
            area2 = MGX(nlev)%x(iipG)-MGX(nlev)%x(iiG )

            sum = diblank_MG(iim,j,k) *area1 *MGX(nlev)%res(iim,j,k) &
                + diblank_MG(ii ,j,k) *area2 *MGX(nlev)%res(ii ,j,k)

            area = diblank_MG(ii ,j,k) *area2  &
                 + diblank_MG(iim,j,k) *area1 

            IF ( area > 0.0_CGREAL ) THEN 
              MGX(nlev+1)%rhs(i,j,k) = sum/area
            ELSE 
              MGX(nlev+1)%rhs(i,j,k) = sum
            END IF ! area
          END DO ! i

!         Treating Odd number of grids, Added by SAMK
!         -------------------------------------------
          i   = MGX(nlev+1)%nxc
          ii  = 2*i
          iim = 2*i-1
          iip = 2*i+1

          iiG  = L2GI_MG(nlev,ii )
          iimG = L2GI_MG(nlev,iim)
          iipG = L2GI_MG(nlev,iip)

          KillForOdd=real(1-(ii-MGX(nlev)%nxc), kind=CGREAL)

!         Prevents ii to exceed the limit (Also keeps an original definition for iiG)
!         ---------------------------------------------------------------------------
          ii = (2*i)*(1-(2*i-MGX(nlev)%nxc)) + MGX(nlev)%nxc*(2*i-MGX(nlev)%nxc)
          
          area1 =  MGX(nlev)%x(iiG) -MGX(nlev)%x(iimG)
          area2 = (MGX(nlev)%x(iipG)-MGX(nlev)%x(iiG ))*KillForOdd

          sum = diblank_MG(iim,j,k) *area1 *MGX(nlev)%res(iim,j,k) &
              + diblank_MG(ii ,j,k) *area2 *MGX(nlev)%res(ii ,j,k)

          area = diblank_MG(ii ,j,k) *area2  &
               + diblank_MG(iim,j,k) *area1 

          IF ( area > 0.0_CGREAL ) THEN 
            MGX(nlev+1)%rhs(i,j,k) = sum/area
          ELSE 
            MGX(nlev+1)%rhs(i,j,k) = sum
          END IF ! area
        END DO ! j
      END DO ! k

!   Restrict in y-direction
!   -----------------------
    CASE (JCOORD)
      DO k=1,nzc
      DO j=1,MGY(nlev)%nyc
      DO i=1,nxc
        diblank_MG(i,j,k) = REAL(1-iblank_MG(i,j,k),KIND=CGREAL) 
      END DO ! i
      END DO ! j
      END DO ! k

      DO k=1,nzc
        DO j=1,MGY(nlev+1)%nyc-1
          jj  = 2*j
          jjm = 2*j-1
          jjp = 2*j+1

          jjG  = L2GJ_MG(nlev,jj )
          jjmG = L2GJ_MG(nlev,jjm)
          jjpG = L2GJ_MG(nlev,jjp)

          area1 = MGY(nlev)%y(jjG )-MGY(nlev)%y(jjmG)
          area2 = MGY(nlev)%y(jjpG)-MGY(nlev)%y(jjG )

          DO i=1,nxc
            sum = diblank_MG(i,jjm,k) *area1 *MGY(nlev)%res(i,jjm,k) &
                + diblank_MG(i,jj ,k) *area2 *MGY(nlev)%res(i,jj ,k)

            area = diblank_MG(i,jj ,k) *area2  &
                 + diblank_MG(i,jjm,k) *area1 

            IF ( area > 0.0_CGREAL ) THEN 
              MGY(nlev+1)%rhs(i,j,k) = sum/area
            ELSE 
              MGY(nlev+1)%rhs(i,j,k) = sum
            END IF ! area
          END DO ! i
        END DO ! j

!       Treating Odd number of grids, Added by SAMK
!       -------------------------------------------
        j   = MGY(nlev+1)%nyc
        jj  = 2*j
        jjm = 2*j-1
        jjp = 2*j+1

        jjG  = L2GJ_MG(nlev,jj )
        jjmG = L2GJ_MG(nlev,jjm)
        jjpG = L2GJ_MG(nlev,jjp)

        KillForOdd=real(1-(jj-MGY(nlev)%nyc), kind=CGREAL)

!       Prevents jj to exceed the limit (Also keeps an original definition for jjG)
!       ---------------------------------------------------------------------------        
        jj = (2*j)*(1-(2*j-MGY(nlev)%nyc)) + MGY(nlev)%nyc*(2*j-MGY(nlev)%nyc)

        area1 = MGY(nlev)%y(jjG )-MGY(nlev)%y(jjmG)
        area2 = MGY(nlev)%y(jjpG)-MGY(nlev)%y(jjG )*KillForOdd

        DO i=1,nxc
          sum = diblank_MG(i,jjm,k) *area1 *MGY(nlev)%res(i,jjm,k) &
              + diblank_MG(i,jj ,k) *area2 *MGY(nlev)%res(i,jj ,k)

          area = diblank_MG(i,jj ,k) *area2  &
               + diblank_MG(i,jjm,k) *area1 

          IF ( area > 0.0_CGREAL ) THEN 
            MGY(nlev+1)%rhs(i,j,k) = sum/area
          ELSE 
            MGY(nlev+1)%rhs(i,j,k) = sum
          END IF ! area
        END DO ! i
      END DO ! k

!   Restrict in z-direction
!   -----------------------
    CASE (KCOORD)
      DO k = 1, MGZ(nlev)%nzc
      DO j = 1, nyc
      DO i = 1, nxc
        diblank_MG(i,j,k) = REAL(1-iblank_MG(i,j,k),KIND=CGREAL)
      END DO ! i
      END DO ! j
      END DO ! k

      DO k = 1, MGZ(nlev+1)%nzc-1
        kk  = 2*k
        kkm = 2*k-1
        kkp = 2*k+1
        
        kkG  = kk
        kkmG = kkm
        kkpG = kkp

        area1 = MGZ(nlev)%z(kkG )-MGZ(nlev)%z(kkmG)
        area2 = MGZ(nlev)%z(kkpG)-MGZ(nlev)%z(kkG )

        DO j = 1, nyc
          DO i = 1, nxc
            sum = diblank_MG(i,j,kkm) *area1 *MGZ(nlev)%res(i,j,kkm) &
                + diblank_MG(i,j,kk ) *area2 *MGZ(nlev)%res(i,j,kk )

            area = diblank_MG(i,j,kkm) *area1 + &
                 + diblank_MG(i,j,kk ) *area2 

            IF ( area > 0.0_CGREAL ) THEN    
              MGZ(nlev+1)%rhs(i,j,k) = sum/area
            ELSE 
              MGZ(nlev+1)%rhs(i,j,k) = sum
            END IF ! area
          END DO ! i
        END DO ! j
      END DO ! k

      k   = MGZ(nlev+1)%nzc
      kk  = 2*k
      kkm = 2*k-1
      kkp = 2*k+1
      
      kkG  = kk
      kkmG = kkm
      kkpG = kkp

      KillForOdd=real(1-(kk-MGZ(nlev)%nzc), kind=CGREAL)

!     Prevents kk to exceed the limit (Also keeps an original definition for kkG)
!     ---------------------------------------------------------------------------
      kk = (2*k)*(1-(2*k-MGZ(nlev)%nzc)) + MGZ(nlev)%nzc*(2*k-MGZ(nlev)%nzc)
          
      area1 =  MGZ(nlev)%z(kkG )-MGZ(nlev)%z(kkmG)
      area2 = (MGZ(nlev)%z(kkpG)-MGZ(nlev)%z(kkG ))*KillForOdd

      DO j = 1, nyc
        DO i = 1, nxc
          sum = diblank_MG(i,j,kkm) *area1 *MGZ(nlev)%res(i,j,kkm) &
              + diblank_MG(i,j,kk ) *area2 *MGZ(nlev)%res(i,j,kk )

          area = diblank_MG(i,j,kkm) *area1 + &
               + diblank_MG(i,j,kk ) *area2 

          IF ( area > 0.0_CGREAL ) THEN    
            MGZ(nlev+1)%rhs(i,j,k) = sum/area
          ELSE 
            MGZ(nlev+1)%rhs(i,j,k) = sum
          END IF ! area
        END DO ! i
      END DO ! j

    END SELECT ! NC_dim

END SUBROUTINE MG_Restrict
!---------------------------------------------------------------------



SUBROUTINE MG_Prolong(nlev,NC_dim)

! ----------------------------------------------------------------------
!  Purpose: Prolongation step of MG method.
!
!           A two-point Lagrange interpolation method is used in this
!           subroutine. (SAMK)
!
!  Input:  correc: errors at coarser level 
!          nLev  : Current level (finer level)
!          nC_dim: direction (1:x, 2:y, 3:z)
!  
!  Output:  phi: errors at finer level 
! ----------------------------------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE MG_arrays

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: nlev, NC_dim
 
!   Local variables
!   ---------------
    INTEGER :: i , icp , icm
    INTEGER :: iG, icpG, icmG

    INTEGER :: j , jcp , jcm
    INTEGER :: jG, jcpG, jcmG

    INTEGER :: k , kcp , kcm

    integer :: iErr

    REAL(KIND=CGREAL) :: delta1, delta2, delta
    REAL(KIND=CGREAL) :: phi


    
!   Start Prolongation
!   ------------------
    SELECT CASE( NC_dim )

!   Prolongate in x-direction
!   -------------------------
    CASE (ICOORD)

      IF (myCoords(1)==0) THEN
        IF (pbcx1 == PBC_DIRICHLET) THEN 

!         Zero error condition
!         --------------------
          MGX(nlev+1)%phi(0,:,:)=0.0_CGREAL

        ELSE

!         Zero error-variation condition
!         ------------------------------
          MGX(nlev+1)%phi(0,:,:)=MGX(nlev+1)%phi(1,:,:)
        END IF
      END IF

      IF (myCoords(1)==Np(1)-1) THEN
        IF (pbcx2 == PBC_DIRICHLET) THEN 

!         Zero error condition
!         --------------------
          MGX(nlev+1)%phi(MGX(nlev+1)%nxc+1,:,:)=0.0_CGREAL

        ELSE

!         Zero error-variation condition
!         ------------------------------
          MGX(nlev+1)%phi(MGX(nlev+1)%nxc+1,:,:)= MGX(nlev+1)%phi(MGX(nlev+1)%nxc,:,:)
        END IF
      END IF

      DO k=1,nzc
      DO j=1,nyc
      DO i=1,MGX(nlev)%nxc
        icp = (i+2)/2
        icm = icp-1
        
        iG   = L2GI_MG(nlev  ,i  )
        icpG = L2GI_MG(nlev+1,icp)
        icmG = L2GI_MG(nlev+1,icm)
        
!       prolongation
!       ------------
        delta1 = (MGX(nlev+1)%xc(icpG) -MGX(nlev  )%xc(iG  )) &
               * REAL(1-MGX(nlev+1)%iblank(icm,j,k),KIND=CGREAL)

        delta2 = (MGX(nlev  )%xc(iG  ) -MGX(nlev+1)%xc(icmG)) &
               * REAL(1-MGX(nlev+1)%iblank(icp,j,k),KIND=CGREAL)

        delta  = delta1+delta2

        phi = delta1*MGX(nlev+1)%phi(icm,j,k) &
            + delta2*MGX(nlev+1)%phi(icp,j,k)

        IF ( delta > 0.0_CGREAL ) THEN
          phi = phi/delta
        END IF ! delta
        
        MGX(nlev)%phi(i,j,k) = (MGX(nlev)%phi(i,j,k) + phi) &
                             * REAL(1-MGX(nlev)%iblank(i,j,k),KIND=CGREAL)
      END DO ! ic
      END DO ! j
      END DO ! k

!   Prolongate in y-direction
!   -------------------------
    CASE (JCOORD)

      IF (myCoords(2)==0) THEN 
        IF (pbcy1 == PBC_DIRICHLET) THEN 

!         Zero error condition
!         --------------------
          MGY(nlev+1)%phi(:,0,:)=0.0_CGREAL

        ELSE

!         Zero error-variation condition
!         ------------------------------
          MGY(nlev+1)%phi(:,0,:)=MGY(nlev+1)%phi(:,1,:)
        END IF
      END IF

      IF (myCoords(2)==Np(2)-1) THEN 
        IF (pbcy2 == PBC_DIRICHLET) THEN 

!         Zero error condition
!         --------------------
          MGY(nlev+1)%phi(:,MGY(nlev+1)%nyc+1,:)=0.0_CGREAL

        ELSE

!         Zero error-variation condition
!         ------------------------------
          MGY(nlev+1)%phi(:,MGY(nlev+1)%nyc+1,:)=MGY(nlev+1)%phi(:,MGY(nlev+1)%nyc,:)
        END IF
      END IF

      DO k = 1, nzc
      DO i = 1, nxc
      DO j = 1, MGY(nlev)%nyc
        jcp = (j+2)/2
        jcm = jcp-1

        jG   = L2GJ_MG(nlev  ,j  )
        jcpG = L2GJ_MG(nlev+1,jcp)
        JcmG = L2GJ_MG(nlev+1,jcm)
        
!       prolongation
!       ------------
        delta1 = (MGY(nlev+1)%yc(jcpG) -MGY(nlev  )%yc(jG  )) &
               * REAL(1-MGY(nlev+1)%iblank(i,jcm,k),KIND=CGREAL)

        delta2 = (MGY(nlev  )%yc(jG  ) -MGY(nlev+1)%yc(jcmG)) &
               * REAL(1-MGY(nlev+1)%iblank(i,jcp,k),KIND=CGREAL)

        delta=delta1+delta2

        phi = delta1*MGY(nlev+1)%phi(i,jcm,k) &
            + delta2*MGY(nlev+1)%phi(i,jcp,k)

        IF ( delta > 0.0_CGREAL ) THEN
          phi =  phi/delta
        END IF ! delta

        MGY(nlev)%phi(i,j,k) = (MGY(nlev)%phi(i,j,k) + phi) &
                             * REAL(1-MGY(nlev)%iblank(i,j,k),KIND=CGREAL)
      END DO ! j
      END DO ! i  
      END DO ! k

!   Prolongate in z-direction
!   -------------------------
    CASE (KCOORD) 

      IF ( pbcz1 == PBC_DIRICHLET ) THEN 

!       Zero error condition
!       --------------------
        MGZ(nlev+1)%phi(:,:,0)=0.0_CGREAL

      ELSE

!       Zero error-variation condition
!       ------------------------------
        MGZ(nlev+1)%phi(:,:,0)= MGZ(nlev+1)%phi(:,:,1)
        
      END IF

      IF ( pbcz2 == PBC_DIRICHLET ) THEN 

!       Zero error condition
!       --------------------
        MGZ(nlev+1)%phi(:,:,MGZ(nlev+1)%nzc+1)=0.0_CGREAL

      ELSE

!       Zero error-variation condition
!       ------------------------------
        MGZ(nlev+1)%phi(:,:,MGZ(nlev+1)%nzc+1)= MGZ(nlev+1)%phi(:,:,MGZ(nlev+1)%nzc)
      END IF

      DO j = 1, nyc
      DO i = 1, nxc
      DO k = 1, MGZ(nlev)%nzc
        kcp = (k+2)/2
        kcm = kcp-1

!       prolongation
!       ------------
        delta1 = (MGZ(nlev+1)%zc(kcp) -MGZ(nlev  )%zc(k  )) &
               * REAL(1-MGZ(nlev+1)%iblank(i,j,kcm),KIND=CGREAL)

        delta2 = (MGZ(nlev  )%zc(k  ) -MGZ(nlev+1)%zc(kcm)) &
               * REAL(1-MGZ(nlev+1)%iblank(i,j,kcp),KIND=CGREAL)

        delta=delta1+delta2

        phi = delta1*MGZ(nlev+1)%phi(i,j,kcm) &
            + delta2*MGZ(nlev+1)%phi(i,j,kcp)


        IF ( delta > 0.0_CGREAL ) THEN
          phi =  phi/delta
        END IF ! delta
        
        MGZ(nlev)%phi(i,j,k) = (MGZ(nlev)%phi(i,j,k) + phi) &
                             * REAL(1-MGZ(nlev)%iblank(i,j,k),KIND=CGREAL)
      END DO ! k
      END DO ! i  
      END DO ! j
    
    END SELECT ! NC_dim 

END SUBROUTINE MG_Prolong
!---------------------------------------------------------------------
