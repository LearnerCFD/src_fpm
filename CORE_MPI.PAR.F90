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
!    Haibo Dong
!    Haoxiang Luo
!    Meliha Bozkurttas
!    Rupesh Babu K. A.
!    Xudong Zheng 
!    Reza Ghias 
!    S. A. Mohsen Karimian
!
!  Filename: CORE_MPI.PAR.F90
!  Latest Modification: Dec, 29 2010 (PAT 2.1.0)
!  by Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following modules:
!     par_init_proc
!     par_domain_decomp
!     par_init_domain
!     par_init_MG_domain
!     par_comm_outermost_int
!     par_comm_outermost_real
!     par_comm_POISSON
!     par_comm_MGX
!     par_comm_MGY
!     par_comm_MGZ
!     par_comm_var
!     par_decomp(nicG,dir,nic,ni,Is,Ie)
!     par_getMaxReal(myReal)
!     par_getSumReal(myReal)
!     par_getSumRealArray(myReal,n)
!     par_getMaxInteger(myInteger)
!     par_getSumInteger(myInteger)
!     par_bcast_marker
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)       myIs+i-1
# define L2GJ(j)       myJs+j-1



subroutine par_init_proc

! -----------------------------------------------------------
!  This subroutine initializes and configures the processes.
! -----------------------------------------------------------

    USE flow_parameters
    use par_parameters
    use mpi

    IMPLICIT NONE

!   Parallel constants
!   ------------------
    integer, parameter :: Xdir=0, Ydir=1
    logical, parameter :: isperiodic(2)=.false., reorder=.true.
    
!   Local parameters
!   ----------------
    INTEGER :: ierr


!   Creating a Cartisian domain
!   ---------------------------
    CALL mpi_cart_create(MPI_COMM_WORLD,2,Np,isperiodic,reorder,commCart, ierr)
    IF(ierr /= ERR_NONE) THEN
     WRITE(STDOUT, *)'Communication error for mpi_cart_create in the par_init_proc subroutine STOP'
     CALL flow_stop
     STOP
    ENDIF

    CALL mpi_comm_rank(commCart,myCartRank, ierr)
    IF(ierr /= ERR_NONE) THEN
     WRITE(STDOUT, *)'Communication error for mpi_comm_rank in the par_init_proc subroutine STOP'
     CALL flow_stop
     STOP
    ENDIF  

!   Determining my coordinates in the parallel grid
!   -----------------------------------------------
    CALL mpi_cart_coords(commCart,myCartRank,2,myCoords, ierr)
    IF(ierr /= ERR_NONE) THEN
     WRITE(STDOUT, *)'Communication error for mpi_cart_coords in the par_init_proc subroutine STOP'
     CALL flow_stop
     STOP
    ENDIF

    CALL mpi_cart_shift(commCart,Xdir,1,myWest,myEast, ierr)
    IF(ierr /= ERR_NONE) THEN
     WRITE(STDOUT, *)'Communication error for mpi_cart_shift in the par_init_proc subroutine STOP'
     CALL flow_stop
     STOP
    ENDIF

    CALL mpi_cart_shift(commCart,Ydir,1,mySouth,myNorth, ierr)
    IF(ierr /= ERR_NONE) THEN
     WRITE(STDOUT, *)'Communication error for mpi_cart_shif in the par_init_proc subroutine STOP'
     CALL flow_stop
     STOP
    ENDIF


end subroutine par_init_proc
!---------------------------------------------------------------------



SUBROUTINE par_domain_decomp

! ----------------------------------------------------------------
!  This subroutine performs domain decomposition among processes.
! ----------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays

    IMPLICIT NONE



!   2D domain decomposition
!   -----------------------
    IF (it_solver_type==IT_SOLVER_TYPE_MG) THEN

!     Project the first level of MG to the main domain
!     ------------------------------------------------
      nxc =MGX(1)%nxc
      nx  =MGX(1)%nx
      myIs=MGX(1)%myIs
      myIe=MGX(1)%myIe
      
      nyc =MGY(1)%nyc
      ny  =MGY(1)%ny
      myJs=MGY(1)%myJs
      myJe=MGY(1)%myJe
    ELSE
      CALL par_decomp(nxc_GLBL,ICOORD,nxc,nx,myIs,myIe)
      CALL par_decomp(nyc_GLBL,JCOORD,nyc,ny,myJs,myJe)

    END IF
    
!   Setting the iblank range and the iblank search range for each subdomain
!   -----------------------------------------------------------------------
    IF (myCoords(1)==0) THEN
      myILL=1
      myImin=MyILL
    ELSE
      myILL=1-Ngl
      myImin=myILL+1
    END IF

    IF (myCoords(1)==Np(1)-1) THEN
      myIUL=nxc
      myImax=myIUL
    ELSE
      myIUL=nxc+Ngl
      myImax=myIUL-1
    END IF

    IF (myCoords(2)==0) THEN
      myJLL=1
      myJmin=myJLL
    ELSE
      myJLL=1-Ngl
      myJmin=myJLL+1
    END IF

    IF (myCoords(2)==Np(2)-1) THEN
      myJUL=nyc
      myJmax=myJUL
    ELSE
      myJUL=nyc+Ngl
      myJmax=myJUL-1
    END IF

END SUBROUTINE par_domain_decomp
!---------------------------------------------------------------------



subroutine par_init_domain

! ----------------------------------------------------------------------------------
!  This subroutine initializes the MPI derived communication vectors in subdomains.
! ----------------------------------------------------------------------------------

    USE flow_parameters
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Local parameters
!   ----------------
    integer :: ierr



!   Communication vectors for main variables
!   ----------------------------------------
    CALL MPI_TYPE_VECTOR(nzc*(nyc+2*Ngl),Ngl,nxc+2*Ngl,mpi_double_precision,parVecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_domain subroutine STOP'
      CALL flow_stop
      STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(ParVecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_domain subroutine STOP'
      CALL flow_stop
      STOP
     ENDIF

    
    CALL MPI_TYPE_VECTOR(nzc,(nxc+2*Ngl)*Ngl,(nyc+2*Ngl)*(nxc+2*Ngl),mpi_double_precision,parVecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_domain subroutine STOP'
      CALL flow_stop
      STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(ParVecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_domain subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF
   
 
end subroutine par_init_domain
!---------------------------------------------------------------------



subroutine par_init_MG_domain

! ----------------------------------------------------------------------------------------------------
!  This subroutine initializes the MPI derived communication vectors for each MG level in subdomains.
! ----------------------------------------------------------------------------------------------------

    USE flow_parameters
    USE mg_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Local parameters
!   ----------------
    integer :: ilevel, glN, ierr



    DO ilevel=1,mgLevels_X

      SELECT CASE(ilevel)
      CASE(1)
        glN=Ngl
      CASE DEFAULT
        glN=1
      END SELECT 
      
      CALL MPI_TYPE_VECTOR(nzc*(nyc+2*glN),glN,MGX(ilevel)%nxc+2*glN,mpi_double_precision,MGX(ilevel)%parVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGX(ilevel)%ParVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
 
      CALL MPI_TYPE_VECTOR(nzc,(MGX(ilevel)%nxc+2*glN)*glN,(nyc+2*glN)*(MGX(ilevel)%nxc+2*glN),mpi_double_precision,MGX(ilevel)%parVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGX(ilevel)%ParVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      
      CALL MPI_TYPE_VECTOR(nzc*(nyc+2*glN),glN,MGX(ilevel)%nxc+2*glN,mpi_integer,MGX(ilevel)%parIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGX(ilevel)%ParIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
        WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
        CALL flow_stop
        STOP
       ENDIF

      CALL MPI_TYPE_VECTOR(nzc,(MGX(ilevel)%nxc+2*glN)*glN,(nyc+2*glN)*(MGX(ilevel)%nxc+2*glN),mpi_integer,MGX(ilevel)%parIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGX(ilevel)%ParIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

    END DO
    
    DO ilevel=1,mgLevels_Y

      SELECT CASE(ilevel)
      CASE(1)
        glN=Ngl
      CASE DEFAULT
        glN=1
      END SELECT
      
      CALL MPI_TYPE_VECTOR(nzc*(MGY(ilevel)%nyc+2*glN),glN,nxc+2*glN,mpi_double_precision,MGY(ilevel)%parVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGY(ilevel)%ParVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
    
      CALL MPI_TYPE_VECTOR(nzc,(nxc+2*glN)*glN,(MGY(ilevel)%nyc+2*glN)*(nxc+2*glN),mpi_double_precision,MGY(ilevel)%parVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGY(ilevel)%ParVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
      
      CALL MPI_TYPE_VECTOR(nzc*(MGY(ilevel)%nyc+2*glN),glN,nxc+2*glN,mpi_integer,MGY(ilevel)%parIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGY(ilevel)%ParIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
   
      CALL MPI_TYPE_VECTOR(nzc,(nxc+2*glN)*glN,(MGY(ilevel)%nyc+2*glN)*(nxc+2*glN),mpi_integer,MGY(ilevel)%parIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGY(ilevel)%ParIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
      
    END DO
    
    DO ilevel=1,mgLevels_Z

      SELECT CASE(ilevel)
      CASE(1)
        glN=Ngl
      CASE DEFAULT
        glN=1
      END SELECT
      
      CALL MPI_TYPE_VECTOR(MGZ(ilevel)%nzc*(nyc+2*glN),glN,nxc+2*glN,mpi_double_precision,MGZ(ilevel)%parVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGZ(ilevel)%ParVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
    
      CALL MPI_TYPE_VECTOR(MGZ(ilevel)%nzc,(nxc+2*glN)*glN,(nyc+2*glN)*(nxc+2*glN),mpi_double_precision,MGZ(ilevel)%parVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGZ(ilevel)%ParVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
      
      CALL MPI_TYPE_VECTOR(MGZ(ilevel)%nzc*(nyc+2*glN),glN,nxc+2*glN,mpi_integer,MGZ(ilevel)%parIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGZ(ilevel)%ParIVecWE, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
    
      CALL MPI_TYPE_VECTOR(MGZ(ilevel)%nzc,(nxc+2*glN)*glN,(nyc+2*glN)*(nxc+2*glN),mpi_integer,MGZ(ilevel)%parIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF

      CALL MPI_TYPE_COMMIT(MGZ(ilevel)%ParIVecSN, ierr)
      IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_init_MG_domain subroutine STOP'
       CALL flow_stop
       STOP
      ENDIF
      
    END DO
    
end subroutine par_init_MG_domain
!---------------------------------------------------------------------



subroutine par_comm_outermost_int(var,myTime)

! -------------------------------------------------------------------------
!  This subroutine communicates the outmost layer of subdomain for iblank.
! -------------------------------------------------------------------------

    USE flow_parameters
    USE boundary_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

    INTEGER, INTENT(INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,1:nzc)

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: VecSN, VecWE
    INTEGER :: ierr



    startTime = MPI_WTIME()

!   ---------------------------------------------------
!    Initialize the MPI derived communication vectors.
!   ---------------------------------------------------
    CALL MPI_TYPE_VECTOR(nzc,(nxc+2*Ngl),(nxc+2*Ngl)*(nyc+2*Ngl),mpi_integer,VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_comm_outermost_int STOP'
       CALL flow_stop
       STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF


    CALL MPI_TYPE_VECTOR(nzc*(nyc+2*Ngl),1,nxc+2*Ngl,mpi_integer,VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_comm_outermost_int STOP'
       CALL flow_stop
       STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(var(1-Ngl,nyc+1-Ngl,1),1,VecSN,myNorth,0,  &
                      var(1-Ngl,    1-Ngl,1),1,VecSN,mySouth,0,  &
                      commCart,status,ierr)

   IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for south to north in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

    
!   North to South
!   --------------
    CALL MPI_SENDRECV(var(1-Ngl,    Ngl,1),1,VecSN,mySouth,0,  &
                      var(1-Ngl,nyc+Ngl,1),1,VecSN,myNorth,0,  &
                      commCart,status,ierr)
   
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for north to south in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF
 
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(var(nxc+1-Ngl,1-Ngl,1),1,VecWE,myEast,0,  &
                      var(    1-Ngl,1-Ngl,1),1,VecWE,myWest,0,  &
                      commCart,status,ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for west to east in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

    
!   East to West
!   ------------
    CALL MPI_SENDRECV(var(    Ngl,1-Ngl,1),1,VecWE,myWest,0,  &
                      var(nxc+Ngl,1-Ngl,1),1,VecWE,myEast,0,  &
                      commCart,status,ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for east to west in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF


!   -----------------------------------------
!    Free MPI derived communication vectors.
!   -----------------------------------------

    CALL MPI_TYPE_FREE(VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_free in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

    CALL MPI_TYPE_FREE(VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_free in the par_comm_outermost_int subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

    endTime = MPI_WTIME()
    myTime=endTime-startTime

end subroutine par_comm_outermost_int
!---------------------------------------------------------------------



subroutine par_comm_outermost_real(var,myTime)

! -------------------------------------------------------------------------
!  This subroutine communicates the outmost layer of subdomain for iblank.
! -------------------------------------------------------------------------

    USE flow_parameters
    USE boundary_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

    REAL(KIND=CGREAL), INTENT(INOUT) :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,1:nzc)

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: VecSN, VecWE
    INTEGER :: ierr



    startTime = MPI_WTIME()

!   ---------------------------------------------------
!    Initialize the MPI derived communication vectors.
!   ---------------------------------------------------
    CALL MPI_TYPE_VECTOR(nzc,(nxc+2*Ngl),(nxc+2*Ngl)*(nyc+2*Ngl),mpi_double_precision,VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_comm_outermost_real STOP'
       CALL flow_stop
       STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

    CALL MPI_TYPE_VECTOR(nzc*(nyc+2*Ngl),1,nxc+2*Ngl,mpi_double_precision,VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN 
       WRITE(STDOUT, *)'Communication error for mpi_type_vector in the par_comm_outermost_real STOP'
       CALL flow_stop 
       STOP
     ENDIF

    CALL MPI_TYPE_COMMIT(VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_commit in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF

!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(var(1-Ngl,nyc+1-Ngl,1),1,VecSN,myNorth,0,  &
                      var(1-Ngl,    1-Ngl,1),1,VecSN,mySouth,0,  &
                      commCart,status,ierr)
   IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for south to north in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF
 
!   North to South
!   --------------
    CALL MPI_SENDRECV(var(1-Ngl,    Ngl,1),1,VecSN,mySouth,0,  &
                      var(1-Ngl,nyc+Ngl,1),1,VecSN,myNorth,0,  &
                      commCart,status,ierr)
   
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for north to south in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(var(nxc+1-Ngl,1-Ngl,1),1,VecWE,myEast,0,  &
                      var(    1-Ngl,1-Ngl,1),1,VecWE,myWest,0,  &
                      commCart,status,ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for west to east in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF   
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(var(    Ngl,1-Ngl,1),1,VecWE,myWest,0,  &
                      var(nxc+Ngl,1-Ngl,1),1,VecWE,myEast,0,  &
                      commCart,status,ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_sendrecv for east to west in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
    ENDIF  
!   -----------------------------------------
!    Free MPI derived communication vectors.
!   -----------------------------------------

    CALL MPI_TYPE_FREE(VecSN, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_free in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
     ENDIF

    CALL MPI_TYPE_FREE(VecWE, ierr)
    IF(ierr /= ERR_NONE) THEN
      WRITE(STDOUT, *)'Communication error for mpi_type_free in the par_comm_outermost_real subroutine STOP'
      CALL flow_stop
      STOP
     ENDIF

    endTime = MPI_WTIME()
    myTime=endTime-startTime

end subroutine par_comm_outermost_real
!---------------------------------------------------------------------



subroutine par_comm_POISSON(myTime)

! ---------------------------------------------------------------
!  This subroutine communicates pPrime after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE pressure_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: ierr



    startTime = MPI_WTIME()

!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(pPrime(1-Ngl,nyc+1-Ngl,1),1,ParVecSN,myNorth,0,  &
                      pPrime(1-Ngl,    1-Ngl,1),1,ParVecSN,mySouth,0,  &
                      commCart,status,ierr)
    
!   North to South
!   --------------
    CALL MPI_SENDRECV(pPrime(1-Ngl,    1,1),1,ParVecSN,mySouth,0,  &
                      pPrime(1-Ngl,nyc+1,1),1,ParVecSN,myNorth,0,  &
                      commCart,status,ierr)
    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(pPrime(nxc+1-Ngl,1-Ngl,1),1,ParVecWE,myEast,0,  &
                      pPrime(    1-Ngl,1-Ngl,1),1,ParVecWE,myWest,0,  &
                      commCart,status,ierr)
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(pPrime(    1,1-Ngl,1),1,ParVecWE,myWest,0,  &
                      pPrime(nxc+1,1-Ngl,1),1,ParVecWE,myEast,0,  &
                      commCart,status,ierr)

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_POISSON
!---------------------------------------------------------------------



subroutine par_comm_MGX(ilevel,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates pPrime after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER          , INTENT(IN ) :: ilevel
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr



    startTime = MPI_WTIME()

    SELECT CASE(ilevel)
    CASE(1)
      glN=Ngl
    CASE DEFAULT
      glN=1
    END SELECT
      
!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(MGX(ilevel)%phi(1-glN,nyc+1-glN,1),1,MGX(ilevel)%ParVecSN,myNorth,0,  &
                      MGX(ilevel)%phi(1-glN,    1-glN,1),1,MGX(ilevel)%ParVecSN,mySouth,0,  &
                      commCart,status,ierr)
    
!   North to South
!   --------------
    CALL MPI_SENDRECV(MGX(ilevel)%phi(1-glN,    1,1),1,MGX(ilevel)%ParVecSN,mySouth,0,  &
                      MGX(ilevel)%phi(1-glN,nyc+1,1),1,MGX(ilevel)%ParVecSN,myNorth,0,  &
                      commCart,status,ierr)
    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(MGX(ilevel)%phi(MGX(ilevel)%nxc+1-glN,1-glN,1),1,MGX(ilevel)%ParVecWE,myEast,0,  &
                      MGX(ilevel)%phi(                1-glN,1-glN,1),1,MGX(ilevel)%ParVecWE,myWest,0,  &
                      commCart,status,ierr)
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(MGX(ilevel)%phi(                1,1-glN,1),1,MGX(ilevel)%ParVecWE,myWest,0,  &
                      MGX(ilevel)%phi(MGX(ilevel)%nxc+1,1-glN,1),1,MGX(ilevel)%ParVecWE,myEast,0,  &
                      commCart,status,ierr)

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_MGX
!---------------------------------------------------------------------



subroutine par_comm_MGY(ilevel,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates pPrime after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER          , INTENT(IN ) :: ilevel
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr



    startTime = MPI_WTIME()

    SELECT CASE(ilevel)
    CASE(1)
      glN=Ngl
    CASE DEFAULT
      glN=1
    END SELECT
      
!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(MGY(ilevel)%phi(1-glN,MGY(ilevel)%nyc+1-glN,1),1,MGY(ilevel)%ParVecSN,myNorth,0,  &
                      MGY(ilevel)%phi(1-glN,                1-glN,1),1,MGY(ilevel)%ParVecSN,mySouth,0,  &
                      commCart,status,ierr)
    
!   North to South
!   --------------
    CALL MPI_SENDRECV(MGY(ilevel)%phi(1-glN,                1,1),1,MGY(ilevel)%ParVecSN,mySouth,0,  &
                      MGY(ilevel)%phi(1-glN,MGY(ilevel)%nyc+1,1),1,MGY(ilevel)%ParVecSN,myNorth,0,  &
                      commCart,status,ierr)
    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(MGY(ilevel)%phi(nxc+1-glN,1-glN,1),1,MGY(ilevel)%ParVecWE,myEast,0,  &
                      MGY(ilevel)%phi(    1-glN,1-glN,1),1,MGY(ilevel)%ParVecWE,myWest,0,  &
                      commCart,status,ierr)
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(MGY(ilevel)%phi(    1,1-glN,1),1,MGY(ilevel)%ParVecWE,myWest,0,  &
                      MGY(ilevel)%phi(nxc+1,1-glN,1),1,MGY(ilevel)%ParVecWE,myEast,0,  &
                      commCart,status,ierr)

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_MGY
!---------------------------------------------------------------------



subroutine par_comm_MGZ(ilevel,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates pPrime after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER          , INTENT(IN ) :: ilevel
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr



    startTime = MPI_WTIME()

    SELECT CASE(ilevel)
    CASE(1)
      glN=Ngl
    CASE DEFAULT
      glN=1
    END SELECT
      
!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------

!   South to North
!   --------------
    CALL MPI_SENDRECV(MGZ(ilevel)%phi(1-glN,nyc+1-glN,1),1,MGZ(ilevel)%ParVecSN,myNorth,0,  &
                      MGZ(ilevel)%phi(1-glN,    1-glN,1),1,MGZ(ilevel)%ParVecSN,mySouth,0,  &
                      commCart,status,ierr)
    
!   North to South
!   --------------
    CALL MPI_SENDRECV(MGZ(ilevel)%phi(1-glN,    1,1),1,MGZ(ilevel)%ParVecSN,mySouth,0,  &
                      MGZ(ilevel)%phi(1-glN,nyc+1,1),1,MGZ(ilevel)%ParVecSN,myNorth,0,  &
                      commCart,status,ierr)
    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(MGZ(ilevel)%phi(nxc+1-glN,1-glN,1),1,MGZ(ilevel)%ParVecWE,myEast,0,  &
                      MGZ(ilevel)%phi(    1-glN,1-glN,1),1,MGZ(ilevel)%ParVecWE,myWest,0,  &
                      commCart,status,ierr)
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(MGZ(ilevel)%phi(    1,1-glN,1),1,MGZ(ilevel)%ParVecWE,myWest,0,  &
                      MGZ(ilevel)%phi(nxc+1,1-glN,1),1,MGZ(ilevel)%ParVecWE,myEast,0,  &
                      commCart,status,ierr)

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_MGZ
!---------------------------------------------------------------------



subroutine par_comm_iblankX(mgLevels_X,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates iblank after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

    INTEGER, INTENT(IN) :: mgLevels_X

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr, ilevel



    startTime = MPI_WTIME()

    DO ilevel = 2,mgLevels_X
!     ----------------------------------
!      Communicating with the neighbors
!     ----------------------------------

!     South to North
!     --------------
!!      CALL MPI_SENDRECV(MGX(ilevel)%iblank(0,nyc,1),1,MGX(ilevel)%ParIVecSN,myNorth,0,  &
!!                        MGX(ilevel)%iblank(0,  0,1),1,MGX(ilevel)%ParIVecSN,mySouth,0,  &
!!                        commCart,status,ierr)
    
!     North to South
!     --------------
!!      CALL MPI_SENDRECV(MGX(ilevel)%iblank(0,    1,1),1,MGX(ilevel)%ParIVecSN,mySouth,0,  &
!!                        MGX(ilevel)%iblank(0,nyc+1,1),1,MGX(ilevel)%ParIVecSN,myNorth,0,  &
!!                        commCart,status,ierr)
    
!     --------------------------------------------------------------------------------------------------
!      NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                  the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!     --------------------------------------------------------------------------------------------------    

!     West to East
!     ------------
      CALL MPI_SENDRECV(MGX(ilevel)%iblank(MGX(ilevel)%nxc,0,1),1,MGX(ilevel)%ParIVecWE,myEast,0,  &
                        MGX(ilevel)%iblank(              0,0,1),1,MGX(ilevel)%ParIVecWE,myWest,0,  &
                        commCart,status,ierr)
    
!     East to West
!     ------------
      CALL MPI_SENDRECV(MGX(ilevel)%iblank(                1,0,1),1,MGX(ilevel)%ParIVecWE,myWest,0,  &
                        MGX(ilevel)%iblank(MGX(ilevel)%nxc+1,0,1),1,MGX(ilevel)%ParIVecWE,myEast,0,  &
                        commCart,status,ierr)
    ENDDO

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_iblankX
!---------------------------------------------------------------------



subroutine par_comm_iblankY(mgLevels_Y,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates iblank after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

    INTEGER, INTENT(IN) :: mgLevels_Y

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr, ilevel



    startTime = MPI_WTIME()

    DO ilevel = 2,mgLevels_Y
!     ----------------------------------
!      Communicating with the neighbors
!     ----------------------------------

!     South to North
!     --------------
      CALL MPI_SENDRECV(MGY(ilevel)%iblank(0,MGY(ilevel)%nyc,1),1,MGY(ilevel)%ParIVecSN,myNorth,0,  &
                        MGY(ilevel)%iblank(0,              0,1),1,MGY(ilevel)%ParIVecSN,mySouth,0,  &
                        commCart,status,ierr)
    
!     North to South
!     --------------
      CALL MPI_SENDRECV(MGY(ilevel)%iblank(0,                1,1),1,MGY(ilevel)%ParIVecSN,mySouth,0,  &
                        MGY(ilevel)%iblank(0,MGY(ilevel)%nyc+1,1),1,MGY(ilevel)%ParIVecSN,myNorth,0,  &
                        commCart,status,ierr)

!     --------------------------------------------------------------------------------------------------
!      NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                  the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!     --------------------------------------------------------------------------------------------------    

!     West to East
!     ------------
!!      CALL MPI_SENDRECV(MGY(ilevel)%iblank(nxc,0,1),1,MGY(ilevel)%ParIVecWE,myEast,0,  &
!!                        MGY(ilevel)%iblank(  0,0,1),1,MGY(ilevel)%ParIVecWE,myWest,0,  &
!!                        commCart,status,ierr)
    
!     East to West
!     ------------
!!      CALL MPI_SENDRECV(MGY(ilevel)%iblank(    1,0,1),1,MGY(ilevel)%ParIVecWE,myWest,0,  &
!!                        MGY(ilevel)%iblank(nxc+1,0,1),1,MGY(ilevel)%ParIVecWE,myEast,0,  &
!!                        commCart,status,ierr)
    ENDDO

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_iblankY
!---------------------------------------------------------------------



subroutine par_comm_iblankZ(mgLevels_Z,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates iblank after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE mg_arrays
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: myTime

    INTEGER, INTENT(IN) :: mgLevels_Z

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: glN, ierr, ilevel



    startTime = MPI_WTIME()

    DO ilevel = 2,mgLevels_Z
!     ----------------------------------
!      Communicating with the neighbors
!     ----------------------------------

!     South to North
!     --------------
      CALL MPI_SENDRECV(MGZ(ilevel)%iblank(0,nyc,1),1,MGZ(ilevel)%ParIVecSN,myNorth,0,  &
                        MGZ(ilevel)%iblank(0,  0,1),1,MGZ(ilevel)%ParIVecSN,mySouth,0,  &
                        commCart,status,ierr)
    
!     North to South
!     --------------
      CALL MPI_SENDRECV(MGZ(ilevel)%iblank(0,    1,1),1,MGZ(ilevel)%ParIVecSN,mySouth,0,  &
                        MGZ(ilevel)%iblank(0,nyc+1,1),1,MGZ(ilevel)%ParIVecSN,myNorth,0,  &
                        commCart,status,ierr)
    
!     --------------------------------------------------------------------------------------------------
!      NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                  the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!     --------------------------------------------------------------------------------------------------    

!     West to East
!     ------------
      CALL MPI_SENDRECV(MGZ(ilevel)%iblank(nxc,0,1),1,MGZ(ilevel)%ParIVecWE,myEast,0,  &
                        MGZ(ilevel)%iblank(  0,0,1),1,MGZ(ilevel)%ParIVecWE,myWest,0,  &
                        commCart,status,ierr)
    
!     East to West
!     ------------
      CALL MPI_SENDRECV(MGZ(ilevel)%iblank(    1,0,1),1,MGZ(ilevel)%ParIVecWE,myWest,0,  &
                        MGZ(ilevel)%iblank(nxc+1,0,1),1,MGZ(ilevel)%ParIVecWE,myEast,0,  &
                        commCart,status,ierr)
    ENDDO

    endTime = MPI_WTIME()
    myTime=endTime-startTime
    
end subroutine par_comm_iblankZ
!---------------------------------------------------------------------



subroutine par_comm_var(var,IL,JL,KL,glN,myTime)

! ---------------------------------------------------------------
!  This subroutine communicates pPrime after one LSOR iteration.
! ---------------------------------------------------------------

    USE flow_parameters
    USE par_parameters
    USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER          , INTENT(IN   ) :: IL, JL, KL, glN
    REAL(KIND=CGREAL), INTENT(INOUT) :: var(1-glN:IL+glN,1-glN:JL+glN,0:KL+1)
    REAL(KIND=CGREAL), INTENT(  OUT) :: myTime

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: startTime, endTime

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: VecSN, VecWE
    INTEGER :: ierr
!DEBUG
 DOUBLE PRECISION :: vart
 INTEGER :: ANY_TAG
!ENDDEBUG


    startTime = MPI_WTIME()

!   ---------------------------------------------------
!    Initialize the MPI derived communication vectors.
!   ---------------------------------------------------

    CALL MPI_TYPE_VECTOR(KL,(IL+2*glN)*glN,(IL+2*glN)*(JL+2*glN),mpi_double_precision,VecSN, ierr)
    CALL MPI_TYPE_COMMIT(VecSN, ierr)

    CALL MPI_TYPE_VECTOR(KL*(JL+2*glN),glN,IL+2*glN,mpi_double_precision,VecWE, ierr)
    CALL MPI_TYPE_COMMIT(VecWE, ierr)

!   ----------------------------------
!    Communicating with the neighbors
!   ----------------------------------
!!!!!!!!!!!!!!!!!DEBUG~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 vart = 1.0d0
! IF(myrank == 0) THEN
!  CALL MPI_SEND(vart,1,MPI_DOUBLE_PRECISION,1,0,MPI_COMM_WORLD,ierr)
! ELSE
!  CALL  MPI_RECV(vart,1,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
! ENDIF

! IF(myrank == 1) THEN
!  CALL MPI_SEND(vart, 1, MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, ierr)
! ELSE
!  CALL MPI_RECV(vart, 1, MPI_DOUBLE_PRECISION,1,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
! ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   South to North
!   --------------
    CALL MPI_SENDRECV(var(1-glN,JL+1-glN,1),1,VecSN,myNorth,0,  &
                      var(1-glN,   1-glN,1),1,VecSN,mySouth,0,  &
                      commCart,status,ierr)
    
!   North to South
!   --------------
    CALL MPI_SENDRECV(var(1-glN,   1,1),1,VecSN,mySouth,0,  &
                      var(1-glN,JL+1,1),1,VecSN,myNorth,0,  &
                      commCart,status,ierr)
    
!   --------------------------------------------------------------------------------------------------
!    NOTE(SAMK): At this point, the West and East neighbors have an up-to-date ghost values including
!                the corners. In this way I do not have to communicate with my NE,NW,SE,SW! COOL!!!
!   --------------------------------------------------------------------------------------------------    

!   West to East
!   ------------
    CALL MPI_SENDRECV(var(IL+1-glN,1-glN,1),1,VecWE,myEast,0,  &
                      var(   1-glN,1-glN,1),1,VecWE,myWest,0,  &
                      commCart,status,ierr)
    
!   East to West
!   ------------
    CALL MPI_SENDRECV(var(   1,1-glN,1),1,VecWE,myWest,0,  &
                      var(IL+1,1-glN,1),1,VecWE,myEast,0,  &
                      commCart,status,ierr)

!   -----------------------------------------
!    Free MPI derived communication vectors.
!   -----------------------------------------

    CALL MPI_TYPE_FREE(VecSN, ierr)
    CALL MPI_TYPE_FREE(VecWE, ierr)

    endTime = MPI_WTIME()
    myTime=endTime-startTime

end subroutine par_comm_var
!---------------------------------------------------------------------
subroutine par_comm_test
  USE flow_parameters
  USE par_parameters
  USE mpi

    IMPLICIT NONE

!   Parameters
!   ----------
!   Local variables
!   ---------------

    INTEGER :: status(MPI_STATUS_SIZE)
    INTEGER :: ierr
    DOUBLE PRECISION :: vart



  vart = 1.0d0
  IF(myrank == 0) THEN
    CALL MPI_SEND(vart,1,MPI_DOUBLE_PRECISION,1,0,MPI_COMM_WORLD,ierr)
  ELSE
    CALL  MPI_RECV(vart,1,MPI_DOUBLE_PRECISION,0,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
  ENDIF


  IF(myrank == 1) THEN
    CALL MPI_SEND(vart, 1, MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD, ierr)
  ELSE
    CALL MPI_RECV(vart, 1, MPI_DOUBLE_PRECISION,1,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
  ENDIF

endsubroutine par_comm_test
!---------------------------------------------------------------------

subroutine par_decomp(nicG,dir,nic,ni,Is,Ie)

! ------------------------------------------------------------------------
!  This subroutine determines the subdomain size for the current process.
! ------------------------------------------------------------------------

    USE flow_parameters
    use par_parameters
!!    use mg_parameters

    implicit none

!   Parameters
!   ----------
    integer, intent(in)  :: nicG, dir
    integer, intent(out) :: nic, ni
    integer, intent(out) :: Is, Ie

!   Local variables
!   ---------------
    integer :: res



    nic=nicG/Np(dir)

    res=mod((nicG),Np(dir))

    if (myCoords(dir)<res) then
      nic=nic+1
      Is=myCoords(dir)*nic+1
    else
      Is=myCoords(dir)*nic+res+1
    end if
    Ie=Is+nic-1

    ni=nic+1

end subroutine par_decomp
!---------------------------------------------------------------------



subroutine par_getMaxReal(myReal)

! -----------------------------------------------------------------------
!  This subroutine synchronizes a local real value with other processes.
! -----------------------------------------------------------------------

    USE par_parameters
    USE mpi

    implicit none

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(INOUT) :: myReal

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: maxReal
    INTEGER :: ierr



    CALL MPI_ALLREDUCE(myReal,maxReal,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD, ierr)
    
    myReal=maxReal

end subroutine par_getMaxReal



subroutine par_getSumReal(myReal)

! -----------------------------------------------------------------------
!  This subroutine synchronizes a local real value with other processes.
! -----------------------------------------------------------------------

    USE par_parameters
    USE mpi

    implicit none

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(INOUT) :: myReal

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: sumReal
    INTEGER :: ierr



    CALL MPI_ALLREDUCE(myReal,sumReal,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
    
    myReal=sumReal

end subroutine par_getSumReal



subroutine par_getSumRealArray(myReal,n)

! -----------------------------------------------------------------------
!  This subroutine synchronizes a local real value with other processes.
! -----------------------------------------------------------------------

    USE par_parameters
    USE mpi

    implicit none

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(INOUT) :: myReal(n)
    INTEGER :: n

!   Local variables
!   ---------------
    REAL(KIND=CGREAL) :: sumReal(n)
    INTEGER :: ierr



    CALL MPI_ALLREDUCE(myReal(1),sumReal(1),n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, ierr)
    
    myReal=sumReal

end subroutine par_getSumRealArray



subroutine par_getMaxInteger(myInteger)

! --------------------------------------------------------------------------
!  This subroutine synchronizes a local integer value with other processes.
! --------------------------------------------------------------------------

    USE par_parameters
    USE mpi

    implicit none

!   Parameters
!   ----------
    INTEGER, INTENT(INOUT) :: myInteger

!   Local variables
!   ---------------
    INTEGER :: maxInteger
    INTEGER :: ierr



    CALL MPI_ALLREDUCE(myInteger,maxInteger,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD, ierr)
    
    myInteger=maxInteger

end subroutine par_getMaxInteger



subroutine par_getSumInteger(myInteger)

! --------------------------------------------------------------------------
!  This subroutine synchronizes a local integer value with other processes.
! --------------------------------------------------------------------------

    USE par_parameters
    USE mpi

    implicit none

!   Parameters
!   ----------
    INTEGER, INTENT(INOUT) :: myInteger

!   Local variables
!   ---------------
    INTEGER :: sumInteger
    INTEGER :: ierr



    CALL MPI_ALLREDUCE(myInteger,sumInteger,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
    
    myInteger=sumInteger

end subroutine par_getSumInteger


subroutine par_bcast_marker(var, ibody, nb, isend, npt)
 
 USE global_parameters
 USE mpi
 USE flow_parameters
 IMPLICIT NONE

 INTEGER :: ibody, nb, isend, npt
 REAL(KIND=CGREAL) :: var(nB, nPtsMax)

 INTEGER :: VecM, ierr
 
 CALL MPI_TYPE_VECTOR(npt,1, nb,mpi_double_precision,VecM, ierr)
 CALL MPI_TYPE_COMMIT(VecM, ierr)

 CALL MPI_BCAST(var(ibody, 1), 1, VecM, isend, MPI_COMM_WORLD, ierr)
 
 CALL MPI_TYPE_FREE(VecM, ierr) 
end subroutine par_bcast_marker


! Added by Rajneesh for implicit coupling 

!-----------------------------------------
SUBROUTINE par_bcast_value(var)
 
 USE global_parameters
 USE mpi
 USE flow_parameters
 IMPLICIT NONE

  REAL(KIND=CGREAL) :: var

 INTEGER :: VecM, ierr
 
 CALL MPI_BCAST(var, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 
END SUBROUTINE par_bcast_value

!-----------------------------------------


