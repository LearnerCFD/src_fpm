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
!  Filename: CORE_SET_OUTER_GHOST.PAR.F90
!  Latest Modification: October 21, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     set_outer_ghost_vel()
!     set_outer_ghost_pres(pres)
! --------------------------------------------------------------------

# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1


SUBROUTINE set_outer_ghost_vel()

! --------------------------------------------------------------------------
!  This subroutine sets outer ghost boundary conditions for velocity field.
! --------------------------------------------------------------------------
    USE global_parameters
    USE flow_parameters
    USE flow_arrays

    IMPLICIT NONE

    INTEGER  :: i,j,k



!   left boundary
!   -------------
    IF (myCoords(1)==0) THEN
      i = 0

      u(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxu(i+1,1:nyc,1:nzc) - u(i+1,1:nyc,1:nzc)
      uGhost(i,1:nyc,1:nzc) = u(i,1:nyc,1:nzc)

      v(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxv(i+1,1:nyc,1:nzc) - v(i+1,1:nyc,1:nzc)
      vGhost(i,1:nyc,1:nzc) = v(i,1:nyc,1:nzc)

      w(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxw(i+1,1:nyc,1:nzc) - w(i+1,1:nyc,1:nzc)
      wGhost(i,1:nyc,1:nzc) = w(i,1:nyc,1:nzc)

    END IF

!   right boundary
!   --------------
    IF (myCoords(1)==Np(1)-1) THEN
      i = nxc+1

      u(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxu(i-1,1:nyc,1:nzc) - u(i-1,1:nyc,1:nzc)
      uGhost(i,1:nyc,1:nzc) = u(i,1:nyc,1:nzc)

      v(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxv(i-1,1:nyc,1:nzc) - v(i-1,1:nyc,1:nzc)
      vGhost(i,1:nyc,1:nzc) = v(i,1:nyc,1:nzc)

      wGhost(i,1:nyc,1:nzc) = w(i,1:nyc,1:nzc)
      w(i,1:nyc,1:nzc) = 2.0_CGREAL*bcxw(i-1,1:nyc,1:nzc) - w(i-1,1:nyc,1:nzc)

    END IF

!   bottom boundary
!   ---------------
    IF (myCoords(2)==0) THEN
      j = 0

      u(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyu(1:nxc,j+1,1:nzc) - u(1:nxc,j+1,1:nzc)
      uGhost(1:nxc,j,1:nzc) = u(1:nxc,j,1:nzc)

      v(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyv(1:nxc,j+1,1:nzc) - v(1:nxc,j+1,1:nzc)
      vGhost(1:nxc,j,1:nzc) = v(1:nxc,j,1:nzc)

      wGhost(1:nxc,j,1:nzc) = w(1:nxc,j,1:nzc)
      w(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyw(1:nxc,j+1,1:nzc) - w(1:nxc,j+1,1:nzc)

    END IF

!   top boundary
!   ------------
    IF (myCoords(2)==Np(2)-1) THEN
      j = nyc+1

      u(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyu(1:nxc,j-1,1:nzc) - u(1:nxc,j-1,1:nzc)
      uGhost(1:nxc,j,1:nzc) = u(1:nxc,j,1:nzc)

      v(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyv(1:nxc,j-1,1:nzc) - v(1:nxc,j-1,1:nzc)
      vGhost(1:nxc,j,1:nzc) = v(1:nxc,j,1:nzc)

      w(1:nxc,j,1:nzc) = 2.0_CGREAL*bcyw(1:nxc,j-1,1:nzc) - w(1:nxc,j-1,1:nzc)
      wGhost(1:nxc,j,1:nzc) = w(1:nxc,j,1:nzc)

    END IF

!   back boundary
!   -------------
    k = 0

    u(1:nxc,1:nyc,k) = 2.0_CGREAL*bczu(1:nxc,1:nyc,k+1) - u(1:nxc,1:nyc,k+1)
    uGhost(1:nxc,1:nyc,k) = u(1:nxc,1:nyc,k)

    v(1:nxc,1:nyc,k) = 2.0_CGREAL*bczv(1:nxc,1:nyc,k+1) - v(1:nxc,1:nyc,k+1)
    vGhost(1:nxc,1:nyc,k) = v(1:nxc,1:nyc,k)

    w(1:nxc,1:nyc,k) = 2.0_CGREAL*bczw(1:nxc,1:nyc,k+1) - w(1:nxc,1:nyc,k+1)
    wGhost(1:nxc,1:nyc,k) = w(1:nxc,1:nyc,k)

!   front boundary
!   --------------
    k = nzc+1

    u(1:nxc,1:nyc,k) = 2.0_CGREAL*bczu(1:nxc,1:nyc,k-1) - u(1:nxc,1:nyc,k-1)
    uGhost(1:nxc,1:nyc,k) = u(1:nxc,1:nyc,k)

    v(1:nxc,1:nyc,k) = 2.0_CGREAL*bczv(1:nxc,1:nyc,k-1) - v(1:nxc,1:nyc,k-1)
    vGhost(1:nxc,1:nyc,k) = v(1:nxc,1:nyc,k)

    w(1:nxc,1:nyc,k) = 2.0_CGREAL*bczw(1:nxc,1:nyc,k-1) - w(1:nxc,1:nyc,k-1)
    wGhost(1:nxc,1:nyc,k) = w(1:nxc,1:nyc,k)

END SUBROUTINE set_outer_ghost_vel
!---------------------------------------------------------------------



SUBROUTINE set_outer_ghost_pres(pres)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), INTENT(INOUT) :: pres(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    INTEGER             :: i,j,k
    INTEGER             :: iG, jG


!   outer ghost boundary conditions
!   -------------------------------
!!   IF (bcx1 == BC_TYPE_PERIODIC .AND. &
!!       bcx2 == BC_TYPE_PERIODIC) THEN
!!     DO k=0,nz
!!     DO j=0,ny
!!        pres(0,j,k)  = pres(nx-1,j,k)
!!        pres(nx,j,k) = pres(1,j,k)
!!     ENDDO ! j
!!     ENDDO ! k
!!   END IF

    IF (pbcx1 == PBC_NEUMANN .AND. myCoords(1)==0) THEN
      iG = L2GI(1)
      DO k=0,nzc+1
      DO j=myJLL,myJUL
        pres(0,j,k)   = pres(1,j,k) - pgradx1(j,k)*dxc(iG)
        pGhost(0,j,k) = pres(0,j,k) 
      ENDDO ! j
      ENDDO ! k
    END IF
 
    IF (pbcx2 == PBC_NEUMANN .AND. myCoords(1)==Np(1)-1) THEN
      iG = L2GI(nxc+1)
      DO k=0    ,nzc+1
      DO j=myJLL,myJUL
        pres(nxc+1,j,k)   = pres(nxc  ,j,k) + pgradx2(j,k)*dxc(iG)
        pGhost(nxc+1,j,k) = pres(nxc+1,j,k)
      ENDDO ! j
      ENDDO ! k
    END IF

!!   IF (bcy1 .EQ. BC_TYPE_PERIODIC .AND. &
!!       bcy2 .EQ. BC_TYPE_PERIODIC) THEN
!!     DO k=0,nz
!!     DO i=0,nx
!!        pres(i,0,k)  = pres(i,ny-1,k)
!!        pres(i,ny,k) = pres(i,1,k)
!!     ENDDO ! i
!!     ENDDO ! k
!!   END IF

    IF (pbcy1 == PBC_NEUMANN .AND. myCoords(2)==0) THEN
      jG = L2GJ(1)
      DO k=0    ,nzc+1
      DO i=myILL,myIUL
        pres(i,0,k)    = pres(i,1,k) - pgrady1(i,k)*dyc(jG)
        pGhost(i,0,k)  = pres(i,0,k)
      ENDDO ! i
      ENDDO ! k
    END IF
 
    IF (pbcy2 == PBC_NEUMANN .AND. myCoords(2)==Np(2)-1) THEN
      jG = L2GJ(nyc+1)
      DO k=0,nz
      DO i=myILL,myIUL
        pres(i,nyc+1,k)   = pres(i,nyc  ,k) + pgrady2(i,k)*dyc(jG)
        pGhost(i,nyc+1,k) = pres(i,nyc+1,k)
      ENDDO ! i
      ENDDO ! k
    END IF

    IF (ndim == DIM_3D) THEN                            !Ehsan added for periodic BC
      IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &            
          bcz2 .EQ. BC_TYPE_PERIODIC) THEN            
       DO j=myJLL,myJUL                            
       DO i=myILL,myIUL                               
          pres(i,j,0)  = pres(i,j,nzc)          
          pres(i,j,nzc+1) = pres(i,j,1)        
       ENDDO ! i                              
       ENDDO ! j                
     END IF
   ENDIF


!    IF (ndim == DIM_3D) THEN
!!      IF (bcz1 .EQ. BC_TYPE_PERIODIC .AND. &
!!          bcz2 .EQ. BC_TYPE_PERIODIC) THEN
!!       DO j=0,ny
!!       DO i=0,nx
!!          pres(i,j,0)  = pres(i,j,nz-1)
!!          pres(i,j,nz) = pres(i,j,1)
!!       ENDDO ! i
!!       ENDDO ! j
!!     ELSE
!!       DO j=0,ny
!!       DO i=0,nx
!!          pres(i,j,0)  = pres(i,j,1)
!!          pres(i,j,nz) = pres(i,j,nz-1)
!!       ENDDO ! i
!!       ENDDO ! j
!!     END IF
!   ENDIF

      IF (pbcz1 == PBC_NEUMANN ) THEN
        DO j=myJLL,myJUL
        DO i=myILL,myIUL
          pres(i,j,0)   = pres(i,j,1) - pgradz1(i,j)*dzc(1)
          pGhost(i,j,0) = pres(i,j,0) 
        ENDDO ! i
        ENDDO ! j
      END IF
 
      IF (pbcz2 == PBC_NEUMANN ) THEN
        DO j=myJLL,myJUL
        DO i=myILL,myIUL
          pres(i,j,nzc+1)   = pres(i,j,nzc) + pgradz2(i,j)*dzc(nz)
          pGhost(i,j,nzc+1) = pres(i,j,nz)
        ENDDO ! i
        ENDDO ! j
      END IF

END SUBROUTINE set_outer_ghost_pres
!-------------------------------------------------------------------------
