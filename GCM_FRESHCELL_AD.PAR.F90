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
!  Filename: GCM_FRESHCELL_AD.PAR.F90
!  Latest Modification: October 21, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_correct_res_ad()
!     GCM_correct_rhs_ad(mDirection,iFr,jFr,kFr,nman,rhs,var)
!------------------------------------------------------------------------------



SUBROUTINE GCM_correct_res_ad(iFr,jFr,kFr,var,res)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER,           INTENT(IN)    :: iFr, jFr, kFr
    REAL(KIND=CGREAL), INTENT(INOUT) :: res
    REAL(KIND=CGREAL), INTENT(IN)    :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

!   Loop Variables
!   --------------
    INTEGER :: n,iRow

!   Local Variables
!   ---------------
    INTEGER :: i,j,k


!   Off diagonal elements contribution
!   First select appropriate fresh cell from list
!   ---------------------------------------------
    IF (nDim == DIM_2D) THEN
      k = kFr
    ENDIF ! nDim

    DO n = 1, nFresh
      IF ( iFr == iFresh(n) .AND. jFr == jFresh(n) .AND. kFr == kFresh(n) ) THEN
        DO iRow = 1,iRowMax
          i  = iFreshCellIndex(n) + incI(iRow)
          j  = jFreshCellIndex(n) + incJ(iRow)
          k  = kFreshCellIndex(n) + incK(iRow)

          IF ( (i==iFr-1 .AND. j==jFr-1) .OR. &     ! these are the stencil elements
               (i==iFr-1 .AND. j==jFr+1) .OR. &     !
               (i==iFr+1 .AND. j==jFr-1) .OR. &     !
               (i==iFr+1 .AND. j==jFr+1) .OR. &     ! that do not fit
               (i==iFr-1 .AND. k==kFr-1) .OR. &     !
               (i==iFr-1 .AND. k==kFr+1) .OR. &     !
               (i==iFr+1 .AND. k==kFr-1) .OR. &     ! within the normal 7-point stencil
               (i==iFr+1 .AND. k==kFr+1) .OR. &     !
               (j==jFr-1 .AND. k==kFr-1) .OR. &     !
               (j==jFr-1 .AND. k==kFr+1) .OR. &     !
               (j==jFr+1 .AND. k==kFr-1) .OR. &     !
               (j==jFr+1 .AND. k==kFr+1)      ) THEN

            res = res + coeffGCMFreshD(iRow,n)*var(i,j,k)

          ENDIF ! i
        ENDDO ! iRow
      ENDIF ! iFr
    ENDDO ! n

   END SUBROUTINE GCM_correct_res_ad
!---------------------------------------------------------------------



SUBROUTINE GCM_correct_rhs_ad(mDirection,iFr,jFr,kFr,nmax,rhs,var)

    USE global_parameters
    USE   flow_parameters
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: mDirection, iFr, jFr, kFr, nmax
    REAL(KIND=CGREAL), INTENT(INOUT) :: rhs(1:nmax)
    REAL(KIND=CGREAL), INTENT(IN)    :: var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

!   Loop Variables
!   --------------
    INTEGER :: n,iRow

!   Local Variables
!   ---------------
    INTEGER :: i,j,k



!   Off diagonal elements contribution: First select appropriate fresh cell from list.
!   ----------------------------------------------------------------------------------
    DO n = 1, nFresh
      IF ( iFr == iFresh(n) .AND. jFr == jFresh(n) .AND. kFr == kFresh(n) ) THEN
        DO iRow = 1,iRowMax
          i  = iFreshCellIndex(n) + incI(iRow)
          j  = jFreshCellIndex(n) + incJ(iRow)
          k  = kFreshCellIndex(n) + incK(iRow)

          IF ( (i==iFr-1 .AND. j==jFr-1) .OR. &     ! these are the stencil elements
               (i==iFr-1 .AND. j==jFr+1) .OR. &     !
               (i==iFr+1 .AND. j==jFr-1) .OR. &     !
               (i==iFr+1 .AND. j==jFr+1) .OR. &     ! that do not fit
               (i==iFr-1 .AND. k==kFr-1) .OR. &     !
               (i==iFr-1 .AND. k==kFr+1) .OR. &     !
               (i==iFr+1 .AND. k==kFr-1) .OR. &     ! within the normal 7-point stencil
               (i==iFr+1 .AND. k==kFr+1) .OR. &     !
               (j==jFr-1 .AND. k==kFr-1) .OR. &     !
               (j==jFr-1 .AND. k==kFr+1) .OR. &     !
               (j==jFr+1 .AND. k==kFr-1) .OR. &     !
               (j==jFr+1 .AND. k==kFr+1)      ) THEN

            SELECT CASE(mDirection)

            CASE(ICOORD)
              rhs(iFr) = rhs(iFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

            CASE(JCOORD)
              rhs(jFr) = rhs(jFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

            CASE(KCOORD)
              rhs(kFr) = rhs(kFr)+coeffGCMFreshD(iRow,n)*var(i,j,k)

            END SELECT ! mDirection

          ENDIF ! i
        ENDDO ! iRow
      ENDIF ! iFr
    ENDDO ! n


END SUBROUTINE GCM_correct_rhs_ad
!---------------------------------------------------------------------
