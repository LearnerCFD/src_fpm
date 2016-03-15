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
!  Filename: CORE_VEL_ADJUST2D.PAR.F90
!  Latest Modification: August 14, 2008 (ver. P1.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     vel_adjust2D()
! --------------------------------------------------------------------



SUBROUTINE vel_adjust2D() 

! ---------------------------------------------------------
!  This subroutine sets velocity field for 2D calculations
! ---------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays

    IMPLICIT NONE

    INTEGER:: i,j,k



!   copy k=1 plane to other planes
!   ------------------------------
    DO k = 2,nzc
    DO j = myJLL,myJUL
    DO i = myILL,myIUL
      u(i,j,k) = u(i,j,1)
      v(i,j,k) = v(i,j,1)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

    DO k = 2,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      face_uw(i,j,k) = face_uw(i,j,1)
      face_vs(i,j,k) = face_vs(i,j,1)
      face_ue(i,j,k) = face_ue(i,j,1)
      face_vn(i,j,k) = face_vn(i,j,1)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    
!   zero w-component
!   ----------------
    w(:,:,:)       = 0.0_CGREAL
    face_wb(:,:,:) = 0.0_CGREAL
    face_wf(:,:,:) = 0.0_CGREAL

END SUBROUTINE  vel_adjust2D
!---------------------------------------------------------------------
