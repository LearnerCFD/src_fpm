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
!     Qian Xue
!     Rupesh Babu K. A.
!     Xudong Zheng 
!     Reza Ghias 
!     S. A. Mohsen Karimian
!
!  Filename: BOUNDARY_FRESH_CELL.PAR.F90
!  Latest Modification: September 10, 2008 (ver. P1.2.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     FreshCell_CalcExpWeight()
!     FreshCell_UpdateRhs()
! --------------------------------------------------------------------



SUBROUTINE FreshCell_CalcExpWeight()

! --------------------------------------------------------
!  Fresh Cells are created in Move_Boundary.
!  Adjust weight for NL term evaluation
!  AB2 for normal cell
!  FE  for cell that was fresh in the previous time step.
! --------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE
    
    INTEGER :: i, j, k



    DO k=1,nzc
    DO j=1,nyc
    DO i=1,nxc
      exp_weight(i,j,k) = 0.5_CGREAL * REAL(3-fresh_cell(i,j,k),KIND=CGREAL)
    END DO
    END DO
    END DO
       
END SUBROUTINE FreshCell_CalcExpWeight
!---------------------------------------------------------------------



SUBROUTINE FreshCell_UpdateRhs() 

! -----------------------------------------------
!  Update Advection-Diffusion RHS for fresh cell
! -----------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      nlu(i,j,k) = nlu(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)
      nlv(i,j,k) = nlv(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)
      nlw(i,j,k) = nlw(i,j,k)* REAL(1-fresh_cell(i,j,k),KIND=CGREAL)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

   END SUBROUTINE  FreshCell_UpdateRhs
!-------------------------------------------------------------------------------
