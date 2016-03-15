! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  Johns Hopkins University
!
!  VICAR3D, a viscous, Cartesian, 3D flow solver.
!
!  This is a contineously developing project.
!
!  Starting Developers:
!  Fady Najjar
!  Rajat Mittal
!
!  Other contributing programmers:
!     S. A. Mohsen Karimian
!
!  Filename: CORE_BC_PERIODIC.F90
!  Created: July 30, 2007 (ver. 12.2.4)
!  Last Modification: September 05, 2007 (ver. 12.2.4)
!  By F. Najjar
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     apply_periodic_bcvel()
!     apply_periodic_pres()
! --------------------------------------------------------------------

!****************************************************************************************
!***************Ehsan added this subrountine to Parallel code for Periodic BC************
!****************************************************************************************


SUBROUTINE apply_periodic_bcvel

! --------------------------------------------------
!  Purpose: Apply periodic conditions on bcVel.
! --------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,n

!   z-direction 

    IF ( ndim == DIM_3D           .AND. &
         bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN

      DO j=1,nyc
      DO i=1,nxc
        bczu(i,j,1)    = 0.5_CGREAL*( u(i,j,1) + u(i,j,nzc) )
        bczv(i,j,1)    = 0.5_CGREAL*( v(i,j,1) + v(i,j,nzc) )
        bczw(i,j,1)    = 0.5_CGREAL*( w(i,j,1) + w(i,j,nzc) )

        bczu(i,j,nzc) = bczu(i,j,1)
        bczv(i,j,nzc) = bczv(i,j,1)
        bczw(i,j,nzc) = bczw(i,j,1)
      ENDDO ! i
      ENDDO ! j
    
    END IF ! bcz1
                
END SUBROUTINE apply_periodic_bcvel
! -------------------------------------------------------------------------




SUBROUTINE apply_periodic_pres(pres,mx,my,mz)
          
! -------------------------------------------------------------------------
!  Purpose: Apply periodic pressure conditions for LSOR and MG methods.
!    
!  Input: pressure values
!                          
!  Output: enforced pressure values.     
! -------------------------------------------------------------------------
          
    USE global_parameters
    USE flow_parameters   
    USE pressure_arrays  

    IMPLICIT NONE

     
    INTEGER,INTENT(IN):: mx,my,mz
    REAL(KIND=CGREAL), DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT(INOUT) :: pres

!   Local variables     

    INTEGER :: i,j,k,n

!   z-direction  

    IF ( ndim == DIM_3D           .AND. &
         bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN

      DO j=0, my 
      DO i=0, mx 
        pres(i,j,0) = pres(i,j,mz-1)
        pres(i,j,mz) = pres(i,j,1) 
      END DO !i
      END DO !j 

    END IF ! bcz1
      
END SUBROUTINE apply_periodic_pres
!-------------------------------------------------------------------------



