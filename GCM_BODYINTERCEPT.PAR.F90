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
!  Filename: GCM_BODYINTERCEPT.PAR.F90
!  Latest Modification: February 13, 2008 (ver. 0.3.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_SetBodyInterceptValues()
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  Purpose: Compute velocity components and normal pressure gradient 
!           at internal body intercept points. These are needed as BC
!           in solvers.
!
!  Input: [u,v,w]BodyMarker  = velocity components of BM points
!         closestMarker      = closest Marker for Ghost point,
!         closestMarkerRatio = closest Marker ratio for Ghost point
!
!  Output: [u,v,w]BodyIntercept = velocity components of BI point,
!          dpdnBodyIntercept    = normal pressure gradient of BI point,
!          dpdtBodyIntercept    = tangential pressure gradient of BI point.
!
!  Notes: Currently dpdnBodyIntercept and dpdtBodyIntercept are set to zero.
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



  SUBROUTINE GCM_SetBodyInterceptValues()

! ---------------------------------------------------------------------------
!  This subroutine computes velocity components and normal pressure gradient
!  at internal body intercept points. These are needed as BC in solvers.
! ---------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE pressure_arrays
    USE grid_arrays
    USE GCM_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

!   Loop variables
!   --------------
    INTEGER :: n, k

!   Local variables
!   ---------------
    INTEGER :: iG, jG, kG

    REAL(KIND=CGREAL) :: xGC, yGC, zGC



!   Loop over all ghost points, compute velocity field
!   --------------------------------------------------
    DO n=1,nGhost
      iG      = iGhost(n)
      jG      = jGhost(n)
      kG      = kGhost(n)
      xGC     = xc(L2GI(iG))
      yGC     = yc(L2GJ(jG))
      zGC     = zc(kG)
      
      CALL GCM_calc_BIVelocity_Unstruc( iG, jG,  kG,        &
                                       xGC, yGC, zGC,       &
                                       closestElementGC(n), &
                                       uBodyIntercept(n),   &
                                       vBodyIntercept(n),   &
                                       wBodyIntercept(n)    )
    ENDDO ! n

!    IF (nGhost/=0) THEN
    IF (MOD(ntime,nmonitor) == 0 ) THEN
      WRITE(ifuParLog,'(A,2D14.5)') 'Min-Max uBI-Actual=     ', MINVAL(uBodyIntercept(1:nGhost)), MAXVAL(uBodyIntercept(1:nGhost))
      WRITE(ifuParLog,'(A,2D14.5)') 'Min-Max vBI-Actual=     ', MINVAL(vBodyIntercept(1:nGhost)), MAXVAL(vBodyIntercept(1:nGhost))
      WRITE(ifuParLog,'(A,2D14.5)') 'Min-Max wBI-Actual=     ', MINVAL(wBodyIntercept(1:nGhost)), MAXVAL(wBodyIntercept(1:nGhost))

      WRITE(ifuParLog,'(A,2D14.5)') 'Min-Max xBINorm-Actual= ', MINVAL(xBodyInterceptNorm(1:nGhost)), MAXVAL(xBodyInterceptNorm(1:nGhost))
      WRITE(ifuParLog,'(A,2D14.5)') 'Min-Max yBINorm-Actual= ', MINVAL(yBodyInterceptNorm(1:nGhost)), MAXVAL(yBodyInterceptNorm(1:nGhost))
    END IF
 
!   ---------------------------------------------------
!    compute body intercept values for pressure field
!    
!       p   - p   = (dP/dn)  *  (total probe length)/2
!        gp    bi          bi           
!   ---------------------------------------------------
    DO n = 1,nGhost
      iG = iGhost(n)
      jG = jGhost(n)

      DO k = 1,nzc

!       IF (pressure_bc == INVISCID) THEN
!         dpdnBodyIntercept(n) = probeLength(n) *    
!                            ( -axBodyIntercept(n)*xBodyInterceptNorm(n) &  
!                              -ayBodyIntercept(n)*yBodyInterceptNorm(n) &
!                              -azBodyIntercept(n)*zBodyInterceptNorm(n) )
!       ELSE

        dpdnBodyIntercept(n) = 0.0_CGREAL

!       ENDIF ! pressure_bc 

! compute dp/dtau (tangential component)
!           dpdtBodyIntercept(n) = 0.0_CGREAL

      ENDDO ! k 
    ENDDO ! n
    
END SUBROUTINE GCM_SetBodyInterceptValues      
!---------------------------------------------------------------------
