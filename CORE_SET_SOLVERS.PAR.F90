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
!  Filename: CORE_SET_SOLVERS.PAR.F90
!  Latest Modification: September 02, 2008 (ver. P1.2.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     set_solve_ad()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE set_solve_ad()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays

    IMPLICIT NONE

!   Loop variables
!   --------------    
    INTEGER :: i, j, k, n, iRow
    
!   Local variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: iFr,jFr,kFr

    REAL(KIND=CGREAL) :: rFreshCell
    REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuF,nuB



!   Initialize coefficients
!   -----------------------
    amx_ad = 0.0_CGREAL
    apx_ad = 0.0_CGREAL
    
    amy_ad = 0.0_CGREAL
    apy_ad = 0.0_CGREAL
    
    amz_ad = 0.0_CGREAL
    apz_ad = 0.0_CGREAL
    
    DO k=1,nzc
    DO j=1,nyc
      jG=L2GJ(j)
      DO i=1,nxc
        iG=L2GI(i)

        nuE = (             fx(iG+1)   *viscTot(i+1,j,k)                            &
            +  ( 1.0_CGREAL-fx(iG+1) ) *viscTot(i,j,k)   )*(1.0_CGREAL-iup(i,j,k))  &
            + bcxvisc(i,j,k)*iup(i,j,k)

        nuW = (             fx(iG)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fx(iG)   ) *viscTot(i-1,j,k) )*(1.0_CGREAL-ium(i,j,k))  &
            + bcxvisc(i,j,k)*ium(i,j,k)

        nuN = (             fy(jG+1)   *viscTot(i,j+1,k)                            &
            +  ( 1.0_CGREAL-fy(jG+1) ) *viscTot(i,j,k)   )*(1.0_CGREAL-jup(i,j,k))  &
            + bcyvisc(i,j,k)*jup(i,j,k)

        nuS = (             fy(jG)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fy(jG)   ) *viscTot(i,j-1,k) )*(1.0_CGREAL-jum(i,j,k))  &
            + bcyvisc(i,j,k)*jum(i,j,k)
        
        nuF = (             fz(k+1)   *viscTot(i,j,k+1)                            &
            +  ( 1.0_CGREAL-fz(k+1) ) *viscTot(i,j,k)   )*(1.0_CGREAL-kup(i,j,k))  &
            + bczvisc(i,j,k)*kup(i,j,k)
 
        nuB = (             fz(k)     *viscTot(i,j,k)                              &
            +  ( 1.0_CGREAL-fz(k)   ) *viscTot(i,j,k-1) )*(1.0_CGREAL-kum(i,j,k))  &
            + bczvisc(i,j,k)*kum(i,j,k)

        amx_ad(i,j,k) =  dxcinv(iG)  * dxinv(iG)
        apx_ad(i,j,k) =  dxcinv(iG+1)* dxinv(iG)

        amx_ad(i,j,k) =- (0.50_CGREAL*dt*nuW)*amx_ad(i,j,k)
        apx_ad(i,j,k) =- (0.50_CGREAL*dt*nuE)*apx_ad(i,j,k)
       
        amy_ad(i,j,k) =  dycinv(jG)  * dyinv(jG)
        apy_ad(i,j,k) =  dycinv(jG+1)* dyinv(jG)

        amy_ad(i,j,k) =- (0.50_CGREAL*dt*nuS)*amy_ad(i,j,k)
        apy_ad(i,j,k) =- (0.50_CGREAL*dt*nuN)*apy_ad(i,j,k)
      
        amz_ad(i,j,k) =  dzcinv(k)  * dzinv(k)
        apz_ad(i,j,k) =  dzcinv(k+1)* dzinv(k)

        amz_ad(i,j,k) =- (0.50_CGREAL*dt*nuB)*amz_ad(i,j,k)*KillFor2D
        apz_ad(i,j,k) =- (0.50_CGREAL*dt*nuF)*apz_ad(i,j,k)*KillFor2D 
      END DO ! i
    END DO ! j
    END DO !k

!   TAKE CARE OF FRESH CELLS
!   ------------------------
    IF ( boundary_motion == MOVING_BOUNDARY ) THEN

      SELECT CASE(boundary_formulation)

      CASE(SSM_METHOD)

!       --------------------------------------------------------------------------    
!        Value of fresh cell is computed through interpolation from six neighbors
!        use Inversed Distance Weighted Intepolation (Shepards Method)
!        we grab bcxu etc values from nearest SSM walls in this procedure
!
!                          1   6  [  -p    ]
!         u              = -  SUM [ h   u  ]   
!          j               H  i=1 [  ij  i ]
!                           j 
!        where 
!                              6  [  -p  ]
!         H              =    SUM [ h    ]   
!          j                  i=1 [  ij  ]
!
!         h  : distance between location i and j
!          ij
!         p  : free parameter (usually taken as 2)
!
!        equation is coded up in the following form: 
!                     6  [  -p    ]
!         H  u    -  SUM [ h   u  ]  = 0
!          j  j      i=1 [  ij  i ]
!       --------------------------------------------------------------------------    

        DO k=1,nzc    
        DO j=1,nyc   
          jG=L2GJ(j)

          DO i=1,nxc
            iG=L2GI(i)

            rFreshCell = REAL(fresh_cell(i,j,k),KIND=CGREAL)

            amx_ad(i,j,k) = amx_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - ( dxcinv(iG) *(1.0_CGREAL + ium(i,j,k) ) )**sidw *rFreshCell

            apx_ad(i,j,k) = apx_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - (dxcinv(iG+1)*(1.0_CGREAL + iup(i,j,k) ) )**sidw *rFreshCell

            amy_ad(i,j,k) = amy_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - ( dycinv(jG) *(1.0_CGREAL + jum(i,j,k) ) )**sidw *rFreshCell

            apy_ad(i,j,k) = apy_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - (dycinv(jG+1)*(1.0_CGREAL + jup(i,j,k) ) )**sidw *rFreshCell
      
            amz_ad(i,j,k) = amz_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - ( dzcinv(k)  *(1.0_CGREAL + kum(i,j,k) ) )**sidw *rFreshCell
      
            apz_ad(i,j,k) = apz_ad(i,j,k)*(1.0_CGREAL - rFreshCell)     &
                          - ( dzcinv(k+1)*(1.0_CGREAL + kup(i,j,k) ) )**sidw *rFreshCell

            amz_ad(i,j,k) = amz_ad(i,j,k)*KillFor2D
            apz_ad(i,j,k) = apz_ad(i,j,k)*KillFor2D

          END DO ! i
        END DO ! j
        END DO ! k

      CASE(GCM_METHOD)

        CALL GCM_set_freshcell()

        DO n=1,nFresh

          iFr = iFresh(n)
          jFr = jFresh(n)
          kFr = kFresh(n)

!         Initialize coefficient to zero for fresh cells
!         ----------------------------------------------
          amx_ad (iFr,jFr,kFr) = 0.0_CGREAL
          apx_ad (iFr,jFr,kFr) = 0.0_CGREAL
          amy_ad (iFr,jFr,kFr) = 0.0_CGREAL
          apy_ad (iFr,jFr,kFr) = 0.0_CGREAL
          amz_ad (iFr,jFr,kFr) = 0.0_CGREAL
          apz_ad (iFr,jFr,kFr) = 0.0_CGREAL

          DO iRow = 1,iRowMax
            i  = iFreshCellIndex(n) + incI(iRow)
            j  = jFreshCellIndex(n) + incJ(iRow)
            k  = kFreshCellIndex(n) + incK(iRow)

            IF ( i==iFr-1 .AND. j==jFr   .AND. k==kFr   ) amx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            IF ( i==iFr+1 .AND. j==jFr   .AND. k==kFr   ) apx_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            IF ( i==iFr   .AND. j==jFr-1 .AND. k==kFr   ) amy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            IF ( i==iFr   .AND. j==jFr+1 .AND. k==kFr   ) apy_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            IF ( i==iFr   .AND. j==jFr   .AND. k==kFr-1 ) amz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
            IF ( i==iFr   .AND. j==jFr   .AND. k==kFr+1 ) apz_ad(iFr,jFr,kFr) =-coeffGCMFreshD(iRow,n)
          ENDDO ! iRow

        ENDDO ! n

      END SELECT 

    ENDIF !  boundary_motion 

END SUBROUTINE set_solve_ad
