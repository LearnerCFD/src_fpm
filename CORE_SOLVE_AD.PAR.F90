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
!  Filename: CORE_SOLVE_AD.F90
!  Latest Modification: October 21, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     solve_ad()
!     calc_residual_ad(var,r,bcxvel,bcyvel,bczvel,ghostVar,resm)
!     itsolv_ad(var,r,bcxvel,bcyvel,bczvel,ghostVar)
!     itsolv_ad_x(var,r)
!     itsolv_ad_y(var,r)
!     itsolv_ad_z(var,r)
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE solve_ad(adsCommTime)

    USE global_parameters
    USE flow_parameters
    USE flow_parameters
    USE flow_arrays
    USE multiuse_arrays
    USE grid_arrays

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: adsCommTime

!   Local Variables
!   ---------------
    INTEGER           :: iter,i,j,k, ierr
    REAL(KIND=CGREAL) :: maxresu,maxresv,maxresw,maxres,restol_ad



    adsCommTime=0.d0

!   Initialize variables
!   --------------------
    restol_ad = restol
    iter      = 0
    maxres    = HUGE(1.0_CGREAL)
    maxresu   = 0.0_CGREAL
    maxresv   = 0.0_CGREAL
    maxresw   = 0.0_CGREAL


!   Compute initial residual
!   ------------------------
    CALL calc_residual_ad(u,nlu,bcxu,bcyu,bczu,uGhost,maxresu)!......................COMPLETE(SAMK)
    CALL calc_residual_ad(v,nlv,bcxv,bcyv,bczv,vGhost,maxresv)!......................COMPLETE(SAMK)
    IF (ndim == DIM_3D) CALL calc_residual_ad(w,nlw,bcxw,bcyw,bczw,wGhost,maxresw)!..COMPLETE(SAMK)

    maxres = MAX(DABS(maxresu),DABS(maxresv),DABS(maxresw))


    IF (monitorON) THEN
      WRITE(STDOUT,'(5X,A)') '--------------------------------------------------------------------------------------------------------------'
      WRITE(STDOUT,'(5X,A,3X,I4,4(1X,1PE19.11))') 'velocity Residuals : ',iter,maxres,maxresu,maxresv,maxresw
    END IF



!   ------------------------------------------------
!    Streamline calls to solver to SSM and GCM 
!    This is possible now with ghostVel formulation
!   ------------------------------------------------
  
  DO WHILE ((iter < itermax) .AND. (maxres > restol_ad))
      CALL itsolv_ad(u,nlu,bcxu,bcyu,bczu,uGhost)!........................................ONLY_PER(SAMK)
      CALL itsolv_ad(v,nlv,bcxv,bcyv,bczv,vGhost)!........................................ONLY_PER(SAMK)
      IF (ndim == DIM_3D) CALL itsolv_ad(w,nlw,bcxw,bcyw,bczw,wGhost)!....................ONLY_PER(SAMK)

      IF ( bcx1 == BC_TYPE_PERIODIC .OR. &                                      !Ehsan added for Periodic BC 
           bcy1 == BC_TYPE_PERIODIC .OR. &
           bcz1 == BC_TYPE_PERIODIC      ) THEN 
         CALL apply_periodic_bcvel
      END IF ! bcx1
 
      CALL set_outer_ghost_vel()!....................................................COMPLETE(SAMK)

      IF (boundary_formulation == GCM_METHOD ) THEN
        CALL GCM_GhostCell_Vel()!....................................................COMPLETE(SAMK)
       END IF ! boundary_formulation


      CALL calc_residual_ad(u,nlu,bcxu,bcyu,bczu,uGhost,maxresu)!....................COMPLETE(SAMK)
      CALL calc_residual_ad(v,nlv,bcxv,bcyv,bczv,vGhost,maxresv)!....................COMPLETE(SAMK)
      IF (ndim == DIM_3D) CALL calc_residual_ad(w,nlw,bcxw,bcyw,bczw,wGhost,maxresw)!COMPLETE(SAMK)

      maxres = MAX(DABS(maxresu),DABS(maxresv),DABS(maxresw))
 
      IF(advec_scheme == CRANK_NICOLSON2) CALL face_vel()!...........................COMPLETE(SAMK)
       
      iter = iter + 1


      IF (monitorON) THEN
        IF ( iter == itermax .AND. maxres > restol_ad ) THEN
          WRITE(STDOUT,'(5X,A,3X,I4,4(1X,1PE19.11))') 'Velocity Residuals : ',iter,maxres,maxresu,maxresv,maxresw
          WRITE(STDOUT,'(5X,A,I4,A,1X,1PE19.11)') '$$$$$$ Velocity did not converge in ',itermax,' iterations. Final residual =',maxres
        ELSE
          WRITE(STDOUT,'(5X,A,3X,I4,4(1X,1PE19.11))') 'Velocity Residuals : ',iter,maxres,maxresu,maxresv,maxresw
        ENDIF ! iter
      END IF
    ENDDO ! do while

    IF (monitorON) WRITE(STDOUT,'(5X,A)') '--------------------------------------------------------------------------------------------------------------'

END SUBROUTINE solve_ad
!---------------------------------------------------------------------



SUBROUTINE calc_residual_ad(var,r,bcxvel,bcyvel,bczvel,ghostVar,resm)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo
 
    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT (IN)    ::      var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL), INTENT (INOUT) :: ghostVar(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT (IN):: r, bcxvel, bcyvel, bczvel

    REAL(KIND=CGREAL), INTENT (OUT) :: resm

!   Loop Variables
!   --------------
    INTEGER :: i, j, k

!   Local Variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL) :: res,killForGCMFresh,killForSSMFresh
    REAL(KIND=CGREAL) :: oned, twod, half, half_dt  ! H. Luo
    REAL(KIND=CGREAL)  :: rhsImplicit
    REAL(KIND=CGREAL)  :: vele, velw, veln, vels, velf, velb

    INTEGER :: iloc, jloc, kloc



!   H. Luo start
!   ------------
    oned = 1.0_CGREAL
    twod = 2.0_CGREAL
    half = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
      half_dt  = 0.5_CGREAL * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
      half_dt  = 0.0_CGREAL
    ENDIF
!   ----------
!   H. Luo end

    resm = 0.0_CGREAL

    DO k=1,nzc
      KP   = k + 1
      KM   = k - 1

    DO j=1,nyc
    
      jG=L2GJ(j)
      JP   = j + 1
      JM   = j - 1

      DO i=1,nxc
        iG=L2GI(i)

        IP   = i + 1 
        IM   = i - 1 
 

!       Added by H. Luo
!       ---------------
        amx(i) = amx_ad(i,j,k)
        apx(i) = apx_ad(i,j,k)
        acx(i) =- ( amx(i) + apx(i) )

        amy(j) = amy_ad(i,j,k)
        apy(j) = apy_ad(i,j,k)
        acy(j) =- ( amy(j) + apy(j) )

        amz(k) = amz_ad(i,j,k)
        apz(k) = apz_ad(i,j,k)
        acz(k) =- ( amz(k) + apz(k) )

!       TO DO-FMN: Set coefficients in z-directions to zero for 2-D computations
!       ----------
!       End adding

      killForGCMFresh = (1.0_CGREAL-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
      killForSSMFresh = (1.0_CGREAL-(1.0_CGREAL-gcmFlag)*REAL(fresh_cell(i,j,k),KIND=CGREAL))

!     Compute velocity in ghostcells
!------------------------------------------------------------------------------      
      ghostVar(i+1,j,k) = gcmFlag*ghostVar(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
             *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*var(i,j,k) )/fx(iG+1)

      ghostVar(i-1,j,k) = gcmFlag*ghostVar(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                   *( bcxvel(i,j,k) - fx(iG)*var(i,j,k) )/(1.0_CGREAL-fx(iG))

      ghostVar(i,j+1,k) = gcmFlag*ghostVar(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
             *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*var(i,j,k) )/fy(jG+1)

      ghostVar(i,j-1,k) = gcmFlag*ghostVar(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
                   *( bcyvel(i,j,k) - fy(jG)*var(i,j,k) )/(1.0_CGREAL-fy(jG))

      ghostVar(i,j,k+1) = gcmFlag*ghostVar(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
             *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*var(i,j,k) )/fz(k+1)

      ghostVar(i,j,k-1) = gcmFlag*ghostVar(i,j,k-1)  + (1.0_CGREAL-gcmFlag) &
                   *( bczvel(i,j,k) - fz(k)*var(i,j,k) )/(1.0_CGREAL-fz(k))

      amx(i) = amx(i) *killForSSMFresh &
             + amx(i) *(1.0_CGREAL-ium(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
      apx(i) = apx(i) *killForSSMFresh &
             + apx(i) *(1.0_CGREAL-iup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
      amy(j) = amy(j) *killForSSMFresh &
             + amy(j) *(1.0_CGREAL-jum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
      apy(j) = apy(j) *killForSSMFresh &
             + apy(j) *(1.0_CGREAL-jup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
      amz(k) = amz(k) *killForSSMFresh &
             + amz(k) *(1.0_CGREAL-kum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
      apz(k) = apz(k) *killForSSMFresh &
             + apz(k) *(1.0_CGREAL-kup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

      acx(i) = 1.0_CGREAL*killForSSMFresh                         &
             + ( acx(i) + acy(j) + acz(k) )*killForGCMFresh

      res    = r(i,j,k) - var(i,j,k)*acx(i)                       &
                     - ( var(IM,j,k)*(1.0_CGREAL-ium(i,j,k))      &
                        + ghostVar(IM,j,k)*ium(i,j,k) )*amx(i)    &
                     - ( var(IP,j,k)*(1.0_CGREAL-iup(i,j,k))      &
                        + ghostVar(IP,j,k)*iup(i,j,k) )*apx(i)    &
                     - ( var(i,JM,k)*(1.0_CGREAL-jum(i,j,k))      &
                        + ghostVar(i,JM,k)*jum(i,j,k) )*amy(j)    &
                     - ( var(i,JP,k)*(1.0_CGREAL-jup(i,j,k))      &
                        + ghostVar(i,JP,k)*jup(i,j,k) )*apy(j)    &
                     - ( var(i,j,KM)*(1.0_CGREAL-kum(i,j,k))      &
                        + ghostVar(i,j,KM)*kum(i,j,k) )*amz(k)    &
                     - ( var(i,j,KP)*(1.0_CGREAL-kup(i,j,k))      &
                        + ghostVar(i,j,KP)*kup(i,j,k) )*apz(k)

!     Compute contribution of implicit CN term.
!      Currently rhsImplicit is run iteratively till convergence is reached.
!     Note: Algorithm does not take into account UPWINDING scheme
!     -----------------------------------------------------------------------------

      IF( advec_scheme == CRANK_NICOLSON1 .OR. &
          advec_scheme == CRANK_NICOLSON2      ) THEN

        vele =        fx(iG+1) * (      var(iP,j,k) *(oned -iup(i,j,k))   &
                               + ghostVar(iP,j,k)        *iup(i,j,k))   &
             + (oned -fx(iG+1))        *var(i ,j,k)

        velw =        fx(iG)          *var(i ,j,k)                       &
             + (oned -fx(iG)) * (      var(iM,j,k) *(oned -ium(i,j,k))   &
                               + ghostVar(iM,j,k)        *ium(i,j,k) )

        veln =        fy(jG+1) * (      var(i,jP,k) *(oned -jup(i,j,k))   &
                               + ghostVar(i,jP,k)        *jup(i,j,k) )  &
             + (oned -fy(jG+1))        *var(i,j ,k)

        vels =        fy(jG)         *var(i,j ,k)                        &
             + (oned -fy(jG)) *(      var(i,jM,k) *(oned -jum(i,j,k))    &
                              + ghostVar(i,jM,k)        *jum(i,j,k) )

        velf =        fz(k+1) *(      var(i,j,kP) *(oned -kup(i,j,k))    &
                              + ghostVar(i,j,kP)        *kup(i,j,k) )   &
             + (oned -fz(k+1))       *var(i,j,k )

        velb =        fz(k)         *var(i,j,k )                        &
             + (oned -fz(k)) *(      var(i,j,kM) *(oned -kum(i,j,k))    &
                              + ghostVar(i,j,kM)        *kum(i,j,k) )

        rhsImplicit = ( vele *face_ue(i,j,k) -velw *face_uw(i,j,k) ) *dxinv(iG) &
                    + ( veln *face_vn(i,j,k) -vels *face_vs(i,j,k) ) *dyinv(jG) &
                    + ( velf *face_wf(i,j,k) -velb *face_wb(i,j,k) ) *dzinv(k)

        rhsImplicit = -half_dt *rhsImplicit

        res = res +rhsImplicit
      ENDIF ! advec_scheme

!  Update Fresh cells for GCM approach
!   Note: This approach should result in improved performance
!         instead of an IF statement in triple do-loop
!------------------------------------------------------------------------------

      IF ( boundary_formulation == GCM_METHOD .AND. &
         fresh_cell(i,j,k) == 1                   ) THEN
        iFr = i; jFr = j; kFr = k;
        CALL GCM_correct_res_ad(iFr,jFr,kFr,var,res)!............COMPLETE(SAMK)
      ENDIF ! fresh_cell
       
      res = DABS(res)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      IF(res > resm) THEN
       iloc = i
       jloc = j
       kloc = k
      ENDIF
      resm = MAX(res,resm)
    ENDDO
    ENDDO
    ENDDO
    
#ifdef MPI
    CALL par_getMaxReal(resm)
#endif

END SUBROUTINE calc_residual_ad  
!----------------------------------------------------------



SUBROUTINE itsolv_ad(var,r,bcxvel,bcyvel,bczvel,ghostVar)

! -----------------------------------------------------------
!  currently coded as Line SOR with Gauss Siedel as smoother
! -----------------------------------------------------------

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT (INOUT) :: Var,ghostVar
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc)                        , INTENT (IN)    :: r,bcxvel,bcyvel,bczvel

!   Local Variables
!   ---------------
    INTEGER :: i,j,k
    
    REAL(KIND=CGREAL) :: myTime



!   Line solve in the x-direction
!   -----------------------------
    CALL itsolv_ad_x(var,r,bcxvel,bcyvel,bczvel,ghostVar)!................................ONLY_PER(SAMK)
#ifdef MPI
    CALL par_comm_var(var,nxc,nyc,nzc,Ngl,myTime)!...................................COMPLETE(SAMK)
#endif

!   Line solve in the y-direction
!   -----------------------------
    CALL itsolv_ad_y(var,r,bcxvel,bcyvel,bczvel,ghostVar)!................................ONLY_PER(SAMK)
#ifdef MPI
    CALL par_comm_var(var,nxc,nyc,nzc,Ngl,myTime)!...................................COMPLETE(SAMK)
#endif

!   Line solve in the z-direction    
!   -----------------------------
    IF (nDim == DIM_3D) THEN
      CALL itsolv_ad_z(var,r,bcxvel,bcyvel,bczvel,ghostVar)!................................ONLY_PER(SAMK)
#ifdef MPI
      CALL par_comm_var(var,nxc,nyc,nzc,Ngl,myTime)!...................................COMPLETE(SAMK)
#endif
    END IF

END SUBROUTINE itsolv_ad
!----------------------------------------------



SUBROUTINE itsolv_ad_x(var,r,bcxvel,bcyvel,bczvel,ghostVar)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT (INOUT) :: var,ghostVar
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc)                        , INTENT (IN)    :: r,bcxvel,bcyvel,bczvel

!   Loop Variables
!   --------------
    INTEGER :: i,j,k

!   Local Variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: iDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL) :: omega_ad, oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: killForGCMFresh,killForSSMFresh
    REAL(KIND=CGREAL) :: rhsImplicit
    REAL(KIND=CGREAL)  :: vele, velw, veln, vels, velf, velb
    REAL(KIND=CGREAL) :: ium_sd, iup_sd

    REAL(KIND=CGREAL) :: riblank(nxc,nyc,nzc)



    iDirection = 1
    omega_ad   = 1.0_CGREAL

!   H. Luo start
!   ------------
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = 0.0_CGREAL
    ENDIF
!   ----------
!   H. Luo end

    riblank(1:nxc,1:nyc,1:nzc) = REAL(1-iblank(1:nxc,1:nyc,1:nzc),KIND=CGREAL)

!   this loop is to transfer to RHS all the weights that do not fit on the x-tdma
!   Line solve in the x-direction
!   -----------------------------------------------------------------------------
    DO k=1,nzc
    DO j=1,nyc

      jG = L2GJ(j)

      JP = j+1 
      JM = j-1 
 
      KP = k+1 
      KM = k-1
 
      DO i=1,nxc
       
       IP = i+1
       IM = i-1
      
!       -------------------------------------------------------------------------------
!        Take care of left __subdomain__ BC:
!        If ium=1 then BC is already imposed, otherwise it should be imposed manually.
!       -------------------------------------------------------------------------------

!       SAMK: Lets BCR8TV
!       -----------------
        ium_sd=REAL(1/i,KIND=CGREAL)                    ! This is an INTEGER operation
        ium_sd=ium_sd*(1._CGREAL-ium(i,j,k))

!       -------------------------------------------------------------------------------
!        Take care of right __subdomain__ BC
!        If iup=1 then BC is already imposed, otherwise it should be imposed manually.
!       -------------------------------------------------------------------------------

!       SAMK: Lets BCR8TV
!       -----------------
        iup_sd=REAL(i/nxc,KIND=CGREAL)                  ! This is an INTEGER operation
        iup_sd=iup_sd*(1._CGREAL-iup(i,j,k))

        iG = L2GI(i)

        amx(i) = amx_ad(i,j,k)
        apx(i) = apx_ad(i,j,k)
        acx(i) = -( amx(i) + apx(i) )      

        amy(j) = amy_ad(i,j,k)
        apy(j) = apy_ad(i,j,k)
        acy(j) = -( amy(j) + apy(j) )     

        amz(k) = amz_ad(i,j,k)
        apz(k) = apz_ad(i,j,k)
        acz(k) = -( amz(k) + apz(k) )    

!       TO DO-FMN: Set coefficients in z-directions to zero for 2-D computations
!       ------------------------------------------------------------------------

!       killForGCMFresh = 0 if GCM = true and frescell== TRUE
!       killForGCMFresh = 1  otherwise
!       -----------------------------------------------------
        killForGCMFresh = (1.0_CGREAL-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))

!       killForSSMFresh = 0 if SSM = true and frescell== TRUE
!       killForSSMFresh = 1  otherwise
!       -----------------------------------------------------
        killForSSMFresh = (1.0_CGREAL-(1.0_CGREAL-gcmFlag)*REAL(fresh_cell(i,j,k),KIND=CGREAL))


!      Compute velocity in ghostcells
!------------------------------------------------------------------------------
        ghostVar(i+1,j,k) = gcmFlag*ghostVar(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*var(i,j,k) )/fx(iG+1)

        ghostVar(i-1,j,k) = gcmFlag*ghostVar(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - fx(iG)*var(i,j,k) )/(1.0_CGREAL-fx(iG))

        ghostVar(i,j+1,k) = gcmFlag*ghostVar(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*var(i,j,k) )/fy(jG+1)

        ghostVar(i,j-1,k) = gcmFlag*ghostVar(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - fy(jG)*var(i,j,k) )/(1.0_CGREAL-fy(jG))

        ghostVar(i,j,k+1) = gcmFlag*ghostVar(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*var(i,j,k) )/fz(k+1)

        ghostVar(i,j,k-1) = gcmFlag*ghostVar(i,j,k-1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - fz(k)*var(i,j,k) )/(1.0_CGREAL-fz(k))

        amx(i) = amx(i) *killForSSMFresh &
               + amx(i) *(1.0_CGREAL-ium(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apx(i) = apx(i) *killForSSMFresh & 
               + apx(i) *(1.0_CGREAL-iup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
        amy(j) = amy(j) *killForSSMFresh &
               + amy(j) *(1.0_CGREAL-jum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apy(j) = apy(j) *killForSSMFresh &
               + apy(j) *(1.0_CGREAL-jup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
        amz(k) = amz(k) *killForSSMFresh &
               + amz(k) *(1.0_CGREAL-kum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apz(k) = apz(k) *killForSSMFresh &
               + apz(k) *(1.0_CGREAL-kup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        rhs(i) = r(i,j,k) - ( ( var(i,JM,k)*(1.0_CGREAL-jum(i,j,k))  &
                               +ghostVar(i,JM,k)*jum(i,j,k) )*amy(j) & 
                             +( var(i,JP,k)*(1.0_CGREAL-jup(i,j,k))  &
                               +ghostVar(i,JP,k)*jup(i,j,k) )*apy(j) &
                             +( var(i,j,KM)*(1.0_CGREAL-kum(i,j,k))  &
                               +ghostVar(i,j,KM)*kum(i,j,k) )*amz(k) &
                             +( var(i,j,KP)*(1.0_CGREAL-kup(i,j,k))  &
                               +ghostVar(i,j,KP)*kup(i,j,k) )*apz(k))

        rhs(i) = rhs(i) - amx(i)*ghostVar(i-1,j,k)*ium(i,j,k) &
                        - amx(i)*     Var(i-1,j,k)*ium_sd     &
                        - apx(i)*ghostVar(i+1,j,k)*iup(i,j,k) &
                        - apx(i)*     Var(i+1,j,k)*iup_sd

!       Compute contribution of implicit CN term.
!        Currently rhsImplicit is run iteratively till convergence is reached.
!     Note: Algorithm does not take into account UPWINDING scheme
!-----------------------------------------------------------------------------

        IF( advec_scheme == CRANK_NICOLSON1 .OR. &
            advec_scheme == CRANK_NICOLSON2      ) THEN
          vele =        fx(iG+1) * (      var(iP,j,k) *(oned -iup(i,j,k))   &
                                 + ghostVar(iP,j,k)        *iup(i,j,k))   &
               + (oned -fx(iG+1))        *var(i ,j,k)

          velw =        fx(iG)          *var(i ,j,k)                       &
               + (oned -fx(iG)) * (      var(iM,j,k) *(oned -ium(i,j,k))   &
                                 + ghostVar(iM,j,k)        *ium(i,j,k) )

          veln =        fy(jG+1) * (      var(i,jP,k) *(oned -jup(i,j,k))   &
                                 + ghostVar(i,jP,k)        *jup(i,j,k) )  &
               + (oned -fy(jG+1))        *var(i,j ,k)

          vels =        fy(jG)         *var(i,j ,k)                        &
               + (oned -fy(jG)) *(      var(i,jM,k) *(oned -jum(i,j,k))    &
                                + ghostVar(i,jM,k)        *jum(i,j,k) )

          velf =        fz(k+1) *(      var(i,j,kP) *(oned -kup(i,j,k))    &
                                + ghostVar(i,j,kP)        *kup(i,j,k) )   &
               + (oned -fz(k+1))       *var(i,j,k )

          velb =        fz(k)         *var(i,j,k )                        &
               + (oned -fz(k)) *(      var(i,j,kM) *(oned -kum(i,j,k))    &
                                + ghostVar(i,j,kM)        *kum(i,j,k) )

          rhsImplicit = ( vele *face_ue(i,j,k) -velw *face_uw(i,j,k) ) *dxinv(iG) &
                      + ( veln *face_vn(i,j,k) -vels *face_vs(i,j,k) ) *dyinv(jG) &
                      + ( velf *face_wf(i,j,k) -velb *face_wb(i,j,k) ) *dzinv(k)

          rhs(i) = rhs(i) -half_dt *rhsImplicit
        ENDIF ! advec_scheme

!       Take into consideration iblank and Freshcells
!-----------------------------------------------------------------------------


        amx(i) = amx(i)*riblank(i,j,k)*(1.0_CGREAL-ium(i,j,k))*(1.0_CGREAL-ium_sd)
        apx(i) = apx(i)*riblank(i,j,k)*(1.0_CGREAL-iup(i,j,k))*(1.0_CGREAL-iup_sd)

        acx(i) = 1.0_CGREAL *killForSSMFresh &
               + ( acx(i)+acy(j)+acz(k) ) *riblank(i,j,k) * killForGCMFresh

        rhs(i) = rhs(i)*riblank(i,j,k)

!       Adding off-pentadiagonal elements for fresh cell to RHS
!       -------------------------------------------------------
        IF ( boundary_formulation == GCM_METHOD .AND. &
            fresh_cell(i,j,k) == 1                   ) THEN
          iFr = i; jFr = j; kFr = k;
          CALL GCM_correct_rhs_ad(iDirection,iFr,jFr,kFr,nxc,rhs(1:nxc),var)
        END iF
       
      ENDDO ! i

!     Calling TDMA Solver
!     -------------------
!!      IF ( bcx1 == BC_TYPE_PERIODIC .AND. &
!!           bcx2 == BC_TYPE_PERIODIC       ) THEN
!!        rhs(1)    = rhs(1)    - var(nx-1,j,k)*amx(1)
!!        rhs(nx-1) = rhs(nx-1) - var(1,j,k)*apx(nx-1)
!!      END IF ! bcx1


      CALL tdma(amx,acx,apx,rhs(1:nxc),dummy(1:nxc),1,nxc)


!     Loading var from dummy
!     ----------------------
      DO i=1,nxc
        var(i,j,k) = var(i,j,k) + omega_ad*(dummy(i)-var(i,j,k))
      ENDDO

    ENDDO ! j
    ENDDO ! k


END SUBROUTINE itsolv_ad_x
!----------------------------------------------



SUBROUTINE itsolv_ad_y(var,r,bcxvel,bcyvel,bczvel,ghostVar)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE GCM_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT (INOUT) :: var,ghostVar
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc)                        , INTENT (IN)    :: r,bcxvel,bcyvel,bczvel

!   Loop Variables
!   --------------
    INTEGER :: i,j,k

!   Local Variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: jDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL) :: omega_ad, oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: killForGCMFresh,killForSSMFresh
    REAL(KIND=CGREAL) :: rhsImplicit
    REAL(KIND=CGREAL)  :: vele, velw, veln, vels, velf, velb
    REAL(KIND=CGREAL) :: jum_sd, jup_sd

    REAL(KIND=CGREAL) :: riblank(nxc,nyc,nzc)



    jDirection = 2
    omega_ad = 1.0_CGREAL

!   H. Luo start
!   ------------
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = 0.0_CGREAL
    ENDIF
!   ----------
!   H. Luo end

    riblank(1:nxc,1:nyc,1:nzc) = REAL(1-iblank(1:nxc,1:nyc,1:nzc),KIND=CGREAL)

    DO k=1,nzc
    DO i=1,nxc

      iG = L2GI(i)
       
      IP = i + 1 
      IM = i - 1 
 
      KP = k + 1 
      KM = k - 1 

      DO j=1,nyc
      
       JP   = j + 1
       JM   = j - 1

      
!       -------------------------------------------------------------------------------
!        Take care of bottom __subdomain__ BC:
!        If jum=1 then BC is already imposed, otherwise it should be imposed manually.
!       -------------------------------------------------------------------------------

!       SAMK: Lets BCR8TV
!       -----------------
        jum_sd=REAL(1/j,KIND=CGREAL)                    ! This is an INTEGER operation
        jum_sd=jum_sd*(1._CGREAL-jum(i,j,k))

!       -------------------------------------------------------------------------------
!        Take care of top __subdomain__ BC
!        If jup=1 then BC is already imposed, otherwise it should be imposed manually.
!       -------------------------------------------------------------------------------

!       SAMK: Lets BCR8TV
!       -----------------
        jup_sd=REAL(j/nyc,KIND=CGREAL)                  ! This is an INTEGER operation
        jup_sd=jup_sd*(1._CGREAL-jup(i,j,k))

        jG=L2GJ(j)

        amx(i) = amx_ad(i,j,k)
        apx(i) = apx_ad(i,j,k)
        acx(i) =- ( amx(i) + apx(i) )                           

        amy(j) = amy_ad(i,j,k)
        apy(j) = apy_ad(i,j,k)
        acy(j) =- ( amy(j) + apy(j) )                          

        amz(k) = amz_ad(i,j,k)
        apz(k) = apz_ad(i,j,k)
        acz(k) =- ( amz(k) + apz(k) )   


!       TO DO-FMN: Set coefficients in z-directions to zero for 2-D computations
!       ------------------------------------------------------------------------

        killForGCMFresh = (1.0_CGREAL-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
        killForSSMFresh = (1.0_CGREAL-(1.0_CGREAL-gcmFlag)*REAL(fresh_cell(i,j,k),KIND=CGREAL))

!       Compute velocity in ghostcells
!------------------------------------------------------------------------------

        ghostVar(i+1,j,k) = gcmFlag*ghostVar(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*var(i,j,k) )/fx(iG+1)

        ghostVar(i-1,j,k) = gcmFlag*ghostVar(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - fx(iG)*var(i,j,k) )/(1.0_CGREAL-fx(iG))

        ghostVar(i,j+1,k) = gcmFlag*ghostVar(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*var(i,j,k) )/fy(jG+1)

        ghostVar(i,j-1,k) = gcmFlag*ghostVar(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - fy(jG)*var(i,j,k) )/(1.0_CGREAL-fy(jG))

        ghostVar(i,j,k+1) = gcmFlag*ghostVar(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*var(i,j,k) )/fz(k+1)

        ghostVar(i,j,k-1) = gcmFlag*ghostVar(i,j,k-1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - fz(k)*var(i,j,k) )/(1.0_CGREAL-fz(k))

        amx(i) = amx(i) *killForSSMFresh &
               + amx(i) *(1.0_CGREAL-ium(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apx(i) = apx(i) *killForSSMFresh & 
               + apx(i) *(1.0_CGREAL-iup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
        amy(j) = amy(j) *killForSSMFresh &
               + amy(j) *(1.0_CGREAL-jum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apy(j) = apy(j) *killForSSMFresh &
               + apy(j) *(1.0_CGREAL-jup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)
       
        amz(k) = amz(k) *killForSSMFresh &
               + amz(k) *(1.0_CGREAL-kum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apz(k) = apz(k) *killForSSMFresh &
               + apz(k) *(1.0_CGREAL-kup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        rhs(j) = r(i,j,k) -( (  var(IM,j,k)*(1.0_CGREAL-ium(i,j,k))  &
                               +ghostVar(IM,j,k)*ium(i,j,k) )*amx(i) &
                             +( var(IP,j,k)*(1.0_CGREAL-iup(i,j,k))  &
                               +ghostVar(IP,j,k)*iup(i,j,k) )*apx(i) &
                             +( var(i,j,KM)*(1.0_CGREAL-kum(i,j,k))  &
                              +ghostVar(i,j,KM)*kum(i,j,k) )*amz(k) &
                             +( var(i,j,KP)*(1.0_CGREAL-kup(i,j,k))  & 
                               +ghostVar(i,j,KP)*kup(i,j,k) )*apz(k)   )  

        rhs(j) = rhs(j) - amy(j)*ghostVar(i,j-1,k)*jum(i,j,k) &
                        - amy(j)*     Var(i,j-1,k)*jum_sd     &
                        - apy(j)*ghostVar(i,j+1,k)*jup(i,j,k) &
                        - apy(j)*     Var(i,j+1,k)*jup_sd

!       Compute contribution of implicit CN term.
!        Currently rhsImplicit is run iteratively till convergence is reached.
!     Note: Algorithm does not take into account UPWINDING scheme
!-----------------------------------------------------------------------------

        IF( advec_scheme == CRANK_NICOLSON1 .OR. &
            advec_scheme == CRANK_NICOLSON2      ) THEN
          vele =        fx(iG+1) * (      var(iP,j,k) *(oned -iup(i,j,k))   &
                                 + ghostVar(iP,j,k)        *iup(i,j,k))   &
               + (oned -fx(iG+1))        *var(i ,j,k)

          velw =        fx(iG)          *var(i ,j,k)                       &
               + (oned -fx(iG)) * (      var(iM,j,k) *(oned -ium(i,j,k))   &
                                 + ghostVar(iM,j,k)        *ium(i,j,k) )

          veln =        fy(jG+1) * (      var(i,jP,k) *(oned -jup(i,j,k))   &
                                 + ghostVar(i,jP,k)        *jup(i,j,k) )  &
               + (oned -fy(jG+1))        *var(i,j ,k)

          vels =        fy(jG)         *var(i,j ,k)                        &
               + (oned -fy(jG)) *(      var(i,jM,k) *(oned -jum(i,j,k))    &
                                + ghostVar(i,jM,k)        *jum(i,j,k) )

          velf =        fz(k+1) *(      var(i,j,kP) *(oned -kup(i,j,k))    &
                                + ghostVar(i,j,kP)        *kup(i,j,k) )   &
               + (oned -fz(k+1))       *var(i,j,k)

          velb =        fz(k)         *var(i,j,k )                        &
               + (oned -fz(k)) *(      var(i,j,kM) *(oned -kum(i,j,k))    &
                                + ghostVar(i,j,kM)        *kum(i,j,k) )

          rhsImplicit = ( vele *face_ue(i,j,k) -velw *face_uw(i,j,k) ) *dxinv(iG) &
                      + ( veln *face_vn(i,j,k) -vels *face_vs(i,j,k) ) *dyinv(jG) &
                      + ( velf *face_wf(i,j,k) -velb *face_wb(i,j,k) ) *dzinv(k)

          rhs(j) = rhs(j) -half_dt *rhsImplicit
        ENDIF ! advec_scheme

!       Take into consideration iblank and Freshcells
!-----------------------------------------------------------------------------

        amy(j) = amy(j)*riblank(i,j,k)*(1.0_CGREAL-jum(i,j,k))*(1.0_CGREAL-jum_sd)
        apy(j) = apy(j)*riblank(i,j,k)*(1.0_CGREAL-jup(i,j,k))*(1.0_CGREAL-jup_sd)

        acy(j) = 1.0_CGREAL *killForSSMFresh &
               + ( acx(i)+acy(j)+acz(k) ) * riblank(i,j,k) *killForGCMFresh 

        rhs(j) = rhs(j)*riblank(i,j,k) 
        
!       Adding off-pentadiagonal elements for fresh cell to RHS
!       -------------------------------------------------------
        IF ( boundary_formulation == GCM_METHOD .AND. &
             fresh_cell(i,j,k) == 1                   ) THEN
          iFr = i; jFr = j; kFr = k;
          CALL GCM_correct_rhs_ad(jDirection,iFr,jFr,kFr,nyc,rhs(1:nyc),var)
        ENDIF ! fresh_cell
       
      ENDDO ! j

!     Calling TDMA solver
!     -------------------
!!      IF ( bcy1 == BC_TYPE_PERIODIC .AND. &
!!           bcy2 == BC_TYPE_PERIODIC       ) THEN
!!        rhs(1) = rhs(1) - var(i,ny-1,k)*amy(1)
!!        rhs(ny-1) = rhs(ny-1) - var(i,1,k)*apy(ny-1)
!!      END IF ! bcy1


      CALL tdma(amy,acy,apy,rhs(1:nyc),dummy(1:nyc),1,nyc)
        

!     Loading var from dummy
!     ----------------------
      DO j=1,nyc
        var(i,j,k) = var(i,j,k) + omega_ad*(dummy(j)-var(i,j,k))
      ENDDO ! j

    ENDDO ! i
    ENDDO ! k
END SUBROUTINE itsolv_ad_y
!----------------------------------------------------------



SUBROUTINE itsolv_ad_z(var,r,bcxvel,bcyvel,bczvel,ghostVar,iCyclic)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE flow_arrays  ! H. Luo

    IMPLICIT NONE

!   Parameters
!   ----------
    INTEGER, INTENT(IN) :: iCyclic

    REAL(KIND=CGREAL), DIMENSION(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1), INTENT (INOUT) :: var,ghostVar
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc)                        , INTENT (IN)    :: r,bcxvel,bcyvel,bczvel

!   Loop Variables
!   --------------
    INTEGER :: i,j,k

!   Local Variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: kDirection, iFr, jFr, kFr
    INTEGER :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL) :: omega_ad, oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: killForGCMFresh, killForSSMFresh
    REAL(KIND=CGREAL) :: rhsImplicit
    REAL(KIND=CGREAL) :: vele, velw, veln, vels, velf, velb
    REAL(KIND=CGREAL) :: riblank(nxc,nyc,nzc)



    kDirection = 3
    omega_ad = 1.0_CGREAL

!   H. Luo start
!   ------------
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON1 .or. advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE IF (advec_scheme == ADAMS_BASHFORTH2) THEN
       half_dt  = 0.0_CGREAL
    ENDIF
!   ----------
!   H. Luo end

    riblank(1:nxc,1:nyc,1:nzc) = REAL(1-iblank(1:nxc,1:nyc,1:nzc),KIND=CGREAL)
    
    DO i=1,nxc
    DO j=1,nyc   
      
      iG=L2GI(i)
      jG=L2GJ(j)

      IP   = i + 1 
      IM   = i - 1 
 
      JP   = j + 1 
      JM   = j - 1 

      DO k=1,nzc
      
        KP   = k + 1
        KM   = k - 1

        amx(i) = amx_ad(i,j,k)
        apx(i) = apx_ad(i,j,k)
        acx(i) =- ( amx(i) + apx(i) )                           

        amy(j) = amy_ad(i,j,k)
        apy(j) = apy_ad(i,j,k)
        acy(j) =- ( amy(j) + apy(j) )                          

        amz(k) = amz_ad(i,j,k)
        apz(k) = apz_ad(i,j,k)
        acz(k) =- ( amz(k) + apz(k) )                         


!       TO DO-FMN: Set coefficients in z-directions to zero for 2-D computations
!       ------------------------------------------------------------------------

        killForGCMFresh = (1.0_CGREAL-gcmFlag*REAL(fresh_cell(i,j,k),KIND=CGREAL))
        killForSSMFresh = (1.0_CGREAL-(1.0_CGREAL-gcmFlag)*REAL(fresh_cell(i,j,k),KIND=CGREAL))

!       Compute velocity in ghostcells
!------------------------------------------------------------------------------

        ghostVar(i+1,j,k) = gcmFlag*ghostVar(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*var(i,j,k) )/fx(iG+1)

        ghostVar(i-1,j,k) = gcmFlag*ghostVar(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcxvel(i,j,k) - fx(iG)*var(i,j,k) )/(1.0_CGREAL-fx(iG))

        ghostVar(i,j+1,k) = gcmFlag*ghostVar(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*var(i,j,k) )/fy(jG+1)

        ghostVar(i,j-1,k) = gcmFlag*ghostVar(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
                 *( bcyvel(i,j,k) - fy(jG)*var(i,j,k) )/(1.0_CGREAL-fy(jG))

        ghostVar(i,j,k+1) = gcmFlag*ghostVar(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*var(i,j,k) )/fz(k+1)

        ghostVar(i,j,k-1) = gcmFlag*ghostVar(i,j,k-1)  + (1.0_CGREAL-gcmFlag) &
                 *( bczvel(i,j,k) - fz(k)*var(i,j,k) )/(1.0_CGREAL-fz(k))

        amx(i) = amx(i) *killForSSMFresh &
               + amx(i) *(1.0_CGREAL-ium(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apx(i) = apx(i) *killForSSMFresh &
               + apx(i) *(1.0_CGREAL-iup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        amy(j) = amy(j) *killForSSMFresh &
               + amy(j) *(1.0_CGREAL-jum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apy(j) = apy(j) *killForSSMFresh &
               + apy(j) *(1.0_CGREAL-jup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        amz(k) = amz(k) *killForSSMFresh &
               + amz(k) *(1.0_CGREAL-kum(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        apz(k) = apz(k) *killForSSMFresh &
               + apz(k) *(1.0_CGREAL-kup(i,j,k)) *(1.0_CGREAL-killForSSMFresh)

        rhs(k) = r(i,j,k) - ( ( var(i,JM,k)*(1.0_CGREAL-jum(i,j,k))  &
                               +ghostVar(i,JM,k)*jum(i,j,k) )*amy(j) &
                             +( var(i,JP,k)*(1.0_CGREAL-jup(i,j,k))  &
                               +ghostVar(i,JP,k)*jup(i,j,k) )*apy(j) &
                             +( var(IM,j,k)*(1.0_CGREAL-ium(i,j,k))  &
                               +ghostVar(IM,j,k)*ium(i,j,k) )*amx(i) &
                             +( var(IP,j,k)*(1.0_CGREAL-iup(i,j,k))  &
                               +ghostVar(IP,j,k)*iup(i,j,k) )*apx(i))

        rhs(k) = rhs(k) - amz(k)*ghostVar(i,j,k-1)*kum(i,j,k) &
                        - apz(k)*ghostVar(i,j,k+1)*kup(i,j,k)

!       Compute contribution of implicit CN term.
!        Currently rhsImplicit is run iteratively till convergence is reached.
!     Note: Algorithm does not take into account UPWINDING scheme
!-----------------------------------------------------------------------------

        IF( advec_scheme == CRANK_NICOLSON1 .OR. &
            advec_scheme == CRANK_NICOLSON2      ) THEN
          vele =        fx(iG+1) * (      var(iP,j,k) *(oned -iup(i,j,k))   &
                                 + ghostVar(iP,j,k)        *iup(i,j,k))   &
               + (oned -fx(iG+1))        *var(i ,j,k)

          velw =        fx(iG)          *var(i ,j,k)                       &
               + (oned -fx(iG)) * (      var(iM,j,k) *(oned -ium(i,j,k))   &
                                 + ghostVar(iM,j,k)        *ium(i,j,k) )

          veln =        fy(jG+1) * (      var(i,jP,k) *(oned -jup(i,j,k))   &
                                 + ghostVar(i,jP,k)        *jup(i,j,k) )  &
               + (oned -fy(jG+1))        *var(i,j ,k)

          vels =        fy(jG)         *var(i,j ,k)                        &
               + (oned -fy(jG)) *(      var(i,jM,k) *(oned -jum(i,j,k))    &
                                + ghostVar(i,jM,k)        *jum(i,j,k) )

          velf =        fz(k+1) *(      var(i,j,kP) *(oned -kup(i,j,k))    &
                                + ghostVar(i,j,kP)        *kup(i,j,k) )   &
               + (oned -fz(k+1))       *var(i,j,k )

          velb =        fz(k)         *var(i,j,k )                        &
               + (oned -fz(k)) *(      var(i,j,kM) *(oned -kum(i,j,k))    &
                                + ghostVar(i,j,kM)        *kum(i,j,k) )

          rhsImplicit = ( vele *face_ue(i,j,k) -velw *face_uw(i,j,k) ) *dxinv(iG) &
                      + ( veln *face_vn(i,j,k) -vels *face_vs(i,j,k) ) *dyinv(jG) &
                      + ( velf *face_wf(i,j,k) -velb *face_wb(i,j,k) ) *dzinv(k)

          rhs(k) = rhs(k) -half_dt *rhsImplicit
        ENDIF ! advec_scheme

!       Take into consideration iblank and Freshcells
!-----------------------------------------------------------------------------

        amz(k) = amz(k)*riblank(i,j,k)*(1.0_CGREAL-kum(i,j,k))
        apz(k) = apz(k)*riblank(i,j,k)*(1.0_CGREAL-kup(i,j,k))

        acz(k) = 1.0_CGREAL *killForSSMFresh &
               + ( acx(i)+acy(j)+acz(k) ) * riblank(i,j,k) *killForGCMFresh 
 
        rhs(k) = rhs(k)*riblank(i,j,k) 
       
       
        IF ( boundary_formulation == GCM_METHOD .AND. &
             fresh_cell(i,j,k) == 1                   ) THEN
          iFr = i; jFr = j; kFr = k;
          CALL GCM_correct_rhs_ad(kDirection,iFr,jFr,kFr,nzc,rhs(1:nzc),var)
        ENDIF ! fresh_cell
       
      ENDDO ! k

!     Calling TDMA solver
!     -------------------
!!      IF ( bcz1 == BC_TYPE_PERIODIC .AND. &
!!           bcz2 == BC_TYPE_PERIODIC       ) THEN
!!        rhs(1) = rhs(1) - var(i,j,nz-1)*amz(1)
!!        rhs(nz-1) = rhs(nz-1) - var(i,j,1)*apz(nz-1)
!!      END IF ! bcz1

      IF ( bcz1 == BC_TYPE_PERIODIC .AND. &                      !Ehsan added for Periodic BC
           bcz2 == BC_TYPE_PERIODIC       ) THEN                 
        rhs(1) = rhs(1) - var(i,j,nzc)*amz(1)                
        rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)         
      END IF ! bcz1                                            

      CALL tdma(amz,acz,apz,rhs(1:nzc),dummy(1:nzc),1,nzc)

!     Loading var from dummy
!     ----------------------
      DO k=1,nzc
        var(i,j,k) = var(i,j,k) + omega_ad*(dummy(k)-var(i,j,k))
      ENDDO !k

    ENDDO !j 
    ENDDO !i
    
END SUBROUTINE itsolv_ad_z
