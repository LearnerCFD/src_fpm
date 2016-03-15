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
!  Filename: CORE_RHS_AD.PAR.F90
!  Latest Modification: September 010, 2008 (ver. P1.2.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     rhs_advec_diff() 
!       rhs_advec()
!       rhs_loadAB(nlvel,nlvelOld)
!       rhs_diff()
!       rhs_vankan()
!     rhs_adjust_fresh()
!     rhs_adjust2D()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE rhs_advec_diff() 

!-------------------------------------------------------------------------------
!          n       {       n      n-1}             {    n         n       n  }
!  RHS  = u   - dt { 3/2 NL -1/2 NL  } + dt(1/2Re) {am u    + ac u  + ap u   }
!                  {       i       i }             {  i i-1     i i     i i+1} 
!-------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE nlold_arrays
    USE pressure_arrays
    USE solver_ad_arrays
    USE turb_parameters
    USE turb_global_parameters
    USE turb_arrays

    IMPLICIT NONE
    
    integer::i,j,k,ierr




!   Non-linear convection term
!   --------------------------
    CALL rhs_advec(u,nlu,bcxu,bcyu,bczu,uGhost)   
    CALL rhs_advec(v,nlv,bcxv,bcyv,bczv,vGhost)    
    IF ( nDim == DIM_3D ) THEN
      CALL rhs_advec(w,nlw,bcxw,bcyw,bczw,wGhost)
    END IF ! nDim

!   Load Adams-Bashforth terms
!   --------------------------
    CALL rhs_loadAB(nlu,nluOld)
    CALL rhs_loadAB(nlv,nlvOld)
    IF ( nDim == DIM_3D ) THEN
      CALL rhs_loadAB(nlw,nlwOld)
    END IF ! nDim

!   Diffusion terms
!   ---------------
    CALL rhs_diff(u,nlu,bcxu,bcyu,bczu,uGhost)    
    CALL rhs_diff(v,nlv,bcxv,bcyv,bczv,vGhost)    
    IF ( nDim == DIM_3D ) THEN
      CALL rhs_diff(w,nlw,bcxw,bcyw,bczw,wGhost)
    END IF ! nDim

!   Cross LES-based Diffusion terms
!   --------------------------------
    IF ( nDim == DIM_3D       .AND. &
         turbActive == ACTIVE .AND. &
         turbModel  > TURB_NOMODEL  ) THEN
          CALL turb_rhs_diff(nlu,nlv,nlw)
    END IF ! nDim

!   Invoke Van-kan formalism
!   ------------------------
    IF (frac_step_type == VAN_KAN) CALL rhs_vankan

!---------------------------------------------------------------------

CONTAINS



  SUBROUTINE rhs_advec(vel,nlvel,bcxvel,bcyvel,bczvel,ghostVel) 
    


    IMPLICIT NONE

!   parameters
!   ----------   
    REAL(KIND=CGREAL), INTENT(IN)      ::      vel(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL), INTENT(INOUT)   :: ghostVel(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT(IN)  :: bcxvel, bcyvel, bczvel
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT(OUT) :: nlvel

!   loop  variables
!   ---------------
    INTEGER           :: i,j,k
    INTEGER           :: imm,jmm,kmm         ! Added by Rupesh
    INTEGER           :: ipp,jpp,kpp         ! used in 2nd Upwinding

!   local variables
!   ---------------
    INTEGER :: iG, jG
    INTEGER :: immG, jmmG
    INTEGER :: ippG, jppG

    REAL(KIND=CGREAL) :: edxWeightm1, edxWeightm2, edxWeightp1, edxWeightp2
    REAL(KIND=CGREAL) :: wdxWeightm1, wdxWeightm2, wdxWeightp1, wdxWeightp2

    REAL(KIND=CGREAL) :: ndxWeightm1, ndxWeightm2, ndxWeightp1, ndxWeightp2
    REAL(KIND=CGREAL) :: sdxWeightm1, sdxWeightm2, sdxWeightp1, sdxWeightp2

    REAL(KIND=CGREAL) :: bdxWeightm1, bdxWeightm2, bdxWeightp1, bdxWeightp2
    REAL(KIND=CGREAL) :: fdxWeightm1, fdxWeightm2, fdxWeightp1, fdxWeightp2

    REAL(KIND=CGREAL) :: Usign, UsignP, Vsign, VsignP, Wsign, WsignP

    REAL(KIND=CGREAL) :: vele,velw,veln,vels,velf,velb
    REAL(KIND=CGREAL) :: vele_Up,velw_Up,veln_Up,vels_Up,velf_Up,velb_Up ! Added by Rupesh 



!---------------------------------------------------------------------------------------------------
! Convective terms:
! nlvel = d(vel*U)/dx + d(vel*V)/dy + d(vel*W)/dz
!---------------------------------------------------------------------------------------------------
!   The following decription was added by Rupesh and pertains to 2nd Upwind scheme for 
!   convected face velocities.  
!    ___________________________________________________
!   |          |         |         |         |          |
!   |          |         |         |         |          |
!   |    o     |    o   w+    o    +e   o    |    o     |
!   |   WW     |    W    |    P    |    E    |    EE    |
!   |__________|_________|_________|_________|__________|
!
!   |<--S_uu-->|<--S_u-->|<--S_c-->|<--S_d-->|<--S_dd-->|  1.Flow: Left --> Right
!                                           
!   |<--S_dd-->|<--S_d-->|<--S_c-->|<--S_u-->|<--S_uu-->|  2.Flow: Left <-- Right 
!
!   Subscripts 'u' and 'd' in S refer to upstream and downstream resply.
!
!   LOGIC: 1. If the flow is from left to right across a face (say, face e), then
!
!                 |     { S_u + 2S_c}      {    S_c    }
!              u_e|   = {-----------}u_P - {-----------}u_W
!                 |up   { S_u + S_c }      { S_u + S_c } 
!
!          2. If the flow is from right to left across a face (say, face e), then
!
!                 |     { S_uu + 2S_u}      {     S_u    }
!              u_e|   = {------------}u_E - {------------}u_EE
!                 |up   { S_uu + S_u }      { S_uu + S_u } 
!        
!   - It should be noted that for u_w, the above formulae are still valid, provided the stencil  
!     is offset by one cell.
!   - These formulae are derived from:  u_face = u + (Grad u)(dot)(dS), where 
!
!     'u' and 'Grad u' are cellcentered value and its gradient in the upstream cell resply, 
!      and dS is the displacement vector from the upstream cell-centroid to the face centroid.
!      'Grad u' is approximated by upwind differencing based on the direction of the wind.
!
!----------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! x-direction 
!------------------------------------------------------------------------------
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc

      imm = MAX(i-2,myImin-1)     
      ipp = MIN(i+2,myImax+1)

      iG  =L2GI(i)
      immG=L2GI(imm)
      ippG=L2GI(ipp)
        
      USign  = SIGN(1.0_CGREAL,face_uw(i,j,k)) 
      USignP = SIGN(1.0_CGREAL,face_ue(i,j,k))

      ghostVel(i+1,j,k) = gcmFlag*ghostVel(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
             *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*vel(i,j,k) )/fx(iG+1)
      ghostVel(i-1,j,k) = gcmFlag*ghostVel(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                   *( bcxvel(i,j,k) - fx(iG)*vel(i,j,k) )/(1.0_CGREAL-fx(iG))

      vele = fx(iG+1)*( vel(i+1,j,k)*(1.0_CGREAL-iup(i,j,k))                &
                     + ghostVel(i+1,j,k)*iup(i,j,k))                        &
                     + ( 1.0_CGREAL-fx(iG+1) ) *vel(i,j,k)   

      velw = fx(iG)*vel(i,j,k)                                              &
            + (1.0_CGREAL-fx(iG))*( vel(i-1,j,k)*(1.0_CGREAL-ium(i,j,k))    &
            + ghostVel(i-1,j,k)*ium(i,j,k) ) 

!     2nd Upwind differencing added by Rupesh 
!     ----------------------------------------
      edxWeightm1 = ( dxinv(iG-1)+2.0_CGREAL*dxinv(iG) ) /   &
                                ( dxinv(iG-1)+dxinv(iG) )
      edxWeightm2 = dxinv(iG) / ( dxinv(iG-1)+dxinv(iG) )
      edxWeightp1 = ( dxinv(ippG)+2.0_CGREAL*dxinv(iG+1) ) / &
                                ( dxinv(ippG)+dxinv(iG+1) )
      edxWeightp2 = dxinv(iG+1) / ( dxinv(ippG)+dxinv(iG+1) )

      wdxWeightm1 = ( dxinv(immG)+2.0_CGREAL*dxinv(iG-1) ) / &
                                ( dxinv(immG)+dxinv(iG-1) )
      wdxWeightm2 = dxinv(iG-1) / ( dxinv(immG)+dxinv(iG-1) )
      wdxWeightp1 = ( dxinv(iG+1)+2.0_CGREAL*dxinv(iG) ) /   &
                              ( dxinv(iG+1)+dxinv(iG) )
      wdxWeightp2 = dxinv(iG) / ( dxinv(iG+1)+dxinv(iG) )

      vele_Up = ( ( 1.0_CGREAL+UsignP )                                       &
                   *( iumm(i,j,k)*vele                                        &
                            + (1.0_CGREAL-iumm(i,j,k) )*                      &
                  ( edxWeightm1*vel(i,j,k)   - edxWeightm2*vel(i-1,j,k) ) )   &
               + ( 1.0_CGREAL-UsignP )                                        &
                 *( iupp(i,j,k)*vele                                          &
                            + (1.0_CGREAL-iupp(i,j,k) )*                      &
                  ( edxWeightp1*vel(i+1,j,k) - edxWeightp1*vel(ipp,j,k) ) )   &
               )*0.5_CGREAL*(1.0_CGREAL-iup(i,j,k)) + bcxvel(i,j,k)*iup(i,j,k) 

      velw_Up = ( ( 1.0_CGREAL+Usign )                                        &
                 *( iumm(i,j,k)*velw                                          &
                            + (1.0_CGREAL-iumm(i,j,k) )*                      &
                  ( wdxWeightm1*vel(i-1,j,k) - wdxWeightm2*vel(imm,j,k) ) )   &
               + ( 1.0_CGREAL-Usign )                                         &
                 *( iupp(i,j,k)*velw                                          &
                            + (1.0_CGREAL-iupp(i,j,k) )*                      &
                  ( wdxWeightp1*vel(i,j,k)   - wdxWeightp2*vel(i+1,j,k) ) )   &
               )*0.5_CGREAL*(1.0_CGREAL-ium(i,j,k)) + bcxvel(i,j,k)*ium(i,j,k)  

!     ----------------------------------------------------------------------
!      nlvel(i,j,k) = ( vele*face_u(i+1,j,k) -velw*face_u(i,j,k) )*dxinv(i)      
!      Original CDS definition commented by Rupesh
!     ----------------------------------------------------------------------

      nlvel(i,j,k) = (  ( (1.0_CGREAL - alfa)*vele + alfa*vele_Up )*       &
                                                  face_ue(i,j,k)           & 
                      - ( (1.0_CGREAL - alfa)*velw + alfa*velw_Up )*       &
                                                  face_uw(i,j,k)           &
                    )*dxinv(iG)           ! Hybrid definition added by Rupesh
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! y-direction
!------------------------------------------------------------------------------
     
    DO k = 1,nzc
    DO j = 1,nyc

    jmm = MAX(j-2,myJmin-1)
    jpp = MIN(j+2,myJmax+1)
      
    jG  =L2GJ(j)
    jmmG=L2GJ(jmm)
    jppG=L2GJ(jpp)

    DO i = 1,nxc
      VSign  = SIGN(1.0_CGREAL,face_vs(i,j,k)) 
      VSignP = SIGN(1.0_CGREAL,face_vn(i,j,k))
    
      ghostVel(i,j+1,k) = gcmFlag*ghostVel(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
              *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*vel(i,j,k) )/fy(jG+1)

      ghostVel(i,j-1,k) = gcmFlag*ghostVel(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
              *( bcyvel(i,j,k) - fy(jG)*vel(i,j,k) )/(1.0_CGREAL-fy(jG))

      veln =     fy(jG+1)*( vel(i,j+1,k)*(1.0_CGREAL-jup(i,j,k))  &
                          + ghostVel(i,j+1,k)*jup(i,j,k)        ) &
          +  ( 1.0_CGREAL-fy(jG+1) ) *vel(i,j,k)   

      vels =  fy(jG)*vel(i,j,k)                              &
          +  (1.0_CGREAL-fy(jG))*( vel(i,j-1,k)*(1.0_CGREAL-jum(i,j,k))   &
                                  + ghostVel(i,j-1,k)*jum(i,j,k) ) 

!     2nd Upwind differencing added by Rupesh
!     ---------------------------------------
      ndxWeightm1 = ( dyinv(jG-1)+2.0_CGREAL*dyinv(jG) ) /   &
                              ( dyinv(jG-1)+dyinv(jG) )
      ndxWeightm2 = dyinv(jG) / ( dyinv(jG-1)+dyinv(jG) )
      ndxWeightp1 = ( dyinv(jppG)+2.0_CGREAL*dyinv(jG+1) ) / &
                                ( dyinv(jppG)+dyinv(jG+1) )
      ndxWeightp2 = dyinv(jG+1) / ( dyinv(jppG)+dyinv(jG+1) )

      sdxWeightm1 = ( dyinv(jmmG)+2.0_CGREAL*dyinv(jG-1) ) / &
                                ( dyinv(jmmG)+dyinv(jG-1) ) 
      sdxWeightm2 = dyinv(jG-1) / ( dyinv(jmmG)+dyinv(jG-1) ) 
      sdxWeightp1 = ( dyinv(jG+1)+2.0_CGREAL*dyinv(jG) ) /   &
                                ( dyinv(jG+1)+dyinv(jG) )
      sdxWeightp2 = dyinv(jG) / ( dyinv(jG+1)+dyinv(jG) )

      veln_Up = ( ( 1.0_CGREAL+VsignP )                                        &
                 *( jumm(i,j,k)*veln                                          &
                            + (1.0_CGREAL-jumm(i,j,k) )*                      &
                  (  ndxWeightm1*vel(i,j,k)   - ndxWeightm2*vel(i,j-1,k) ) )  &
               + ( 1.0_CGREAL-VsignP )                                        &
                 *( jupp(i,j,k)*veln                                          & 
                            + (1.0_CGREAL-jupp(i,j,k) )*                      &
                  (  ndxWeightp1*vel(i,j+1,k) - ndxWeightp2*vel(i,jpp,k) ) )  &
               )*0.5_CGREAL*(1.0_CGREAL-jup(i,j,k)) + bcyvel(i,j,k)*jup(i,j,k)  

      vels_Up = ( ( 1.0_CGREAL+Vsign )                                         &
                 *( jumm(i,j,k)*vels                                          &
                            + (1.0_CGREAL-jumm(i,j,k) )*                      &
                  ( sdxWeightm1*vel(i,j-1,k) - sdxWeightm2*vel(i,jmm,k) ) )   &
               + ( 1.0_CGREAL-Vsign )                                         &
                 *( jupp(i,j,k)*vels                                          &
                            + (1.0_CGREAL-jupp(i,j,k) )*                      &
                  ( sdxWeightp1*vel(i,j,k)   - sdxWeightp2*vel(i,j+1,k) ) )   &
               )*0.5_CGREAL*(1.0_CGREAL-jum(i,j,k)) + bcyvel(i,j,k)*jum(i,j,k) 

!   -----------------------------------------------------------------------------
!    nlvel(i,j,k) = nlvel(i,j,k)                                             & 
!                     + ( veln*face_v(i,j+1,k) -vels*face_v(i,j,k) )*dyinv(j)  &
!      Original CDS definition commented by Rupesh
!   -----------------------------------------------------------------------------

      nlvel(i,j,k) = nlvel(i,j,k) +                                           &
                   (   ( (1.0_CGREAL - alfa)*veln + alfa*veln_Up )*          &
                                                   face_vn(i,j,k)           &
                     - ( (1.0_CGREAL - alfa)*vels + alfa*vels_Up )*          &
                                                   face_vs(i,j,k)             &
                   )*dyinv(jG)            ! Hybrid definition added by Rupesh
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

!------------------------------------------------------------------------------
! z-direction
!------------------------------------------------------------------------------
    IF ( nDim == DIM_3D ) THEN     
      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc  
        WSign  = SIGN(1.0_CGREAL,face_wb(i,j,k)) 
        WSignP = SIGN(1.0_CGREAL,face_wf(i,j,k))

        kmm = MAX(k-2,0)
        kpp = MIN(k+2,nzc+1)

        ghostVel(i,j,k+1) = gcmFlag*ghostVel(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
              *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*vel(i,j,k) )/fz(k+1)

        ghostVel(i,j,k-1) = gcmFlag*ghostVel(i,j,k-1)  + (1.0_CGREAL-gcmFlag) &
              *( bczvel(i,j,k) - fz(k)*vel(i,j,k) )/(1.0_CGREAL-fz(k))

        velf =   fz(k+1)*( vel(i,j,k+1)*(1.0_CGREAL-kup(i,j,k))  &
                       + ghostVel(i,j,k+1)*kup(i,j,k)        ) &
          +  ( 1.0_CGREAL-fz(k+1) ) *vel(i,j,k)   

         velb =  fz(k)*vel(i,j,k)                              &
          +  (1.0_CGREAL-fz(k))*( vel(i,j,k-1)*(1.0_CGREAL-kum(i,j,k))   &
                                  + ghostVel(i,j,k-1)*kum(i,j,k) ) 

!       2nd Upwind differencing added by Rupesh
!       ---------------------------------------
        fdxWeightm1 = ( dzinv(k-1)+2.0_CGREAL*dzinv(k) ) /   &
                              ( dzinv(k-1)+dzinv(k) )
        fdxWeightm2 = dzinv(k) / ( dzinv(k-1)+dzinv(k) )
        fdxWeightp1 = ( dzinv(kpp)+2.0_CGREAL*dzinv(k+1) ) / &
                              ( dzinv(kpp)+dzinv(k+1) )
        fdxWeightp2 = dzinv(k+1) / ( dzinv(kpp)+dzinv(k+1) )

        bdxWeightm1 = ( dzinv(kmm)+2.0_CGREAL*dzinv(k-1) ) / &
                              ( dzinv(kmm)+dzinv(k-1) ) 
        bdxWeightm2 = dzinv(k-1) / ( dzinv(kmm)+dzinv(k-1) )
        bdxWeightp1 = ( dzinv(k+1)+2.0_CGREAL*dzinv(k) ) /   &
                              ( dzinv(k+1)+dzinv(k) )
        bdxWeightp2 = dzinv(k) / ( dzinv(k+1)+dzinv(k) )

        velf_Up = ( ( 1.0_CGREAL+WsignP )                                        &
                 *( kumm(i,j,k)*velf                                          &
                            + (1.0_CGREAL-kumm(i,j,k) )*                      &
                  ( fdxWeightm1*vel(i,j,k)   - fdxWeightm2*vel(i,j,k-1) ) )   &
               + ( 1.0_CGREAL-WsignP )                                        &
                 *( kupp(i,j,k)*velf                                          &
                            + (1.0_CGREAL-kupp(i,j,k) )*                      &
                  ( fdxWeightp1*vel(i,j,k+1) - fdxWeightp2*vel(i,j,kpp) ) )   &
               )*0.5_CGREAL*(1.0_CGREAL-kup(i,j,k)) + bczvel(i,j,k)*kup(i,j,k)

        velb_Up = ( ( 1.0_CGREAL+Wsign )                                         &
                 *( kumm(i,j,k)*velb                                          &
                            + (1.0_CGREAL-kumm(i,j,k) )*                      &
                  ( bdxWeightm1*vel(i,j,k-1) - bdxWeightm2*vel(i,j,kmm) ) )   &   
               + ( 1.0_CGREAL-Wsign )                                         &
                 *( kupp(i,j,k)*velb                                          &
                            + (1.0_CGREAL-kupp(i,j,k) )*                      &
                  ( bdxWeightp1*vel(i,j,k)   - bdxWeightp2*vel(i,j,k+1) ) )   &
               )*0.5_CGREAL*(1.0_CGREAL-kum(i,j,k)) + bczvel(i,j,k)*kum(i,j,k) 

!     ---------------------------------------------------------------------------
!      nlvel(i,j,k) = nlvel(i,j,k)                                             &                                          
!                   + ( velf*face_w(i,j,k+1) -velb*face_w(i,j,k) )*dzinv(k)    &
!      Original CDS definition commented by Rupesh
!     ---------------------------------------------------------------------------

        nlvel(i,j,k) = nlvel(i,j,k) +                                           &
                    (   ( (1.0_CGREAL - alfa)*velf + alfa*velf_Up )*         &
                                                    face_wf(i,j,k)          &
                      - ( (1.0_CGREAL - alfa)*velb + alfa*velb_Up )*         &
                                                    face_wb(i,j,k)            &
                    )*dzinv(k)           ! Hybrid definition added by Rupesh     
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k
    END IF ! nDim

  END SUBROUTINE rhs_advec
!------------------------------------------------------------------------------



  SUBROUTINE rhs_loadAB(nlvel,nlvelOld) 
    
    IMPLICIT NONE

!   parameters
!   ----------
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT(INOUT) :: nlvel,nlvelOld

!   loop  variables
!   ---------------
    INTEGER :: i,j,k

!   local variables
!   ---------------
    REAL(KIND=CGREAL) :: tempvel    



!   ------------------------------------------------------------
!    The CN scheme for the advection terms were added by H. Luo
!   ------------------------------------------------------------

    IF (advec_scheme == ADAMS_BASHFORTH2) THEN    
      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc      
        tempvel =   exp_weight(i,j,k)              *   nlvel(i,j,k) &
                 -( exp_weight(i,j,k)-1.0_CGREAL ) *nlvelOld(i,j,k)

        nlvelOld(i,j,k) = nlvel(i,j,k)
        nlvel(i,j,k)    = tempvel
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k    

    ELSEIF (advec_scheme == CRANK_NICOLSON1 .OR. &
            advec_scheme == CRANK_NICOLSON2) THEN
      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc      
        nlvel(i,j,k) = nlvel(i,j,k) * 0.5_CGREAL
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k  
    ENDIF

  END SUBROUTINE rhs_loadAB
!------------------------------------------------------------------------------



  SUBROUTINE rhs_diff(vel,nlvel,bcxvel,bcyvel,bczvel,ghostVel) 
    
    IMPLICIT NONE

!   parameters
!   ----------   
    REAL(KIND=CGREAL), INTENT(IN)    ::      vel(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL), INTENT(INOUT) :: ghostVel(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT(IN)    :: bcxvel, bcyvel, bczvel
    REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc), INTENT(INOUT) :: nlvel

!   loop  variables
!   ---------------
    INTEGER :: i,j,k

!   local variables
!   ---------------
    INTEGER :: iG, jG

    REAL(KIND=CGREAL) :: acx,amx,apx,acy,amy,apy,acz,amz,apz
    REAL(KIND=CGREAL) :: diffvel
    REAL(KIND=CGREAL) :: nuE,nuW,nuS,nuN,nuF,nuB

!******************************************************************************

!------------------------------------------------------------------------------
! Diffusive terms
! nlvel = d/dx_j [(1/Re+nut) d(vel)/dx_j ]
!
! Note: Since amx_ad coefficients have already a negative sign in them (-1/2 dt)
!       diffvel needs only to be subtracted in computing nlvel term
!------------------------------------------------------------------------------

    DO k = 1,nzc
    DO j = 1,nyc
    
    jG=L2GJ(j)

    DO i = 1,nxc
    
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

      amx = dxcinv(iG)*dxinv(iG) 
      apx = dxcinv(iG+1)*dxinv(iG) 
      amx = amx*nuW
      apx = apx*nuE
      acx = - ( amx + apx )

      ghostVel(i+1,j,k) = gcmFlag*ghostVel(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
              *( bcxvel(i,j,k) - (1.0_CGREAL - fx(iG+1))*vel(i,j,k) )/fx(iG+1)

      ghostVel(i-1,j,k) = gcmFlag*ghostVel(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
              *( bcxvel(i,j,k) - fx(iG)*vel(i,j,k) )/(1.0_CGREAL-fx(iG))

      diffvel  = amx*(      vel(i-1,j,k)*(1.0_CGREAL-ium(i,j,k)) &
                     + ghostVel(i-1,j,k)*            ium(i,j,k)) &
               + acx*       vel(i,j,k)                           &
               + apx*(      vel(i+1,j,k)*(1.0_CGREAL-iup(i,j,k)) &
                     + ghostVel(i+1,j,k)*            iup(i,j,k))  

      amy =  dycinv(jG)*dyinv(jG) 
      apy =  dycinv(jG+1)*dyinv(jG) 
      amy = amy * nuS
      apy = apy * nuN
      acy = - ( amy + apy )

      ghostVel(i,j+1,k) = gcmFlag*ghostVel(i,j+1,k)  + (1.0_CGREAL-gcmFlag) &
              *( bcyvel(i,j,k) - (1.0_CGREAL - fy(jG+1))*vel(i,j,k) )/fy(jG+1)
     
      ghostVel(i,j-1,k) = gcmFlag*ghostVel(i,j-1,k)  + (1.0_CGREAL-gcmFlag) & 
              *( bcyvel(i,j,k) - fy(jG)*vel(i,j,k) )/(1.0_CGREAL-fy(jG))         

      diffvel  = diffvel                                         &
               + amy*(      vel(i,j-1,k)*(1.0_CGREAL-jum(i,j,k)) &
                     + ghostVel(i,j-1,k)*            jum(i,j,k)) &
               + acy*       vel(i,j,k)                           &
               + apy*(      vel(i,j+1,k)*(1.0_CGREAL-jup(i,j,k)) &
                     + ghostVel(i,j+1,k)*            jup(i,j,k))  

      amz =  dzcinv(k)*dzinv(k)*KillFor2D
      apz =  dzcinv(k+1)*dzinv(k)*KillFor2D
      amz = amz * nuB
      apz = apz * nuF
      acz = - ( amz + apz )

      ghostVel(i,j,k+1) = gcmFlag*ghostVel(i,j,k+1)  + (1.0_CGREAL-gcmFlag) & 
              *( bczvel(i,j,k) - (1.0_CGREAL - fz(k+1))*vel(i,j,k) )/fz(k+1)

      ghostVel(i,j,k-1) = gcmFlag*ghostVel(i,j,k-1)  + (1.0_CGREAL-gcmFlag) & 
              *( bczvel(i,j,k) - fz(k)*vel(i,j,k) )/(1.0_CGREAL-fz(k))        
 
      diffvel  = diffvel                                         &
               + amz*(      vel(i,j,k-1)*(1.0_CGREAL-kum(i,j,k)) &
                     + ghostVel(i,j,k-1)*            kum(i,j,k)) &
               + acz*         vel(i,j,k)                         &
               + apz*(      vel(i,j,k+1)*(1.0_CGREAL-kup(i,j,k)) &
                     + ghostVel(i,j,k+1)*            kup(i,j,k))

      nlvel(i,j,k) = vel(i,j,k) - dt*nlvel(i,j,k) + 0.5_CGREAL*dt*diffvel 
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

  END SUBROUTINE rhs_diff
!------------------------------------------------------------------------------



  SUBROUTINE rhs_vankan
    
    IMPLICIT NONE

!   loop variables
!   --------------
    INTEGER :: i,j,k

!   local variables
!   ---------------
    INTEGER :: iG, jG

    REAL(KIND=CGREAL) :: pe,pw,pn,ps,pf,pb,pgx,pgy,pgz

!******************************************************************************

!------------------------------------------------------------------------------
! add pressure gradient
!------------------------------------------------------------------------------

    DO k = 1,nzc
    DO j = 1,nyc
    jG=L2GJ(j)

    DO i = 1,nxc
      iG=L2GI(i)

      pGhost(i+1,j,k) = gcmFlag*pGhost(i+1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
      pGhost(i-1,j,k) = gcmFlag*pGhost(i-1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

      pGhost(i,j+1,k) = gcmFlag*pGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
      pGhost(i,j-1,k) = gcmFlag*pGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

      pGhost(i,j,k+1) = gcmFlag*pGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
      pGhost(i,j,k-1) = gcmFlag*pGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

      pe =              fx(iG+1) *(      p(i+1,j,k)*(1.0_CGREAL-iup(i,j,k)) &
                                  + pGhost(i+1,j,k)*            iup(i,j,k)) &
          + (1.0_CGREAL-fx(iG+1))*       p(i,j,k)  

      pw =  fx(iG)*p(i,j,k)   &
              + (1.0_CGREAL-fx(iG))*( p(i-1,j,k)*(1.0_CGREAL-ium(i,j,k)) &
                                    + pGhost(i-1,j,k)*ium(i,j,k) )

      pn =   fy(jG+1)*( p(i,j+1,k)*(1.0_CGREAL-jup(i,j,k)) &
                      + jup(i,j,k)*pGhost(i,j+1,k)      ) &
              + (1.0_CGREAL-fy(jG+1))*p(i,j,k)  

      ps =  fy(jG)*p(i,j,k)   &
              + (1.0_CGREAL-fy(jG))*( p(i,j-1,k)*(1.0_CGREAL-jum(i,j,k)) &
                                    + pGhost(i,j-1,k)*jum(i,j,k) )

      pf =   fz(k+1)*( p(i,j,k+1)*(1.0_CGREAL-kup(i,j,k)) &
                      + kup(i,j,k)*pGhost(i,j,k+1)      ) &
              + (1.0_CGREAL-fz(k+1))*p(i,j,k)  

      pb =  fz(k)*p(i,j,k)   &
              + (1.0_CGREAL-fz(k))*( p(i,j,k-1)*(1.0_CGREAL-kum(i,j,k)) &
                                    + pGhost(i,j,k-1)*kum(i,j,k) )

      pgx= (pe-pw)*dxinv(iG)
      pgy= (pn-ps)*dyinv(jG)
      pgz= (pf-pb)*dzinv(k)

      nlu(i,j,k) = nlu(i,j,k) - dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)
      nlv(i,j,k) = nlv(i,j,k) - dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)
      nlw(i,j,k) = nlw(i,j,k) - dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)

    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE  rhs_vankan
!------------------------------------------------------------------------------

END SUBROUTINE  rhs_advec_diff 
!-------------------------------------------------------------------------------



SUBROUTINE rhs_adjust_fresh()

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE boundary_arrays
  USE grid_arrays
  USE multiuse_arrays
  USE GCM_arrays

  IMPLICIT NONE

  INTEGER :: i,j,k,iG,jG,n,iFr,jFr,kFr,iRow

  REAL(KIND=CGREAL) :: amxd,apxd,acxd
  REAL(KIND=CGREAL) :: amyd,apyd,acyd
  REAL(KIND=CGREAL) :: amzd,apzd,aczd
  REAL(KIND=CGREAL) :: rFreshCell



  IF ( boundary_motion      == MOVING_BOUNDARY .AND.  &
       boundary_formulation == SSM_METHOD   ) THEN

    DO k=1,nzc
    DO j=1,nyc
      jG=L2GJ(j)
      DO i=1,nxc

        iG=L2GI(i)

!       Take care of fresh cells in SSM, only applies to fresh cells
!       ------------------------------------------------------------

        rFreshCell = REAL(fresh_cell(i,j,k),KIND=CGREAL)

        IF ( fresh_cell(i,j,k) == 1) THEN

          amxd = - (2.0_CGREAL*dxcinv(iG)  )**sidw 
          apxd = - (2.0_CGREAL*dxcinv(iG+1))**sidw   

          amyd = - (2.0_CGREAL*dycinv(jG)  )**sidw 
          apyd = - (2.0_CGREAL*dycinv(jG+1))**sidw 

          amzd = -((2.0_CGREAL*dzcinv(k )  )**sidw)*KillFor2D
          apzd = -((2.0_CGREAL*dzcinv(k+1 ))**sidw)*KillFor2D

          nlu(i,j,k) = - amxd*bcxu(i,j,k)*ium(i,j,k) &
                       - apxd*bcxu(i,j,k)*iup(i,j,k) &
                       - amyd*bcyu(i,j,k)*jum(i,j,k) &
                       - apyd*bcyu(i,j,k)*jup(i,j,k) &
                       - amzd*bczu(i,j,k)*kum(i,j,k) &
                       - apzd*bczu(i,j,k)*kup(i,j,k)
 
          nlv(i,j,k) = - amxd*bcxv(i,j,k)*ium(i,j,k) &
                       - apxd*bcxv(i,j,k)*iup(i,j,k) &
                       - amyd*bcyv(i,j,k)*jum(i,j,k) &
                       - apyd*bcyv(i,j,k)*jup(i,j,k) &
                       - amzd*bczv(i,j,k)*kum(i,j,k) &
                       - apzd*bczv(i,j,k)*kup(i,j,k)  
 
          nlw(i,j,k) = - amxd*bcxw(i,j,k)*ium(i,j,k) &
                       - apxd*bcxw(i,j,k)*iup(i,j,k) &
                       - amyd*bcyw(i,j,k)*jum(i,j,k) &
                       - apyd*bcyw(i,j,k)*jup(i,j,k) &
                       - amzd*bczw(i,j,k)*kum(i,j,k) &
                       - apzd*bczw(i,j,k)*kup(i,j,k)

        ENDIF ! fresh_cell
      ENDDO
    ENDDO
    ENDDO

  ENDIF 

! Take care of fresh cells in GCM
! -------------------------------
  IF ( boundary_motion      == MOVING_BOUNDARY .AND.  &
       boundary_formulation == GCM_METHOD .AND. nFresh > 0   ) THEN

    DO n=1,nFresh
      iFr = iFresh(n)
      jFr = jFresh(n)
      kFr = kFresh(n)
      DO iRow = 1,iRowMax
        i  = iFreshCellIndex(n) + incI(iRow)
        j  = jFreshCellIndex(n) + incJ(iRow)
        k  = kFreshCellIndex(n) + incK(iRow)
        IF ( i==iFr .AND. j==jFr .AND. k==kFr) THEN
          nlu(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*uBodyInterceptFresh(n)
          nlv(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*vBodyInterceptFresh(n)
          nlw(iFr,jFr,kFr) = coeffGCMFreshD(iRow,n)*wBodyInterceptFresh(n)
        ENDIF
      ENDDO
    ENDDO
  ENDIF


END SUBROUTINE rhs_adjust_fresh
!------------------------------------------------------------------------------- 



SUBROUTINE rhs_adjust2D() 

  USE global_parameters
  USE flow_parameters
  USE multiuse_arrays

  IMPLICIT NONE

  INTEGER :: i,j,k



! -----------------------------------------------------
!  Sets RHS for momentum equations for 2D calculations
! -----------------------------------------------------

! copy k=1 plane to other planes
! ------------------------------
  DO k = 2,nzc
  DO j = 1,nyc
  DO i = 1,nxc
    nlu(i,j,k) = nlu(i,j,1)
    nlv(i,j,k) = nlv(i,j,1)
  ENDDO ! i
  ENDDO ! j
  ENDDO ! k
    
! zero w-component
! ----------------
  nlw(:,:,:) = 0.0_CGREAL

END SUBROUTINE  rhs_adjust2D
!---------------------------------------------------------------------
