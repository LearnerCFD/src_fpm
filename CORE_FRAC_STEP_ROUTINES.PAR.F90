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
!  Filename: CORE_FRAC_STEP_ROUTINES.PAR.F90
!  Latest Modification: Oct, 20 2010 (ver. PAT 1.2.0)
!  by JHSeo
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     face_vel()
!     correct_vel()
!     update_pressure()
!     set_face_vel_body()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE face_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE pressure_arrays
    USE solver_arrays

    IMPLICIT NONE

    INTEGER :: i, j, k
    INTEGER :: iG, jG

    REAL(KIND=CGREAL) :: pe, pw, pn, ps, pf, pb
    REAL(KIND=CGREAL) :: pgx, pgy, pgz
    REAL(KIND=CGREAL) :: pgxe, pgxw, pgyn, pgys, pgzf, pgzb
    REAL(KIND=CGREAL) :: myTime
    
!DEBUG    
!    CHARACTER*20              :: fname1
!DEBUG    



!DEBUG    
!    WRITE(fname1,"('q.',I3.3,'.',I7.7)") myRank,ntime
!    OPEN(UNIT=70,FILE=fname1,FORM='UNFORMATTED')
!    WRITE(70) nxc,myIs,myIe,nyc,MyJs,myJe,nzc
!DEBUG    

!   compute face velocities
!   -----------------------
!    face_ue = 0.0_CGREAL
!    face_uw = 0.0_CGREAL
!    face_vn = 0.0_CGREAL
!    face_vs = 0.0_CGREAL
!    face_wf = 0.0_CGREAL
!    face_wb = 0.0_CGREAL

    IF (frac_step_type == VAN_KAN) THEN

      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc
      
        iG=L2GI(i)
        jG=L2GJ(j)

        pGhost(i+1,j,k) = gcmFlag*pGhost(i+1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i-1,j,k) = gcmFlag*pGhost(i-1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pGhost(i,j+1,k) = gcmFlag*pGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i,j-1,k) = gcmFlag*pGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pGhost(i,j,k+1) = gcmFlag*pGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i,j,k-1) = gcmFlag*pGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pe = fx(iG+1)*(p(i+1,j,k)*(1.0_CGREAL-iup(i,j,k)) + iup(i,j,k)*pGhost(i+1,j,k)) &
            +(1.0_CGREAL-fx(iG+1))*p(i,j,k)

        pw = fx(iG)*p(i,j,k)   &
            +(1.0_CGREAL-fx(iG))*(p(i-1,j,k)*(1.0_CGREAL-ium(i,j,k)) + pGhost(i-1,j,k)*ium(i,j,k))

        pn = fy(jG+1)*(p(i,j+1,k)*(1.0_CGREAL-jup(i,j,k)) + jup(i,j,k)*pGhost(i,j+1,k)) &
            +(1.0_CGREAL-fy(jG+1))*p(i,j,k)

        ps = fy(jG)*p(i,j,k)   &
            +(1.0_CGREAL-fy(jG))*(p(i,j-1,k)*(1.0_CGREAL-jum(i,j,k)) + pGhost(i,j-1,k)*jum(i,j,k))

        pf = fz(k+1)*(p(i,j,k+1)*(1.0_CGREAL-kup(i,j,k)) + kup(i,j,k)*pGhost(i,j,k+1)) &
            +(1.0_CGREAL-fz(k+1))*p(i,j,k)

        pb = fz(k)*p(i,j,k)   &
            +(1.0_CGREAL-fz(k))*(p(i,j,k-1)*(1.0_CGREAL-kum(i,j,k)) + pGhost(i,j,k-1)*kum(i,j,k))

        pgx= (pe-pw)*dxinv(iG)
        pgy= (pn-ps)*dyinv(jG)
        pgz= (pf-pb)*dzinv(k)

        uTilde(i,j,k) = u(i,j,k) + dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)
        vTilde(i,j,k) = v(i,j,k) + dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)
        wTilde(i,j,k) = w(i,j,k) + dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)

      ENDDO
      ENDDO
      ENDDO

!     Updating Parallel Ghost Layer
!     -----------------------------
#     ifdef MPI
        CALL par_comm_var(uTilde,nxc,nyc,nzc,Ngl,myTime)
        CALL par_comm_var(vTilde,nxc,nyc,nzc,Ngl,myTime)
        CALL par_comm_var(wTilde,nxc,nyc,nzc,Ngl,myTime)
#     endif

    ELSE

	  DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc

       uTilde(i,j,k) = u(i,j,k)*(1.-iblank_solid(i,j,k)) + uGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
       vTilde(i,j,k) = v(i,j,k)*(1.-iblank_solid(i,j,k)) + vGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
       wTilde(i,j,k) = w(i,j,k)*(1.-iblank_solid(i,j,k)) + wGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
	  
	  ENDDO
	  ENDDO
	  ENDDO
	  
#     ifdef MPI
        CALL par_comm_var(uTilde,nxc,nyc,nzc,Ngl,myTime)
        CALL par_comm_var(vTilde,nxc,nyc,nzc,Ngl,myTime)
        CALL par_comm_var(wTilde,nxc,nyc,nzc,Ngl,myTime)
#     endif	  

!     ------------------------------------------
!     Parallel Ghost Layer is already up-to-date
!     ------------------------------------------

    ENDIF! frac_step_type

    DO k = 1,nzc
    DO j = 1,nyc
      jG=L2GJ(j)

      DO i = 1,nxc
        iG=L2GI(i)

        uGhost(i+1,j,k) = gcmFlag*uGhost(i+1,j,k)  + (1.0_CGREAL-gcmFlag) &
                *( bcxu(i,j,k) - (1.0_CGREAL - fx(iG+1))*uTilde(i,j,k) )/fx(iG+1)

        uGhost(i-1,j,k) = gcmFlag*uGhost(i-1,j,k)  + (1.0_CGREAL-gcmFlag) &
                *( bcxu(i,j,k) - fx(iG)*uTilde(i,j,k) )/(1.0_CGREAL-fx(iG))

        vGhost(i,j+1,k) = gcmFlag*vGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag) & 
                *( bcyv(i,j,k) - (1.0_CGREAL - fy(jG+1))*vTilde(i,j,k) )/fy(jG+1)  
                 
        vGhost(i,j-1,k) = gcmFlag*vGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag) &
                *( bcyv(i,j,k) - fy(jG)*vTilde(i,j,k) )/(1.0_CGREAL-fy(jG))        

        wGhost(i,j,k+1) = gcmFlag*wGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag) &
                *( bczw(i,j,k) - (1.0_CGREAL - fz(k+1))*wTilde(i,j,k) )/fz(k+1)
       
        wGhost(i,j,k-1) = gcmFlag*wGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag) & 
                *( bczw(i,j,k) - fz(k)*wTilde(i,j,k) )/(1.0_CGREAL-fz(k))        

!!        IF(iblank_solid(i, j, k) == 0) THEN                     
         
         face_ue(i,j,k) = fx(iG+1) *(uTilde(i+1,j,k)*(1.0_CGREAL-iup(i,j,k))+uGhost(i+1,j,k)*iup(i,j,k)) &
                        + (1.0_CGREAL-fx(iG+1))* uTilde(i  ,j,k)

         face_uw(i,j,k) = fx(iG) * uTilde(i  ,j,k)                                                       &
                       + (1.0_CGREAL-fx(iG))*(uTilde(i-1,j,k)*(1.0_CGREAL-ium(i,j,k))+uGhost(i-1,j,k)*ium(i,j,k))

         face_vn(i,j,k) = fy(jG+1) *(vTilde(i,j+1,k)*(1.0_CGREAL-jup(i,j,k))+vGhost(i,j+1,k)*jup(i,j,k)) &
                       + (1.0_CGREAL-fy(jG+1))* vTilde(i,j  ,k) 

         face_vs(i,j,k) = fy(jG) * vTilde(i,j  ,k)                                                       &
                       + (1.0_CGREAL-fy(jG))*(vTilde(i,j-1,k)*(1.0_CGREAL-jum(i,j,k))+vGhost(i,j-1,k)*jum(i,j,k))

         face_wf(i,j,k) = fz(k+1) *(wTilde(i,j,k+1)*(1.0_CGREAL-kup(i,j,k))+wGhost(i,j,k+1)*kup(i,j,k))  &
                       + (1.0_CGREAL-fz(k+1))* wTilde(i,j,k)

         face_wb(i,j,k) = fz(k) * wTilde(i,j,k)                                                          &
                       + (1.0_CGREAL-fz(k))*(wTilde(i,j,k-1)*(1.0_CGREAL-kum(i,j,k))+wGhost(i,j,k-1)*kum(i,j,k))
 !!      ENDIF !iblank_solid

!DEBUG
!        WRITE(70) u(i,j,k), v(i,j,k), w(i,j,k), face_ue(i,j,k), face_uw(i,j,k), face_vn(i,j,k), face_vs(i,j,k), face_wf(i,j,k), face_wb(i,j,k), pprime(i,j,k)
!DEBUG
      ENDDO
    ENDDO
    ENDDO
!DEBUG
! close(70)
! call flow_stop
! stop
!DEBUG

    IF (frac_step_type == VAN_KAN) THEN

!     add pressure gradient to face velocities
!     ----------------------------------------
      DO k = 1,nzc
      DO j = 1,nyc
      DO i = 1,nxc
        iG=L2GI(i)

        pGhost(i+1,j,k) = gcmFlag*pGhost(i+1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i-1,j,k) = gcmFlag*pGhost(i-1,j,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pgxw       = ( p(i,j,k)  - ( p(i-1,j,k)*(1.0_CGREAL-ium(i,j,k))  &
                                   + pGhost(i-1,j,k)*ium(i,j,k) )      )*dxcinv(iG)
        pgxe       = ( ( p(i+1,j,k)*(1.0_CGREAL-iup(i,j,k)) &
                       + pGhost(i+1,j,k)*iup(i,j,k) )   - p(i,j,k)  )*dxcinv(iG+1)

        face1(i)   = face_uw(i,j,k) - dt*pgxw
        face2(i)   = face_ue(i,j,k) - dt*pgxe
      ENDDO
      DO i = 1,nxc
        face_uw(i,j,k)   = face1(i)
        face_ue(i,j,k)   = face2(i)
      ENDDO
      ENDDO
      ENDDO

      DO k = 1,nzc
      DO i = 1,nxc
      DO j = 1,nyc
        jG=L2GJ(j)
          
        pGhost(i,j+1,k) = gcmFlag*pGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i,j-1,k) = gcmFlag*pGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pgys      = ( p(i,j,k)  - ( p(i,j-1,k)*(1.0_CGREAL-jum(i,j,k)) &
                                  + pGhost(i,j-1,k)*jum(i,j,k) )  )*dycinv(jG)
        pgyn      = ( ( p(i,j+1,k)*(1.0_CGREAL-jup(i,j,k)) &
                      + pGhost(i,j+1,k)*jup(i,j,k) ) - p(i,j,k)  )*dycinv(jG+1)

        face1(j)   = face_vs(i,j,k) - dt*pgys
        face2(j)   = face_vn(i,j,k) - dt*pgyn
      ENDDO
      DO j = 1,nyc
        face_vs(i,j,k) = face1(j)
        face_vn(i,j,k) = face2(j)
      ENDDO
      ENDDO
      ENDDO

      DO j = 1,nyc
      DO i = 1,nxc
      DO k = 1,nzc
        pGhost(i,j,k+1) = gcmFlag*pGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)
        pGhost(i,j,k-1) = gcmFlag*pGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag)*p(i,j,k)

        pgzb     = ( p(i,j,k)  - ( p(i,j,k-1)*(1.0_CGREAL-kum(i,j,k)) &
                                 + pGhost(i,j,k-1)*kum(i,j,k) ) )*dzcinv(k)
        pgzf     = ( ( p(i,j,k+1)*(1.0_CGREAL-kup(i,j,k))  &
                     + pGhost(i,j,k+1)*kup(i,j,k) ) - p(i,j,k)  )*dzcinv(k+1)

        face1(k) = face_wb(i,j,k) - dt*pgzb
        face2(k) = face_wf(i,j,k) - dt*pgzf
      ENDDO
      DO k = 1,nzc
        face_wb(i,j,k) = face1(k)
        face_wf(i,j,k) = face2(k)
      ENDDO
      ENDDO
      ENDDO

    ENDIF ! frac_step_type

!   Zeroing out face velocities for solid cells
!   Also for "dead" faces of ghost cells
!   -------------------------------------------
    CALL set_face_vel_body()

END SUBROUTINE face_vel
!-------------------------------------------------------------------------------



SUBROUTINE correct_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE grid_arrays
    USE solver_arrays
	  USE Cutcell_arrays
	  USE unstructured_surface_arrays	

    IMPLICIT NONE
    
    INTEGER :: i, j, k
    INTEGER :: iG, jG, cElem, iBody

    REAL(KIND=CGREAL) :: pe, pw, pn, ps, pf, pb, pc
    REAL(KIND=CGREAL) :: pgx, pgy, pgz
    REAL(KIND=CGREAL) :: pgxe, pgxw, pgyn, pgys, pgzf, pgzb
    REAL(KIND=CGREAL) :: gcmflag2, myTime
 	  REAL(KIND=CGREAL) :: pxn, pyn, pzn, dpdn, dist, killPBC
	  REAL(KIND=CGREAL) :: xBI,yBI,zBI,xBIN,yBIN,zBIN,xp,yp,zp, killw, kille, killn, kills, killf, killb
!DEBUG
!  CHARACTER*20 :: fname1
!DEBUG

    killPBC = 1.0
    
    IF ( pbcx1 == PBC_DIRICHLET .OR. &                                 ! SER_TO_PAR. QX. CH23
         pbcx2 == PBC_DIRICHLET .OR. &
         pbcy1 == PBC_DIRICHLET .OR. &
         pbcy2 == PBC_DIRICHLET .OR. &
         pbcz1 == PBC_DIRICHLET .OR. &
         pbcz2 == PBC_DIRICHLET       )  killPBC = 0.0    


!DEBUG
!  WRITE(fname1,"('q.',I3.3,'.',I7.7)") myRank,ntime
!  OPEN(UNIT=70,FILE=fname1,FORM='UNFORMATTED')
!  WRITE(70) nxc,myIs,myIe,nyc,MyJs,myJe,nzc
!DEBUG


!!    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &
!!        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &
!!        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
!!       CALL apply_periodic_pres(pPrime)
!!    END IF

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &             !Ehsan added for Periodic BC
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. &             
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN             
          do i=0,nxc+1                         
            do j=0,nyc+1                             
             pPrime (i,j,0) = pPrime(i,j,nzc)       
              pPrime(i,j,nzc+1) = pPrime(i,j,1)     
            enddo                                
          enddo                             
    END IF                               

      IF (pbcx1 == PBC_DIRICHLET .OR. &   !
        pbcy1 == PBC_DIRICHLET .OR. &   ! Added by H. Luo
        pbcz1 == PBC_DIRICHLET .OR. &   !
        pbcx2 == PBC_DIRICHLET .OR. &   ! Modified by SAMK
        pbcy2 == PBC_DIRICHLET .OR. &   !
        pbcz2 == PBC_DIRICHLET) THEN    !

      CALL set_pressure_dirichlet_bc()!..........................COMPLETE(SAMK)
    END IF


!   correct nodal velocities
!   ------------------------
    DO k = 1,nzc
    DO j = 1,nyc
      jG = L2GJ(j)

      DO i = 1,nxc
        iG = L2GI(i)
		
		kille = 1.0-ghostcellmark(i+1,j,k)*iSSMP*iCC
		killw = 1.0-ghostcellmark(i-1,j,k)*iSSMP*iCC
		killn = 1.0-ghostcellmark(i,j+1,k)*iSSMP*iCC
		kills = 1.0-ghostcellmark(i,j-1,k)*iSSMP*iCC
		killf = 1.0-ghostcellmark(i,j,k+1)*iSSMP*iCC
		killb = 1.0-ghostcellmark(i,j,k-1)*iSSMP*iCC
		

        pGhost(i+1,j,k) = gcmFlag*kille*pGhost(i+1,j,k)  + (1.0_CGREAL-gcmFlag*kille)*pPrime(i,j,k)
        pGhost(i-1,j,k) = gcmFlag*killw*pGhost(i-1,j,k)  + (1.0_CGREAL-gcmFlag*killw)*pPrime(i,j,k)

        pGhost(i,j+1,k) = gcmFlag*killn*pGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag*killn)*pPrime(i,j,k)
        pGhost(i,j-1,k) = gcmFlag*kills*pGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag*kills)*pPrime(i,j,k)

        pGhost(i,j,k+1) = gcmFlag*killf*pGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag*killf)*pPrime(i,j,k)
        pGhost(i,j,k-1) = gcmFlag*killb*pGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag*killb)*pPrime(i,j,k)

        pe =  fx(iG+1)*( pPrime(i+1,j,k)*(1.0_CGREAL-iup(i,j,k)) &
                        + iup(i,j,k)*pGhost(i+1,j,k)      ) &
                + (1.0_CGREAL-fx(iG+1))*pPrime(i,j,k)
!DEBUG
!        WRITE(70) fx(iG+1), pPrime(i+1,j,k), iup(i,j,k), pGhost(i+1,j,k), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
!DEBUG

        pw =  fx(iG)*pPrime(i,j,k)   &
                + (1.0_CGREAL-fx(iG))*( pPrime(i-1,j,k)*(1.0_CGREAL-ium(i,j,k)) &
                                      + pGhost(i-1,j,k)*ium(i,j,k) )

        pn =  fy(jG+1)*( pPrime(i,j+1,k)*(1.0_CGREAL-jup(i,j,k)) &
                        + jup(i,j,k)*pGhost(i,j+1,k)      ) &
                + (1.0_CGREAL-fy(jG+1))*pPrime(i,j,k)

        ps =  fy(jG)*pPrime(i,j,k)   &
                + (1.0_CGREAL-fy(jG))*( pPrime(i,j-1,k)*(1.0_CGREAL-jum(i,j,k)) &
                                      + pGhost(i,j-1,k)*jum(i,j,k) )

        pf =  fz(k+1)*( pPrime(i,j,k+1)*(1.0_CGREAL-kup(i,j,k)) &
                        + kup(i,j,k)*pGhost(i,j,k+1)      ) &
                + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)

        pb =  fz(k)*pPrime(i,j,k)   &
                + (1.0_CGREAL-fz(k))*( pPrime(i,j,k-1)*(1.0_CGREAL-kum(i,j,k)) &
                                      + pGhost(i,j,k-1)*kum(i,j,k) )
	
	    IF(ghostcellsolid(i,j,k) == 1 .and. iCC.eq.1) THEN
		 xp = xc(iG) !cent_x(i,j,k) !xc(iG)
		 yp = yc(jG) !cent_y(i,j,k) !yc(jG)
		 zp = zc(k)  !cent_z(i,j,k) !zc(k)
		 
		 iBody = BodyNum(i,j,k)
		 IF(iBody == 0) THEN
          iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
          IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
         ENDIF
		 
		 CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI, yBI, zBI, cElem ) 

		 xBIN = triElemNormx(iBody,cElem)
	     yBIN = triElemNormy(iBody,cElem)
         zBIN = triElemNormz(iBody,cElem)  
		 
		 !dist = sqrt( (xp-xBI)**2 + (yp-yBI)**2 + (zp-zBI)**2 )
		 !CALL cutcell_dpdn(i,j,k,dpdn)
		 
		 pxn = pPrime(i,j,k)*xBIN !(pPrime(i,j,k)+dpdn*dist)*xBIN
		 pyn = pPrime(i,j,k)*yBIN !(pPrime(i,j,k)+dpdn*dist)*yBIN
		 pzn = pPrime(i,j,k)*zBIN !(pPrime(i,j,k)+dpdn*dist)*zBIN
		 
		ELSE
		 
		 pxn = 0.0
		 pyn = 0.0
		 pzn = 0.0
		 
		ENDIF ! ghostcellsolid .and. iCC
		
		pgx= (pe*AREA_E(i,j,k)-pw*AREA_W(i,j,k))*dxinv(iG) + pxn*Cut_Area(i,j,k)*dxinv(iG)*dyinv(jG)*dzinv(k)
        pgy= (pn*AREA_N(i,j,k)-ps*AREA_S(i,j,k))*dyinv(jG) + pyn*Cut_Area(i,j,k)*dxinv(iG)*dyinv(jG)*dzinv(k)
        pgz= (pf*AREA_F(i,j,k)-pb*AREA_B(i,j,k))*dzinv(k ) + pzn*Cut_Area(i,j,k)*dxinv(iG)*dyinv(jG)*dzinv(k)
!!DEBUG
!        WRITE(70) pe, pw, pn, ps, pf, pb, pgx, pgy, pgz, 0.0d0
!!DEBUG

        if(iblank_solid(i,j,k) .eq. 0) THEN
        
		u(i,j,k) = u(i,j,k) - dt*pgx*REAL(1-iblank(i,j,k),KIND=CGREAL)/VOL_CELL(i,j,k)
        v(i,j,k) = v(i,j,k) - dt*pgy*REAL(1-iblank(i,j,k),KIND=CGREAL)/VOL_CELL(i,j,k)
        w(i,j,k) = w(i,j,k) - dt*pgz*REAL(1-iblank(i,j,k),KIND=CGREAL)/VOL_CELL(i,j,k)
        
		endif

      ENDDO
    ENDDO
    ENDDO
	
!DEBUG
!  CLOSE(70)
!  CALL flow_stop
!  STOP
!DEBUG

    IF (bcx1 .EQ. BC_TYPE_PERIODIC .OR. &                       !Ehsan added for Periodic BC 
        bcy1 .EQ. BC_TYPE_PERIODIC .OR. & 
        bcz1 .EQ. BC_TYPE_PERIODIC) THEN
      CALL apply_periodic_bcvel
    END IF

!   update parallel ghost layer
!   ---------------------------
#   ifdef MPI
      CALL par_comm_var(u,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(v,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(w,nxc,nyc,nzc,Ngl,myTime)
#   endif

!   update value of velocity at ghost points through interpolation
!   --------------------------------------------------------------
    IF (boundary_formulation == GCM_METHOD) CALL GCM_vel_set_bc_internal()

!   correct face velocities
!   -----------------------
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      iG=L2GI(i)

      gcmFlag2 = gcmFlag*(1.0_CGREAL - iBound(i)*killPBC)
      pGhost(i+1,j,k) = gcmFlag2*pGhost(i+1,j,k) + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)
      pGhost(i-1,j,k) = gcmFlag2*pGhost(i-1,j,k) + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)

	  pc = pPrime(i,j,k)*(1.-iblank_solid(i,j,k))     + pGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
	  pw = pPrime(i-1,j,k)*(1.-iblank_solid(i-1,j,k)) + pGhost(i-1,j,k)*ghostcellmark(i-1,j,k)*iblank_solid(i-1,j,k)
	  pe = pPrime(i+1,j,k)*(1.-iblank_solid(i+1,j,k)) + pGhost(i+1,j,k)*ghostcellmark(i+1,j,k)*iblank_solid(i+1,j,k)

      pgxw       = ( pc  - ( pw*(1.0_CGREAL-ium(i,j,k))  &
                                + pGhost(i-1,j,k)*ium(i,j,k) )      )*dxcinv(iG)
      pgxe       = ( ( pe*(1.0_CGREAL-iup(i,j,k)) &
                     + pGhost(i+1,j,k)*iup(i,j,k) )   - pc  )*dxcinv(iG+1)
      face1(i)   = face_uw(i,j,k) - dt*pgxw*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)
      face2(i)   = face_ue(i,j,k) - dt*pgxe*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)
    ENDDO
    
    DO i = 1,nxc
      face_uw(i,j,k)  = face_uw(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)      &
                           +face1(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      face_ue(i,j,k)  = face_ue(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)      &
                           +face2(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)
    ENDDO
    ENDDO
    ENDDO

    DO k = 1,nzc
    DO i = 1,nxc
    DO j = 1,nyc
      jG=L2GJ(j)
      
      gcmFlag2 = gcmFlag*(1.0_CGREAL - jBound(j)*killPBC)
      pGhost(i,j+1,k) = gcmFlag2*pGhost(i,j+1,k)  + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)
      pGhost(i,j-1,k) = gcmFlag2*pGhost(i,j-1,k)  + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)

	  pc = pPrime(i,j,k)*(1.-iblank_solid(i,j,k)) + pGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
	  ps = pPrime(i,j-1,k)*(1.-iblank_solid(i,j-1,k)) + pGhost(i,j-1,k)*ghostcellmark(i,j-1,k)*iblank_solid(i,j-1,k)
	  pn = pPrime(i,j+1,k)*(1.-iblank_solid(i,j+1,k)) + pGhost(i,j+1,k)*ghostcellmark(i,j+1,k)*iblank_solid(i,j+1,k)	  

      pgys      = ( pc  - ( ps*(1.0_CGREAL-jum(i,j,k)) &
                                      + pGhost(i,j-1,k)*jum(i,j,k) )  )*dycinv(jG)
      pgyn      = ( ( pn*(1.0_CGREAL-jup(i,j,k)) &
                            + pGhost(i,j+1,k)*jup(i,j,k) ) - pc  )*dycinv(jG+1)
      face1(j)   = face_vs(i,j,k) - dt*pgys*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)
      face2(j)   = face_vn(i,j,k) - dt*pgyn*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)
    ENDDO

    DO j = 1,nyc
      face_vs(i,j,k) = face_vs(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)      &
                          +face1(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      face_vn(i,j,k) = face_vn(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)      &
                          +face2(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
    ENDDO
    ENDDO
    ENDDO

    DO j = 1,nyc
    DO i = 1,nxc
    DO k = 1,nzc
            gcmFlag2 = gcmFlag*(1.0_CGREAL - kBound(k)*killPBC)
            pGhost(i,j,k+1) = gcmFlag2*pGhost(i,j,k+1)  + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)
            pGhost(i,j,k-1) = gcmFlag2*pGhost(i,j,k-1)  + (1.0_CGREAL-gcmFlag2)*pPrime(i,j,k)

	  pc = pPrime(i,j,k)*(1.-iblank_solid(i,j,k)) + pGhost(i,j,k)*ghostcellmark(i,j,k)*iblank_solid(i,j,k)
	  pb = pPrime(i,j,k-1)*(1.-iblank_solid(i,j,k-1)) + pGhost(i,j,k-1)*ghostcellmark(i,j,k-1)*iblank_solid(i,j,k-1)
	  pf = pPrime(i,j,k+1)*(1.-iblank_solid(i,j,k+1)) + pGhost(i,j,k+1)*ghostcellmark(i,j,k+1)*iblank_solid(i,j,k+1)			

      pgzb     = ( pc  - ( pb*(1.0_CGREAL-kum(i,j,k)) &
                                    + pGhost(i,j,k-1)*kum(i,j,k) ) )*dzcinv(k)
      pgzf     = ( ( pf*(1.0_CGREAL-kup(i,j,k))  &
                          + pGhost(i,j,k+1)*kup(i,j,k) ) - pc )*dzcinv(k+1)
      face1(k) = face_wb(i,j,k) - dt*pgzb*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)
      face2(k) = face_wf(i,j,k) - dt*pgzf*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC)
    ENDDO
    
    DO k = 1,nzc
      face_wb(i,j,k) = face_wb(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)     &
                          +face1(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
      face_wf(i,j,k) = face_wf(i,j,k)*REAL(iblank(i,j,k),KIND=CGREAL)     &
                          +face2(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)
    ENDDO
    ENDDO
    ENDDO

!   set face vel. to zero for solid cells. also for ghost-ghost cell faces
!   ----------------------------------------------------------------------
    CALL set_face_vel_body
  
    IF (pbcx1 == PBC_DIRICHLET .OR. &   !
        pbcy1 == PBC_DIRICHLET .OR. &   ! Added by H. Luo
        pbcz1 == PBC_DIRICHLET .OR. &   !
        pbcx2 == PBC_DIRICHLET .OR. &   ! Modified by SAMK
        pbcy2 == PBC_DIRICHLET .OR. &   !
        pbcz2 == PBC_DIRICHLET) THEN    !
            
      CALL unset_pressure_dirichlet_bc()!........................COMPLETE(SAMK)
    END IF

END SUBROUTINE correct_vel
!---------------------------------------------------------------------



SUBROUTINE update_pressure()

    USE global_parameters
    USE flow_parameters
    USE pressure_arrays

    IMPLICIT NONE



    IF (frac_step_type == VAN_KAN) THEN
      p(:,:,:) = p(:,:,:) + pPrime(:,:,:)
    ELSE
      p(:,:,:) = pPrime(:,:,:)
    ENDIF

!   Update pressure for ghost cells
!   -------------------------------  
    IF (boundary_formulation == GCM_METHOD) CALL GCM_p_set_bc_internal(p)!......COMPLETE(SAMK)
    
!   ---------------------------------------------------------------
!    NOTE(SAMK): There is no need of updating parallel ghost layer
!                as it is already calculated.
!   ---------------------------------------------------------------

END SUBROUTINE update_pressure
!---------------------------------------------------------------------



SUBROUTINE set_face_vel_body() 

! ----------------------------------------------------------------------
!  Zeroing out face velocities for solid cells.
!  Also for "dead" faces of Ghost Cells
! ----------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
	USE cutcell_arrays

    IMPLICIT NONE

!   Loop variables
!   --------------
    INTEGER :: i,j,k

    IF(iCC) THEN
	
	DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      IF ( isolid(i,j,k) == 1 ) THEN
		face_ue(i,j,k) = 0.0_CGREAL
        face_uw(i,j,k) = 0.0_CGREAL
        face_vn(i,j,k) = 0.0_CGREAL
        face_vs(i,j,k) = 0.0_CGREAL
        face_wf(i,j,k) = 0.0_CGREAL
        face_wb(i,j,k) = 0.0_CGREAL
      ENDIF ! iblank
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
	
	ELSE


    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      IF ( iblank_solid(i,j,k) == 1 ) THEN
        IF (iblank_solid(i+1,j,k) == 1 ) face_ue(i,j,k) = 0.0_CGREAL
        IF (iblank_solid(i-1,j,k) == 1 ) face_uw(i,j,k) = 0.0_CGREAL
        IF (iblank_solid(i,j+1,k) == 1 ) face_vn(i,j,k) = 0.0_CGREAL
        IF (iblank_solid(i,j-1,k) == 1 ) face_vs(i,j,k) = 0.0_CGREAL
        IF (iblank_solid(i,j,k+1) == 1 ) face_wf(i,j,k) = 0.0_CGREAL
        IF (iblank_solid(i,j,k-1) == 1 ) face_wb(i,j,k) = 0.0_CGREAL
      ENDIF ! iblank
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
	
	ENDIF ! iCC

  END SUBROUTINE set_face_vel_body    
!----------------------------------------------------------------------
