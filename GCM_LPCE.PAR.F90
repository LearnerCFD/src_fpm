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
!  Filename: GCM_LPCE.PAR.F90
!  Latest Modification: Sep, 11 2010 (ver. P2.0.0)
!  Made by Jung-Hee Seo
! --------------------------------------------------------------------
SUBROUTINE GCM_LPCE_vel_set_bc_norm(utemp,vtemp,wtemp)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE grid_arrays
    USE unstructured_surface_arrays
	USE mpi

    IMPLICIT NONE	

	REAL(KIND=CGREAL) :: utemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: vtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: wtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    
!   Loop variables
!   --------------
    INTEGER :: iBody, iRow, k, n

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iGh, jGh, kGh, iterGC
    INTEGER :: ii,jj,kk
	   
    REAL(KIND=CGREAL)  :: uIP, vIP, wIP,        &
                          resVelMax, resVel,    &
                          uPrev, vPrev, wPrev,  &
                          sameSideGhost,        &
						  dslength, Velnorm,    &
						  xBI,  yBI,  zBI,      &
						  xBIn, yBIn, zBIn,     &
						  xGC,  yGC,  zGC



!   Compute body intercept velocity field
!   -------------------------------------
!   CALL GCM_SetBodyInterceptValues()
   
!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resVelMax = 1.0D10
      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC<itermax) .AND. (resVelMax>restol))

      resVelMax = 0.0_CGREAL

      DO n = 1,nGhost
        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
             CALL flow_stop
              STOP
         ENDIF

!       Initialize values
!       -----------------
        uIP = 0.0_CGREAL
        vIP = 0.0_CGREAL
        wIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
           ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + inck(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) .EQV. .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            uIP = uIP + coeffGCMD(iRow,n)*(    utemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                        +utemp(ii,jj,kk)*            sameSideGhost)

            vIP = vIP + coeffGCMD(iRow,n)*(    vtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +vtemp(ii,jj,kk)*            sameSideGhost)

            wIP = wIP + coeffGCMD(iRow,n)*(    wtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +wtemp(ii,jj,kk)*            sameSideGhost)
          ! ELSE
              !uIP = uIP !+ coeffGCMD(iRow,n)* 0.0 !uBodyIntercept(n)
              !vIP = vIP !+ coeffGCMD(iRow,n)* 0.0 !vBodyIntercept(n) 
              !wIP = wIP !+ coeffGCMD(iRow,n)* 0.0 !wBodyIntercept(n)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        uPrev = utemp(iGh,jGh,kGh)
        vPrev = vtemp(iGh,jGh,kGh)
        wPrev = wtemp(iGh,jGh,kGh)
 
!       Apply Dirichlet conditions on Ghost Nodes
!       -----------------------------------------
        utemp(iGh,jGh,kGh) =  - uIP*imagePointWeight
        vtemp(iGh,jGh,kGh) =  - vIP*imagePointWeight
        wtemp(iGh,jGh,kGh) =  - wIP*imagePointWeight

!        utemp(iGh,jGh,kGh) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
!        vtemp(iGh,jGh,kGh) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
!        wtemp(iGh,jGh,kGh) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight
 
!       Compute residual
!       ----------------
        resVel = ABS( utemp(iGh,jGh,kGh)-uPrev ) &
               + ABS( vtemp(iGh,jGh,kGh)-vPrev ) &
               + ABS( wtemp(iGh,jGh,kGh)-wPrev )

        resVelMax=MAX(resVelMax,resVel)
 
      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resVelMax)
      CALL par_comm_var(utemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(vtemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(wtemp,nxc,nyc,nzc,Ngl,myTime)

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax .AND. resVelMax .GT. restol ) THEN
      PRINT*,'[LPCE] GCM_vel_set_bc for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell Velocity convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF
	
	
	Do n=1, nGhost
    
     iGh = iGhost(n)
     jGh = jGhost(n)
     kGh = kGhost(n)
     
     xGC = xc(iGh)
     yGC = yc(jGh)
     zGC = zc(kGh)
     
     xBI = xBodyIntercept(n)
     yBI = yBodyIntercept(n)
     zBI = zBodyIntercept(n)
     
     dslength = sqrt( (xBI-xGC)**2 + (yBI-yGC)**2 + (zBI-zGC)**2 )
     
     xBIn = (xBI-xGC)/dslength
     yBIn = (yBI-yGC)/dslength
     zBIn = (zBI-zGC)/dslength
     
     Velnorm = utemp(iGh,jGh,kGh)*xBIn + vtemp(iGh,jGh,kGh)*yBIn + wtemp(iGh,jGh,kGh)*zBIn
     
     utemp(iGh,jGh,kGh) = Velnorm*xBIn
     vtemp(iGh,jGh,kGh) = Velnorm*yBIn
     wtemp(iGh,jGh,kGh) = Velnorm*zBIn
     
    ENDDO ! nGhost 


END SUBROUTINE GCM_LPCE_vel_set_bc_norm
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE GCM_LPCE_vel_set_bc_noslip(utemp,vtemp,wtemp)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE mpi

    IMPLICIT NONE	

	REAL(KIND=CGREAL) :: utemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: vtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: wtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    
!   Loop variables
!   --------------
    INTEGER :: iBody, iRow, k, n

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iGh, jGh, kGh, iterGC
    INTEGER :: ii,jj,kk
	   
	REAL(KIND=CGREAL)  :: uIP, vIP, wIP,        &
                          resVelMax, resVel,    &
                          uPrev, vPrev, wPrev,  &
                          sameSideGhost



!   Compute body intercept velocity field
!   -------------------------------------
!   CALL GCM_SetBodyInterceptValues()
   
!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resVelMax = 1.0D10
      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC<itermax) .AND. (resVelMax>restol))

      resVelMax = 0.0_CGREAL

      DO n = 1,nGhost
        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
             CALL flow_stop
              STOP
         ENDIF

!       Initialize values
!       -----------------
        uIP = 0.0_CGREAL
        vIP = 0.0_CGREAL
        wIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
           ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + inck(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) .EQV. .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            uIP = uIP + coeffGCMD(iRow,n)*(    utemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                        +utemp(ii,jj,kk)*            sameSideGhost)

            vIP = vIP + coeffGCMD(iRow,n)*(    vtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +vtemp(ii,jj,kk)*            sameSideGhost)

            wIP = wIP + coeffGCMD(iRow,n)*(    wtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +wtemp(ii,jj,kk)*            sameSideGhost)
          ! ELSE
              !uIP = uIP !+ coeffGCMD(iRow,n)* 0.0 !uBodyIntercept(n)
              !vIP = vIP !+ coeffGCMD(iRow,n)* 0.0 !vBodyIntercept(n) 
              !wIP = wIP !+ coeffGCMD(iRow,n)* 0.0 !wBodyIntercept(n)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        uPrev = utemp(iGh,jGh,kGh)
        vPrev = vtemp(iGh,jGh,kGh)
        wPrev = wtemp(iGh,jGh,kGh)
 
!       Apply Dirichlet conditions on Ghost Nodes
!       -----------------------------------------
        utemp(iGh,jGh,kGh) =  - uIP*imagePointWeight
        vtemp(iGh,jGh,kGh) =  - vIP*imagePointWeight
        wtemp(iGh,jGh,kGh) =  - wIP*imagePointWeight

!        utemp(iGh,jGh,kGh) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
!        vtemp(iGh,jGh,kGh) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
!        wtemp(iGh,jGh,kGh) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight
 
!       Compute residual
!       ----------------
        resVel = ABS( utemp(iGh,jGh,kGh)-uPrev ) &
               + ABS( vtemp(iGh,jGh,kGh)-vPrev ) &
               + ABS( wtemp(iGh,jGh,kGh)-wPrev )

        resVelMax=MAX(resVelMax,resVel)
 
      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resVelMax)
      CALL par_comm_var(utemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(vtemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(wtemp,nxc,nyc,nzc,Ngl,myTime)

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax .AND. resVelMax .GT. restol ) THEN
      PRINT*,'[LPCE] GCM_vel_set_bc for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell Velocity convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF


END SUBROUTINE GCM_LPCE_vel_set_bc_noslip
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE GCM_LPCE_vel_set_bc_SLIP(utemp,vtemp,wtemp)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE grid_arrays
    USE unstructured_surface_arrays
	USE mpi
	
    IMPLICIT NONE

	REAL(KIND=CGREAL) :: utemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: vtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: wtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    
!   Loop variables
!   --------------
    INTEGER :: iBody, iRow, k, n

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iGh, jGh, kGh, iterGC
    INTEGER :: ii,jj,kk
	   
    REAL(KIND=CGREAL)  :: uIP, vIP, wIP,        &
                          resVelMax, resVel,    &
                          uPrev, vPrev, wPrev,  &
                          sameSideGhost,        &
						  dslength, Velnorm,    &
						  xBI,  yBI,  zBI,      &
						  xBIn, yBIn, zBIn,     &
						  xGC,  yGC,  zGC



!   Compute body intercept velocity field
!   -------------------------------------
!   CALL GCM_SetBodyInterceptValues()
   
!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resVelMax = 1.0D10
      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC<itermax) .AND. (resVelMax>restol))

      resVelMax = 0.0_CGREAL

      DO n = 1,nGhost
        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
             CALL flow_stop
              STOP
         ENDIF

!       Initialize values
!       -----------------
        uIP = 0.0_CGREAL
        vIP = 0.0_CGREAL
        wIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
           ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + inck(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) .EQV. .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            uIP = uIP + coeffGCMN(iRow,n)*(    utemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                        +utemp(ii,jj,kk)*            sameSideGhost)

            vIP = vIP + coeffGCMN(iRow,n)*(    vtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +vtemp(ii,jj,kk)*            sameSideGhost)

            wIP = wIP + coeffGCMN(iRow,n)*(    wtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +wtemp(ii,jj,kk)*            sameSideGhost)
          ! ELSE
              !uIP = uIP !+ coeffGCMD(iRow,n)* 0.0 !uBodyIntercept(n)
              !vIP = vIP !+ coeffGCMD(iRow,n)* 0.0 !vBodyIntercept(n) 
              !wIP = wIP !+ coeffGCMD(iRow,n)* 0.0 !wBodyIntercept(n)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        uPrev = utemp(iGh,jGh,kGh)
        vPrev = vtemp(iGh,jGh,kGh)
        wPrev = wtemp(iGh,jGh,kGh)
 
!       Apply Dirichlet conditions on Ghost Nodes
!       -----------------------------------------
        utemp(iGh,jGh,kGh) =   uIP
        vtemp(iGh,jGh,kGh) =   vIP
        wtemp(iGh,jGh,kGh) =   wIP

!        utemp(iGh,jGh,kGh) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
!        vtemp(iGh,jGh,kGh) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
!        wtemp(iGh,jGh,kGh) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight
 
!       Compute residual
!       ----------------
        resVel = ABS( utemp(iGh,jGh,kGh)-uPrev ) &
               + ABS( vtemp(iGh,jGh,kGh)-vPrev ) &
               + ABS( wtemp(iGh,jGh,kGh)-wPrev )

        resVelMax=MAX(resVelMax,resVel)
 
      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resVelMax)
      CALL par_comm_var(utemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(vtemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(wtemp,nxc,nyc,nzc,Ngl,myTime)

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax .AND. resVelMax .GT. restol ) THEN
      PRINT*,'[LPCE] GCM_vel_set_bc for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell (SLIP) Velocity convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF
	
	! Correct Normal Velocity
	
	Do n=1, nGhost
    
     iGh = iGhost(n)
     jGh = jGhost(n)
     kGh = kGhost(n)
     
     xGC = xc(iGh)
     yGC = yc(jGh)
     zGC = zc(kGh)
     
     xBI = xBodyIntercept(n)
     yBI = yBodyIntercept(n)
     zBI = zBodyIntercept(n)
     
     dslength = sqrt( (xBI-xGC)**2 + (yBI-yGC)**2 + (zBI-zGC)**2 )
     
     xBIn = (xBI-xGC)/dslength
     yBIn = (yBI-yGC)/dslength
     zBIn = (zBI-zGC)/dslength
     
     Velnorm = utemp(iGh,jGh,kGh)*xBIn + vtemp(iGh,jGh,kGh)*yBIn + wtemp(iGh,jGh,kGh)*zBIn
     
     utemp(iGh,jGh,kGh) = utemp(iGh,jGh,kGh) - 2.*Velnorm*xBIn
     vtemp(iGh,jGh,kGh) = vtemp(iGh,jGh,kGh) - 2.*Velnorm*yBIn
     wtemp(iGh,jGh,kGh) = wtemp(iGh,jGh,kGh) - 2.*Velnorm*zBIn
     
    ENDDO ! nGhost 


END SUBROUTINE GCM_LPCE_vel_set_bc_SLIP
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE GCM_LPCE_p_set_bc(pres)

! ------------------------------------------------------------------------
!  Compute Ghost point values at internal boundary for the pressure field
!     
!   P  =   P   
!    gp     ip 
!
!
! and
!        (          ) 
!   P  = (coeffGCMN ) P
!    ip  (         i)  i
! ------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE mpi
    
    IMPLICIT NONE

    REAL(KIND=CGREAL) :: pres(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

!   Loop variables
!   --------------
    INTEGER :: n, iBody, iRow

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iterGC
    INTEGER :: iGh, jGh, kGh
    INTEGER :: ii, jj, kk 
    INTEGER :: i, j, k

    REAL(KIND=CGREAL) :: resPresMax, resPres, &
                         pIP, pPrev, sameSideGhost



!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resPresMax = 1.0E10_CGREAL

      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC < itermax) .AND. (resPresMax > restol))

      resPresMax = 0.0_CGREAL

      DO n = 1,nGhost

        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
          WRITE(STDOUT,'(4(A,I3))') 'Invalid body number for ghost cell (',i,',',j,',',k,'), reported by process ',myRank
          CALL flow_stop
          STOP
        END IF

!       Initialize values
!       -----------------    
        pIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
          ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + incK(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) == .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            pIP = pIP + coeffGCMN(iRow,n)*( pres(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +pres(ii,jj,kk)*            sameSideGhost)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        pPrev = pres(iGh,jGh,kGh)
 
!       Apply the Neumann condition on Ghost Nodes
!       ------------------------------------------
        pres(iGh,jGh,kGh) = pIP
 
!       Compute residual
!       ----------------
        resPres = ABS( pres(iGh,jGh,kGh) - pPrev )

        resPresMax = MAX(resPres,resPresMax)

      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resPresMax)
      CALL par_comm_var(pres,nxc,nyc,nzc,Ngl,myTime)
    
    ENDDO ! iterGC

    IF (resPresMax>restol) THEN
      PRINT*,'GCM_p_set_bc_internal for iBody :', iBody
      PRINT*,'GhostCell Pressure did not converge in ',itermax,' iterations'
      PRINT*,'Final residual = ',resPresMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell Pressure convergence : ',iterGC,resPresMax
      END IF
    ENDIF

!   Ensure that pGhost is zero for all non-ghost cells.
!   This is used as a test in update_pressure_freshcells
!   ----------------------------------------------------
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      pres(i,j,k) = pres(i,j,k)*REAL(1-iblank_solid(i,j,k)*(1-ghostCellMark(i,j,k)),KIND=CGREAL)
    ENDDO
    ENDDO
    ENDDO


END SUBROUTINE GCM_LPCE_p_set_bc     
!----------------------------------------------------------------------


SUBROUTINE update_new_presGC(var)

! -------------------------------------------------------------------------------------------------
!  A fresh-cell was a ghost-cell in the previous time step and its value is stored in pGhost array
!  At current time-step we need to transfer this value to pPrime so it can be correctly used
!  as initial guess for current time-step. POC RM
! -------------------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE
    
    REAL(KIND=CGREAL), INTENT (INOUT) ::var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    REAL(KIND=CGREAL) :: denom, myTime

    INTEGER :: i,j,k,ii,jj,kk,n,iBody
    INTEGER :: iG,jG,kG,ib,ibG,ibF
    INTEGER :: iBeg,iEnd,jBeg,jEnd,kBeg,kEnd

!     Now obtain estimate of pressure for newly created ghost cell.
!     -------------------------------------------------------------
      DO n=1,nGhost
        iG = iGhost(n)
        jG = jGhost(n)
        kG = kGhost(n)
        
        IF (iG<1 .OR. iG>nxc .OR. &
            jG<1 .OR. jG>nyc .OR. &
            kG<1 .OR. kG>nzc) THEN
          IF (ImTheBOSS) WRITE(*,*) 'Error in ghost cell number ',n,'. Location out of domain. Terminating...'
          CALL flow_stop
          STOP
        END IF

        iBody = bodyNum(iG,jG,kG)

        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
          IF (ImTheBOSS) WRITE(*,*) 'A ghost cell is assigned to a non-solid, non-membrane body! Terminating ...'
          CALL flow_stop
          stop
        END IF

!       this identifies newly created ghost-cell:
!       -----------------------------------------
        IF ( ABS(var(iG,jG,kG)) <= 1.0E-10 ) THEN   
          
          denom = 0.0_CGREAL
          
          iBeg = MAX(myILL,iG-1)
          iEnd = MIN(myIUL,iG+1)
          jBeg = MAX(myJLL,jG-1)
          jEnd = MIN(myJUL,jG+1)
          kBeg = MAX(1    ,kG-1)
          kEnd = MIN(nzc  ,kG+1)

          DO kk=kBeg,kEnd
          DO jj=jBeg,jEnd
          DO ii=iBeg,iEnd
            IF ( unstruc_surface_type(iBody) == SOLID_BODY) THEN
              ibG = iblank_solid(iG,jG,kG) 
              ib  = iblank_solid(ii,jj,kk) 
            ELSE
              IF (ALLOCATED(iblank_memb)) THEN
                ibG = iblank_memb(iG,jG,kG,iBody-nBody_solid)
                ib  = iblank_memb(ii,jj,kk,iBody-nBody_solid)
              ENDIF
            ENDIF
            IF (ibG /= ib) THEN  

!             pGhost is obtained as mean of adjacent fluid nodal values
!             ---------------------------------------------------------
              var(iG,jG,kG) =  var(iG,jG,kG) + var(ii,jj,kk) 
              denom            = denom + 1.0_CGREAL
            ENDIF
          ENDDO
          ENDDO
          ENDDO

          var(iG,jG,kG) =  var(iG,jG,kG)/denom
   
        ENDIF ! pGhost

      ENDDO ! nGhost
    

END SUBROUTINE update_new_presGC
!---------------------------------------------------------------------


SUBROUTINE GCM_vel0_set_bc(utemp,vtemp,wtemp)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)
!    gp                       ip      b
!
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE mpi

    IMPLICIT NONE	

	REAL(KIND=CGREAL) :: utemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: vtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
	REAL(KIND=CGREAL) :: wtemp(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    
!   Loop variables
!   --------------
    INTEGER :: iBody, iRow, k, n

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iGh, jGh, kGh, iterGC
    INTEGER :: ii,jj,kk
	   
	REAL(KIND=CGREAL)  :: uIP, vIP, wIP,        &
                          resVelMax, resVel,    &
                          uPrev, vPrev, wPrev,  &
                          sameSideGhost



!   Compute body intercept velocity field
!   -------------------------------------
!   CALL GCM_SetBodyInterceptValues()
   
!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resVelMax = 1.0D10
	
! get guess for ghost cell value from flow data...
	
	DO n = 1,nGhost
      
	  iGh = iGhost(n)
      jGh = jGhost(n)
      kGh = kGhost(n)
	  
	  utemp(iGh,jGh,kGh) = uGhost(iGh,jGh,kGh)
	  vtemp(iGh,jGh,kGh) = vGhost(iGh,jGh,kGh)
	  wtemp(iGh,jGh,kGh) = wGhost(iGh,jGh,kGh)
	  
	ENDDO
	  
		
      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC<itermax) .AND. (resVelMax>restol))

      resVelMax = 0.0_CGREAL

      DO n = 1,nGhost
        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
             CALL flow_stop
              STOP
         ENDIF

!       Initialize values
!       -----------------
        uIP = 0.0_CGREAL
        vIP = 0.0_CGREAL
        wIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
           ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + inck(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) .EQV. .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            uIP = uIP + coeffGCMD(iRow,n)*(    utemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                        +utemp(ii,jj,kk)*            sameSideGhost)

            vIP = vIP + coeffGCMD(iRow,n)*(    vtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +vtemp(ii,jj,kk)*            sameSideGhost)

            wIP = wIP + coeffGCMD(iRow,n)*(    wtemp(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +wtemp(ii,jj,kk)*            sameSideGhost)
          ELSE
              uIP = uIP + coeffGCMD(iRow,n)* uBodyIntercept(n)
              vIP = vIP + coeffGCMD(iRow,n)* vBodyIntercept(n) 
              wIP = wIP + coeffGCMD(iRow,n)* wBodyIntercept(n)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        uPrev = utemp(iGh,jGh,kGh)
        vPrev = vtemp(iGh,jGh,kGh)
        wPrev = wtemp(iGh,jGh,kGh)
 
!       Apply Dirichlet conditions on Ghost Nodes
!       -----------------------------------------
        utemp(iGh,jGh,kGh) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
        vtemp(iGh,jGh,kGh) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
        wtemp(iGh,jGh,kGh) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight
 
!       Compute residual
!       ----------------
        resVel = ABS( utemp(iGh,jGh,kGh)-uPrev ) &
               + ABS( vtemp(iGh,jGh,kGh)-vPrev ) &
               + ABS( wtemp(iGh,jGh,kGh)-wPrev )

        resVelMax=MAX(resVelMax,resVel)
 
      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resVelMax)
      CALL par_comm_var(utemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(vtemp,nxc,nyc,nzc,Ngl,myTime)
      CALL par_comm_var(wtemp,nxc,nyc,nzc,Ngl,myTime)

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax .AND. resVelMax .GT. restol ) THEN
      PRINT*,'[LPCE] GCM_vel_set_bc for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell Vel0 convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF


END SUBROUTINE GCM_vel0_set_bc
!-----------------------------------------------------------------------------------------------------------------------



SUBROUTINE GCM_p0_set_bc(pres)

! ------------------------------------------------------------------------
!  Compute Ghost point values at internal boundary for the pressure field
!     
!   P  =   P   
!    gp     ip 
!
!
! and
!        (          ) 
!   P  = (coeffGCMN ) P
!    ip  (         i)  i
! ------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
	USE mpi
    
    IMPLICIT NONE

    REAL(KIND=CGREAL) :: pres(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

!   Loop variables
!   --------------
    INTEGER :: n, iBody, iRow

!   Local variables
!   ---------------
    INTEGER :: myTime
    INTEGER :: iterGC
    INTEGER :: iGh, jGh, kGh
    INTEGER :: ii, jj, kk 
    INTEGER :: i, j, k

    REAL(KIND=CGREAL) :: resPresMax, resPres, &
                         pIP, pPrev, sameSideGhost



!   Loop over all immersed bodies
!   -----------------------------
    iterGC    = 0
    resPresMax = 1.0E10_CGREAL
	
	DO n = 1,nGhost

      iGh = iGhost(n)
      jGh = jGhost(n)
      kGh = kGhost(n)
	  
	  pres(iGh,jGh,kGh) = pGhost(iGh,jGh,kGh)
	  
	ENDDO

      
!   Iterate to correct interdependent ghost points
!   ----------------------------------------------
    DO WHILE ((iterGC < itermax) .AND. (resPresMax > restol))

      resPresMax = 0.0_CGREAL

      DO n = 1,nGhost

        iGh = iGhost(n)
        jGh = jGhost(n)
        kGh = kGhost(n)

        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
          WRITE(STDOUT,'(4(A,I3))') 'Invalid body number for ghost cell (',i,',',j,',',k,'), reported by process ',myRank
          CALL flow_stop
          STOP
        END IF

!       Initialize values
!       -----------------    
        pIP = 0.0_CGREAL

!       Compute velocity field at Image points
!       --------------------------------------
        DO iRow = 1, iRowMax
          ii = iCellIndex(n) + incI(iRow) 
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + incK(iRow)

          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /=kGh ) THEN

            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) == .TRUE. ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            pIP = pIP + coeffGCMN(iRow,n)*( pres(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +pres(ii,jj,kk)*            sameSideGhost)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        pPrev = pres(iGh,jGh,kGh)
 
!       Apply the Neumann condition on Ghost Nodes
!       ------------------------------------------
        pres(iGh,jGh,kGh) = pIP
 
!       Compute residual
!       ----------------
        resPres = ABS( pres(iGh,jGh,kGh) - pPrev )

        resPresMax = MAX(resPres,resPresMax)

      ENDDO ! n    

      iterGC = iterGC + 1
	  
	  CALL par_getMaxReal(resPresMax)
      CALL par_comm_var(pres,nxc,nyc,nzc,Ngl,myTime)
    
    ENDDO ! iterGC

    IF (resPresMax>restol) THEN
      PRINT*,'GCM_p_set_bc_internal for iBody :', iBody
      PRINT*,'GhostCell Pressure did not converge in ',itermax,' iterations'
      PRINT*,'Final residual = ',resPresMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') '[LPCE] GhostCell P0 convergence : ',iterGC,resPresMax
      END IF
    ENDIF

!   Ensure that pGhost is zero for all non-ghost cells.
!   This is used as a test in update_pressure_freshcells
!   ----------------------------------------------------
!    DO k = 1,nzc
!    DO j = 1,nyc
!    DO i = 1,nxc
!      pres(i,j,k) = pres(i,j,k)*REAL(1-iblank_solid(i,j,k)*(1-ghostCellMark(i,j,k)),KIND=CGREAL)
!    ENDDO
!    ENDDO
!    ENDDO


END SUBROUTINE GCM_p0_set_bc     
!----------------------------------------------------------------------
