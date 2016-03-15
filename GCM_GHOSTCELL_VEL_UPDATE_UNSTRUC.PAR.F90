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
!  Filename: GCM_GHOSTCELL_VEL_UPDATE_UNSTRUC.PAR.F90
!  Latest Modification: October 22, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_GhostCell_Vel()
!     GCM_vel_set_bc_internal()
!     GCM_p_set_bc_internal(pres)
! --------------------------------------------------------------------



SUBROUTINE GCM_GhostCell_Vel() 

! -------------------------------------------------------------------------
!  update u*, v* and w* values for ghost cell
!
!
!  Ghost cell velocity satisfies the following equations
!
!  Integral [ DEL.(u*) dv ] =- Uo.nDS  with coeff of "dead" faces = 0
!
! [                                                               ]
! [ U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)   ] . tau
! [  gp                       ip      b                           ]   ---
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i
! -------------------------------------------------------------------------


    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE flow_arrays
    USE GCM_arrays
    USE mpi

    IMPLICIT NONE

!   Loop variables
!   --------------
    INTEGER :: iterGC, n, iBody, iRow
    
!   Local variables
!   ---------------
    INTEGER :: iGh, jGh, kGh
    INTEGER  :: ii, jj, kk

    REAL(KIND=CGREAL) :: resVelMax, resVel
    REAL(KIND=CGREAL) :: uIP, vIP, wIP
    REAL(KIND=CGREAL) :: uPrev, vPrev, wPrev
    REAL(KIND=CGREAL) :: sameSideGhost
    REAL(KIND=CGREAL) :: myTime
    LOGICAL :: flag

!DEBGU
  INTEGER :: ierr
!ENDDEBUG

! ------------------------------------------------
!  Iterate to correct interdependent ghost points
!  Initialize values
! ------------------------------------------------
    iterGC    = 0
    resVelMax = 1.0E10_CGREAL
    DO WHILE ((iterGC < itermax) .AND. (resVelMax > restol))

      resVelMax = 0.0_CGREAL

      DO n = 1, nGhost

        iGh=iGhost(n)
        jGh=jGhost(n)
        kGh=kGhost(n)
        
        iBody = bodyNum(iGh,jGh,kGh)
        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
            
            WRITE(*,*) 'GCM_GhostCell_Vel(): Bad body number!'
            CALL flow_stop
            STOP
        END IF

        uIP = 0.0_CGREAL
        vIP = 0.0_CGREAL
        wIP = 0.0_CGREAL
      DO iRow = 1, iRowMax

          ii = iCellIndex(n) + incI(iRow)
          jj = jCellIndex(n) + incJ(iRow)
          kk = kCellIndex(n) + incK(iRow)
  
          IF ( ii /= iGh .OR. jj /= jGh .OR. kk /= kGh) THEN
                  
            sameSideGhost = ghostCellMark(ii,jj,kk)

            IF (unstruc_surface_type(iBody) == MEMBRANE) THEN
              IF ( ALLOCATED(iblank_memb) ) THEN
                IF (   iblank_memb(iGh,jGh,kGh,iBody-nBody_solid)  &
                    /= iblank_memb(ii,jj,kk,iBody-nBody_solid) ) sameSideGhost = 0.0_CGREAL
              ENDIF ! ALLOCATED(iblank_memb) 
            ENDIF ! unstruc_surface_type

            uIP = uIP + coeffGCMD(iRow,n)*(   u(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                       +   uGhost(ii,jj,kk)*            sameSideGhost)
            vIP = vIP + coeffGCMD(iRow,n)*(   v(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                       +  vGhost(ii,jj,kk)*            sameSideGhost)
            wIP = wIP + coeffGCMD(iRow,n)*(   w(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                       +  wGhost(ii,jj,kk)*            sameSideGhost)

           ELSE
            uIP = uIP + coeffGCMD(iRow,n)* uBodyIntercept(n)
            vIP = vIP + coeffGCMD(iRow,n)* vBodyIntercept(n)
            wIP = wIP + coeffGCMD(iRow,n)* wBodyIntercept(n)
          ENDIF ! ii

        ENDDO ! iRow
        uPrev = uGhost(iGh,jGh,kGh)
        vPrev = vGhost(iGh,jGh,kGh)
        wPrev = wGhost(iGh,jGh,kGh)
        

        uGhost(iGh,jGh,kGh) =  uBodyIntercept(n)*bodyInterceptWeight - uIP*imagePointWeight
        vGhost(iGh,jGh,kGh) =  vBodyIntercept(n)*bodyInterceptWeight - vIP*imagePointWeight
        wGhost(iGh,jGh,kGh) =  wBodyIntercept(n)*bodyInterceptWeight - wIP*imagePointWeight
!       Compute residual
!       ----------------
        resVel = ABS( uGhost(iGh,jGh,kGh)-uPrev ) &
               + ABS( vGhost(iGh,jGh,kGh)-vPrev ) &
               + ABS( wGhost(iGh,jGh,kGh)-wPrev )

        resVelMax = MAX(resVel,resVelMax)
      ENDDO ! 
#ifdef MPI
    CALL par_getMaxReal(resVelMax)
    CALL par_comm_var(uGhost,nxc,nyc,nzc,Ngl,myTime)
    CALL par_comm_var(vGhost,nxc,nyc,nzc,Ngl,myTime)
    CALL par_comm_var(wGhost,nxc,nyc,nzc,Ngl,myTime)
#endif
      IF ( monitorON ) WRITE(STDOUT,'(7X,A,I10,1PE15.7)') ' Ghostcell Velocity Convergence: ',iterGC,resVelMax

      iterGC = iterGC + 1

    ENDDO ! iterGC


    IF ( iterGC == itermax .AND. resVelMax > restol .AND. monitorON) THEN
      PRINT*,'GCM_vel_set_bc_internal for iBody :', iBody
      PRINT*,'   GhostCell did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
      
!      CALL flow_stop
!      STOP
!    ELSE
!      IF (monitorON) THEN
!        WRITE(STDOUT,'(5X,A)') 'GhostCell u* convergence : k=',k,iterGC,resVelMax
!      END IF ! ntime
    ENDIF

    CALL GCM_enforce_global_mass_consv()

END SUBROUTINE  GCM_GhostCell_Vel
!----------------------------------------------------------------------



SUBROUTINE GCM_vel_set_bc_internal()

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
    
    IMPLICIT NONE

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
    CALL GCM_SetBodyInterceptValues()
   
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

            uIP = uIP + coeffGCMD(iRow,n)*(    u(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                        +uGhost(ii,jj,kk)*            sameSideGhost)

            vIP = vIP + coeffGCMD(iRow,n)*(    v(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +vGhost(ii,jj,kk)*            sameSideGhost)

            wIP = wIP + coeffGCMD(iRow,n)*(    w(ii,jj,kk)*(1.0_CGREAL-sameSideGhost) &
                                         +wGhost(ii,jj,kk)*            sameSideGhost)
          ELSE
              uIP = uIP + coeffGCMD(iRow,n)* uBodyIntercept(n)
              vIP = vIP + coeffGCMD(iRow,n)* vBodyIntercept(n) 
              wIP = wIP + coeffGCMD(iRow,n)* wBodyIntercept(n)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        uPrev = uGhost(iGh,jGh,kGh)
        vPrev = vGhost(iGh,jGh,kGh)
        wPrev = wGhost(iGh,jGh,kGh)
 
!       Apply Dirichlet conditions on Ghost Nodes
!       -----------------------------------------
        uGhost(iGh,jGh,kGh) = uBodyIntercept(n) *bodyInterceptWeight - uIP*imagePointWeight
        vGhost(iGh,jGh,kGh) = vBodyIntercept(n) *bodyInterceptWeight - vIP*imagePointWeight
        wGhost(iGh,jGh,kGh) = wBodyIntercept(n) *bodyInterceptWeight - wIP*imagePointWeight
 
!       Compute residual
!       ----------------
        resVel = ABS( uGhost(iGh,jGh,kGh)-uPrev ) &
               + ABS( vGhost(iGh,jGh,kGh)-vPrev ) &
               + ABS( wGhost(iGh,jGh,kGh)-wPrev )

        resVelMax=MAX(resVelMax,resVel)
 
      ENDDO ! n    

      iterGC = iterGC + 1

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax .AND. resVelMax .GT. restol ) THEN
      PRINT*,'GCM_vel_set_bc_internal for iBody :', iBody
      PRINT*,'   GhostCell Velocity did not converge in ',itermax,' iterations'
      PRINT*,'   Final residual = ',resVelMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') 'GhostCell Velocity convergence : ',iterGC,resVelMax
      END IF ! ntime
    ENDIF

#ifdef MPI
    CALL par_comm_var(uGhost,nxc,nyc,nzc,Ngl,myTime)
    CALL par_comm_var(vGhost,nxc,nyc,nzc,Ngl,myTime)
    CALL par_comm_var(wGhost,nxc,nyc,nzc,Ngl,myTime)
#endif

END SUBROUTINE GCM_vel_set_bc_internal



SUBROUTINE GCM_p_set_bc_internal(pres)

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
    
    IMPLICIT NONE

    REAL(KIND=CGREAL), INTENT(IN) :: pres(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

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
                                         +pGhost(ii,jj,kk)*            sameSideGhost)
          ENDIF ! ii          

        ENDDO ! iRow

!       Load temporary values
!       ---------------------
        pPrev = pGhost(iGh,jGh,kGh)
 
!       Apply the Neumann condition on Ghost Nodes
!       ------------------------------------------
        pGhost(iGh,jGh,kGh) = pIP
 
!       Compute residual
!       ----------------
        resPres = ABS( pGhost(iGh,jGh,kGh) - pPrev )

        resPresMax = MAX(resPres,resPresMax)

      ENDDO ! n    

      iterGC = iterGC + 1

    ENDDO ! iterGC

    IF (resPresMax>restol) THEN
      PRINT*,'GCM_p_set_bc_internal for iBody :', iBody
      PRINT*,'GhostCell Pressure did not converge in ',itermax,' iterations'
      PRINT*,'Final residual = ',resPresMax
    ELSE
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,I5,D15.5)') 'GhostCell Pressure convergence : ',iterGC,resPresMax
      END IF
    ENDIF

!   Ensure that pGhost is zero for all non-ghost cells.
!   This is used as a test in update_pressure_freshcells
!   ----------------------------------------------------
    DO k = 1,nzc
    DO j = 1,nyc
    DO i = 1,nxc
      pGhost(i,j,k) = pGhost(i,j,k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)
    ENDDO
    ENDDO
    ENDDO

#ifdef MPI
    CALL par_comm_var(pGhost,nxc,nyc,nzc,Ngl,myTime)
#endif

END SUBROUTINE GCM_p_set_bc_internal      
!----------------------------------------------------------------------
