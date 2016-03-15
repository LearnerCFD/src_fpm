! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  The George Washington University
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
!  Filename: CORE_TIME_STEP.PAR.F90
!  Latest Modification: Dec, 29 2010 (PAT 2.1.0)
!  by Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     time_step_viscous()
! --------------------------------------------------------------------


SUBROUTINE time_step_viscous()

  USE global_parameters
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
  USE boundary_arrays  
  USE usr_module  
  USE derivative 
  USE grid_arrays
  USE cutcell_arrays
  USE scalar 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
  USE flow_parameters
  USE flow_arrays
  USE pressure_arrays
  USE multiuse_arrays
  USE stat_arrays
  USE finite_element_parameters           ! SER_TO_PAR. QX. CH24
  USE finite_element_arrays
  USE unstructured_surface_arrays
  USE turb_parameters
  USE turb_arrays
  USE turbulent_channel
  USE tahoe_parameters             ! Added by Rajneesh
  USE implicit_coupling_parameters  ! Added by Rajneesh
  


# ifdef MPI
    use mpi
# endif

  IMPLICIT NONE

! Local variables
! ---------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
  INTEGER :: iGroup        
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
  INTEGER :: iBody  ,i,j,k,i1,i2,i3
  INTEGER :: clock1, clock2, clock_rate,mytime
  character text  
  
  REAL(KIND=CGREAL) :: sum, time0
  REAL(KIND=CGREAL) :: startTime, endTime, ppsTime, ppsCommTime, adsTime, adsCommTime
  
! Profiling Variables
! -------------------
# ifdef MPI
    REAL(KIND=CGREAL) :: ThisTimeEnd, elapsedTime
# else
    INTEGER           :: ThisTimeEnd, clock_rate
    REAL(KIND=CGREAL) :: elapsedTime
# endif

  INTEGER :: pwks,pdys,phrs,pmts,psec
  INTEGER :: ierr, status

  REAL(KIND=CGREAL) :: vartmp
  
  REAL(KIND=CGREAL) :: dPdXMean

      
! Start time stepping loop
! ------------------------

DO ntime = ntime_start+1,ntime_start+no_tsteps        ! move solution from n --> n+1
!dddddddddddddddddddddddddddddddddddddddddd
uold = u
vold = v
wold = w
!dddddddddddddddddddddddddddddddddddddddddd


    
    time = time + dt
    time0= time - ninit*dt

    
    IF (MOD(ntime,nmonitor) == 0) THEN
      monitorIT = .TRUE.
      monitorON = ImtheBOSS
    ELSE
      monitorIT = .FALSE.
      monitorON = .FALSE.
    END IF
    
    dumpIT = .FALSE.
    IF (time0>0.d0) THEN
      IF (MOD(ntime,ndump) == 0) dumpIT = .TRUE.
    ELSE
      IF (MOD(ntime,qdump) == 0) dumpIT = .TRUE.
    END IF

    IF (monitorON) THEN  
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(2X,A)')'================================================================'
      WRITE(STDOUT,'(3X,A,I6,A,F15.7,A,I6)') 'NTIME = ',ntime,',  Time = ',time,',  NTIME TO GO = ',no_tsteps+ntime_start-ntime
      WRITE(STDOUT,'(3X,A)')'--------------------------------------------------------------'
    ENDIF

!   Let us have a pressure reference: ZERO!
!   ---------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    IF(fea_boundary_flag_tahoe /= 1 .AND. implicit_coupling_flag_combined_motion /= 1) CALL subtract_mean_pressure()!...............................COMPLETE(SAMK)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
!   ---------------------------------------------------
!    Compute turbulent viscosity
!    Compute AD coefficients for non-moving boundaries
!    Set BC for viscosity
!    Add molecular viscosity
!    Note: nuTot is computed at n-time step level 
!   ---------------------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    IF(VDVactive == 1) CALL advection_term   
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    CALL CALCULATE_DIFFUSION(u,ughost,1,diff_prev)      
    CALL CALCULATE_DIFFUSION(v,vghost,2,diff_prev)
    IF(ndim == 3)THEN
    CALL CALCULATE_DIFFUSION(w,wghost,3,diff_prev)
    ENDIF
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

    IF ( turbActive == ACTIVE ) THEN
      IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Entering TURB_CalcVisc '
      CALL TURB_Calcvisc()

      IF(monitorON) WRITE(*,*) 'Entering TURB_Visc_set_bc '  
      CALL TURB_Visc_set_bc()

     IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Entering TURB_Visc_SetBoundCells '
     CALL TURB_Visc_SetBoundCells()

 #ifdef MPI
       CALL par_comm_var(visTurb,nxc,nyc,nzc,Ngl,myTime)
 #endif

      DO k =0, nzc+1
      DO j =0, nyc+1
      DO i =0, nxc+1
        viscTot(i,j,k) = visTurb(i,j,k) +reInv
        bcxvisc(i,j,k) = bcxvisc(i,j,k) +reInv
        bcyvisc(i,j,k) = bcyvisc(i,j,k) +reInv
        bczvisc(i,j,k) = bczvisc(i,j,k) +reInv
      END DO ! i
      END DO ! j
      END DO ! k


      IF ( boundary_motion /= MOVING_BOUNDARY ) THEN
        IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Entering set_solve_ad()'
        CALL set_solve_ad()
      END IF ! boundary_motion 
    END IF   ! turbActive

!   ------------------------------------------------------------------
!    Compute advection-diffusion terms
!    Note: NL should be at time step level n, hence they are computed 
!          BEFORE moving the body boundary
!   ------------------------------------------------------------------

    IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Computing advection-diffusion RHS'
    CALL rhs_advec_diff()!.......................................COMPLETE(SAMK)

    ! Specific for turbulent channel flow
    IF (turbchannel==1) THEN

     CALL  TURB_Streamwise_Source(dPdXMean)

     DO k = 1, nzc
     DO j = 1, nyc
     DO i = 1, nxc
       nlu(i,j,k) = nlu(i,j,k) + dt*dPdXMean
     ENDDO ! i
     ENDDO ! j
     ENDDO ! k

    ENDIF
    
   
    IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Calculating Adam-Bashford Coefficients'
    CALL FreshCell_CalcExpWeight()!..............................COMPLETE(SAMK)
        
!   Move Body Boundary to (n+1) and compute coefficients
!   ----------------------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    k_derivative = 0
259 IF(nread == 0) k_derivative = k_derivative + 1    !1---u derivative 2---v derivative 3---q derivative 4---original

    IF(derivative_flag == 1 .AND. ntime > nstart_FIM)THEN
       CALL FIM_STORE_RESTORE()    
    ENDIF 
    
    IF(derivative_flag == 1 .AND. ntime > nstart_FIM)THEN
      CALL calculate_arclength_norm_ds()
    !  CALL FIM_calculate_new_iblank               
    !  CALL rhs_advec_diff() !Computing advection-diffusion RHS
    !  CALL set_solve_ad() ! advection eqs coefficents
    ENDIF     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

    IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime > 1) THEN
      IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Moving the body to the next time step'     
      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 	      
      IF(implicit_coupling_flag_combined_motion /= 1) CALL move_boundary()!......................................COMPLETE(SAMK)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)            
  
	  IF(iCC .eq. 1) THEN 
	  IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Calc. CutCell fractions'	  
	  CALL init_CutCell	
      ENDIF ! iCC  


      IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Entering set_solve_ad()'
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 	      
      IF(implicit_coupling_flag_combined_motion /= 1) CALL set_solve_ad()!.......................................COMPLETE(SAMK)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 	
	
      IF ( turbActive == ACTIVE ) THEN
        IF(monitorON) WRITE(*,*) 'Entering TURB_Visc_SetBoundCells '  
        CALL TURB_Visc_SetBoundCells()
      END IF ! turbActive 

    END IF ! boundary_motion
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)     
!===Loop for implicit coupling of combined motion(START) 
	kimplicit_FIM = 0
257 kimplicit_FIM = kimplicit_FIM + 1
    IF(ntime > max_ntime_implicit_FIM)THEN
     implicit_coupling_flag_combined_motion = 0
    ENDIF
    IF(implicit_coupling_flag_combined_motion == 1 .AND. ntime > nstart_FIM)THEN
	IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Entering implicit coupling'
	CALL FIM_implicit_coupling
	END IF 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  

!   Added by Rajneesh for implicit coupling
# ifdef TAHOE
	kimplicit = 0
256 kimplicit = kimplicit + 1
    IF(implicit_coupling_flag == 1 .AND. ntime>=nstart_tahoe)THEN
	IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Entering implicit coupling'
	CALL implicit_coupling
	END IF 
# endif
	
    DO iBody = 1,nBody
      IF (WALL_TYPE(iBody) == POROUS_OR_SLIP) THEN  
        CALL wall_velocity(iBody)!...............................COMPLETE(SAMK)
      ENDIF                                        
    END DO

    IF(monitorON .AND. boundary_formulation == SSM_METHOD) THEN
      WRITE(STDOUT,'(5X,A,2(1X,1PE15.7))') 'AREAX1; AREAX2 = ',areax1,areax2
      WRITE(STDOUT,'(5X,A,2(1X,1PE15.7))') 'AREAY1; AREAY2 = ',areay1,areay2
      WRITE(STDOUT,'(5X,A,2(1X,1PE15.7))') 'AREAZ1; AREAZ2 = ',areaz1,areaz2
      WRITE(STDOUT,*)
    END IF

!   Set viscosity for fresh cells
!   -----------------------------
    IF ( turbActive == ACTIVE ) THEN
      IF(monitorON) WRITE(*,*) 'Entering TURB_Visc_SetFreshCell '  
      CALL TURB_Visc_SetFreshCell()
    END IF ! turbActive

!   Compute weights and update RHS for Fresh cells
!   ----------------------------------------------    
    IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Nullifying RHS for Fresh Cells'
    CALL FreshCell_UpdateRhs()!..................................COMPLETE(SAMK)


!   Fill boundary arrays with (u)^n+1
!   ---------------------------------

    IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Setting Up Boundary Arrays'
    CALL set_bc()!...............................................COMPLETE(SAMK)

    IF(boundary_motion == MOVING_BOUNDARY) then
      IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Adjusting RHS for Fresh Cells'
      CALL rhs_adjust_fresh()!...................................COMPLETE(SAMK)
    end if

!   Adjust RHS for 2D computations
!   ------------------------------
    IF( nDim == DIM_2D) CALL rhs_adjust2D!.......................COMPLETE(SAMK)

!   Solve intermediate velocity field using Thomas algorithm
!   --------------------------------------------------------
#   ifdef MPI
      startTime = MPI_WTIME()
#   else
      call system_clock(clock1)
#   endif

    IF(monitorON) THEN
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(3X,A)') 'Solving Advection-Diffusion Equations'
    END IF

    CALL solve_ad(adsCommTime)!.......................................ONLY_PER(SAMK)
!ddddddddddddddddddddddddddddddddddddddddddddddddd    
    CALL CALCULATE_DIFFUSION(u,ughost,1,diff_star)
    CALL CALCULATE_DIFFUSION(v,vghost,2,diff_star)
    IF(ndim == 3)THEN
    CALL CALCULATE_DIFFUSION(w,wghost,3,diff_star)
    ENDIF
!ddddddddddddddddddddddddddddddddddddddddddddddddd

!   Time Profiling
!   --------------
#   ifdef MPI
      endTime = MPI_WTIME()
      adsTime = endTime-startTime
      call par_getMaxReal(adsCommTime)
#   else
      call system_clock(clock2, clock_rate)
      adsTime= REAL(clock2-clock1)/REAL(clock_rate)
#   endif

  IF (monitorON) WRITE(STDOUT,'(5X,A,F8.3,A,F6.1,A)') 'Ad-Diff. Eq. Time: ', adsTime, ' sec. (',100.*adsCommTime/adsTime,'% Comm.)'

   CALL GCM_enforce_global_mass_consv()   ! added by JHSeo

!   -----------------------
    IF(monitorON) THEN
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(3X,A)') 'Calculating Face Velocities '
    END IF
    CALL face_vel()!.............................................COMPLETE(SAMK)

!   Compute RHS for the Poisson Pressure Equation (PPE)
!   ---------------------------------------------------
    IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Forming RHS for Poisson Eq.'
    CALL rhs_poisson(sum)!.......................................COMPLETE(SAMK)

    IF(monitorOn) THEN
      WRITE(STDOUT,'(5X,A,1X,1PE21.12)') 'Sum Of Poisson RHS =',sum
    END IF

!   Solve the Poisson Pressure Equation (PPE)
!   -----------------------------------------
#   ifdef MPI
      startTime = MPI_WTIME()
#   else
      CALL system_clock(clock1)
#   endif

    IF (monitorON) THEN
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(3X,A)') 'Solving the Poisson Eq.'
    END IF
    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    POSSION_INDEX = index_pressure
    CALL solve_poisson(ppsCommTime)


    IF(VDVactive == 1 .AND. MOD(ntime,n_visual) == 0)THEN 
       POSSION_INDEX = index_scalar   
       CALL solve_laplace
       CALL calculate_scalar_gradien
    ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

#   ifdef MPI
      endTime = MPI_WTIME()
      ppsTime = endTime-startTime
      CALL par_getMaxReal(ppsCommTime)
#   else
      CALL system_clock(clock2, clock_rate)
      ppsTime= REAL(clock2-clock1)/REAL(clock_rate)
#   endif

    IF (monitorON) WRITE(STDOUT,'(5X,A,F8.3,A,F6.1,A)') 'Poisson Eq. Time: ', ppsTime, ' sec. (',100.*ppsCommTime/ppsTime,'% Comm.)'


!   Correct velocity field and update pressure
!   ------------------------------------------
    IF (monitorON) THEN
      WRITE(STDOUT,*)
      WRITE(STDOUT,'(3X,A)') 'Updating the velocity field'
    END IF
    CALL correct_vel()!..........................................COMPLETE(SAMK)

    IF (monitorON) THEN
      WRITE(STDOUT,'(3X,A)') 'Updating the pressure field'
    END IF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
	IF(fea_boundary_flag_tahoe == 1 .OR. implicit_coupling_flag_combined_motion == 1)THEN
	
	 CALL subtract_mean_pressure()!.... Added by Rajneesh
	
	ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
	
    CALL update_pressure()!......................................COMPLETE(SAMK)

!   Adjust velocity field for 2D computations
!   -----------------------------------------
    IF( nDim == DIM_2D) CALL vel_adjust2D !.......................COMPLETE(SAMK)
!DEBUG
! IF (monitorON) WRITE(STDOUT,'(3X,A)') 'Writing out dump file'
! CALL write_dump_debug()
! return
!DEBUG

     IF(iacou) then
       IF(monitorON) WRITE(STDOUT, *) 'Enter LPCE Solver'
	   CALL LPCE_solve
	ENDIF

!========  Solve TAHOE (Added by Rajneesh)======================= 

# ifdef TAHOE

CALL fea_body_marker_pressure !Gather interpolated pressures on the immersed bodies

IF(fea_boundary_flag_tahoe == 1 .AND. IMtheBOSS .AND. ntime>=nstart_tahoe) THEN
		   IF(implicit_coupling_flag == 1)THEN
               DO iBody = 1,nbody 
			   SELECT CASE (boundary_motion_type(iBody))
			   CASE(FEA_TAHOE) 
			   CALL update_tahoe_variables(iBody) !old to new
			   END SELECT 
		       END DO
		   END IF 
                         	
            DO iBody=1,nBody	
               SELECT CASE (boundary_motion_type(iBody))
			   CASE(FEA_TAHOE) 
			   IF(monitorON) WRITE(STDOUT, *) 'Enter TAHOE at ntime',ntime, &
			                       'and kimplicit',kimplicit,'for body',ibody
			   IF(dtratio_tahoe>1 .AND. ntime==nstart_tahoe)iflag_restart_tahoe=0
					 DO itime_tahoe = 1, dtratio_tahoe 		 				   
					 CALL tahoe_solver(iBody)
					 END DO
			   END SELECT 
            END DO
			
END IF !fea_boundary_flag_tahoe

IF(implicit_coupling_flag == 1 .AND. ntime>=nstart_tahoe)THEN	
   CALL convergence_implicit
	IF(ImplicitResidual > ImplicitResidualMax)THEN
	   IF(ImtheBoss) CALL WRITE_CONVERGENCE_DATA   
		IF(kimplicit < kimplicitMax)THEN	
		GO TO 256   !going back to implicit loop
		ELSE
		IF(monitorON) WRITE(STDOUT, *) 'Code did not converge in implicit scheme'
   		STOP
		END IF
    ELSE
	IF(ImtheBoss)CALL WRITE_CONVERGENCE_DATA
	END IF !ImplicitResidual 
END IF !implicit_coupling_flag 
 
# endif 

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    IF(PFactive == 1 .AND. MOD(ntime,n_potential) == 0)THEN     
       POSSION_INDEX = index_potential          
       CALL solve_laplace_potential        
       CALL calculate_potential_vorticity_U      
    ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 


!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
IF(derivative_flag == 1 .AND. ntime > nstart_FIM)THEN

   CALL FIM_CALC_DERIV()
   IF(k_derivative /= 4)  GOTO 259
    
ENDIF !derivative_flag == 1 .AND. ntime > nstart_FIM
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
!===Loop for implicit coupling of combined motion(END) 
IF(implicit_coupling_flag_combined_motion == 1 .AND. ntime > nstart_FIM)THEN

   IF(iDragLift == 1) THEN
     SELECT CASE(boundary_formulation)
     CASE (SSM_METHOD)
       CALL SSM_drag_lift()
     CASE (GCM_METHOD)
       CALL GCM_drag_lift() 
     END SELECT ! boundary_formulation
   ENDIF ! iDragLift  
	

   CALL move_boundary()   
   CALL FIM_convergence_implicit
!---Update geometrical information for all bodies  
   DO iBody = 1,nBody
     CALL MofI_CofG_par(iBody)
   ENDDO   
   
!---Update geometrical information for all combined groups  
IF(nGroup_Combined /= 0)THEN 
   DO iGroup = 1,nGroup_Combined
     CALL MofI_CofG_combined(iGroup) 
   ENDDO
ENDIF
      
	IF(ImplicitResidual_FIM > ImplicitResidualMax_FIM)THEN	
	   IF(ImtheBoss) CALL FIM_WRITE_CONVERGENCE_DATA   
		IF(kimplicit_FIM < kimplicitMax_FIM)THEN	
		GO TO 257   !going back to implicit loop	
		ELSE
		IF(monitorON) WRITE(imconverg, *) 'Code did not converge in FIM implicit scheme'
   		!STOP
		END IF
    ELSE
	IF(ImtheBoss)CALL FIM_WRITE_CONVERGENCE_DATA
	END IF !ImplicitResidual_FIM 
	
	
END IF !implicit_coupling_flag for combined motion
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
    
   
!=========================Solve FEA===========================================

    IF(fea_boundary_flag == 1) THEN                        ! SER_TO_PAR. QX. CH14
 
     IF(monitorON) WRITE(STDOUT, *) 'Enter Finite Element Structure Solver'
!compute the pressure value at marker points      
       CALL fea_body_marker_pressure

!Here in order to compaitable with multicore cpu cluster, only the master solve the FEA and then brocasting the results to save the memroy      
      IF(IMtheBOSS) THEN 
       DO fea_itime = 1, fea_dtratio
        DO ibody = 1, fea_nbody
        CALL fea_dynamicsolver(ibody, ntime)
        END DO !ibody
       ENDDO

       IF ( monitorON .AND. MOD(ntime,nmonitor_probe_liftdrag) == 0 ) THEN
        write (STDOUT, *) 'Write fea_probe_out'
        CALL fea_write_probe_files()
       ENDIF
     ENDIF
     
     CALL mpi_barrier(MPI_COMM_WORLD, ierr)
     
     IF(fea_boundary_flag == 1) CALL output_flux()                                      ! SER_TO_PAR. QX. CH15

    END IF

!===============================DYE========================================
    IF(is_scalar_on == SCALAR_ON) THEN 
     CALL  rhs_advec_diff_T
     CALL  solve_dif_trans
    ENDIF
!============================FEA===================================================
!   Monitor output
!   --------------
    IF (monitorIT) THEN
      CALL write_monitor()!......................................COMPLETE(SAMK)
     ! IF (turbActive == ACTIVE) THEN
     !    if (MonitorON) CALL TURB_write_monitor()
     !  ENDIF
    ENDIF
    


    IF ( nStat > STATS_NONE .AND. MOD(ntime,nstat) == 0 .AND. ntime > 1) THEN
        CALL calc_statistics!.................................COMPLETE(SAMK)
!         CALL calc_statistics_vorticity(0)
        statCtr = statCtr + 1
!        statCtrv = statCtrv + 1
    ENDIF   
    
    
     
    IF ( nmonitor_probe_liftdrag > STATS_NONE .AND. &
         MOD(ntime,nmonitor_probe_liftdrag) == 0 )  THEN



      CALL write_probe_files()!..................................COMPLETE(SAMK)


      IF(iDragLift == 1) THEN     
      SELECT CASE(boundary_formulation)
      CASE (SSM_METHOD)
        CALL SSM_drag_lift()
      CASE (GCM_METHOD)           
       CALL GCM_drag_lift()            
      END SELECT ! boundary_formulation       
     
    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
      IF(VDVactive == 1 .AND. MOD(ntime,n_visual) == 0)     CALL surface_kinetic_energy_difference
      IF(VDVactive == 1 .AND. MOD(ntime,n_visual) == 0)     CALL Gradien_KineticEnergy
      IF(PFactive  == 1 .AND. MOD(ntime,n_potential) == 0)  CALL Gradien_potential_KineticEnergy
    
      DO iBody = 1,nBody
        CALL CALC_FORCE_MOMENT(iBody)
      ENDDO
   
      DO iGroup = 1,nGroup_Combined
        CALL CALC_FORCE_MOMENT_combined(iGroup) 
      ENDDO     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)	  
	  ENDIF ! iDragLift
    END IF

     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
!dddddddddddddddddddddddddddddddddddddddddd
    CALL Calculate_gradU2
    CALL Calculate_unsteady   
    CALL Calculate_gradP
    CALL outer_surface_unsteady
!dddddddddddddddddddddddddddddddddddddddddd

    IF (mod(ntime, nFIM_dump) == 0) THEN
      IF(VDVactive == 1 .AND. MOD(ntime,n_visual) == 0)  call vorticity()
      IF(Imtheboss) CALL write_FIM_dump()
    ENDIF     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
  
    IF (dumpIT) THEN
     
      CALL write_dump()
             
      IF(iacou) call outdata_A(ntime)   ! LPCE output

      IF(fea_boundary_flag == 1 .AND. monitorON) THEN                      ! SER_TO_PAR. QX. CH16
         WRITE(STDOUT,*) 'Output Finite Element Solution'
         IF(ImtheBOSS) CALL fea_dynamicsolutionoutput(fea_nbody, ntime)                  ! added for FEA
      END IF
      
      CALL mpi_barrier(MPI_COMM_WORLD, ierr)     
    
   ENDIF
 

 
    IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN
     CALL write_restart()

     IF ( monitorON .and. fea_boundary_flag == 1) THEN
       IF(ImtheBOSS)  CALL fea_writerestart()                               ! added for FEA ! SER_TO_PAR. QX. CH17
     ENDIF

      CALL mpi_barrier(MPI_COMM_WORLD, ierr)
     IF ( nStat > STATS_NONE ) CALL write_statistics
!!!   IF ( nStat > STATS_NONE ) CALL calc_statistics_vorticity(1)
    ENDIF
   

#   ifdef MPI
      ThisTimeEnd=MPI_WTIME()
      elapsedTime=ThisTimeEnd-MainTimeStart
#   else
      CALL system_clock(ThisTimeEnd, clock_rate)
      elapsedTime=REAL(ThisTimeEnd-MainTimeStart)/REAL(clock_rate)
#   endif


!Added by Rajneesh
# ifdef TAHOE
     IF (fea_boundary_flag_tahoe == 1) THEN
     CALL output_flux()
     IF(ImtheBOSS.AND.iTahoeProbe==1) CALL WRITE_DATA_FILES_TAHOE
     IF(ImtheBOSS.AND.iTahoeOutput==1) CALL WRITE_OUTPUT_TAHOE
	 IF(ImtheBOSS .AND. iTahoeRestart==1)  CALL write_restart_files_tahoe
    END IF
# endif  
   

! Here to have a clean exit, when the wall time is about to be ran out, we force the code to write the restart file and exit.
    IF(iSAFE_EXIT.eq.1 .and. system_wall_time - elapsedTime < exit_safe_time) THEN
     CALL write_restart()
     IF ( fea_boundary_flag == 1) THEN
       IF(ImtheBOSS)  CALL fea_writerestart()                               ! added for FEA ! SER_TO_PAR. QX. CH17
     ENDIF

     CALL mpi_barrier(MPI_COMM_WORLD, ierr)
     IF ( nStat > STATS_NONE ) CALL write_statistics
      
     IF(ImtheBOSS) THEN
       WRITE(STDOUT,*)'TOTAL TIME IS', elapsedTime
       psec=mod(NINT(elapsedTime),60)
       elapsedTime=(NINT(elapsedTime)-psec)/60
       pmts=mod(NINT(elapsedTime),60)
       elapsedTime=(NINT(elapsedTime)-pmts)/60
       phrs=mod(NINT(elapsedTime),24)
       elapsedTime=(NINT(elapsedTime)-phrs)/24
       pdys=mod(NINT(elapsedTime),7)
       elapsedTime=(NINT(elapsedTime)-pdys)/7
       pwks=NINT(elapsedTime)
       WRITE(STDOUT,'(3X,A)') '--------------------------------------------------------------------'
       WRITE(STDOUT,'(3X,A,I2,A,I1,A,I2,A,I2,A,I2,A)') 'Elapsed Time: ',pwks,' weeks, ',pdys,' days, ',phrs,' hours, ',pmts,' minutes and ',psec,' seconds.'
       WRITE(STDOUT,'(3X,A)') '--------------------------------------------------------------------'
     ENDIF 	
     exit_safe = .TRUE.
     RETURN
    ENDIF
    
    IF (monitorON ) THEN
      WRITE(STDOUT,*)'TOTAL TIME IS', elapsedTime
      psec=mod(NINT(elapsedTime),60)
      elapsedTime=(NINT(elapsedTime)-psec)/60
      pmts=mod(NINT(elapsedTime),60)
      elapsedTime=(NINT(elapsedTime)-pmts)/60
      phrs=mod(NINT(elapsedTime),24)
      elapsedTime=(NINT(elapsedTime)-phrs)/24
      pdys=mod(NINT(elapsedTime),7)
      elapsedTime=(NINT(elapsedTime)-pdys)/7
      pwks=NINT(elapsedTime)
    
      WRITE(STDOUT,'(3X,A)') '--------------------------------------------------------------------'
      WRITE(STDOUT,'(3X,A,I2,A,I1,A,I2,A,I2,A,I2,A)') 'Elapsed Time: ',pwks,' weeks, ',pdys,' days, ',phrs,' hours, ',pmts,' minutes and ',psec,' seconds.'
      WRITE(STDOUT,'(3X,A)') '--------------------------------------------------------------------'
    END IF

  ENDDO ! ntime
       
END SUBROUTINE time_step_viscous
