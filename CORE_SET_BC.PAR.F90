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
!  Filename: CORE_SET_BC.PAR.F90
!  Latest Modification: Jan 20, 2008 (ver. P1.5.5)
!  Made by X. Zheng
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     set_bc()
!     SSM_enforce_global_mass_consv()
!     remove_up_um()
!     set_pressure_dirichlet_bc()
!     update_pressure_dirichlet_bc()
!     enforce_p_periodic(pres)
!     add_up_um
!     unset_pressure_dirichlet_bc()
!     enforce_u_periodic()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE set_bc()

! -------------------------------------------------
!  Tis subroutine set the boundary conditions.
!
!  bcflag = 1  => dirichlet bc
!  bcflag = 2  => neumann bc
!  bcflag = 3  => pulsatile bc
!
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |       
!     |   o   |   o   |   o   |   o   |   o   
!     |       |  bcy  |  bcy  |       |       
!     |-------|---+---|-------|-------|-------
!     |       |*******|*******|       |       
!     |   obcx+*******|*******|bcxo   |   o   
!     |       |*******|*******|       |       
!     |-------|-------|-------|-------|-------
!
!
! -------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays


    IMPLICIT NONE

    INTEGER :: i,j,k
    
    REAL (KIND=CGREAL) :: Ly, Lz, yp, zp



!   Internal Boundary Conditions
!   Need to set iup for GCM to satisfy mass conservation on Stair-step boundary
!   ---------------------------------------------------------------------------
    IF ( boundary_formulation /= NO_INTERNAL_BOUNDARY ) THEN
      IF (boundary_formulation == SSM_METHOD) THEN
        CALL SSM_set_bc_internal()
      ELSE IF (boundary_formulation == GCM_METHOD) THEN
        CALL GCM_SetBodyInterceptValues()
        CALL GCM_SetBodyInterceptValuesFresh()
      ENDIF
    ENDIF

!   Outer boundary conditions
!   -------------------------
    IF (bcx1==BC_TYPE_USER_SPECIFIED) THEN
      IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
      CALL flow_stop
      STOP

!!!      CALL blasius_velocity
    END IF

!   left boundary
!   -------------
    i = 1

    IF (myCoords(1)==0) THEN          ! Computational Domain Boundary
      DO k=1,nzc
      DO j=1,nyc

        SELECT CASE (bcx1)

        CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc
          bcxu(i,j,k) = ux1 
          bcxv(i,j,k) = vx1 
          bcxw(i,j,k) = wx1
        CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)
          bcxu(i,j,k) = u(i,j,k)
          bcxv(i,j,k) = v(i,j,k)
          bcxw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcxu(i,j,k) = ux1*sin(2.0_CGREAL*pi*freq_ux1*time)
          bcxv(i,j,k) = vx1*sin(2.0_CGREAL*pi*freq_vx1*time)
          bcxw(i,j,k) = wx1*sin(2.0_CGREAL*pi*freq_wx1*time)            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
          bcxu(i,j,k) = 0.0_CGREAL
          bcxv(i,j,k) = v(i,j,k)
          bcxw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PERIODIC)       ! periodic bc  (explicit & dirty implementation)
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled'
          CALL flow_stop
          STOP

!!!          bcxu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(nx-1,j,k) )
!!!          bcxv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(nx-1,j,k) )
!!!          bcxw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(nx-1,j,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          IF (nDim==DIM_3D) THEN
            Ly=yout/2.0_CGREAL
            Lz=Zout/2.0_CGREAL

            yp=(yc(L2GJ(j))-Ly)/Ly
            zp=(zc(     k )-Lz)/Lz

            bcxu(i,j,k)=2.25d0*ux1*(1.0_CGREAL-(yp*yp))*(1.0_CGREAL-(zp*zp))
          ELSE
            Ly=yout/2.0_CGREAL

            yp=(yc(L2GJ(j))-Ly)/Ly

            bcxu(i,j,k)=1.5d0*ux1*(1.0_CGREAL-(yp*yp))
          END IF
          bcxv(i,j,k) = vx1 
          bcxw(i,j,k) = wx1
        CASE (BC_TYPE_SHEAR)       ! shear bc 
          bcxu(i,j,k) = ux1*y(L2GJ(j))/yout
          bcxv(i,j,k) = vx1 
          bcxw(i,j,k) = wx1
        END SELECT
      END DO
      END DO
    ELSE                          ! Sub-Domain Boundary
!       No need to setup subdomain BC
    END IF

!   right boundary
!   --------------
    i = nxc

    IF (myCoords(1)==Np(1)-1) THEN          ! Computational Domain Boundary
      DO k=1,nzc
      DO j=1,nyc

        SELECT CASE (bcx2)

        CASE (BC_TYPE_DIRICHLET)      ! dirichlet bc
          bcxu(i,j,k) = ux2 
          bcxv(i,j,k) = vx2 
          bcxw(i,j,k) = wx2
        CASE (BC_TYPE_ZERO_GRADIENT)  ! outflow bc ( zero gradient ; explicit)
          bcxu(i,j,k) = u(i,j,k)
          bcxv(i,j,k) = v(i,j,k)
          bcxw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcxu(i,j,k) = ux2*sin(2.0_CGREAL*pi*freq_ux2*time)
          bcxv(i,j,k) = vx2*sin(2.0_CGREAL*pi*freq_vx2*time)
          bcxw(i,j,k) = wx2*sin(2.0_CGREAL*pi*freq_wx2*time)            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
          bcxu(i,j,k) = 0.0_CGREAL
          bcxv(i,j,k) = v(i,j,k)
          bcxw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PERIODIC)       ! periodic bc
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled' 
          CALL flow_stop
          STOP

!!!        bcxu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(1,j,k) )
!!!        bcxv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(1,j,k) )
!!!        bcxw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(1,j,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
          CALL flow_stop
          STOP

!!!        CALL blasius_velocity
!!!        ux2 = u_blasius(j)
!!!        bcxu(i,j,k) = ux2
!!!        bcxv(i,j,k) = vx2
!!!        bcxw(i,j,k) = wx2
        CASE (BC_TYPE_SHEAR)          ! shear bc
          bcxu(i,j,k) = ux2*y(L2GJ(j))/yout 
          bcxv(i,j,k) = vx2 
          bcxw(i,j,k) = wx2
        END SELECT 
      ENDDO ! j
      ENDDO ! k
    ELSE                              ! Sub-Domain Boundary
!       No need to setup subdomain BC
    END IF

!   bottom boundary
!   ---------------
    j = 1

    IF (myCoords(2)==0) THEN
      DO k=1,nzc
      DO i=1,nxc

        SELECT CASE (bcy1)

        CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
          bcyu(i,j,k) = uy1 
          bcyv(i,j,k) = vy1 
          bcyw(i,j,k) = wy1
        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
          bcyu(i,j,k) = u(i,j,k)
          bcyv(i,j,k) = v(i,j,k)
          bcyw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcyu(i,j,k) = uy1*sin(2.0_CGREAL*pi*freq_uy1*time)
          bcyv(i,j,k) = vy1*sin(2.0_CGREAL*pi*freq_vy1*time)
          bcyw(i,j,k) = wy1*sin(2.0_CGREAL*pi*freq_wy1*time)            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
          bcyu(i,j,k) = u(i,j,k)
          bcyv(i,j,k) = 0.0_CGREAL
          bcyw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled' 
          CALL flow_stop
          STOP

!!!        bcyu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(i,ny-1,k) )
!!!        bcyv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(i,ny-1,k) )
!!!        bcyw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(i,ny-1,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
          CALL flow_stop
          STOP

!!!        CALL blasius_velocity
!!!        uy1 = u_blasius(j)
!!!        bcyu(i,j,k) = uy1
!!!        bcyv(i,j,k) = vy1
!!!        bcyw(i,j,k) = wy1
        CASE (BC_TYPE_SHEAR)                ! shear bc
          bcyu(i,j,k) = uy1
          bcyv(i,j,k) = vy1*x(L2GI(i))/xout
          bcyw(i,j,k) = wy1
        END SELECT
      END DO
      END DO
    ELSE
!       No need to setup subdomain BC
    END IF

!   top boundary
!   ------------
    j = nyc

    IF (myCoords(2)==Np(2)-1) THEN
      DO k=1,nzc
      DO i=1,nxc

        SELECT CASE (bcy2)

        CASE (BC_TYPE_DIRICHLET)             ! dirichlet bc
          bcyu(i,j,k) = uy2 
          bcyv(i,j,k) = vy2 
          bcyw(i,j,k) = wy2
        CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
          bcyu(i,j,k) = u(i,j,k)
          bcyv(i,j,k) = v(i,j,k)
          bcyw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PULSATILE_INFLOW)
          bcyu(i,j,k) = uy2*sin(2.0_CGREAL*pi*freq_uy2*time)
          bcyv(i,j,k) = vy2*sin(2.0_CGREAL*pi*freq_vy2*time)
          bcyw(i,j,k) = wy2*sin(2.0_CGREAL*pi*freq_wy2*time)            
        CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
          bcyu(i,j,k) = u(i,j,k)
          bcyv(i,j,k) = 0.0_CGREAL
          bcyw(i,j,k) = w(i,j,k)
        CASE (BC_TYPE_PERIODIC)       ! periodic bc 
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled' 
          CALL flow_stop
          STOP

!!!        bcyu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(i,1,k) )
!!!        bcyv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(i,1,k) )
!!!        bcyw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(i,1,k) )
        CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
          IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
          CALL flow_stop
          STOP

!!!        CALL blasius_velocity
!!!        uy2 = u_blasius(j)
!!!        bcyu(i,j,k) = uy2
!!!        bcyv(i,j,k) = vy2
!!!        bcyw(i,j,k) = wy2
        CASE (BC_TYPE_SHEAR)             ! shear bc
          bcyu(i,j,k) = uy2
          bcyv(i,j,k) = vy2*x(L2GI(i))/xout 
          bcyw(i,j,k) = wy2
!dddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!Zero Vorticity BC
        CASE(18)
          bcyu(i,j,k) = u(i,j,k)
         IF( xc(L2GI(i))>=0.4464  .AND. xc(L2GI(i))<=0.6474  )THEN 
          bcyv(i,j,k) = 0.8_CGREAL*SIN( 2.0_CGREAL*pi*( xc(L2GI(i))-0.4464 )/0.201_CGREAL )          
         ELSE
          bcyv(i,j,k) = 0.0_CGREAL
         ENDIF  
          bcyw(i,j,k) = w(i,j,k)            
!dddddddddddddddddddddddddddddddddddddddddddddddddddddddd          
        END SELECT
      ENDDO ! i
      ENDDO ! k
    ELSE
!       No need to setup subdomain BC
    END IF

    DO j=1,nyc
    DO i=1,nxc

!     front boundary
!     --------------
      k = 1

      SELECT CASE (bcz1)

      CASE (BC_TYPE_DIRICHLET)             ! diriclet bc
        bczu(i,j,k) = uz1 
        bczv(i,j,k) = vz1 
        bczw(i,j,k) = wz1
      CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
        bczu(i,j,k) = u(i,j,k)
        bczv(i,j,k) = v(i,j,k)
        bczw(i,j,k) = w(i,j,k)
      CASE (BC_TYPE_PULSATILE_INFLOW)
        bczu(i,j,k) = uz1*sin(2.0_CGREAL*pi*freq_uz1*time)
        bczv(i,j,k) = vz1*sin(2.0_CGREAL*pi*freq_vz1*time)
        bczw(i,j,k) = wz1*sin(2.0_CGREAL*pi*freq_wz1*time)            
      CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
        bczu(i,j,k) = u(i,j,k)
        bczv(i,j,k) = v(i,j,k)
        bczw(i,j,k) = 0.0
      CASE (BC_TYPE_PERIODIC)       ! periodic bc 
!        IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled'        !Ehsan removed for Periodic BC 
!        CALL flow_stop                                                            !Ehsan removed for Periodic BC
!        STOP                                                                      !Ehsan removed for Periodic BC

        bczu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(i,j,nzc) )                        !Ehsan added for Periodic BC
        bczv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(i,j,nzc) )                        !Ehsan added for Periodic BC
        bczw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(i,j,nzc) )                        !Ehsan added for Periodic BC

      CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
        IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
        CALL flow_stop
        STOP

!!!        CALL blasius_velocity
!!!        uz1 = u_blasius(j)
!!!        bczu(i,j,k) = uz1
!!!        bczv(i,j,k) = vz1
!!!        bczw(i,j,k) = wz1
      CASE (BC_TYPE_SHEAR)          ! shear bc
        bczu(i,j,k) = uz1
        bczv(i,j,k) = vz1 
        bczw(i,j,k) = wz1*y(L2GJ(j))/yout
      END SELECT 

!     back boundary
!     -------------
      k = nzc

      SELECT CASE (bcz2)

      CASE (BC_TYPE_DIRICHLET)             ! diriclet bc
        bczu(i,j,k) = uz2 
        bczv(i,j,k) = vz2 
        bczw(i,j,k) = wz2
      CASE (BC_TYPE_ZERO_GRADIENT)         ! outflow bc ( zero gradient ; explicit)
        bczu(i,j,k) = u(i,j,k)
        bczv(i,j,k) = v(i,j,k)
        bczw(i,j,k) = w(i,j,k)
      CASE (BC_TYPE_PULSATILE_INFLOW)
        bczu(i,j,k) = uz2*sin(2.0_CGREAL*pi*freq_uz2*time)
        bczv(i,j,k) = vz2*sin(2.0_CGREAL*pi*freq_vz2*time)
        bczw(i,j,k) = wz2*sin(2.0_CGREAL*pi*freq_wz2*time)            
      CASE (BC_TYPE_SYMMETRY)       ! symmery bc (explicit)
        bczu(i,j,k) = u(i,j,k)
        bczv(i,j,k) = v(i,j,k)
        bczw(i,j,k) = 0.0_CGREAL
      CASE (BC_TYPE_PERIODIC)       ! periodic bc 
!        IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_PERIODIC is currently disabled'              !Ehsan removed for Periodic BC 
!        CALL flow_stop                                                                  !Ehsan removed for Periodic BC
!        STOP                                                                            !Ehsan removed for Periodic BC

        bczu(i,j,k) = 0.5_CGREAL*( u(i,j,k) + u(i,j,1) )                                 !Ehsan added for Periodic BC
        bczv(i,j,k) = 0.5_CGREAL*( v(i,j,k) + v(i,j,1) )                                 !Ehsan added for Periodic BC
        bczw(i,j,k) = 0.5_CGREAL*( w(i,j,k) + w(i,j,1) )                                 !Ehsan added for Periodic BC

      CASE (BC_TYPE_USER_SPECIFIED)       !  user specified
        IF (ImtheBOSS) WRITE(*,*) 'BC_TYPE_USER_SPECIFIED is currently disabled' 
        CALL flow_stop
        STOP

!!!        CALL blasius_velocity
!!!        uz2 = u_blasius(j)
!!!        bczu(i,j,k) = uz2
!!!        bczv(i,j,k) = vz2
!!!        bczw(i,j,k) = wz2
      CASE (BC_TYPE_SHEAR)          ! shear bc
        bczu(i,j,k) = uz2
        bczv(i,j,k) = vz2 
        bczw(i,j,k) = wz2*y(L2GJ(j))/yout
      END SELECT 

    ENDDO ! i
    ENDDO ! j

!   adjust BC to satisfy global mass conservation for SSM
!   -----------------------------------------------------
    IF (boundary_formulation == NO_INTERNAL_BOUNDARY .OR. &
        boundary_formulation == SSM_METHOD                 ) CALL SSM_enforce_global_mass_consv()

!   set values for outer ghost points
!   ---------------------------------
    CALL set_outer_ghost_vel()

END SUBROUTINE set_bc
!---------------------------------------------------------------------



SUBROUTINE SSM_enforce_global_mass_consv()

! -------------------------------------------------------------------------
!  This subroutine computes mass flux at all boundaries and adjust outflow 
!  BC so as to satisfy global mass conservation.
! -------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays

    IMPLICIT NONE

    INTEGER  :: i ,j ,k
    INTEGER  :: iG,jG,n

    REAL(KIND=CGREAL) :: massflux,correction_vel



!   Added by SAMK: To calculate left wall flow rate
!   -----------------------------------------------
    massflux = 0.0_CGREAL

    IF (myCoords(1)==0) THEN
      DO k=1,nzc
      DO j=1,nyc
        massflux = massflux + bcxu(1,j,k)*ium(1,j,k)*dy(j)*dz(k) &
                 * REAL(1-iblank(1,j,k),KIND=CGREAL)
      ENDDO ! j
      ENDDO ! k
    END IF

#   ifdef MPI
      CALL par_getSumReal(massflux)
#   endif

    IF (monitorON) THEN
      WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: inflow rate              = ',massflux
      WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: inflow area              = ',inflow_area
      WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: inflow perimeter         = ',prim_left
      WRITE(STDOUT,*)
      !WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: charactersitic velocity  = ',massflux/inflow_area
      !WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: hydraulic diameter       = ',4.d0*inflow_area/prim_left
      !WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC: physical Reynolds number = ',4.d0*massflux*Re/prim_left
      WRITE(STDOUT,*)
      END IF



    massflux = 0.0_CGREAL

    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

      DO k=1,nzc
      DO j=1,nyc
        jG=L2GJ(j)
        DO i=1,nxc
          iG=L2GI(i)
        
          massflux = massflux +                              &
                   ( -bcxu(i,j,k)*ium(i,j,k)*dy(jG)*dz(k)    &
                     +bcxu(i,j,k)*iup(i,j,k)*dy(jG)*dz(k)    &
                     -bcyv(i,j,k)*jum(i,j,k)*dx(iG)*dz(k)    &
                     +bcyv(i,j,k)*jup(i,j,k)*dx(iG)*dz(k)    &
                     -bczw(i,j,k)*kum(i,j,k)*dx(iG)*dy(jG)   &
                     +bczw(i,j,k)*kup(i,j,k)*dx(iG)*dy(jG) ) &
                  *REAL(1-iblank(i,j,k),KIND=CGREAL)
        ENDDO ! i
      ENDDO ! j
      ENDDO ! k

#     ifdef MPI
        CALL par_getSumReal(massflux)
#     endif

      correction_vel =-massflux/outflow_area  

      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC:massflux       = ',massflux
        WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC:outflow_area   = ',outflow_area
        WRITE(STDOUT,'(5X,A,1PE12.5)') 'SET_BC:correction_vel = ',correction_vel
      END IF ! ntime

      IF (bcx1 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(1)==0) &
        bcxu(1,:,:) = bcxu(1,:,:) - correction_vel

      IF (bcx2 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(1)==Np(1)-1) &
        bcxu(nxc,:,:) = bcxu(nxc,:,:) + correction_vel

      IF (bcy1 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(2)==0) &
        bcyv(:,1,:) = bcyv(:,1,:) - correction_vel
      
      IF (bcy2 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(2)==Np(2)-1) &
        bcyv(:,nyc,:) = bcyv(:,nyc,:) + correction_vel

      IF (bcz1 == BC_TYPE_ZERO_GRADIENT) bczw(:,:,1)   = bczw(:,:,1)   - correction_vel
      IF (bcz2 == BC_TYPE_ZERO_GRADIENT) bczw(:,:,nzc) = bczw(:,:,nzc) + correction_vel

    ENDIF ! bcx1

END SUBROUTINE SSM_enforce_global_mass_consv
!---------------------------------------------------------------------

SUBROUTINE remove_up_um                                                 !Ehsan added for periodic bc

! ----------------------------------------------------------------------------------
!  Purpose: Removing IUP, IUM, JUP, JUM, KUP, KUM for pressure periodic conditions.
!  Written by H. Dong
! ----------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

    IF (bcz1 == BC_TYPE_PERIODIC .AND. &
      bcz2 == BC_TYPE_PERIODIC) THEN
        kum(:,:,1) = 0.0_CGREAL
        kup(:,:,nzc) = 0.0_CGREAL
    END IF

END SUBROUTINE remove_up_um
!---------------------------------------------------------------------

SUBROUTINE add_up_um                                            !Ehsan added for periodic bc

! -------------------------------------------------------------------------------------
!  Purpose: Adding IUP, IUM, JUP, JUM, KUP, KUM back for pressure periodic conditions.
!
!  Input: No.
!
!  Output: No.
!  Written by H. Dong
! -------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k,m

     IF (bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC) THEN
         DO i=1,nxc
         DO j=1,nyc
           kum(i,j,1)    = 1.0_CGREAL
           kup(i,j,nzc) = 1.0_CGREAL
         ENDDO
         ENDDO
     END IF

END SUBROUTINE add_up_um
!---------------------------------------------------------------------

SUBROUTINE set_pressure_dirichlet_bc
! ---------------------------------------------------
!  Purpose: Setting IUP, IUM, JUP, JUM, KUP, KUM  
!  to be zero at the outer boundary for Dirichlet
!  pressure condition.                           
!                                               
!  Written by H. Luo
!
!  NOTE: The boundary condition at the subdomain
!        borders is either Dirichlet or Neumann
!        REGARDLESS of the outer boundary condition.
!        It depends on the neighbour cell. (SAMK)
! ---------------------------------------------------
                
    USE global_parameters
    USE flow_parameters
    USE pressure_arrays
    USE boundary_arrays 
    
    IMPLICIT NONE  
   
       

    IF (pbcx1 == PBC_DIRICHLET .AND. myCoords(1)==0      ) THEN
      ium(1  ,1:nyc,1:nzc) = 0.0_CGREAL
    END IF
 
    IF (pbcx2 == PBC_DIRICHLET .AND. myCoords(1)==Np(1)-1) THEN
      iup(nxc,1:nyc,1:nzc) = 0.0_CGREAL
    END IF
 
    IF (pbcy1 == PBC_DIRICHLET .AND. myCoords(2)==0      ) THEN
      jum(1:nxc,1  ,1:nzc) = 0.0_CGREAL
    END IF
            
    IF (pbcy2 == PBC_DIRICHLET .AND. myCoords(2)==Np(2)-1) THEN
      jup(1:nxc,nyc,1:nzc) = 0.0_CGREAL
    END IF
            
    IF (pbcz1 == PBC_DIRICHLET)                            THEN
      kum(1:nxc,1:nyc,1  ) = 0.0_CGREAL
    END IF
              
    IF (pbcz2 == PBC_DIRICHLET)                            THEN
      kup(1:nxc,1:nyc,nzc) = 0.0_CGREAL
    END IF

END SUBROUTINE set_pressure_dirichlet_bc
!---------------------------------------------------------------------
SUBROUTINE update_pressure_dirichlet_bc
! ---------------------------------------------------
!  Purpose: Setting IUP, IUM, JUP, JUM, KUP, KUM
!  to be zero at the outer boundary for Dirichlet
!  pressure condition.
!
!  Written by H. Luo
!
!  NOTE: The boundary condition at the subdomain
!        borders is either Dirichlet or Neumann
!        REGARDLESS of the outer boundary condition.
!        It depends on the neighbour cell. (SAMK)
! ---------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE pressure_arrays
    USE boundary_arrays

    IMPLICIT NONE



    IF (pbcx1 == PBC_DIRICHLET .AND. myCoords(1)==0      ) THEN
      pprime(0,:,:)= 2.0_CGREAL * pppx1 - pprime(1,:,:) 
    END IF

    IF (pbcx2 == PBC_DIRICHLET .AND. myCoords(1)==Np(1)-1) THEN
      pprime(nxc+1,:,:)=2.0_CGREAL * pppx2 - pprime(nxc,:,:)
    END IF

    IF (pbcy1 == PBC_DIRICHLET .AND. myCoords(2)==0      ) THEN
      pprime(:,0,:)=2.0_CGREAL * pppy1 - pprime(:,1,:)
    END IF

    IF (pbcy2 == PBC_DIRICHLET .AND. myCoords(2)==Np(2)-1) THEN
      pprime(:,nyc+1,:)=2.0_CGREAL * pppy2 - pprime(:,nyc, :)
    END IF

    IF (pbcz1 == PBC_DIRICHLET)                            THEN
      pprime(:,:,0)=2.0_CGREAL * pppz1 - pprime(:,:,1)
    END IF

    IF (pbcz2 == PBC_DIRICHLET)                            THEN
      pprime(:,:,nzc+1)=2.0_CGREAL * pppz2 - pprime(:,:,nzc)
    END IF

END SUBROUTINE update_pressure_dirichlet_bc
!---------------------------------------------------------------------



SUBROUTINE unset_pressure_dirichlet_bc

! ---------------------------------------------
!  Purpose: Reset IUP, IUM, JUP, JUM, KUP, KUM 
!  to be 1 at the out boundary                 
!                                              
!  Written by H. Luo                           
! ---------------------------------------------
 
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
 
    IMPLICIT NONE
          



    IF (pbcx1 == PBC_DIRICHLET .AND. myCoords(1)==0      ) &
        ium(1  ,1:nyc,1:nzc) = 1.0_CGREAL
 
    IF (pbcx2 == PBC_DIRICHLET .AND. myCoords(1)==Np(1)-1) &
        iup(nxc,1:nyc,1:nzc) = 1.0_CGREAL
 
    IF (pbcy1 == PBC_DIRICHLET .AND. myCoords(2)==0      ) &
        jum(1:nxc,1  ,1:nzc) = 1.0_CGREAL
            
    IF (pbcy2 == PBC_DIRICHLET .AND. myCoords(2)==Np(2)-1) &
        jup(1:nxc,nyc,1:nzc) = 1.0_CGREAL
            
    IF (pbcz1 == PBC_DIRICHLET)                            &
        kum(1:nxc,1:nyc,1  ) = 1.0_CGREAL
              
    IF (pbcz2 == PBC_DIRICHLET)                            &
        kup(1:nxc,1:nyc,nzc) = 1.0_CGREAL

END SUBROUTINE unset_pressure_dirichlet_bc
!------------------------------------------------------------------------------
