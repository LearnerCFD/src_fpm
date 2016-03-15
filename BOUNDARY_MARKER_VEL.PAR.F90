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
!     Jung Hee Seo
!     Rajneesh Bhardwaj
!
!  Filename: BOUNDARY_MARKER_VEL.PAR.F90
!  Latest Modification: October 17, 2010 (ver PAT 1.1.0)
!  By Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     wall_velocity(iBody)
!     forced_motion(iBody)
!     flow_induced_motion(iBody)
!     MofI_CofG(iBody)
! --------------------------------------------------------------------


# define L2GI(i)       myIs+i-1   ! Added by Rajneesh
# define L2GJ(j)       myJs+j-1 

SUBROUTINE wall_velocity(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER,INTENT (IN) :: iBody

    INTEGER             :: m, midx ,k ,num, mm                                                   !Ehsan added for SJ (num,mm)

    REAL(KIND=CGREAL)   :: uBodyMarkerRel, vBodyMarkerRel, wBodyMarkerRel, Vx_coeff, Vy_coeff    !Ehsan added for SJ (Vx_coeff, Vy_coeff)



    IF (boundary_motion /= MOVING_BOUNDARY) THEN      
      uBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL
      vBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL
      wBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL
    ENDIF

    OPEN(1217, FILE='vel_amp.dat', FORM='FORMATTED')        !Ehsan added for SJ
    READ(1217,*)num                                          !Ehsan added for SJ   

    DO mm = 1,num                                            !Ehsan added for SJ
      
      READ(1217,*)m, Vx_coeff, Vy_coeff                      !Ehsan added for SJ 

      uBodyMarkerRel = Vx_coeff*ampVelX(iBody)*COS(  2.0_CGREAL*PI*freqVelX(iBody)*time  &
                                          + phaseVelX(iBody)*pi/180.0_CGREAL  )

      vBodyMarkerRel = Vy_coeff*ampVelY(iBody)*COS(  2.0_CGREAL*PI*freqVelY(iBody)*time  &
                                          + phaseVelY(iBody)*pi/180.0_CGREAL  )
!dddddddddddddddddddddddddddddddddddddd
IF(IMTHEBOSS)THEN
OPEN(21212,FILE='UUUUU.DAT',POSITION='APPEND')
WRITE(21212,*) vBodyMarkerRel
ENDIF 
!dddddddddddddddddddddddddddddddddddddd
      wBodyMarkerRel = ampVelZ(iBody)*COS(  2.0_CGREAL*PI*freqVelZ(iBody)*time  &
                                          + phaseVelZ(iBody)*pi/180.0_CGREAL  )

      uBodyMarker(iBody,m) = uBodyMarker(iBody,m) +  uBodyMarkerRel
      vBodyMarker(iBody,m) = vBodyMarker(iBody,m) +  vBodyMarkerRel
      wBodyMarker(iBody,m) = wBodyMarker(iBody,m) +  wBodyMarkerRel
      print*,m,xBodyMarker(iBody,m),yBodyMarker(iBody,m)
      print*,zBodyMarker(iBody,m),uBodyMarker(iBody,m),vBodyMarker(iBody,m)
    ENDDO                                                 ! Rupeshs additions end here

    close(1217)                                             ! Rupeshs additions end here

!   Extending velocity across span for 2D body
!   ------------------------------------------
!    IF ( body_dim(iBody) == BODY_DIM2 .AND. ntime > 0 ) THEN

!      DO k=2,nz
!        DO m=mMinWallVel(iBody),mMaxWallVel(iBody)
!          midx = (k-1)*nPtsBodymarker(iBody)/nz + m
!          uBodyMarker(iBody,midx) = uBodyMarker(iBody,m)
!          vBodyMarker(iBody,midx) = vBodyMarker(iBody,m)
!          wBodyMarker(iBody,midx) = wBodyMarker(iBody,m)
!        ENDDO
!      ENDDO

!    ENDIF

END SUBROUTINE wall_velocity
!---------------------------------------------------------------------



SUBROUTINE forced_motion(bodyNumber)

! -------------------------------------------------------------------------------
!  This subroutine enforces general mean + sinusoidal component on translational
!  and angular velocity
! -------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
    USE boundary_arrays 
    USE usr_module 
    USE implicit_coupling_parameters  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)      

    IMPLICIT NONE

    INTEGER,INTENT (IN) ::bodyNumber 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER  :: i
!dddddddddddddddddddddddddddddddd
    CHARACTER*50 :: fname1  
    REAL(KIND=CGREAL)   :: theta1    
    REAL(KIND=CGREAL)   :: theta
    REAL(KIND=CGREAL)   :: xx,yy,zz    
    REAL(KIND=CGREAL)   :: R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33  
    REAL(KIND=CGREAL)   :: hingvx, hingvy, hingvz    
!dddddddddddddddddddddddddddddddd 
    REAL(KIND=CGREAL)   :: angvx_pre, angvy_pre, angvz_pre
    REAL(KIND=CGREAL)   :: angvx__old, angvy__old, angvz__old    
    REAL(KIND=CGREAL)   :: angvx_rel, angvy_rel, angvz_rel
    REAL(KIND=CGREAL)   :: temp_angvx,temp_angvy,temp_angvz 
    REAL(KIND=CGREAL)   :: vxcent_rel,vycent_rel,vzcent_rel 
    INTEGER            :: iGroup
    iGroup = combined_Group_index(bodyNumber)

   IF(combined_type(bodyNumber) == UNCOMBINED .OR. ntime <= nstart_FIM)THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    vxcent(bodyNumber) = vxcentTrans(bodyNumber)  &
                        + ampx(bodyNumber)*SIN(2.0_CGREAL*PI*freqx(bodyNumber)*time)
                                                                        
    vycent(bodyNumber) = vycentTrans(bodyNumber)  &
                        + ampy(bodyNumber)*SIN(2.0_CGREAL*PI*freqy(bodyNumber)*time)
    vzcent(bodyNumber) = vzcentTrans(bodyNumber)  &
                        + ampz(bodyNumber)*SIN(2.0_CGREAL*PI*freqz(bodyNumber)*time)


    angvx(bodyNumber) = angvxinit(bodyNumber)    &
                        + ampangx(bodyNumber)    &
                        *( SIN(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*cosphase(bodyNumber) &
                          +COS(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*sinphase(bodyNumber) )   
    angvy(bodyNumber) = angvyinit(bodyNumber)    &
                        + ampangy(bodyNumber)    &
                        *( SIN(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*cosphase(bodyNumber) &
                          +COS(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*sinphase(bodyNumber) )
    angvz(bodyNumber) = angvzinit(bodyNumber)    &
                        + ampangz(bodyNumber)    &
                        *( SIN(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*cosphase(bodyNumber) &
                          +COS(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*sinphase(bodyNumber) )
                          
                    
                          
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    angvx(bodyNumber) = 0.d0
    angvy(bodyNumber) = 0.d0  
    
  IF(pendulum_mother_body(bodyNumber) /= 0)THEN
      temp_angvx  = 0.5_CGREAL*(angvx_old(bodyNumber)+angvx(bodyNumber))
      temp_angvy  = 0.5_CGREAL*(angvy_old(bodyNumber)+angvy(bodyNumber))
      temp_angvz  = 0.5_CGREAL*(angvz_old(bodyNumber)+angvz(bodyNumber))     
      
      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      theta = dt * sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2) 
      xx    = temp_angvx/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)
      yy    = temp_angvy/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)
      zz    = temp_angvz/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)      
      
      R_11 = cos(theta) + (1.0_CGREAL-cos(theta))*xx**2.0_CGREAL
      R_12 = (1.0_CGREAL-cos(theta))*xx*yy - sin(theta)*zz
      R_13 = (1.0_CGREAL-cos(theta))*xx*zz + sin(theta)*yy
      R_21 = (1.0_CGREAL-cos(theta))*yy*xx + sin(theta)*zz
      R_22 = cos(theta) + (1.0_CGREAL-cos(theta))*yy**2.0_CGREAL 
      R_23 = (1.0_CGREAL-cos(theta))*yy*zz - sin(theta)*xx
      R_31 = (1.0_CGREAL-cos(theta))*zz*xx - sin(theta)*yy
      R_32 = (1.0_CGREAL-cos(theta))*zz*yy + sin(theta)*xx    
      R_33 = cos(theta) + (1.0_CGREAL-cos(theta))*zz**2.0_CGREAL  
      
      CALL Calc_Single_Marker_Vel(bodyNumber, theta,temp_angvx,temp_angvy,temp_angvz,                        &
                                   hinge_x(BodyNumber),  hinge_y(BodyNumber),  hinge_z(BodyNumber),      &
                                   hingvx, hingvy, hingvz,                                               &                                  
                                   R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33)
                                   
      vxcent(bodyNumber) = vxcent(bodyNumber) - hingvx
      vycent(bodyNumber) = vycent(bodyNumber) - hingvy
      vzcent(bodyNumber) = vzcent(bodyNumber) - hingvz
    ENDIF !pendulum_mother_body(bodyNumber)  
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddd

    angvx__old = angvx(bodyNumber)
    angvy__old = angvy(bodyNumber)
    angvz__old = angvz(bodyNumber)
    
    !IF(nread == 0)THEN
    !theta1 = angvx__old*dt*ntime
    !ELSEIF(nread == 1)THEN
    !theta1 = angvx__old*dt*(ntime-ntime_start)
    !ENDIF
    !angvy(bodyNumber)  = cos(theta1)*angvy__old - sin(theta1)*angvz__old 
    !angvz(bodyNumber)  = sin(theta1)*angvy__old + cos(theta1)*angvz__old 
    


   ELSE !===>combined_type(bodyNumber) == COMBINED .and. ntime > nstart_FIM 
   
    !i = i_fixed(bodyNumber)

    angvx_pre = angvxinit(bodyNumber)    &
                + ampangx(bodyNumber)    &
                *( SIN(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*cosphase(bodyNumber) &
                  +COS(2.0_CGREAL*PI*freqangx(bodyNumber)*time)*sinphase(bodyNumber) )
    angvy_pre = angvyinit(bodyNumber)    &
                + ampangy(bodyNumber)    &
                *( SIN(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*cosphase(bodyNumber) &
                  +COS(2.0_CGREAL*PI*freqangy(bodyNumber)*time)*sinphase(bodyNumber) )
    angvz_pre = angvzinit(bodyNumber)    &
                + ampangz(bodyNumber)    &
                *( SIN(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*cosphase(bodyNumber) &
                  +COS(2.0_CGREAL*PI*freqangz(bodyNumber)*time)*sinphase(bodyNumber) ) 

    angvx_rel = angvx_pre*C11(iGroup) + angvy_pre*C12(iGroup) + angvz_pre*C13(iGroup)
    angvy_rel = angvx_pre*C21(iGroup) + angvy_pre*C22(iGroup) + angvz_pre*C23(iGroup)
    angvz_rel = angvx_pre*C31(iGroup) + angvy_pre*C32(iGroup) + angvz_pre*C33(iGroup) 
      
    angvx(bodyNumber)  = angvx(bodyNumber) + angvx_rel
    angvy(bodyNumber)  = angvy(bodyNumber) + angvy_rel
    angvz(bodyNumber)  = angvz(bodyNumber) + angvz_rel
    
	IF(body_dim(bodyNumber) == 2)THEN
    angvx(bodyNumber)  = 0.0_CGREAL
    angvy(bodyNumber)  = 0.0_CGREAL
    angvz(bodyNumber)  = angvz(bodyNumber)
	ENDIF     
		
     
   
   ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)                          


END SUBROUTINE forced_motion
!---------------------------------------------------------------



!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
SUBROUTINE CALC_FORCE_MOMENT(bodyNumber)

USE global_parameters
USE flow_parameters
USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
USE boundary_arrays, ONLY : iblank
Use usr_module 
USE derivative
IMPLICIT NONE

INTEGER             :: i,j,k
INTEGER,INTENT (IN) ::bodyNumber

    density_ratio = density_solid(bodyNumber) / density_fluid

    lScale = 1.0_CGREAL      !=> Length scale
    vScale = 1.0_CGREAL      !=> Velocity scale
    
    IF ( body_dim(bodyNumber) == 2 ) THEN
    non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**2 )
    ELSE
    non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**3 )    
    ENDIF

    IF ( body_dim(bodyNumber) == 2 ) THEN
 
      force_x(bodyNumber) = 2.0_CGREAL * sCx(bodyNumber) / ( density_fluid * lScale* (vScale)**2*zout )
!-------No buoyance force for membrane    
!      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2*zout ) + &
!                            2.0_CGREAL * ( density_fluid-density_solid(bodyNumber) )*volume(bodyNumber)*zout*Fr / ( density_fluid * lScale* (vScale)**2*zout ) 
    IF(membrane_type(bodyNumber) == OPEN_MEMBRANE .OR. membrane_type(bodyNumber) == CLOSED_MEMBRANE)THEN
      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2*zout ) + &
                            2.0_CGREAL * ( 0.0_CGREAL   -density_solid(bodyNumber) )*volume(bodyNumber)*zout*Fr / ( density_fluid * (lScale)**2*zout ) 
                                            
    ELSE
      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2*zout ) + &
                            2.0_CGREAL * ( density_fluid-density_solid(bodyNumber) )*volume(bodyNumber)*zout*Fr / ( density_fluid * (lScale)**2*zout )                      
    ENDIF
    

      force_z(bodyNumber)  = 2.0_CGREAL * sCz(bodyNumber) / ( density_fluid * lScale* (vScale)**2*zout )
      moment_x(bodyNumber) = 2.0_CGREAL * scmx(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)
      moment_y(bodyNumber) = 2.0_CGREAL * scmy(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)
      moment_z(bodyNumber) = 2.0_CGREAL * scmz(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)*zout)
      
    ELSE !body_dim(bodyNumber) /= 2
      force_x(bodyNumber) = 2.0_CGREAL * sCx(bodyNumber) / ( density_fluid * lScale* (vScale)**2 )
!-------No buoyance force for membrane     
!      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2 ) + &
!                            2.0_CGREAL * ( density_fluid-density_solid(bodyNumber) )*volume(bodyNumber)*Fr / ( density_fluid * lScale* (vScale)**2 )
    IF(membrane_type(bodyNumber) == OPEN_MEMBRANE .OR. membrane_type(bodyNumber) == CLOSED_MEMBRANE)THEN
      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2 ) + &
                            2.0_CGREAL * ( 0.0_CGREAL   -density_solid(bodyNumber) )*volume(bodyNumber)*Fr / ( density_fluid * (lScale)**2 )
                            
    ELSE
      force_y(bodyNumber) = 2.0_CGREAL * sCy(bodyNumber) / ( density_fluid * lScale* (vScale)**2 ) + &
                            2.0_CGREAL * ( density_fluid-density_solid(bodyNumber) )*volume(bodyNumber)*Fr / ( density_fluid * (lScale)**2 )
    ENDIF 
                  
      force_z(bodyNumber)  = 2.0_CGREAL * sCz(bodyNumber) / ( density_fluid * lScale* (vScale)**2 )
      moment_x(bodyNumber) = 2.0_CGREAL * scmx(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3))
      moment_y(bodyNumber) = 2.0_CGREAL * scmy(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3))
      moment_z(bodyNumber) = 2.0_CGREAL * scmz(bodyNumber)/ ( density_fluid*(vScale**2)*(lScale**3)) 
      
    ENDIF !body_dim(bodyNumber) == 2  
         

END SUBROUTINE CALC_FORCE_MOMENT
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)





SUBROUTINE flow_induced_motion(bodyNumber)
! -----------------------------------------------------------
!  subroutine to specify wall velocity for  body whose motion is flow induced
! -----------------------------------------------------------

! Note the following eqns solved here

! Linear Momentum Equ..(Fx i + Fy j + Fy k) = m[(du/dt)i + (dv/dt - g)j + (dw/dt)k]

! Angular Momentum Equ..
!          n+1         n      n+1
!T = d([I]{omega}) /dt => ([I]{w})   - ([I]{w})    = dt * T

!Approximation used:        n    n+1     n            n+1
!T = d([I]{omega}) /dt => [I] ({w}   -  {w}  )    = dt * T
!
!              n        n+1
! rhside = ([I]{w})   + dt * T


    USE global_parameters
    USE flow_parameters
    USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
    USE boundary_arrays, ONLY : iblank
    Use usr_module 
    IMPLICIT NONE

    INTEGER             :: i,j,k
    INTEGER,INTENT (IN) ::bodyNumber

    REAL(KIND=CGREAL)                :: det
    REAL(KIND=CGREAL),DIMENSION(1:3) :: rhside

!    density_solid and density_fluid are defined 
!    in canonical_body_in.dat (Rajneesh)

    density_ratio = density_solid(bodyNumber) / density_fluid

    lScale = 1.0_CGREAL      !=> Length scale
    vScale = 1.0_CGREAL      !=> Velocity scale

    !CALL MofI_CofG(bodyNumber)  ! Moment of Inertia and Center of Gravity

!   Parallel verison, added by Rajneesh

	
	CALL MofI_CofG_par(bodyNumber)  ! Moment of Inertia and Center of Gravity	
	
	CALL CALC_FORCE_MOMENT(bodyNumber) !Calculate force and moment for single body
		
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    IF ( body_dim(bodyNumber) == 2 ) THEN
    non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**2 )
    ELSE
    non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**3 )    
    ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
	
!   bodyforce_moment

    IF ( body_dim(bodyNumber) == 2 ) THEN

      vxcent_prev = vxcent(bodyNumber)
      vycent_prev = vycent(bodyNumber)
      vzcent_prev = 0.0_CGREAL

      angvx_prev  = 0.0_CGREAL
      angvy_prev  = 0.0_CGREAL
      angvz_prev  = angvz(bodyNumber)

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      vxcent(bodyNumber)= xcentConstr(bodyNumber)*( vxcent_prev + dt*force_x(bodyNumber)*0.5/  &
	                      non_dim_mass(bodyNumber) )

      vycent(bodyNumber)= ycentConstr(bodyNumber)*( vycent_prev + dt*force_y(bodyNumber)*0.5/  &
	                      non_dim_mass(bodyNumber) )

      vzcent(bodyNumber)=  0.0_CGREAL

      rhside(1) = 0.0_CGREAL
      rhside(2) = 0.0_CGREAL
      rhside(3) = nonDimM_I(3,3,bodyNumber)*angvz_prev + dt*moment_z(bodyNumber)

      angvx(bodyNumber) = 0.0_CGREAL
      angvy(bodyNumber) = 0.0_CGREAL
      angvz(bodyNumber) = rhside(3)/nonDimM_I(3,3,bodyNumber)           
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

    ELSE

      vxcent_prev = vxcent(bodyNumber)
      vycent_prev = vycent(bodyNumber)
      vzcent_prev = vzcent(bodyNumber)

      angvx_prev  = angvx(bodyNumber)
      angvy_prev  = angvy(bodyNumber)
      angvz_prev  = angvz(bodyNumber)


!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      vxcent(bodyNumber)= xcentConstr(bodyNumber)*( vxcent_prev + dt*force_x(bodyNumber)*0.5/  &
	                      non_dim_mass(bodyNumber) )

      vycent(bodyNumber)= ycentConstr(bodyNumber)*( vycent_prev + dt*force_y(bodyNumber)*0.5/  &
	                      non_dim_mass(bodyNumber) )

      vzcent(bodyNumber)= zcentConstr(bodyNumber)*( vzcent_prev + dt*force_z(bodyNumber)*0.5/  &
	                      non_dim_mass(bodyNumber) )


      rhside(1) = nonDimM_I(1,1,bodyNumber)*angvx_prev + nonDimM_I(1,2,bodyNumber)*angvy_prev + &
	              nonDimM_I(1,3,bodyNumber)*angvz_prev + dt*moment_x(bodyNumber)
      rhside(2) = nonDimM_I(2,1,bodyNumber)*angvx_prev + nonDimM_I(2,2,bodyNumber)*angvy_prev + &
	              nonDimM_I(2,3,bodyNumber)*angvz_prev + dt*moment_y(bodyNumber)
      rhside(3) = nonDimM_I(3,1,bodyNumber)*angvx_prev + nonDimM_I(3,2,bodyNumber)*angvy_prev + &
	              nonDimM_I(3,3,bodyNumber)*angvz_prev + dt*moment_z(bodyNumber)

      det = nonDimM_I(1,1,bodyNumber)*(nonDimM_I(2,2,bodyNumber)*nonDimM_I(3,3,bodyNumber) - &
	        nonDimM_I(2,3,bodyNumber)*nonDimM_I(3,2,bodyNumber)) -                           &
            nonDimM_I(1,2,bodyNumber)*(nonDimM_I(2,1,bodyNumber)*nonDimM_I(3,3,bodyNumber) - &
			nonDimM_I(2,3,bodyNumber)*nonDimM_I(3,1,bodyNumber)) +                           &
            nonDimM_I(1,3,bodyNumber)*(nonDimM_I(2,1,bodyNumber)*nonDimM_I(3,2,bodyNumber) - &
			nonDimM_I(3,1,bodyNumber)*nonDimM_I(2,2,bodyNumber))


      invMI(1,1,bodyNumber) = (nonDimM_I(2,2,bodyNumber)*nonDimM_I(3,3,bodyNumber)-      &
	                           nonDimM_I(2,3,bodyNumber)*nonDimM_I(3,2,bodyNumber))/det
      invMI(1,2,bodyNumber) = (nonDimM_I(1,3,bodyNumber)*nonDimM_I(3,2,bodyNumber)-      &
	                           nonDimM_I(1,2,bodyNumber)*nonDimM_I(3,3,bodyNumber))/det
      invMI(1,3,bodyNumber) = (nonDimM_I(1,2,bodyNumber)*nonDimM_I(2,3,bodyNumber)-      &
	                           nonDimM_I(1,3,bodyNumber)*nonDimM_I(2,2,bodyNumber))/det

      invMI(2,1,bodyNumber) = (nonDimM_I(2,3,bodyNumber)*nonDimM_I(3,1,bodyNumber)-      &
	                           nonDimM_I(2,1,bodyNumber)*nonDimM_I(3,3,bodyNumber))/det
      invMI(2,2,bodyNumber) = (nonDimM_I(1,1,bodyNumber)*nonDimM_I(3,3,bodyNumber)-      &
	                           nonDimM_I(1,3,bodyNumber)*nonDimM_I(3,1,bodyNumber))/det
      invMI(2,3,bodyNumber) = (nonDimM_I(1,3,bodyNumber)*nonDimM_I(2,1,bodyNumber)-      &
	                           nonDimM_I(1,1,bodyNumber)*nonDimM_I(2,3,bodyNumber))/det
 
      invMI(3,1,bodyNumber) = (nonDimM_I(2,1,bodyNumber)*nonDimM_I(3,2,bodyNumber)-      &
	                           nonDimM_I(3,1,bodyNumber)*nonDimM_I(2,2,bodyNumber))/det
      invMI(3,2,bodyNumber) = (nonDimM_I(1,2,bodyNumber)*nonDimM_I(3,1,bodyNumber)-      &
	                           nonDimM_I(1,1,bodyNumber)*nonDimM_I(3,2,bodyNumber))/det
      invMI(3,3,bodyNumber) = (nonDimM_I(1,1,bodyNumber)*nonDimM_I(2,2,bodyNumber)-      &
	                           nonDimM_I(1,2,bodyNumber)*nonDimM_I(2,1,bodyNumber))/det

      angvx(bodyNumber)= invMI(1,1,bodyNumber)*rhside(1) + invMI(1,2,bodyNumber)*rhside(2) + invMI(1,3,bodyNumber)*rhside(3)
      angvy(bodyNumber)= invMI(2,1,bodyNumber)*rhside(1) + invMI(2,2,bodyNumber)*rhside(2) + invMI(2,3,bodyNumber)*rhside(3)
      angvz(bodyNumber)= invMI(3,1,bodyNumber)*rhside(1) + invMI(3,2,bodyNumber)*rhside(2) + invMI(3,3,bodyNumber)*rhside(3)
	      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    ENDIF


    write(1001,*)'--------------------'
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    write(1001,*)'Volume=', volume(bodyNumber)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    write(1001,*)'x,y,z',  xcent(bodyNumber), ycent(bodyNumber), zcent(bodyNumber)
    DO I = 1,3
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      write(1001,*)  (nonDimM_I(I,J,bodyNumber),J=1,3)

    ENDDO
    write(1002,200)time,vxcent(bodyNumber),vycent(bodyNumber),vzcent(bodyNumber),angvx(bodyNumber), &
                   angvy(bodyNumber),angvz(bodyNumber),force_x(bodyNumber),force_y(bodyNumber),force_z(bodyNumber), &
                   moment_x(bodyNumber),moment_y(bodyNumber),moment_z(bodyNumber)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
200  FORMAT(13(1x,E12.5))


END SUBROUTINE flow_induced_motion

!-----------------------------------------------------------------------------

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
SUBROUTINE MofI_CofG_par(bodyNumber) ! Created by Rajneesh

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays, ONLY : xc,yc,zc,dx,dy,dz
    USE derivative
    USE boundary_arrays

    USE usr_module

    USE unstructured_surface_arrays


Implicit None
INTEGER :: i,j,k, bodyNumber,iElem,iG, jG

        

   lScale = 1.0_CGREAL      !=> Length scale
   vScale = 1.0_CGREAL      !=> Velocity scale

   volume(bodyNumber) = 0.0_CGREAL
   I_XX(bodyNumber)   = 0.0_CGREAL
   I_YY(bodyNumber)   = 0.0_CGREAL
   I_ZZ(bodyNumber)   = 0.0_CGREAL
   I_XY(bodyNumber)   = 0.0_CGREAL
   I_YZ(bodyNumber)   = 0.0_CGREAL
   I_XZ(bodyNumber)   = 0.0_CGREAL
        
        
   xcent(bodyNumber)  = 0.0_CGREAL
   ycent(bodyNumber)  = 0.0_CGREAL
   zcent(bodyNumber)  = 0.0_CGREAL
 

! Centroid and Volume

    

 IF(membrane_type(bodyNumber) == OPEN_MEMBRANE .OR. membrane_type(bodyNumber) == CLOSED_MEMBRANE)THEN

       DO iElem=1,totNumTriElem(bodyNumber)

          xcent(bodyNumber) = xcent(bodyNumber) + triElemCentx(bodyNumber,iElem)*triElemArea(bodyNumber,iElem)
          ycent(bodyNumber) = ycent(bodyNumber) + triElemCenty(bodyNumber,iElem)*triElemArea(bodyNumber,iElem)
          zcent(bodyNumber) = zcent(bodyNumber) + triElemCentz(bodyNumber,iElem)*triElemArea(bodyNumber,iElem)
          volume(bodyNumber)= volume(bodyNumber)+ triElemArea(bodyNumber,iElem) !Assuming the thickness of membrane is unit.
	  
	   ENDDO !iElem 
!  -----------------------------Centroid of mass and non-dimension mass for membrane body
       xcent(bodyNumber) = xcent(bodyNumber) / volume(bodyNumber)
       ycent(bodyNumber) = ycent(bodyNumber) / volume(bodyNumber)
       zcent(bodyNumber) = zcent(bodyNumber) / volume(bodyNumber)    
	   non_dim_mass(bodyNumber) = volume(bodyNumber)*density_solid(bodyNumber) / (density_fluid*lScale**3)
	   
!  -----------------------------Moment of inertia for membrane body
       DO iElem=1,totNumTriElem(bodyNumber)

        I_XX(bodyNumber) = I_XX(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                                      
		                   ((triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**2 +   &
						   (triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**2)     

        I_YY(bodyNumber) = I_YY(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                                
		                   ((triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**2 +   &
						   (triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**2)  

        I_ZZ(bodyNumber) = I_ZZ(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                               
		                   ((triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**2 +   &
						   (triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**2)  

        I_XY(bodyNumber) = I_XY(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                                 
		                   ( (triElemCentx(bodyNumber,iElem)-xcent(bodyNumber)) *     &
						   (triElemCenty(bodyNumber,iElem)-ycent(bodyNumber)) )      

        I_YZ(bodyNumber) = I_YZ(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                                 
		                   ( (triElemCentz(bodyNumber,iElem)-zcent(bodyNumber)) *     &
						   (triElemCenty(bodyNumber,iElem)-ycent(bodyNumber)) )      

        I_XZ(bodyNumber) = I_XZ(bodyNumber) +  triElemArea(bodyNumber,iElem) *        &                                                                
		                   ( (triElemCentx(bodyNumber,iElem)-xcent(bodyNumber)) *     &
						   (triElemCentz(bodyNumber,iElem)-zcent(bodyNumber)) )
      
	   ENDDO !iElem	        
	   
      I_XX(bodyNumber) = I_XX(bodyNumber) * density_solid(bodyNumber)
      I_YY(bodyNumber) = I_YY(bodyNumber) * density_solid(bodyNumber)
      I_ZZ(bodyNumber) = I_ZZ(bodyNumber) * density_solid(bodyNumber)
      I_XY(bodyNumber) = I_XY(bodyNumber) * density_solid(bodyNumber)
      I_YZ(bodyNumber) = I_YZ(bodyNumber) * density_solid(bodyNumber)
      I_XZ(bodyNumber) = I_XZ(bodyNumber) * density_solid(bodyNumber)	   

ELSE !==>SOLID_BODY

!-----------------------------------------------SURFACE INTEGRAL METHOD (CALCULATE MOMENT OF INERTIA)
   IF(SURFACE_INTEGRAL(bodyNumber) == 1)THEN

      IF ( body_dim(bodyNumber) == 2 ) THEN   ! Area Moment of Inertia
!---------------Calculate Volume and Center of Mass (2D)    
       DO iElem=1,totNumTriElem(bodyNumber)      
       
       volume(bodyNumber) = volume(bodyNumber) - (triElemNormx(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)/2.0_CGREAL +  &
                                                  triElemNormy(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem) 
                                                                                            
       xcent(bodyNumber)  = xcent(bodyNumber)  - (triElemNormx(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)**2/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem) 
       
       
       ycent(bodyNumber)  = ycent(bodyNumber)  - (triElemNormy(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)**2/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
                                                                                                                
       ENDDO  
      
      xcent(bodyNumber) = xcent(bodyNumber)/volume(bodyNumber)
      ycent(bodyNumber) = ycent(bodyNumber)/volume(bodyNumber)
      volume(bodyNumber) = volume(bodyNumber)/zout       
      non_dim_mass(bodyNumber) = volume(bodyNumber)*density_solid(bodyNumber) / (density_fluid*lScale**2)  
     
                      
      
       DO iElem=1,totNumTriElem(bodyNumber)
       
   
       I_ZZ(bodyNumber) = I_ZZ(bodyNumber) - (triElemNormx(bodyNumber,iElem)*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**3/3.0_CGREAL +  &
                                              triElemNormy(bodyNumber,iElem)*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**3/3.0_CGREAL)*  &
                                              triElemArea(bodyNumber,iElem)/zout
                                            
!       I_ZZ(bodyNumber) = I_ZZ(bodyNumber) - (triElemNormx(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)**2/3.0_CGREAL +  &
!                                              triElemNormy(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)**2/3.0_CGREAL)*  &
!                                              triElemArea(bodyNumber,iElem)      
                                                                                 
       ENDDO 

      ELSE   ! ===> body_dim(bodyNumber) == 3
!---------------Calculate Volume and Center of Mass (3D) 
       DO iElem=1,totNumTriElem(bodyNumber)      
       
       volume(bodyNumber) = volume(bodyNumber) - (triElemNormx(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)/3.0_CGREAL +  &
                                                  triElemNormy(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)/3.0_CGREAL +  & 
                                                  triElemNormz(bodyNumber,iElem)*triElemCentz(bodyNumber,iElem)/3.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem) 
                                                  
                                         
       xcent(bodyNumber)  = xcent(bodyNumber)  - (triElemNormx(bodyNumber,iElem)*triElemCentx(bodyNumber,iElem)**2/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem) 
       
       ycent(bodyNumber)  = ycent(bodyNumber)  - (triElemNormy(bodyNumber,iElem)*triElemCenty(bodyNumber,iElem)**2/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
                                                  
       zcent(bodyNumber)  = zcent(bodyNumber)  - (triElemNormz(bodyNumber,iElem)*triElemCentz(bodyNumber,iElem)**2/2.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
                                                                                                                                                                                                              
       ENDDO 
       
       
      xcent(bodyNumber) = xcent(bodyNumber)/volume(bodyNumber)
      ycent(bodyNumber) = ycent(bodyNumber)/volume(bodyNumber)
      zcent(bodyNumber) = zcent(bodyNumber)/volume(bodyNumber)     
      non_dim_mass(bodyNumber) = volume(bodyNumber)*density_solid(bodyNumber) / (density_fluid*lScale**3)         
       
       
        DO iElem=1,totNumTriElem(bodyNumber)
     
           I_XX(bodyNumber) = I_XX(bodyNumber) - (triElemNormy(bodyNumber,iElem)*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**3/3.0_CGREAL +  &
                                                  triElemNormz(bodyNumber,iElem)*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**3/3.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
          
     
           I_YY(bodyNumber) = I_YY(bodyNumber) - (triElemNormx(bodyNumber,iElem)*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**3/3.0_CGREAL +  &
                                                  triElemNormz(bodyNumber,iElem)*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**3/3.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
     
     
           I_ZZ(bodyNumber) = I_ZZ(bodyNumber) - (triElemNormx(bodyNumber,iElem)*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**3/3.0_CGREAL +  &
                                                  triElemNormy(bodyNumber,iElem)*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**3/3.0_CGREAL)*  &
                                                  triElemArea(bodyNumber,iElem)
                                            
                                            
           I_XY(bodyNumber) = I_XY(bodyNumber) - (triElemNormx(bodyNumber,iElem)*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**2*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))/4.0_CGREAL +   &
                                                  triElemNormy(bodyNumber,iElem)*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**2*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))/4.0_CGREAL )*  &
                                                  triElemArea(bodyNumber,iElem)
	                                 

           I_YZ(bodyNumber) = I_YZ(bodyNumber) - (triElemNormy(bodyNumber,iElem)*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))**2*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))/4.0_CGREAL +   &
                                                  triElemNormz(bodyNumber,iElem)*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**2*(triElemCenty(bodyNumber,iElem)-ycent(bodyNumber))/4.0_CGREAL )*  &
                                                  triElemArea(bodyNumber,iElem)
	                                 
	
           I_XZ(bodyNumber) = I_XZ(bodyNumber) - (triElemNormx(bodyNumber,iElem)*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))**2*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))/4.0_CGREAL +   &
                                                  triElemNormz(bodyNumber,iElem)*(triElemCentz(bodyNumber,iElem)-zcent(bodyNumber))**2*(triElemCentx(bodyNumber,iElem)-xcent(bodyNumber))/4.0_CGREAL )*  &
                                                  triElemArea(bodyNumber,iElem)                                                                          
        ENDDO
        
      ENDIF   ! ===> body_dim(bodyNumber) == 2 or 3  
            
      I_XX(bodyNumber) = I_XX(bodyNumber) * density_solid(bodyNumber)
      I_YY(bodyNumber) = I_YY(bodyNumber) * density_solid(bodyNumber)
      I_ZZ(bodyNumber) = I_ZZ(bodyNumber) * density_solid(bodyNumber)
      I_XY(bodyNumber) = I_XY(bodyNumber) * density_solid(bodyNumber)
      I_YZ(bodyNumber) = I_YZ(bodyNumber) * density_solid(bodyNumber)
      I_XZ(bodyNumber) = I_XZ(bodyNumber) * density_solid(bodyNumber)        
      
   ELSE  !SURFACE_INTEGRAL(bodyNumber) == 0===>Volume Integral
	
     IF ( body_dim(bodyNumber) == 2 ) THEN
	
	
      k = 1
      DO j = 1, nyc
        jG=L2GJ(j)

        DO i = 1, nxc
          iG=L2GI(i)

          IF(bodyNum(i,j,k) == bodyNumber)THEN
          volume(bodyNumber) = volume(bodyNumber) +        dx(iG)*dy(jG)*REAL(iblank(i,j,k),KIND=CGREAL)  	   
	      xcent(bodyNumber)  = xcent(bodyNumber)  + xc(iG)*dx(iG)*dy(jG)*REAL(iblank(i,j,k),KIND=CGREAL) 
          ycent(bodyNumber)  = ycent(bodyNumber)  + yc(jG)*dx(iG)*dy(jG)*REAL(iblank(i,j,k),KIND=CGREAL)
          ENDIF
      END DO
      END DO
 
#ifdef MPI
      CALL par_getSumReal(volume(bodyNumber))         
	  CALL par_getSumReal(xcent(bodyNumber))
	  CALL par_getSumReal(ycent(bodyNumber))
#endif
      
 
      xcent(bodyNumber) = xcent(bodyNumber)/volume(bodyNumber)
      ycent(bodyNumber) = ycent(bodyNumber)/volume(bodyNumber)
      non_dim_mass(bodyNumber) = volume(bodyNumber)*density_solid(bodyNumber) / (density_fluid*lScale**2)            
                                                                                    	  
        k = 1		
		 DO j = 1, nyc
          jG=L2GJ(j)

        DO i = 1, nxc
          iG=L2GI(i)

          IF(bodyNum(i,j,k) == bodyNumber)THEN	  
		  I_ZZ(bodyNumber) = I_ZZ(bodyNumber) + ((yc(jG)-ycent(bodyNumber))**2 + (xc(iG)-xcent(bodyNumber))**2)&
	                                 *dx(iG)*dy(jG)*REAL(iblank(i,j,k),KIND=CGREAL)
	      ENDIF	
 	                                 
        ENDDO
        ENDDO
        	
#ifdef MPI
	  CALL par_getSumReal(I_ZZ(bodyNumber))	  
#endif		

      I_ZZ(bodyNumber) = I_ZZ(bodyNumber) * density_solid(bodyNumber)

	  
      ELSE  ! ===> body_dim(bodyNumber) == 3
	  
	  DO k = 1, nzc
      DO j = 1, nyc
        jG=L2GJ(j)

        DO i = 1, nxc
          iG=L2GI(i)

          IF(bodyNum(i,j,k) == bodyNumber)THEN          
          volume(bodyNumber) = volume(bodyNumber) +        dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL) 	     
	      xcent(bodyNumber)  = xcent(bodyNumber)  + xc(iG)*dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL) 
          ycent(bodyNumber)  = ycent(bodyNumber)  + yc(jG)*dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
          zcent(bodyNumber)  = zcent(bodyNumber)  + zc(k) *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
          ENDIF
 		  
		  
      END DO
      END DO
      END DO
	  
#ifdef MPI
      CALL par_getSumReal(volume(bodyNumber))     
	  CALL par_getSumReal(xcent(bodyNumber))
	  CALL par_getSumReal(ycent(bodyNumber))
	  CALL par_getSumReal(zcent(bodyNumber))
#endif
    
      xcent(bodyNumber) = xcent(bodyNumber)/volume(bodyNumber)
      ycent(bodyNumber) = ycent(bodyNumber)/volume(bodyNumber)
      zcent(bodyNumber) = zcent(bodyNumber)/volume(bodyNumber)
      non_dim_mass(bodyNumber) = volume(bodyNumber)*density_solid(bodyNumber) / (density_fluid*lScale**3)

! Moment of Inertia Tensor

     DO k = 1, nzc
      DO j = 1, nyc
        jG=L2GJ(j)

        DO i = 1, nxc
          iG=L2GI(i)

          IF(bodyNum(i,j,k) == bodyNumber)THEN   
          I_XX(bodyNumber) = I_XX(bodyNumber) + ((yc(jG)-ycent(bodyNumber))**2 + (zc(k)-zcent(bodyNumber))**2)&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
									 
          I_YY(bodyNumber) = I_YY(bodyNumber) + ((xc(iG)-xcent(bodyNumber))**2 + (zc(k)-zcent(bodyNumber))**2)&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
									 
	      I_ZZ(bodyNumber) = I_ZZ(bodyNumber) + ((yc(jG)-ycent(bodyNumber))**2 + (xc(iG)-xcent(bodyNumber))**2)&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)								 
 
          I_XY(bodyNumber) = I_XY(bodyNumber) + ((xc(iG)-xcent(bodyNumber))*(yc(jG)-ycent(bodyNumber)))&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)

          I_YZ(bodyNumber) = I_YZ(bodyNumber) + ((zc(k) -zcent(bodyNumber))*(yc(jG)-ycent(bodyNumber)))&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
	
          I_XZ(bodyNumber) = I_XZ(bodyNumber) + ((xc(iG)-xcent(bodyNumber))*(zc(k)-zcent(bodyNumber)))&
	                                 *dx(iG)*dy(jG)*dz(k)*REAL(iblank(i,j,k),KIND=CGREAL)
	      ENDIF	
	  
      END DO
      END DO
      END DO

#ifdef MPI
      CALL par_getSumReal(I_XX(bodyNumber))
	  CALL par_getSumReal(I_YY(bodyNumber))
	  CALL par_getSumReal(I_ZZ(bodyNumber))
	  CALL par_getSumReal(I_XY(bodyNumber))
	  CALL par_getSumReal(I_YZ(bodyNumber))
	  CALL par_getSumReal(I_XZ(bodyNumber)) 
#endif

   
      I_XX(bodyNumber) = I_XX(bodyNumber) * density_solid(bodyNumber)
      I_YY(bodyNumber) = I_YY(bodyNumber) * density_solid(bodyNumber)
      I_ZZ(bodyNumber) = I_ZZ(bodyNumber) * density_solid(bodyNumber)
      I_XY(bodyNumber) = I_XY(bodyNumber) * density_solid(bodyNumber)
      I_YZ(bodyNumber) = I_YZ(bodyNumber) * density_solid(bodyNumber)
      I_XZ(bodyNumber) = I_XZ(bodyNumber) * density_solid(bodyNumber)

      ENDIF  ! ===> body_dim(bodyNumber) == 2 OR 3

   ENDIF !SURFACE_INTEGRAL 
      
ENDIF !WHETHER MEMBRANE OR NOT
	  
      nonDimM_I(1,1,bodyNumber) =2.0_CGREAL *I_XX(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(1,2,bodyNumber) =2.0_CGREAL *I_XY(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(1,3,bodyNumber) =2.0_CGREAL *I_XZ(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(2,1,bodyNumber) =2.0_CGREAL *I_XY(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(2,2,bodyNumber) =2.0_CGREAL *I_YY(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(2,3,bodyNumber) =2.0_CGREAL *I_YZ(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(3,1,bodyNumber) =2.0_CGREAL *I_XZ(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(3,2,bodyNumber) =2.0_CGREAL *I_YZ(bodyNumber)/(density_fluid*lScale**5)
      nonDimM_I(3,3,bodyNumber) =2.0_CGREAL *I_ZZ(bodyNumber)/(density_fluid*lScale**5)

END SUBROUTINE MofI_CofG_par
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)



! ======================Added for FEA===========================
! This subroutine transfers the velocity and coordinate of structure vortex
! from FEA arrays to Marker arrays

!---------------------------------------------------------------
SUBROUTINE fea_flow_structure(ibody)                            !SER_TO_PAR. QX. CH10

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE finite_element_parameters
    USE finite_element_arrays

    IMPLICIT NONE

    INTEGER,INTENT (IN) ::ibody
    INTEGER :: ifea_body, ipt_fea
    INTEGER :: i, j, k, n
    REAL(KIND=CGREAL) :: contactlocation

 IF(ImtheBOSS) THEN
    DO i = 1, nPtsBodyMarker(ibody)
      j = markertofeapointtable(i, 1, ibody)
      n = markertofeapointtable(i, 2, ibody)
      IF (j<0 .or. n<0) THEN
       uBodyMarker(ibody, i) = 0.0_CGREAL
       vBodyMarker(ibody, i) = 0.0_CGREAL
       wBodyMarker(ibody, i) = 0.0_CGREAL
       xBodyMarker(ibody, i) = xBodyMarker(ibody, i)
       yBodyMarker(ibody, i) = yBodyMarker(ibody, i)
       zBodyMarker(ibody, i) = zBodyMarker(ibody, i)
      ELSE
        IF (fea_nd == 3) THEN
          uBodyMarker(ibody, i) =  fea_v(fea_nd * (j - 1) + 1, n)
          vBodyMarker(ibody, i) =  fea_v(fea_nd * (j - 1) + 2, n)
          wBodyMarker(ibody, i) =  fea_v(fea_nd * (j - 1) + 3, n)
          xBodyMarker(ibody, i) =  fea_cood(j,  1,  n) + fea_d(fea_nd * (j - 1) + 1, n)
          yBodyMarker(ibody, i) =  fea_cood(j,  2,  n) + fea_d(fea_nd * (j - 1) + 2, n)
          zBodyMarker(ibody, i) =  fea_cood(j,  3,  n) + fea_d(fea_nd * (j - 1) + 3, n)
        ELSE IF (fea_nd == 2) THEN
          uBodyMarker(ibody, i) =  fea_v(fea_nd * (j - 1) + 1, n)
          vBodyMarker(ibody, i) =  fea_v(fea_nd * (j - 1) + 2, n)
          wBodyMarker(ibody, i) =  wBodyMarker(ibody, i)
          xBodyMarker(ibody, i) =  fea_cood(j,  1,  n) + fea_d(fea_nd * (j - 1) + 1, n)
          yBodyMarker(ibody, i) =  fea_cood(j,  2,  n) + fea_d(fea_nd * (j - 1) + 2, n)
          zBodyMarker(ibody, i) =  zBodyMarker(ibody, i)
        ENDIF
      ENDIF
   ENDDO ! i
!Gometric contact model    
   IF(fea_contactprobflag==1 .and. fea_contact_model == GEO_CONTACT_MODEL) THEN
    DO ifea_body = 1, fea_nbody
     IF(fea_contactflag(ifea_body) == 1) THEN
      DO i = 1, fea_ncontactpoint(ifea_body)
       ipt_fea = fea_icontactpoint(i, ifea_body)
       j = featomarkerpointtable(ipt_fea, 1, ifea_body)
       n = featomarkerpointtable(ipt_fea, 2, ifea_body)
       IF(n /= ibody) CYCLE      
 
       SELECT CASE(fea_icontactdir(ifea_body))
        CASE(ICOORD)
         contactlocation = x(fea_icontactplane(ifea_body))
         xBodyMarker(n, j) = contactlocation
         uBodyMarker(n, j) = 0.0_CGREAL
         IF(ndim == DIM_2D) THEN
          DO k = 2, nz
            xBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = contactlocation
            uBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = 0.0_CGREAL
          ENDDO
         ENDIF
        CASE(JCOORD)
         contactlocation = y(fea_icontactplane(ifea_body))
         yBodyMarker(n, j) = contactlocation
         vBodyMarker(n, j) = 0.0_CGREAL
         IF(ndim == DIM_2D) THEN
          DO k = 2, nz
            yBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = contactlocation
            vBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = 0.0_CGREAL
          ENDDO
         ENDIF

        CASE(KCOORD)
         contactlocation = z(fea_icontactplane(ifea_body))
         zBodyMarker(n, j) = contactlocation
         wBodyMarker(n, j) = 0.0_CGREAL
         IF(ndim == DIM_2D) THEN
          DO k = 2, nz
            zBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = contactlocation
            wBodyMarker(n, j+nPtsBodyMarker(n)*(k-1)/nz) = 0.0_CGREAL
          ENDDO
         ENDIF

        ENDSELECT
         
      ENDDO !i
     ENDIF !fea_contactflag(ifea_body) 
    ENDDO !ifea_body
   ENDIF !fea_contactproflag & fea_contact_model

 ENDIF

#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(uBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(vBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(wBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
#  endif
END SUBROUTINE fea_flow_structure

!---------------------------------------------------------------



!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
SUBROUTINE CALC_FORCE_MOMENT_combined(groupNumber) 
USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
Use usr_module
Use implicit_coupling_parameters
IMPLICIT NONE

INTEGER,INTENT (IN)  :: groupNumber
INTEGER               :: bodyNumber

lScale = 1.0_CGREAL      !=> Length scale
vScale = 1.0_CGREAL      !=> Velocity scale

force_x_combined(:)  = 0.0_CGREAL 
force_y_combined(:)  = 0.0_CGREAL 
force_z_combined(:)  = 0.0_CGREAL 
moment_x_combined(:) = 0.0_CGREAL 
moment_y_combined(:) = 0.0_CGREAL 
moment_z_combined(:) = 0.0_CGREAL 


    DO bodyNumber = 1, nbody
    
    IF(combined_Group_index(bodyNumber) == groupNumber)THEN
     
     CALL CALC_FORCE_MOMENT(bodyNumber)       
              
     !Right Hand Side Coordinate
     force_x_combined(groupNumber)  = force_x_combined(groupNumber)  + force_x(bodyNumber)
     force_y_combined(groupNumber)  = force_y_combined(groupNumber)  + force_y(bodyNumber)
     force_z_combined(groupNumber)  = force_z_combined(groupNumber)  + force_z(bodyNumber)
     
     moment_x_combined(groupNumber) = moment_x_combined(groupNumber) + moment_x(bodyNumber) -                        &
                                      force_z(bodyNumber)*(ycent_combined(groupNumber)-ycent(bodyNumber))/lScale +   &     
                                      force_y(bodyNumber)*(zcent_combined(groupNumber)-zcent(bodyNumber))/lScale
                                         
     moment_y_combined(groupNumber) = moment_y_combined(groupNumber) + moment_y(bodyNumber) -                        &
                                      force_x(bodyNumber)*(zcent_combined(groupNumber)-zcent(bodyNumber))/lScale +   &
                                      force_z(bodyNumber)*(xcent_combined(groupNumber)-xcent(bodyNumber))/lScale
                                              
     moment_z_combined(groupNumber) = moment_z_combined(groupNumber) + moment_z(bodyNumber) -                        &
                                      force_y(bodyNumber)*(xcent_combined(groupNumber)-xcent(bodyNumber))/lScale +   &     
                                      force_x(bodyNumber)*(ycent_combined(groupNumber)-ycent(bodyNumber))/lScale                                         

     
    ENDIF 
     
    ENDDO !bodyNumber


END SUBROUTINE CALC_FORCE_MOMENT_combined



SUBROUTINE combined_body_motion()

! -----------------------------------------------------------
!  subroutine to calculate wall velocity for combined multi-bodies
! -----------------------------------------------------------
   USE global_parameters
   USE flow_parameters
   USE grid_arrays
   USE boundary_arrays
   USE usr_module
   USE derivative
   USE implicit_coupling_parameters
   IMPLICIT NONE

   INTEGER             :: i,j,k
   INTEGER             :: groupNumber, bodyNumber
   REAL(KIND=CGREAL)                :: det
   REAL(KIND=CGREAL)                :: temp_angvx,temp_angvy,temp_angvz
   REAL(KIND=CGREAL)                :: temp_vxcent_pre,temp_vycent_pre,temp_vzcent_pre  
   REAL(KIND=CGREAL)                :: temp_vcent_magni_pre 
   REAL(KIND=CGREAL)                :: temp_vxcent_new,temp_vycent_new,temp_vzcent_new  
   REAL(KIND=CGREAL)                :: temp_vcent_magni_new
   REAL(KIND=CGREAL)                :: average_vxcent,average_vycent,average_vzcent 
   REAL(KIND=CGREAL)                :: temp_angvx_combined,temp_angvy_combined,temp_angvz_combined        
   REAL(KIND=CGREAL),DIMENSION(1:3) :: rhside

   lScale = 1.0_CGREAL      !=> Length scale
   vScale = 1.0_CGREAL      !=> Velocity scale



   
!  -----------------------------Update velocity for combined group
   DO groupNumber = 1, nGroup_Combined

!  -----------------------------Calculate total mass and moment of inertia for combined bodies(MofI_CofG_combined())
   CALL MofI_CofG_combined(groupNumber)  
   
!  -----------------------------Calculate total force and moment for combined bodies in different group
   CALL CALC_FORCE_MOMENT_combined(groupNumber) 
   
   
       IF ( ndim == 2 ) THEN

         vxcent_combined_prev = vxcent_combined(groupNumber)
         vycent_combined_prev = vycent_combined(groupNumber)
         vzcent_combined_prev = 0.0_CGREAL

         angvx_combined_prev  = 0.0_CGREAL
         angvy_combined_prev  = 0.0_CGREAL
         angvz_combined_prev  = angvz_combined(groupNumber)
         

         vxcent_combined(groupNumber)= vxcent_combined_prev + dt*force_x_combined(groupNumber)*0.5/non_dim_mass_combined(groupNumber) 

         vycent_combined(groupNumber)= vycent_combined_prev + dt*force_y_combined(groupNumber)*0.5/non_dim_mass_combined(groupNumber)

         vzcent_combined(groupNumber)= 0.0_CGREAL        
                 
         
         rhside(1) = 0.0_CGREAL
         rhside(2) = 0.0_CGREAL
         rhside(3) = nonDimM_I_combined(3,3,groupNumber)*angvz_combined_prev + dt*moment_z_combined(groupNumber)

         angvx_combined(groupNumber) = 0.0_CGREAL
         angvy_combined(groupNumber) = 0.0_CGREAL
         angvz_combined(groupNumber) = rhside(3)/nonDimM_I_combined(3,3,groupNumber)               


       ELSE !===>( ndim == 3 )

         vxcent_combined_prev = vxcent_combined(groupNumber)
         vycent_combined_prev = vycent_combined(groupNumber)
         vzcent_combined_prev = vzcent_combined(groupNumber)

         angvx_combined_prev  = angvx_combined(groupNumber)
         angvy_combined_prev  = angvy_combined(groupNumber)
         angvz_combined_prev  = angvz_combined(groupNumber)


         vxcent_combined(groupNumber)=  vxcent_combined_prev + dt*force_x_combined(groupNumber)*0.5/non_dim_mass_combined(groupNumber) 

         vycent_combined(groupNumber)=  vycent_combined_prev + dt*force_y_combined(groupNumber)*0.5/non_dim_mass_combined(groupNumber) 

         vzcent_combined(groupNumber)=  vzcent_combined_prev + dt*force_z_combined(groupNumber)*0.5/non_dim_mass_combined(groupNumber) 
         
         rhside(1) = nonDimM_I_combined(1,1,groupNumber)*angvx_combined_prev + nonDimM_I_combined(1,2,groupNumber)*angvy_combined_prev +   &
	                 nonDimM_I_combined(1,3,groupNumber)*angvz_combined_prev + dt*moment_x_combined(groupNumber)
         rhside(2) = nonDimM_I_combined(2,1,groupNumber)*angvx_combined_prev + nonDimM_I_combined(2,2,groupNumber)*angvy_combined_prev +   &
	                 nonDimM_I_combined(2,3,groupNumber)*angvz_combined_prev + dt*moment_y_combined(groupNumber)
         rhside(3) = nonDimM_I_combined(3,1,groupNumber)*angvx_combined_prev + nonDimM_I_combined(3,2,groupNumber)*angvy_combined_prev +   &
	                 nonDimM_I_combined(3,3,groupNumber)*angvz_combined_prev + dt*moment_z_combined(groupNumber)


         det = nonDimM_I_combined(1,1,groupNumber)*(nonDimM_I_combined(2,2,groupNumber)*nonDimM_I_combined(3,3,groupNumber) -   & 
	           nonDimM_I_combined(2,3,groupNumber)*nonDimM_I_combined(3,2,groupNumber)) -                                       &
               nonDimM_I_combined(1,2,groupNumber)*(nonDimM_I_combined(2,1,groupNumber)*nonDimM_I_combined(3,3,groupNumber) -   &
		       nonDimM_I_combined(2,3,groupNumber)*nonDimM_I_combined(3,1,groupNumber)) +                                       &
               nonDimM_I_combined(1,3,groupNumber)*(nonDimM_I_combined(2,1,groupNumber)*nonDimM_I_combined(3,2,groupNumber) -   &
		       nonDimM_I_combined(3,1,groupNumber)*nonDimM_I_combined(2,2,groupNumber))

         invMI_combined(1,1,groupNumber) = (nonDimM_I_combined(2,2,groupNumber)*nonDimM_I_combined(3,3,groupNumber)-nonDimM_I_combined(2,3,groupNumber)*nonDimM_I_combined(3,2,groupNumber))/det
         invMI_combined(1,2,groupNumber) = (nonDimM_I_combined(1,3,groupNumber)*nonDimM_I_combined(3,2,groupNumber)-nonDimM_I_combined(1,2,groupNumber)*nonDimM_I_combined(3,3,groupNumber))/det
         invMI_combined(1,3,groupNumber) = (nonDimM_I_combined(1,2,groupNumber)*nonDimM_I_combined(2,3,groupNumber)-nonDimM_I_combined(1,3,groupNumber)*nonDimM_I_combined(2,2,groupNumber))/det

         invMI_combined(2,1,groupNumber) = (nonDimM_I_combined(2,3,groupNumber)*nonDimM_I_combined(3,1,groupNumber)-nonDimM_I_combined(2,1,groupNumber)*nonDimM_I_combined(3,3,groupNumber))/det
         invMI_combined(2,2,groupNumber) = (nonDimM_I_combined(1,1,groupNumber)*nonDimM_I_combined(3,3,groupNumber)-nonDimM_I_combined(1,3,groupNumber)*nonDimM_I_combined(3,1,groupNumber))/det
         invMI_combined(2,3,groupNumber) = (nonDimM_I_combined(1,3,groupNumber)*nonDimM_I_combined(2,1,groupNumber)-nonDimM_I_combined(1,1,groupNumber)*nonDimM_I_combined(2,3,groupNumber))/det

         invMI_combined(3,1,groupNumber) = (nonDimM_I_combined(2,1,groupNumber)*nonDimM_I_combined(3,2,groupNumber)-nonDimM_I_combined(3,1,groupNumber)*nonDimM_I_combined(2,2,groupNumber))/det
         invMI_combined(3,2,groupNumber) = (nonDimM_I_combined(1,2,groupNumber)*nonDimM_I_combined(3,1,groupNumber)-nonDimM_I_combined(1,1,groupNumber)*nonDimM_I_combined(3,2,groupNumber))/det
         invMI_combined(3,3,groupNumber) = (nonDimM_I_combined(1,1,groupNumber)*nonDimM_I_combined(2,2,groupNumber)-nonDimM_I_combined(1,2,groupNumber)*nonDimM_I_combined(2,1,groupNumber))/det
     

         angvx_combined(groupNumber)= invMI_combined(1,1,groupNumber)*rhside(1) + invMI_combined(1,2,groupNumber)*rhside(2) + invMI_combined(1,3,groupNumber)*rhside(3)
         angvy_combined(groupNumber)= invMI_combined(2,1,groupNumber)*rhside(1) + invMI_combined(2,2,groupNumber)*rhside(2) + invMI_combined(2,3,groupNumber)*rhside(3)
         angvz_combined(groupNumber)= invMI_combined(3,1,groupNumber)*rhside(1) + invMI_combined(3,2,groupNumber)*rhside(2) + invMI_combined(3,3,groupNumber)*rhside(3)   

!dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!For perturbation analysis
OPEN(365,FILE='perturb_input.dat')
READ(365,*)
READ(365,*) perturb_switch
IF(perturb_switch == 1)THEN
READ(365,*)
READ(365,*) vxcent_combined(groupNumber), vycent_combined(groupNumber), vzcent_combined(groupNumber)
READ(365,*) 
READ(365,*) angvx_combined(groupNumber), angvy_combined(groupNumber), angvz_combined(groupNumber)

OPEN(366,FILE='outout.dat')
WRITE(366,*) vxcent_combined(groupNumber), vycent_combined(groupNumber), vzcent_combined(groupNumber)
WRITE(366,*) angvx_combined(groupNumber), angvy_combined(groupNumber), angvz_combined(groupNumber)
ENDIF

CLOSE(365)
!dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd


       ENDIF ! Judge dimension
	  
   ENDDO !groupNumber
!  -----------------------------Calculate partial difference coefficients
   IF(derivative_flag == 1)THEN
     
     CALL non_inertia_UVQ_increase
     
   ENDIF
!  -----------------------------Calculate vector rotation matrix for different combined group
   DO groupNumber = 1, nGroup_Combined
	  
      CALL Vector_Rotation_Matrix(groupNumber)
	  
   ENDDO 

!  -----------------------------Update velocity for bodies in combined group respectively
   DO groupNumber = 1, nGroup_Combined

    DO bodyNumber = 1, nbody

     IF(combined_Group_index(bodyNumber) == groupNumber)THEN

       IF ( body_dim(bodyNumber) == 2 ) THEN       
        
         angvx(bodyNumber) = 0.0_CGREAL
         angvy(bodyNumber) = 0.0_CGREAL
         angvz(bodyNumber) = angvz_combined(groupNumber)
         
         temp_angvx  = 0.5_CGREAL*(angvx_combined_old(groupNumber)+angvx_combined(groupNumber))
         temp_angvy  = 0.5_CGREAL*(angvy_combined_old(groupNumber)+angvy_combined(groupNumber))
         temp_angvz  = 0.5_CGREAL*(angvz_combined_old(groupNumber)+angvz_combined(groupNumber))     

         !Velocity for different body center is the combination of Velocity of combined group's center plus relative rotation (Right Hand Side Coordinate)
         temp_vxcent_pre=   temp_angvz*(ycent_combined(groupNumber)-ycent(bodyNumber))

         temp_vycent_pre= - temp_angvz*(xcent_combined(groupNumber)-xcent(bodyNumber))

         temp_vzcent_pre= 0.0_CGREAL  

       ELSE !===>( body_dim(bodyNumber) == 3 )

         angvx(bodyNumber)= angvx_combined(groupNumber)
         angvy(bodyNumber)= angvy_combined(groupNumber)
         angvz(bodyNumber)= angvz_combined(groupNumber)
         
         temp_angvx  = 0.5_CGREAL*(angvx_combined_old(groupNumber)+angvx_combined(groupNumber))
         temp_angvy  = 0.5_CGREAL*(angvy_combined_old(groupNumber)+angvy_combined(groupNumber))
         temp_angvz  = 0.5_CGREAL*(angvz_combined_old(groupNumber)+angvz_combined(groupNumber))            
    
               
         !Velocity for different body center is the combination of Velocity of combined group's center plus relative rotation (Right Hand Side Coordinate)

         temp_vxcent_pre = temp_angvz*(ycent_combined(groupNumber)-ycent(bodyNumber))    - &
                           temp_angvy*(zcent_combined(groupNumber)-zcent(bodyNumber)) 

         temp_vycent_pre = temp_angvx*(zcent_combined(groupNumber)-zcent(bodyNumber))    - &
                           temp_angvz*(xcent_combined(groupNumber)-xcent(bodyNumber))   

         temp_vzcent_pre = temp_angvy*(xcent_combined(groupNumber)-xcent(bodyNumber))    - &
                           temp_angvx*(ycent_combined(groupNumber)-ycent(bodyNumber))
                            


       ENDIF ! Judge dimension


         temp_vcent_magni_pre = sqrt(temp_vxcent_pre**2 + temp_vycent_pre**2 + temp_vzcent_pre**2)


         temp_vxcent_new = temp_vxcent_pre*R11(groupNumber) + temp_vycent_pre*R12(groupNumber) + temp_vzcent_pre*R13(groupNumber)

         temp_vycent_new = temp_vxcent_pre*R21(groupNumber) + temp_vycent_pre*R22(groupNumber) + temp_vzcent_pre*R23(groupNumber)
		
         temp_vzcent_new = temp_vxcent_pre*R31(groupNumber) + temp_vycent_pre*R32(groupNumber) + temp_vzcent_pre*R33(groupNumber)

         average_vxcent = temp_vxcent_pre + temp_vxcent_new

         average_vycent = temp_vycent_pre + temp_vycent_new
		
         average_vzcent = temp_vzcent_pre + temp_vzcent_new

		 temp_vcent_magni_new = sqrt(average_vxcent**2 + average_vycent**2 + average_vzcent**2)

      IF(temp_vcent_magni_new > 0.0_CGREAL)THEN
      
        vxcent(bodyNumber) = vxcent_combined(groupNumber) + average_vxcent * temp_vcent_magni_pre * sin(0.5_CGREAL*theta_group(groupNumber))/(0.5_CGREAL*theta_group(groupNumber)) / temp_vcent_magni_new

        vycent(bodyNumber) = vycent_combined(groupNumber) + average_vycent * temp_vcent_magni_pre * sin(0.5_CGREAL*theta_group(groupNumber))/(0.5_CGREAL*theta_group(groupNumber)) / temp_vcent_magni_new
		
        vzcent(bodyNumber) = vzcent_combined(groupNumber) + average_vzcent * temp_vcent_magni_pre * sin(0.5_CGREAL*theta_group(groupNumber))/(0.5_CGREAL*theta_group(groupNumber)) / temp_vcent_magni_new 
        
      ELSE
                         
        vxcent(bodyNumber) = vxcent_combined(groupNumber) 

        vycent(bodyNumber) = vycent_combined(groupNumber) 
		
        vzcent(bodyNumber) = vzcent_combined(groupNumber) 
                                                       
      ENDIF
!hahahaha

                                             
     ENDIF ! group number
	  
    ENDDO !bodyNumber

   ENDDO !groupNumber
   

END SUBROUTINE combined_body_motion






SUBROUTINE MofI_CofG_combined(groupNumber) 

USE global_parameters
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE usr_module

Implicit None

INTEGER :: i,j,k, bodyNumber
INTEGER, INTENT(IN) :: groupNumber

   non_dim_mass_combined = 0.0_CGREAL

   I_XX_COMBINED(groupNumber) = 0.0_CGREAL
   I_YY_COMBINED(groupNumber) = 0.0_CGREAL
   I_ZZ_COMBINED(groupNumber) = 0.0_CGREAL
   I_XY_COMBINED(groupNumber) = 0.0_CGREAL
   I_YZ_COMBINED(groupNumber) = 0.0_CGREAL
   I_XZ_COMBINED(groupNumber) = 0.0_CGREAL

   xcent_combined(groupNumber) = 0.0_CGREAL
   ycent_combined(groupNumber) = 0.0_CGREAL
   zcent_combined(groupNumber) = 0.0_CGREAL

!  -----------------------------Calculate total mass of combined bodies
   DO bodyNumber = 1,  nbody
   


     IF(combined_Group_index(bodyNumber) == groupNumber)THEN

     CALL MofI_CofG_par(bodyNumber)
     
     IF ( body_dim(bodyNumber) == 2 .OR. membrane_type(bodyNumber) == 1 ) THEN
     non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**2*zout )
     ELSE
     non_dim_volume(bodyNumber) = volume(bodyNumber) / ( lScale**3 )    
     ENDIF

     non_dim_mass_combined(groupNumber)  = non_dim_mass_combined(groupNumber)  + non_dim_mass(bodyNumber)

	 ENDIF
     
	                              
   ENDDO !bodyNumber

!  -----------------------------Calculate barycenter of combined bodies
   DO bodyNumber = 1,  nbody

     IF(combined_Group_index(bodyNumber) == groupNumber)THEN

      xcent_combined(groupNumber) = xcent_combined(groupNumber) + xcent(bodyNumber)*non_dim_mass(bodyNumber)
      ycent_combined(groupNumber) = ycent_combined(groupNumber) + ycent(bodyNumber)*non_dim_mass(bodyNumber)
      zcent_combined(groupNumber) = zcent_combined(groupNumber) + zcent(bodyNumber)*non_dim_mass(bodyNumber)

	 ENDIF

   ENDDO !bodyNumber
   
   xcent_combined(groupNumber) = xcent_combined(groupNumber)/non_dim_mass_combined(groupNumber)
   ycent_combined(groupNumber) = ycent_combined(groupNumber)/non_dim_mass_combined(groupNumber)
   zcent_combined(groupNumber) = zcent_combined(groupNumber)/non_dim_mass_combined(groupNumber)


!  -----------------------------Calculate moment of inertia for combined bodies(Parallel Axis  Theorem)
   DO bodyNumber = 1,  nbody

     IF(combined_Group_index(bodyNumber) == groupNumber)THEN

	      IF ( body_dim(bodyNumber) == 3 ) THEN   


	      I_XX_COMBINED(groupNumber) = I_XX_COMBINED(groupNumber) + I_XX(bodyNumber) +                                                &
	                                   density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
					                   ( (ycent(bodyNumber)-ycent_combined(groupNumber))**2 + (zcent(bodyNumber)-zcent_combined(groupNumber))**2 ) 

          I_YY_COMBINED(groupNumber) = I_YY_COMBINED(groupNumber) + I_YY(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
					                   ( (xcent(bodyNumber)-xcent_combined(groupNumber))**2 + (zcent(bodyNumber)-zcent_combined(groupNumber))**2 ) 

          I_ZZ_COMBINED(groupNumber) = I_ZZ_COMBINED(groupNumber) + I_ZZ(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
		                              ( (xcent(bodyNumber)-xcent_combined(groupNumber))**2 + (ycent(bodyNumber)-ycent_combined(groupNumber))**2 )

          I_XY_COMBINED(groupNumber) = I_XY_COMBINED(groupNumber) + I_XY(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
					                   ( (xcent(bodyNumber)-xcent_combined(groupNumber)) * (ycent(bodyNumber)-ycent_combined(groupNumber)) )

          I_YZ_COMBINED(groupNumber) = I_YZ_COMBINED(groupNumber) + I_YZ(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
					                   ( (ycent(bodyNumber)-ycent_combined(groupNumber)) * (zcent(bodyNumber)-zcent_combined(groupNumber)) )

          I_XZ_COMBINED(groupNumber) = I_XZ_COMBINED(groupNumber) + I_XZ(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
					                   ( (xcent(bodyNumber)-xcent_combined(groupNumber)) * (zcent(bodyNumber)-zcent_combined(groupNumber)) )
		  ELSE !==>body_dim(bodyNumber) == 2
          I_ZZ_COMBINED(groupNumber) = I_ZZ_COMBINED(groupNumber) + I_ZZ(bodyNumber) +                                                &
                                       density_solid(bodyNumber)*volume(bodyNumber)*                                                  &
	                                   ( (xcent(bodyNumber)-xcent_combined(groupNumber))**2 + (ycent(bodyNumber)-ycent_combined(groupNumber))**2 )
          ENDIF !dimension=2 or 3 
	 ENDIF
   
   ENDDO !bodyNumber
!-----------------------For 2D case, there is only I_ZZ left
   IF ( body_dim(bodyNumber) == 2 ) THEN   ! Area Moment of Inertia

        I_XX_COMBINED(groupNumber) =0.0_CGREAL
        I_YY_COMBINED(groupNumber) =0.0_CGREAL
   
        I_XY_COMBINED(groupNumber) =0.0_CGREAL
        I_YZ_COMBINED(groupNumber) =0.0_CGREAL
        I_XZ_COMBINED(groupNumber) =0.0_CGREAL

   ENDIF



   nonDimM_I_combined(1,1,groupNumber) = 2.0_CGREAL *I_XX_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(1,2,groupNumber) = 2.0_CGREAL *I_XY_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(1,3,groupNumber) = 2.0_CGREAL *I_XZ_COMBINED(groupNumber)/(density_fluid*lScale**5)

   nonDimM_I_combined(2,1,groupNumber) = 2.0_CGREAL *I_XY_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(2,2,groupNumber) = 2.0_CGREAL *I_YY_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(2,3,groupNumber) = 2.0_CGREAL *I_YZ_COMBINED(groupNumber)/(density_fluid*lScale**5)

   nonDimM_I_combined(3,1,groupNumber) = 2.0_CGREAL *I_XZ_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(3,2,groupNumber) = 2.0_CGREAL *I_YZ_COMBINED(groupNumber)/(density_fluid*lScale**5)
   nonDimM_I_combined(3,3,groupNumber) = 2.0_CGREAL *I_ZZ_COMBINED(groupNumber)/(density_fluid*lScale**5)

END SUBROUTINE MofI_CofG_combined


!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
SUBROUTINE wave_motion(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi

    IMPLICIT NONE

    INTEGER,INTENT (IN) :: iBody

    INTEGER             :: m, midx ,k ,num, mm, ierr                                                   !Ehsan added for SJ (num,mm)



    REAL(KIND=CGREAL)   :: uBodyMarkerRel, vBodyMarkerRel, wBodyMarkerRel
    REAL(KIND=CGREAL),PARAMETER   :: body_length = 2.0_CGREAL
    REAL(KIND=CGREAL),PARAMETER   :: wave_speed = 40.0_CGREAL
    REAL(KIND=CGREAL),PARAMETER   :: amplitude = 0.5_CGREAL

    
      uBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL
      vBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL
      wBodyMarker(iBody,1:nPtsBodyMarker(iBody)) = 0.0_CGREAL

 
    DO m = 1,nPtsBodyMarker(iBody)                                            

      uBodyMarkerRel = 2.0_CGREAL*pi/body_length*wave_speed*amplitude*COS(  2.0_CGREAL*pi/body_length*wave_speed*time  &
                                                                         + (yBodyMarker(iBody,m)-3.d0)*2.0_CGREAL*pi/body_length  )

      vBodyMarkerRel = 0.0_CGREAL

      wBodyMarkerRel = 0.0_CGREAL

      uBodyMarker(iBody,m) = uBodyMarker(iBody,m) +  uBodyMarkerRel
      vBodyMarker(iBody,m) = vBodyMarker(iBody,m) +  vBodyMarkerRel
      wBodyMarker(iBody,m) = wBodyMarker(iBody,m) +  wBodyMarkerRel

    ENDDO                                                 ! Rupeshs additions end here
                                         
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
!   Extending velocity across span for 2D body
!   ------------------------------------------
!    IF ( body_dim(iBody) == BODY_DIM2 .AND. ntime > 0 ) THEN

!      DO k=2,nz
!        DO m=mMinWallVel(iBody),mMaxWallVel(iBody)
!          midx = (k-1)*nPtsBodymarker(iBody)/nz + m
!          uBodyMarker(iBody,midx) = uBodyMarker(iBody,m)
!          vBodyMarker(iBody,midx) = vBodyMarker(iBody,m)
!          wBodyMarker(iBody,midx) = wBodyMarker(iBody,m)
!        ENDDO
!      ENDDO

!    ENDIF

END SUBROUTINE wave_motion
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!---------------------------------------------------------------------

!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
SUBROUTINE fruitfly_motion(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi
    USE scalar

    IMPLICIT NONE

    INTEGER              :: i
    INTEGER,INTENT (IN) :: iBody
    REAL(KIND=CGREAL),DIMENSION(2:2500)    :: pitch1, theta1, alpha1    
    REAL(KIND=CGREAL),DIMENSION(1,2480) :: x0, y0, z0
    REAL(KIND=CGREAL),DIMENSION(1) :: angvalfa, angvphi, angvtheta 
    REAL(KIND=CGREAL),PARAMETER  :: span_length  = 0.00196141743341347535878915024604D0  
    !REAL(KIND=CGREAL),PARAMETER  :: chord_length = 1.6575121714079955912352512063492d-4
    !REAL(KIND=CGREAL),PARAMETER  :: product      = -8.0179364701387917888D-8     
    REAL(KIND=CGREAL),PARAMETER  :: chord_length = 2.9534143282534957375705769771415D-4 
    REAL(KIND=CGREAL),PARAMETER  :: product      = 4.1800035558432875022D-8
    REAL(KIND=CGREAL)  :: X_P_tip, Y_P_tip, Z_P_tip      
    REAL(KIND=CGREAL)  :: X_P_end, Y_P_end, Z_P_end            
    REAL(KIND=CGREAL)  :: dummy   
    REAL(KIND=CGREAL)  :: A, B, C           
    
                                     
!==========PHI  THETA  ALFA  prescribed motion========== 
!---Readin PHI  THETA  ALFA angles
IF(NTIME==2)THEN
x0 = xBodyMarker
y0 = yBodyMarker
z0 = zBodyMarker
ENDIF


OPEN(333,FILE='fruitfly_angles.dat',STATUS='UNKNOWN')
!OPEN(333,FILE='moth_angles.dat',STATUS='UNKNOWN')
READ(333,*)
DO i = 2, 2500
   READ(333,*) dummy, pitch1(i), theta1(i)
   READ(333,*) alpha1(i)   
ENDDO !i

CLOSE(333) 


IF(NTIME == 2)THEN
  DO i=1,nPtsBodyMarker(iBody)
    uBodyMarker(iBody,i) = ( xBodyMarker(iBody,72)  -                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME))*cos(pitch1(NTIME))                                                      +  & 
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME))*sin(pitch1(NTIME)) - cos(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME))) +  &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME))*sin(pitch1(NTIME)) + sin(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME)))  -  &
                          x0(iBody,i) )/DT
    vBodyMarker(iBody,i) = ( yBodyMarker(iBody,72) +                                                                                                         &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME))*cos(pitch1(NTIME)) - sin(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME))*cos(pitch1(NTIME)) + cos(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME))*sin(pitch1(NTIME)) - y0(iBody,i) )/DT
    wBodyMarker(iBody,i) = ( zBodyMarker(iBody,72) +                                             &
                          (-z0(iBody,i)+z0(iBody,72))*cos(alpha1(NTIME))*cos(theta1(NTIME)) -   &
                          (-x0(iBody,i)+x0(iBody,72))*sin(theta1(NTIME)) -                     &
                          (y0(iBody,i)-y0(iBody,72))*cos(theta1(NTIME))*sin(alpha1(NTIME)) - z0(iBody,i) )/DT
  ENDDO !i
ELSE
  DO i=1,nPtsBodyMarker(iBody)
    uBodyMarker(iBody,i) = ( xBodyMarker(iBody,72)  -                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME))*cos(pitch1(NTIME))                                                      +  & 
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME))*sin(pitch1(NTIME)) - cos(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME))) +  &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME))*sin(pitch1(NTIME)) + sin(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME)))  -  &
                           xBodyMarker(iBody,72)  +                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME-1))*cos(pitch1(NTIME-1))                                                      -  & 
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME-1))*sin(pitch1(NTIME-1)) - cos(alpha1(NTIME-1))*cos(pitch1(NTIME-1))*sin(theta1(NTIME-1))) -  &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME-1))*sin(pitch1(NTIME-1)) + sin(alpha1(NTIME-1))*cos(pitch1(NTIME-1))*sin(theta1(NTIME-1))) )/DT                          
    vBodyMarker(iBody,i) = ( yBodyMarker(iBody,72) +                                                                                                         &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME))*cos(pitch1(NTIME)) - sin(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME))*cos(pitch1(NTIME)) + cos(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME))*sin(pitch1(NTIME))  -  &
                          yBodyMarker(iBody,72) -                                                                                                         &
                          (y0(iBody,i)-y0(iBody,72))*(cos(alpha1(NTIME-1))*cos(pitch1(NTIME-1)) - sin(alpha1(NTIME-1))*sin(theta1(NTIME-1))*sin(pitch1(NTIME-1))) -  &
                          (-z0(iBody,i)+z0(iBody,72))*(sin(alpha1(NTIME-1))*cos(pitch1(NTIME-1)) + cos(alpha1(NTIME-1))*sin(theta1(NTIME-1))*sin(pitch1(NTIME-1))) -  &
                          (-x0(iBody,i)+x0(iBody,72))*cos(theta1(NTIME-1))*sin(pitch1(NTIME-1))  )/DT                        
    wBodyMarker(iBody,i) = ( zBodyMarker(iBody,72) +                                             &
                          (-z0(iBody,i)+z0(iBody,72))*cos(alpha1(NTIME))*cos(theta1(NTIME)) -   &
                          (-x0(iBody,i)+x0(iBody,72))*sin(theta1(NTIME)) -                     &
                          (y0(iBody,i)-y0(iBody,72))*cos(theta1(NTIME))*sin(alpha1(NTIME)) -  &
                          zBodyMarker(iBody,72) -                                             &
                          (-z0(iBody,i)+z0(iBody,72))*cos(alpha1(NTIME-1))*cos(theta1(NTIME-1)) +   &
                          (-x0(iBody,i)+x0(iBody,72))*sin(theta1(NTIME-1)) +                     &
                          (y0(iBody,i)-y0(iBody,72))*cos(theta1(NTIME-1))*sin(alpha1(NTIME-1))  )/DT
  ENDDO !i
ENDIF
        
IF(IMTHEBOSS)THEN           
OPEN(666,FILE='Trajectory_Tip.dat',POSITION='APPEND')
OPEN(888,FILE='Trajectory_TE.dat',POSITION='APPEND')
OPEN(999,FILE='fort.77',POSITION='APPEND')
WRITE(666,*) xBodyMarker(iBody,4),yBodyMarker(iBody,4),zBodyMarker(iBody,4)
WRITE(888,*) xBodyMarker(iBody,127),yBodyMarker(iBody,127),zBodyMarker(iBody,127)
WRITE(999,*) dt, time-dt, nPtsBodyMarker(1)
DO i=1,nPtsBodyMarker(iBody)
   IF(NTIME == 2)THEN
    WRITE(999,*) 0.D0, 0.D0, 0.D0
   ELSE
    WRITE(999,*) uBodyMarker(iBody,i), vBodyMarker(iBody,i), wBodyMarker(iBody,i)
   ENDIF
ENDDO       
close(666)
close(888)
close(999)                                 

ENDIF !IMTHEBOSS

!X_P_tip = X_P_fixed - span_length*COS(theta1(ntime))*COS( pitch1(1)-pitch1(ntime) )
!Y_P_tip = Y_P_fixed - span_length*COS(theta1(ntime))*SIN( pitch1(1)-pitch1(ntime) )
!Z_P_tip = Z_P_fixed + span_length*SIN(theta1(ntime))


!Z_P_end = REAL(Z_P_tip - chord_length*COS( alpha1(ntime) ),KIND=CGREAL)
!A = 1.D0 + (Y_P_fixed-Y_P_tip)**2.d0/(X_P_fixed-X_P_tip)**2.d0
!B = -2.D0*(Y_P_fixed-Y_P_tip)*(product- (Z_P_end-Z_P_tip)*(Z_P_fixed-Z_P_tip) )/(X_P_fixed-X_P_tip)**2.d0
!C = (Z_P_end-Z_P_tip)**2.d0 - chord_length**2.d0 + ( product - (Z_P_end-Z_P_tip)*(Z_P_fixed-Z_P_tip) )**2.d0/(X_P_fixed-X_P_tip)**2.d0
!IF( alpha1(ntime)-alpha1(1) > 0.D0 )THEN
!  Y_P_end = Y_P_tip + ( -B + SQRT(B**2-4.D0*A*C))/2.D0/A
!ELSE
!  Y_P_end = Y_P_tip + ( -B - SQRT(B**2-4.D0*A*C))/2.D0/A
!ENDIF
!X_P_end = X_P_tip + ( product - (Z_P_end-Z_P_tip)*(Z_P_fixed-Z_P_tip) )/(X_P_fixed-X_P_tip) -  &
!         (Y_P_end-Y_P_tip)*(Y_P_fixed-Y_P_tip)/(X_P_fixed-X_P_tip)

!IF(IMTHEBOSS)THEN         
!   OPEN(361,FILE='dot.dat',POSITION='APPEND')
!   WRITE(361,*) '========================' 
!   WRITE(361,*) ntime  
!   WRITE(361,*) X_P_tip, Y_P_tip, Z_P_tip
!   WRITE(361,*) X_P_end, Y_P_end, Z_P_end
!   WRITE(361,*) A, B, C, B**2-4.D0*A*C   
!   WRITE(361,*) X_P_fixed, Y_P_fixed, Z_P_fixed   
!ENDIF

!==========PHI  THETA  ALFA  prescribed motion==========    
!-------Sane & Dickinson's Model
!    IF(NTIME <= 400)THEN
!       angvphi(iBody) = pi/2.5D-3
!    ELSEIF(NTIME > 400 .AND. NTIME <= 800)THEN
!       angvphi(iBody) = -pi/2.5D-3
!    ENDIF   
!-------Hedrick's Model
!    angvphi(iBody) = (2.D0*pi**2/1.5D-2)*sin(2.D0*PI*ntime/800)    
          
    
!    angvtheta(iBody) = (60.D0/180.D0)*pi**2*cos(2.D0*pi*time/5.D-3) / 5.D-3  
    
!    IF(NTIME <= 200)THEN
!       angvalfa(iBody) = 0.0_CGREAL
!    ELSEIF(NTIME > 200 .AND. NTIME <= 272)THEN
!       angvalfa(iBody) = -(pi/2.D0)/(0.16*5.0D-3) * (NTIME-200)/72.D0           
!    ELSEIF(NTIME > 272 .AND. NTIME <= 400)THEN  !flip
!       angvalfa(iBody) = -(pi/2.D0)/(0.16*5.0D-3)
!    ELSEIF(NTIME > 400 .AND. NTIME <= 600)THEN
!       angvalfa(iBody) = 0.0_CGREAL 
!    ELSEIF(NTIME > 600 .AND. NTIME <= 672)THEN
!       angvalfa(iBody) = (pi/2.D0)/(0.16*5.0D-3) * (NTIME-600)/72.D0    
!    ELSEIF(NTIME > 672 .AND. NTIME <= 800)THEN  !flip
!       angvalfa(iBody) = (pi/2.D0)/(0.16*5.0D-3)                 
!    ENDIF      
    
!    angle(iBody)= angle(iBody) + dt * angvphi(iBody)
!==========Calculate Angular Velocity in XYZ Frame==========    
!    angvx(iBody)= cos(-angle(iBody))*angvalfa(iBody)  - sin(-angle(iBody))*angvtheta(iBody) 
!    angvy(iBody)= angvphi(iBody)    
!    angvz(iBody)= cos(-angle(iBody))*angvtheta(iBody) + sin(-angle(iBody))*angvalfa(iBody)            

END SUBROUTINE fruitfly_motion


!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
SUBROUTINE Dickinson_motion(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi
    USE scalar

    IMPLICIT NONE

    INTEGER              :: i
    INTEGER,INTENT (IN) :: iBody

          
vxcent(iBody) = 0.272D0
IF(NTIME< 50 )THEN
   angvz(iBody) = (NTIME-0.D0)*1.5/50.D0 
ELSEIF(NTIME> 50 )THEN
   angvz(iBody) = 1.5D0    
ENDIF  

    
                                

     

END SUBROUTINE Dickinson_motion



SUBROUTINE fruitfly_hinge(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi
    USE scalar

    IMPLICIT NONE

    INTEGER,INTENT (IN) :: iBody
    INTEGER              :: i
    REAL(KIND=CGREAL) :: uBodyMarker_fixed, vBodyMarker_fixed, wBodyMarker_fixed 
    
    uBodyMarker_fixed = uBodyMarker(iBody,i_fixed(iBody))
    vBodyMarker_fixed = vBodyMarker(iBody,i_fixed(iBody))
    wBodyMarker_fixed = wBodyMarker(iBody,i_fixed(iBody))
                                                    
    DO i=1,nPtsBodyMarker(iBody)
        uBodyMarker(iBody,i) = uBodyMarker(iBody,i) - uBodyMarker_fixed
        vBodyMarker(iBody,i) = vBodyMarker(iBody,i) - vBodyMarker_fixed
        wBodyMarker(iBody,i) = wBodyMarker(iBody,i) - wBodyMarker_fixed           
    ENDDO ! i    

END SUBROUTINE fruitfly_hinge
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd


!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
SUBROUTINE moth_motion(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi
    USE scalar

    IMPLICIT NONE

    INTEGER              :: i
    INTEGER,INTENT (IN) :: iBody
    REAL(KIND=CGREAL),DIMENSION(2:4200)    :: pitch1, theta1, alpha1    
    REAL(KIND=CGREAL),DIMENSION(1,28861) :: x0, y0, z0
    REAL(KIND=CGREAL),DIMENSION(1) :: angvalfa, angvphi, angvtheta 
    REAL(KIND=CGREAL),PARAMETER  :: span_length  = 0.00196141743341347535878915024604D0  
    !REAL(KIND=CGREAL),PARAMETER  :: chord_length = 1.6575121714079955912352512063492d-4
    !REAL(KIND=CGREAL),PARAMETER  :: product      = -8.0179364701387917888D-8     
    REAL(KIND=CGREAL),PARAMETER  :: chord_length = 2.9534143282534957375705769771415D-4 
    REAL(KIND=CGREAL),PARAMETER  :: product      = 4.1800035558432875022D-8
    REAL(KIND=CGREAL)  :: X_P_tip, Y_P_tip, Z_P_tip      
    REAL(KIND=CGREAL)  :: X_P_end, Y_P_end, Z_P_end            
    REAL(KIND=CGREAL)  :: dummy   
    REAL(KIND=CGREAL)  :: A, B, C           
    
                                     
!==========PHI  THETA  ALFA  prescribed motion========== 
!---Readin PHI  THETA  ALFA angles
IF(NTIME==2)THEN
x0 = xBodyMarker
y0 = yBodyMarker
z0 = zBodyMarker
ENDIF


OPEN(333,FILE='moth_angles.dat',STATUS='UNKNOWN')
!OPEN(333,FILE='fruitfly_angles.dat',STATUS='UNKNOWN')
READ(333,*)
DO i = 2, 4200
   READ(333,*) dummy, pitch1(i), theta1(i)
   READ(333,*) alpha1(i)   
ENDDO !i

CLOSE(333) 


IF(NTIME == 2)THEN
  DO i=1,nPtsBodyMarker(iBody)
    uBodyMarker(iBody,i) = ( xBodyMarker(iBody,2)  -                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME))*cos(pitch1(NTIME))                                                      +  & 
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME))*sin(pitch1(NTIME)) - cos(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME))) +  &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME))*sin(pitch1(NTIME)) + sin(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME)))  -  &
                          x0(iBody,i) )/DT
    vBodyMarker(iBody,i) = ( yBodyMarker(iBody,2) +                                                                                                         &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME))*cos(pitch1(NTIME)) - sin(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME))*cos(pitch1(NTIME)) + cos(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME))*sin(pitch1(NTIME)) - y0(iBody,i) )/DT
    wBodyMarker(iBody,i) = ( zBodyMarker(iBody,2) +                                             &
                          (-z0(iBody,i)+z0(iBody,2))*cos(alpha1(NTIME))*cos(theta1(NTIME)) -   &
                          (-x0(iBody,i)+x0(iBody,2))*sin(theta1(NTIME)) -                     &
                          (y0(iBody,i)-y0(iBody,2))*cos(theta1(NTIME))*sin(alpha1(NTIME)) - z0(iBody,i) )/DT
  ENDDO !i
ELSE
  DO i=1,nPtsBodyMarker(iBody)
    uBodyMarker(iBody,i) = ( xBodyMarker(iBody,2)  -                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME))*cos(pitch1(NTIME))                                                      +  & 
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME))*sin(pitch1(NTIME)) - cos(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME))) +  &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME))*sin(pitch1(NTIME)) + sin(alpha1(NTIME))*cos(pitch1(NTIME))*sin(theta1(NTIME)))  -  &
                           xBodyMarker(iBody,2)  +                                                                                                 &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME-1))*cos(pitch1(NTIME-1))                                                      -  & 
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME-1))*sin(pitch1(NTIME-1)) - cos(alpha1(NTIME-1))*cos(pitch1(NTIME-1))*sin(theta1(NTIME-1))) -  &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME-1))*sin(pitch1(NTIME-1)) + sin(alpha1(NTIME-1))*cos(pitch1(NTIME-1))*sin(theta1(NTIME-1))) )/DT                          
    vBodyMarker(iBody,i) = ( yBodyMarker(iBody,2) +                                                                                                         &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME))*cos(pitch1(NTIME)) - sin(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME))*cos(pitch1(NTIME)) + cos(alpha1(NTIME))*sin(theta1(NTIME))*sin(pitch1(NTIME))) +  &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME))*sin(pitch1(NTIME))  -  &
                          yBodyMarker(iBody,2) -                                                                                                         &
                          (y0(iBody,i)-y0(iBody,2))*(cos(alpha1(NTIME-1))*cos(pitch1(NTIME-1)) - sin(alpha1(NTIME-1))*sin(theta1(NTIME-1))*sin(pitch1(NTIME-1))) -  &
                          (-z0(iBody,i)+z0(iBody,2))*(sin(alpha1(NTIME-1))*cos(pitch1(NTIME-1)) + cos(alpha1(NTIME-1))*sin(theta1(NTIME-1))*sin(pitch1(NTIME-1))) -  &
                          (-x0(iBody,i)+x0(iBody,2))*cos(theta1(NTIME-1))*sin(pitch1(NTIME-1))  )/DT                        
    wBodyMarker(iBody,i) = ( zBodyMarker(iBody,2) +                                             &
                          (-z0(iBody,i)+z0(iBody,2))*cos(alpha1(NTIME))*cos(theta1(NTIME)) -   &
                          (-x0(iBody,i)+x0(iBody,2))*sin(theta1(NTIME)) -                     &
                          (y0(iBody,i)-y0(iBody,2))*cos(theta1(NTIME))*sin(alpha1(NTIME)) -  &
                          zBodyMarker(iBody,2) -                                             &
                          (-z0(iBody,i)+z0(iBody,2))*cos(alpha1(NTIME-1))*cos(theta1(NTIME-1)) +   &
                          (-x0(iBody,i)+x0(iBody,2))*sin(theta1(NTIME-1)) +                     &
                          (y0(iBody,i)-y0(iBody,2))*cos(theta1(NTIME-1))*sin(alpha1(NTIME-1))  )/DT
  ENDDO !i
ENDIF
        
IF(IMTHEBOSS)THEN          
OPEN(666,FILE='Trajectory_Tip.dat',POSITION='APPEND')
OPEN(888,FILE='Trajectory_TE.dat',POSITION='APPEND')
WRITE(666,*) xBodyMarker(iBody,367),yBodyMarker(iBody,367),zBodyMarker(iBody,367)
WRITE(888,*) xBodyMarker(iBody,505),yBodyMarker(iBody,505),zBodyMarker(iBody,505)
OPEN(999,FILE='fort.77',POSITION='APPEND')
WRITE(999,*) dt, time-dt, nPtsBodyMarker(1)
DO i=1,nPtsBodyMarker(iBody)
   IF(NTIME == 2)THEN
    WRITE(999,*) 0.D0, 0.D0, 0.D0
   ELSE
    WRITE(999,*) uBodyMarker(iBody,i), vBodyMarker(iBody,i), wBodyMarker(iBody,i)
   ENDIF
ENDDO     
close(999)  
close(666)
close(888)
ENDIF !IMTHEBOSS     

END SUBROUTINE moth_motion


SUBROUTINE moth_hinge(iBody)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE mpi
    USE scalar

    IMPLICIT NONE

    INTEGER,INTENT (IN) :: iBody
    INTEGER              :: i
    REAL(KIND=CGREAL) :: uBodyMarker_fixed, vBodyMarker_fixed, wBodyMarker_fixed 
    
    uBodyMarker_fixed = uBodyMarker(iBody,i_fixed(iBody))
    vBodyMarker_fixed = vBodyMarker(iBody,i_fixed(iBody))
    wBodyMarker_fixed = wBodyMarker(iBody,i_fixed(iBody))
                                                    
    DO i=1,nPtsBodyMarker(iBody)
        uBodyMarker(iBody,i) = uBodyMarker(iBody,i) - uBodyMarker_fixed
        vBodyMarker(iBody,i) = vBodyMarker(iBody,i) - vBodyMarker_fixed
        wBodyMarker(iBody,i) = wBodyMarker(iBody,i) - wBodyMarker_fixed           
    ENDDO ! i    

END SUBROUTINE moth_hinge
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!---------------------------------------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
