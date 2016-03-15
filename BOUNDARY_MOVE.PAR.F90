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
!  Filename: BOUNDARY_MOVE.PAR.F90
!  Latest Modification: September 11, 2008 (ver. P1.2.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     read_marker_vel_check(iBody,ntime_check)
!     move_boundary()
!     compute_marker_vel(iBody)
!     read_marker_vel(iBody)
! --------------------------------------------------------------------




SUBROUTINE read_marker_vel_check(iBody, ntime_check)
 
    USE global_parameters
    USE flow_parameters
 
    IMPLICIT NONE
 
    INTEGER, PARAMETER  :: ifort=64
    
    INTEGER,INTENT(IN)  :: iBody
    INTEGER,INTENT(OUT) :: ntime_check 

    INTEGER             :: m, ifortCent, nptsbm_temp

    REAL(KIND=CGREAL)   :: temp_time, temp_dt
    REAL(KIND=CGREAL)   :: uBodyMarker_check, vBodyMarker_check, wBodyMarker_check

    INTEGER             :: istat

    CHARACTER*20 :: fname1



    ifortCent = ifort+iBody
    IF (ifortCent > 80) PRINT*,'Might have problems with file numbers'
 
!    WRITE(fname1,"('body_mesh_velocity.',I3.3,'.dat')") iBody
!    OPEN(UNIT=ifortCent,FILE=fname1)

    ntime_check = 0
 
    DO
      READ(ifortCent,*, IOSTAT=istat) temp_dt, temp_time, nptsbm_temp
      IF (istat <0) exit
     
      DO m=1,nptsbm_temp
        READ(ifortCent,*) uBodyMarker_check, vBodyMarker_check, wBodyMarker_check
      END DO ! m 

      ntime_check = ntime_check + 1
    END DO

    IF (MonitorON) WRITE(STDOUT,'(5X,A,I5,A,I5)') 'Number of check points per period for body ', iBody , ' is ', ntime_check  
    
    CLOSE(ifortCent) 
!    OPEN(UNIT=ifortCent,FILE=fname1)

  END SUBROUTINE read_marker_vel_check
!------------------------------------------------------------------------------



SUBROUTINE move_boundary()

! ------------------------------------------------------------------------------------
!  Move body markers and centroids and set new velocity of marker points and centroid 
! ------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
    USE tahoe_parameters
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)       
    USE usr_module 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
    IMPLICIT NONE

    REAL(KIND=CGREAL) :: myCommTime
!dddddddddddddddddddddddddddddddd
    CHARACTER*50 :: fname1
!dddddddddddddddddddddddddddddddd 

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER :: i,ifortCent,iGroup,iBody,iFort   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)


    angvx_old(1:nbody) = angvx(1:nbody)
    angvy_old(1:nbody) = angvy(1:nbody)
    angvz_old(1:nbody) = angvz(1:nbody)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
IF(nGroup_Combined /= 0)THEN 
    angvx_combined_old(1:nGroup_Combined) = angvx_combined(1:nGroup_Combined)
    angvy_combined_old(1:nGroup_Combined) = angvy_combined(1:nGroup_Combined)
    angvz_combined_old(1:nGroup_Combined) = angvz_combined(1:nGroup_Combined) 
ENDIF     

!===Update pendulum force first
    IF( ntime > nstart_FIM )THEN
     DO iBody=1,nBody	
	   IF(boundary_motion_type(iBody) == PENDULUM)THEN
	     CALL ADD_PENDULUM_INTERNALFORCE(iBody)
	   ENDIF
	 ENDDO
	ENDIF
	
	IF(nGroup_Combined == 0)THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)      

!   Update body velocity
!   --------------------
    DO iBody=1,nBody
      SELECT CASE (boundary_motion_type(iBody))
      CASE (FORCED)
        CALL forced_motion(iBody)
        CALL compute_marker_vel(iBody)            
      CASE (FLOW_INDUCED)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)     
	    IF( ntime > nstart_FIM)THEN    	      
         CALL flow_induced_motion(iBody)       
         CALL compute_marker_vel(iBody)
        ENDIF
      CASE (PENDULUM)           
	    IF( ntime > nstart_FIM)THEN      
         CALL flow_induced_motion(iBody)
         CALL compute_marker_vel(iBody)
        ENDIF
!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD                                 
      CASE (18) !work for fruitfly motion model
         CALL fruitfly_motion(iBody)
         CALL fruitfly_hinge(iBody)
      CASE (19) !work for fruitfly motion model
         CALL moth_motion(iBody)
         CALL moth_hinge(iBody)  
      CASE (20) !work for fruitfly motion model
         CALL Dickinson_motion(iBody)   
         CALL compute_marker_vel(iBody)           
!DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
!        WRITE(*,*) 'Flow Induced Motion is not available in this version.'
!        CALL flow_stop
!        STOP
!!        CALL flow_induced_motion(iBody)
!!        CALL compute_marker_vel(iBody)
      CASE (PRESCRIBED)
        CALL read_marker_vel(iBody)      
      CASE(FEA_FLOW_STRUC_INTERACTION)                                   ! added for FEA ! SER_TO_PAR. QX. CH11
        CALL fea_flow_structure(iBody)
	  CASE(FEA_TAHOE)   
        !CALL markers_positions_tahoe(iBody)                                             ! Added by Rajneesh

# ifdef TAHOE
        CALL update_markers_positions_velo_tahoe(iBody)                       ! Updates bodymarkers positions
# endif
      END SELECT
    ENDDO
         
    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ELSE !===> nGoup_Combined /= 0      
     
!   Update body velocity
!   --------------------
    DO iBody=1,nBody
    
    IF(combined_type(iBody) == UNCOMBINED)THEN
	
      SELECT CASE (boundary_motion_type(iBody))
      CASE (FORCED)
        CALL forced_motion(iBody)
        CALL compute_marker_vel(iBody)
	  CASE (FLOW_INDUCED)
	    IF( ntime > nstart_FIM)THEN
        CALL flow_induced_motion(iBody)
        CALL compute_marker_vel(iBody)
        ENDIF
	  CASE (PENDULUM)        
	    IF( ntime > nstart_FIM)THEN    		          
         CALL flow_induced_motion(iBody)
         CALL compute_marker_vel(iBody)      
        ENDIF        
      CASE (PRESCRIBED) 
	    CALL read_marker_vel(iBody)
      CASE(FEA_FLOW_STRUC_INTERACTION)                                   ! added for FEA ! SER_TO_PAR. QX. CH11
        CALL fea_flow_structure(iBody)
      CASE(FEA_TAHOE)   
        !CALL markers_positions_tahoe(iBody)                                             ! Added by Rajneesh

# ifdef TAHOE
        CALL update_markers_positions_velo_tahoe(iBody)                       ! Updates bodymarkers positions
# endif  
      END SELECT
!----------------------------------------------------------COMBINED MOTION 
    ENDIF !combined_type(iBody) == UNCOMBINED

    ENDDO !iBody
    
      
    
!----------------Calculate velocity for different combined group(Update angv)       
    IF(ntime > nstart_FIM)THEN
    
	  CALL combined_body_motion()	
	    
	      
    ENDIF
    
!   --------------------Update body velocity for forced motion & rotate velocity vectors
    DO iBody=1,nBody	
	  IF(boundary_motion_type(iBody) == FORCED)THEN	
	    CALL forced_motion(iBody)
      ENDIF		
    ENDDO !iBody     
    
	
!   --------------------Update body velocity for prescribed motion & rotate velocity vectors
    DO iBody=1,nBody	
	  IF(boundary_motion_type(iBody) == PRESCRIBED)THEN	
	    CALL read_marker_vel(iBody)
      ENDIF		
    ENDDO !iBody      
    
!   --------------------Update velocity for every marker points
    DO iGroup = 1, nGroup_Combined
     DO iBody  = 1, nbody
	   IF(combined_Group_index(iBody) == iGroup)THEN
	     IF( ntime > nstart_FIM)THEN
         CALL compute_marker_vel(iBody)  
         ENDIF  
	   ENDIF
	 ENDDO
    ENDDO

	ENDIF !nGoup_Combined == 0
	
		
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)        

! ============================
!      n+1     n
!    x    -  x                
!    - i     - i        n+1
!   ------------   =  v   
!         dt          - i 
! ============================

!   ------------------------------------
!    Update location of boundary points
!   ------------------------------------

!   Update centroid location
!   ------------------------
    DO iBody=1,nBody

      IF (boundary_motion_type(iBody) /= STATIONARY) THEN

        SELECT CASE (body_dim(iBody))
          CASE (BODY_DIM2)
           IF (boundary_motion_type(iBody)/=FEA_FLOW_STRUC_INTERACTION & ! added for FEA
		       .AND.boundary_motion_type(iBody)/=FEA_TAHOE &
               .AND. boundary_motion_type(iBody)/=EIGEN_MOTION) THEN     ! SER_TO_PAR. QX. CH12



            DO i=1,nPtsBodyMarker(iBody)
              xBodyMarker(iBody,i) = xBodyMarker(iBody,i) + dt*uBodyMarker(iBody,i) 
              yBodyMarker(iBody,i) = yBodyMarker(iBody,i) + dt*vBodyMarker(iBody,i) 
              zBodyMarker(iBody,i) = zBodyMarker(iBody,i)                            
            ENDDO ! i
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)         
!=========Update the hinge location for 2D             
           hinge_x(iBody)  =  hinge_x(iBody) + dt*hinge_vx(iBody)
           hinge_y(iBody)  =  hinge_y(iBody) + dt*hinge_vy(iBody)
           hinge_z(iBody)  =  hinge_z(iBody)
                        

  IF(boundary_motion_type(iBody) == FORCED .AND. nGroup_Combined == 0)THEN
            x_rot_cent(iBody) = x_rot_cent(iBody) + dt*vxcent(iBody)
            y_rot_cent(iBody) = y_rot_cent(iBody) + dt*vycent(iBody)
            z_rot_cent(iBody) = z_rot_cent(iBody)  
            
  ELSE            
            xcent(iBody) = xcent(iBody) + dt*vxcent(iBody)
            ycent(iBody) = ycent(iBody) + dt*vycent(iBody)
            zcent(iBody) = zcent(iBody)            
  ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
!           updates the angle of attack for rotating bodies
!           -----------------------------------------------
            alpha(iBody)    = alpha(iBody) + dt*angvz(iBody)*180.0_CGREAL/PI
            cosalpha(iBody) = COS(alpha(iBody)*PI/180.0_CGREAL)
            sinalpha(iBody) = SIN(alpha(iBody)*PI/180.0_CGREAL)
           ENDIF

          CASE (BODY_DIM3)
           IF (boundary_motion_type(iBody)/=FEA_FLOW_STRUC_INTERACTION & ! added for FEA 
		       .AND.boundary_motion_type(iBody)/=FEA_TAHOE &
               .AND. boundary_motion_type(iBody)/=EIGEN_MOTION) THEN     ! SER_TO_PAR. QX. CH13

            DO i=1,nPtsBodyMarker(iBody)
              xBodyMarker(iBody,i) = xBodyMarker(iBody,i) + dt*uBodyMarker(iBody,i) 
              yBodyMarker(iBody,i) = yBodyMarker(iBody,i) + dt*vBodyMarker(iBody,i) 
              zBodyMarker(iBody,i) = zBodyMarker(iBody,i) + dt*wBodyMarker(iBody,i)              
            ENDDO ! i
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)   
!=========Update the hinge location for 3D      
            hinge_x(iBody)  =  hinge_x(iBody) + dt*hinge_vx(iBody)
            hinge_y(iBody)  =  hinge_y(iBody) + dt*hinge_vy(iBody)
            hinge_z(iBody)  =  hinge_z(iBody) + dt*hinge_vz(iBody)
           

  IF(boundary_motion_type(iBody) == FORCED .AND. nGroup_Combined == 0)THEN
            x_rot_cent(iBody) = x_rot_cent(iBody) + dt*vxcent(iBody)
            y_rot_cent(iBody) = y_rot_cent(iBody) + dt*vycent(iBody)
            z_rot_cent(iBody) = z_rot_cent(iBody) + dt*vzcent(iBody)    
  ELSE   
            xcent(iBody) = xcent(iBody) + dt*vxcent(iBody)
            ycent(iBody) = ycent(iBody) + dt*vycent(iBody)
            zcent(iBody) = zcent(iBody) + dt*vzcent(iBody) 
  ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
!           updates the angle of attack for rotating bodies
!           -----------------------------------------------
            alpha(iBody)    = alpha(iBody) + dt*angvz(iBody)*180.0_CGREAL/PI
            cosalpha(iBody) = COS(alpha(iBody)*PI/180.0_CGREAL)
            sinalpha(iBody) = SIN(alpha(iBody)*PI/180.0_CGREAL)
           ENDIF

        END SELECT ! canonical_body_type
   
      ENDIF ! boundary_motion_type

    ENDDO ! iBody

    IF (monitorON) THEN
      ifort = 264
      DO iBody = 1,nbody
        ifortCent = ifort+iBody
        write(ifortCent,'(i6,10(2x,e14.7))')  ntime,time,xcent(iBody),ycent(iBody),zcent(iBody),&
                                              vxcent(iBody),vycent(iBody),vzcent(iBody),&
                                              alpha(iBody),angvz(iBody)
      END DO ! iBody
    END IF ! ntime

    CALL calculate_arclength_norm_ds()
    CALL set_boundary(myCommTime)   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
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

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
END SUBROUTINE move_boundary
!---------------------------------------------------------------------------



SUBROUTINE compute_marker_vel(iBody)

! -------------------------------------------------------------------
!  this is a second-order accurate algorithm for motion developed by
!  R. Mittal that preserves length during rotation.
! -------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: iBody

    INTEGER            :: i,m,j
    REAL(KIND=CGREAL)   :: uB,vB,wB
    REAL(KIND=CGREAL)   :: temp_angvx, temp_angvy, temp_angvz
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER            :: BodyNumber
    INTEGER            :: iGroup
    REAL(KIND=CGREAL)   :: theta
    REAL(KIND=CGREAL)   :: xx,yy,zz    
    REAL(KIND=CGREAL)   :: R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33 
    REAL(KIND=CGREAL)   :: C_11_old,C_12_old,C_13_old,C_21_old,C_22_old,C_23_old,C_31_old,C_32_old,C_33_old  
    REAL(KIND=CGREAL)   :: uBody_Marker,vBody_Marker,wBody_Marker 
    REAL(KIND=CGREAL)   :: hingvx, hingvy, hingvz       
!dddddddddddddddddddddddddddddddd
    CHARACTER*50 :: fname1
!dddddddddddddddddddddddddddddddd    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

! =================================================================================================                             
!        n+1        n+1     n+1/2
!      v       =  v     +  v
!      - i        - c      - P 
!
! where
!
!    n+1/2     -1[       n+1/2     (  n    n ) ]
!   v       = A  [  omega      X   ( x  - r  ) ]
!   -P           [  ----           ( -    -c ) ]
!
! where
!                                       2  2                           2                   2
!                                 [ 1+ox dt /4         -ozdt/2 - oxoydt /4 -oydt/2 + oxozdt /4 ]
!      -1          1              |               2         2  2                           2   |  
!     A    = ---------------------|ozdt/2 - oxoydt /4   1+oy dt /4         -oxdt/2 - oyozdt /4 |
!                2   2    2    2  |               2                    2        2  2           |
!            1+dt (ox + oy + oz ) |oydt/2 + oxozdt /4   oxdt/2 - oyozdt /4  1+oz dt /4         |
!
! and
!                                        n    n+1         n    n+1      
!                   [   0            -(oz + oz   )/2   (oy + oy   )/2   ]
!          n+1/2    |   n    n+1                            n    n+1    |
!     omega       = |(oz + oz   )/2         0           -(ox + ox   )/2 |
!                   |     n    n+1       n    n+1                       |
!                   [ -(oy + oy   )/2 (ox + ox   )/2          0         ]
! =================================================================================================                             

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    iGroup = combined_Group_index(iBody)
    
    IF(body_dim(bodyNumber) == 2)THEN
      temp_angvx  = 0.0_CGREAL
      temp_angvy  = 0.0_CGREAL
      temp_angvz  = temp_angvz
    ENDIF      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

      temp_angvx  = 0.5_CGREAL*(angvx_old(iBody)+angvx(iBody))
      temp_angvy  = 0.5_CGREAL*(angvy_old(iBody)+angvy(iBody))
      temp_angvz  = 0.5_CGREAL*(angvz_old(iBody)+angvz(iBody))     
      
      
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
      
      C_11_old = C_11(iBody)
      C_12_old = C_12(iBody) 
      C_13_old = C_13(iBody)
      C_21_old = C_21(iBody)
      C_22_old = C_22(iBody) 
      C_23_old = C_23(iBody)   
      C_31_old = C_31(iBody) 
      C_32_old = C_32(iBody)
      C_33_old = C_33(iBody)       
      
      C_11(iBody) = R_11*C_11_old + R_12*C_21_old +R_13*C_31_old
      C_12(iBody) = R_11*C_12_old + R_12*C_22_old +R_13*C_32_old
      C_13(iBody) = R_11*C_13_old + R_12*C_23_old +R_13*C_33_old
      C_21(iBody) = R_21*C_11_old + R_22*C_21_old +R_23*C_31_old
      C_22(iBody) = R_21*C_12_old + R_22*C_22_old +R_23*C_32_old
      C_23(iBody) = R_21*C_13_old + R_22*C_23_old +R_23*C_33_old 
      C_31(iBody) = R_31*C_11_old + R_32*C_21_old +R_33*C_31_old
      C_32(iBody) = R_31*C_12_old + R_32*C_22_old +R_33*C_32_old  
      C_33(iBody) = R_31*C_13_old + R_32*C_23_old +R_33*C_33_old   
      
! -----Calculating Hinges Velocity  
DO BodyNumber = 1, nBody
   IF(boundary_motion_type(BodyNumber) == PENDULUM .AND. BodyNumber /= iBody)THEN
   
       IF(pendulum_mother_body(BodyNumber) == iBody)THEN 
            
       CALL Calc_Single_Marker_Vel(iBody, theta,temp_angvx,temp_angvy,temp_angvz,                        &
                                    hinge_x(BodyNumber),  hinge_y(BodyNumber),  hinge_z(BodyNumber),      &
                                    hingvx, hingvy, hingvz,                                               &                                  
                                    R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33)
                                                
       hinge_vx(BodyNumber) = hingvx 
       hinge_vy(BodyNumber) = hingvy 
       hinge_vz(BodyNumber) = hingvz      
       
       ENDIF
       
    ELSEIF(boundary_motion_type(BodyNumber) == PENDULUM .AND. BodyNumber == iBody)THEN  
       
       hinge_vx(BodyNumber) = 0.0_CGREAL
       hinge_vy(BodyNumber) = 0.0_CGREAL
       hinge_vz(BodyNumber) = 0.0_CGREAL
               
    ENDIF
ENDDO


   
! -----Calculating Markerpoints Velocity                   
DO i=1,nPtsBodyMarker(iBody)


   CALL Calc_Single_Marker_Vel(iBody, theta,temp_angvx,temp_angvy,temp_angvz,                             &
                                xBodyMarker(iBody,i),yBodyMarker(iBody,i),zBodyMarker(iBody,i),     &
                                uBody_Marker,vBody_Marker,wBody_Marker,                             &
                                R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33)
                                
                            
                                
    IF(boundary_motion_type(iBody) == PRESCRIBED)THEN
        uBodyMarker(iBody,i)  =  uBody_Marker + uBodyMarker_rel(ibody,i)
        vBodyMarker(iBody,i)  =  vBody_Marker + vBodyMarker_rel(ibody,i)
        wBodyMarker(iBody,i)  =  wBody_Marker + wBodyMarker_rel(ibody,i)
    ELSE
        uBodyMarker(iBody,i)  =  uBody_Marker
        vBodyMarker(iBody,i)  =  vBody_Marker
        wBodyMarker(iBody,i)  =  wBody_Marker                
    ENDIF                               

ENDDO ! i              
! ===============================For forced rotating motion===============================     
    IF(boundary_motion_type(iBody) == FORCED .AND. combined_type(iBody) == COMBINED)THEN
    
       DO i=1,nPtsBodyMarker(iBody)

	     IF(i /= i_fixed(iBody))THEN

         uBodyMarker(iBody,i) =  uBodyMarker(iBody,i) - uBodyMarker(iBody,i_fixed(iBody))
         
         vBodyMarker(iBody,i) =  vBodyMarker(iBody,i) - vBodyMarker(iBody,i_fixed(iBody))
         
         wBodyMarker(iBody,i) =  wBodyMarker(iBody,i) - wBodyMarker(iBody,i_fixed(iBody))

	     ENDIF	

       ENDDO ! i  
       
! -----Calculate the velocity of joint point fixed on the mother body   
       CALL calc_fixed_point_vel(uBodyMarker(iBody,i_fixed(iBody)), vBodyMarker(iBody,i_fixed(iBody)),   &
                                  wBodyMarker(iBody,i_fixed(iBody)), iBody, i_fixed(iBody), fixed_mother_body(iBody))
                                   
         
! -----Calculate total velocity by summing relative velocity with carrier velocity (Take body 1 to be mother body)
       DO i=1,nPtsBodyMarker(iBody)
       
          IF(i /= i_fixed(iBody))THEN
          
          uBodyMarker(iBody,i) =  uBodyMarker(iBody,i) + uBodyMarker(iBody,i_fixed(iBody))
         
          vBodyMarker(iBody,i) =  vBodyMarker(iBody,i) + vBodyMarker(iBody,i_fixed(iBody))
         
          wBodyMarker(iBody,i) =  wBodyMarker(iBody,i) + wBodyMarker(iBody,i_fixed(iBody))
          
          ENDIF
       
       ENDDO !i
    		 	                    
    ENDIF


      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)  
END SUBROUTINE compute_marker_vel
!---------------------------------------------------------------------------



SUBROUTINE read_marker_vel(ibody)

!   ------------------------------------------------------------------
!    This subroutine sets new velocity of marker points and centroid.
!   ------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE derivative    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

    IMPLICIT NONE

    INTEGER, PARAMETER  :: ifort=64
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER              :: iGroup
    INTEGER,INTENT(IN)  :: ibody
!ddddddddddddddddddddddddddddddddddddddddddddddddddd
    REAL(KIND=CGREAL)     :: amplify
!ddddddddddddddddddddddddddddddddddddddddddddddddddd    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    INTEGER             :: m,ifortCent,nptsbm_temp
    INTEGER             :: i, ntime_check, nptsbm_skip         ! added by Haibo

    REAL(KIND=CGREAL)   :: temp_time,temp_dt

    REAL(KIND=CGREAL)   :: time_skip, dt_skip
    REAL(KIND=CGREAL)   :: uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip  
    
    CHARACTER*20 :: fname1    
  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    iGroup = combined_Group_index(iBody)
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    
    ifortCent = ifort+iBody
    IF (ifortCent > 80) PRINT*,'Might have problems with file numbers'
 
!   --------- new added by Haibo ---------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ntime_check = mod(ntime - 1, ntimePerCycle(ibody))
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    ntime_skip  = mod(ntime_start, ntimePerCycle(ibody))

    IF (ntime > 1 .AND. ntime_check == 0) THEN
      CLOSE (ifortCent)
!      WRITE(fname1,"('body_mesh_velocity.',I3.3,'.dat')") iBody
!      OPEN(UNIT=ifortCent,FILE=fname1)
    END IF

    IF (nread == 1 .AND. ntime == ntime_start+1)THEN
      IF ( ntime_skip /= 0) THEN
        CLOSE (ifortCent)
!        WRITE(fname1,"('body_mesh_velocity.',I3.3,'.dat')") iBody
!        OPEN(UNIT=ifortCent,FILE=fname1)

        DO i = 1, ntime_skip
          READ(ifortCent,*) dt_skip,time_skip,nptsbm_skip
          DO m=1,nPtsBodyMarker(ibody)
            READ(ifortCent,*) uBodyMarker_skip,vBodyMarker_skip,wBodyMarker_skip
          END DO ! m
        END DO
        ntime_skip = 0
      END IF
    ENDIF
!   ----------------------------------------------------
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
   IF(k_derivative == 1 .or. nread == 1)THEN
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)   
    READ(ifortCent,*)temp_dt,temp_time,nptsbm_temp
    !PRINT*,temp_dt,temp_time,nptsbm_temp 
    !PRINT*,time

    IF (ABS(dt-temp_dt) > 1.0E-8) THEN
      IF (ImtheBOSS) THEN
        PRINT*,'DT in code and marker_vel files do not match'
        PRINT*,dt,temp_dt
      END IF
!      PRINT*,'Aborting Run'
!      STOP
    ENDIF

    IF (ABS(time-temp_time) > 1.0E-8) THEN
      IF (ImtheBOSS) THEN
        PRINT*,'Time stamp in code and marker_vel files do not match'
        PRINT*,time,temp_time
      END IF
!     PRINT*,'Aborting Run'
!     STOP
    ENDIF

    IF (nptsbm_temp /= nPtsBodyMarker(ibody)) THEN
      IF (ImtheBOSS) THEN
        PRINT*,'nPtsBodyMarker in marker_vel file does not match code'
        PRINT*,nPtsBodyMarker(iBody),nptsbm_temp
        PRINT*,'Aborting Run'
      END IF
      CALL flow_stop
      STOP
    ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
   ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    IF(combined_type(iBody) == UNCOMBINED .OR. ntime <= nstart_FIM)THEN
         
      DO m=1,nPtsBodyMarker(ibody)
        READ(ifortCent,*) uBodyMarker(ibody,m),vBodyMarker(ibody,m),wBodyMarker(ibody,m)
      END DO ! m
      
!ddddddddddddddddddddddddddddddddddddddddddddddddddd
    OPEN(6688,FILE='amplify_vel.dat',STATUS='UNKNOWN')
    READ(6688,*) amplify
    CLOSE(6688)
    uBodyMarker(ibody,:)=uBodyMarker(ibody,:)*amplify
    vBodyMarker(ibody,:)=vBodyMarker(ibody,:)*amplify
    wBodyMarker(ibody,:)=wBodyMarker(ibody,:)*amplify
!ddddddddddddddddddddddddddddddddddddddddddddddddddd      

!---------------UPDATE VELOCITY FOR PRESCRIBED MOTION 
    ELSE !===>combined_type(iBody) == COMBINED .AND. ntime > nstart_FIM
!---------------The velocity becomes relative velocity in "fort.66&fort.67" 
      IF(k_derivative == 1 .or. nread == 1)THEN
        
      DO m=1,nPtsBodyMarker(ibody)
        READ(ifortCent,*) uBodyMarker_rel_old(ibody,m),vBodyMarker_rel_old(ibody,m),wBodyMarker_rel_old(ibody,m)
      END DO ! m  
      
!ddddddddddddddddddddddddddddddddddddddddddddddddddd
    OPEN(6688,FILE='amplify_vel.dat',STATUS='UNKNOWN')
    READ(6688,*) amplify
    CLOSE(6688)
    uBodyMarker_rel_old(ibody,:)=uBodyMarker_rel_old(ibody,:)*amplify
    vBodyMarker_rel_old(ibody,:)=vBodyMarker_rel_old(ibody,:)*amplify
    wBodyMarker_rel_old(ibody,:)=wBodyMarker_rel_old(ibody,:)*amplify
!ddddddddddddddddddddddddddddddddddddddddddddddddddd       
      
      ENDIF !k_derivative == 1
     
!---------------Rotate the velocity vector in original files 
!     [ u ]      [ C11         C12         C13][ u_old ]
!     |   |      |                            ||       |  
!     | v |   =  | C21         C22         C23|| v_old |
!     |   |      |                            ||       |
!     | w |      | C31         C32         C33|| w_old |      
      DO m=1,nPtsBodyMarker(ibody) 
      
      uBodyMarker_rel(ibody,m) = uBodyMarker_rel_old(ibody,m)*C11(iGroup) + vBodyMarker_rel_old(ibody,m)*C12(iGroup) + wBodyMarker_rel_old(ibody,m)*C13(iGroup)
      vBodyMarker_rel(ibody,m) = uBodyMarker_rel_old(ibody,m)*C21(iGroup) + vBodyMarker_rel_old(ibody,m)*C22(iGroup) + wBodyMarker_rel_old(ibody,m)*C23(iGroup)
      wBodyMarker_rel(ibody,m) = uBodyMarker_rel_old(ibody,m)*C31(iGroup) + vBodyMarker_rel_old(ibody,m)*C32(iGroup) + wBodyMarker_rel_old(ibody,m)*C33(iGroup)

      END DO ! m  
                                    
                     
    ENDIF
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

  END SUBROUTINE read_marker_vel
!---------------------------------------------------------------------


!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
SUBROUTINE Vector_Rotation_Matrix(iGroup)
!   ------------------------------------------------------------------
!    This subroutine calculate new displacement in location of mass center and angle
!   ------------------------------------------------------------------
 USE global_parameters
 USE boundary_arrays
 USE flow_parameters 
 USE derivative
 
 INTEGER, INTENT(IN) :: iGroup
 REAL(KIND=CGREAL)   :: temp_angvx,temp_angvy,temp_angvz
 REAL(KIND=CGREAL)   :: x,y,z
 REAL(KIND=CGREAL)   :: C11_old,C12_old,C13_old,C21_old,C22_old,C23_old,C31_old,C32_old,C33_old 
 REAL(KIND=CGREAL)   :: det
!ddddddddddddddddddddddddddddddddddddddddd
    CHARACTER*50 :: fname1
!ddddddddddddddddddddddddddddddddddddddddd 
! =================================================================================================================================
!               [ cos(theta) + (1-cos(theta))*x^2      (1-cos(theta))*x*y - sin(theta)*z      (1-cos(theta))*x*z + sin(theta)*y]       
!               |                                                                                                              |                                   
!     R(n)   =  | (1-cos(theta))*y*x + sin(theta)*z      cos(theta) + (1-cos(theta))*y^2      (1-cos(theta))*y*z - sin(theta)*x|  
!               |                                                                                                              |
!               | (1-cos(theta))*z*x - sin(theta)*y      (1-cos(theta))*z*y + sin(theta)*x      cos(theta) + (1-cos(theta))*z^2|

    
!             [ C11         C12         C13]       
!             |                            |                                   
!     C    =  | C21         C22         C23|   =  R(1) * R(2) * ... ... * R(n) 
!             |                            |
!             | C31         C32         C33|
! =================================================================================================================================
    temp_angvx  = 0.5_CGREAL*(angvx_combined_old(iGroup)+angvx_combined(iGroup))
    temp_angvy  = 0.5_CGREAL*(angvy_combined_old(iGroup)+angvy_combined(iGroup))
    temp_angvz  = 0.5_CGREAL*(angvz_combined_old(iGroup)+angvz_combined(iGroup))
        

    IF(sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2) > 0.0_CGREAL)THEN
    
    x     = temp_angvx/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)
    y     = temp_angvy/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)
    z     = temp_angvz/sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)
    theta_group(iGroup) = dt * sqrt(temp_angvx**2 + temp_angvy**2 + temp_angvz**2)

    R11(iGroup) = cos(theta_group(iGroup)) + (1.0_CGREAL-cos(theta_group(iGroup)))*x**2.0_CGREAL
    R12(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*x*y - sin(theta_group(iGroup))*z
    R13(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*x*z + sin(theta_group(iGroup))*y
    R21(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*y*x + sin(theta_group(iGroup))*z
    R22(iGroup) = cos(theta_group(iGroup)) + (1.0_CGREAL-cos(theta_group(iGroup)))*y**2.0_CGREAL 
    R23(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*y*z - sin(theta_group(iGroup))*x 
    R31(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*z*x - sin(theta_group(iGroup))*y 
    R32(iGroup) = (1.0_CGREAL-cos(theta_group(iGroup)))*z*y + sin(theta_group(iGroup))*x    
    R33(iGroup) = cos(theta_group(iGroup)) + (1.0_CGREAL-cos(theta_group(iGroup)))*z**2.0_CGREAL
    
    C11_old = C11(iGroup)
    C12_old = C12(iGroup) 
    C13_old = C13(iGroup)
    C21_old = C21(iGroup)
    C22_old = C22(iGroup) 
    C23_old = C23(iGroup)   
    C31_old = C31(iGroup) 
    C32_old = C32(iGroup)
    C33_old = C33(iGroup) 
    
!------Calculate velocities in non-inertia coordinates (Body-fixed Coordinates)
    non_vxcent(iGroup) = C11_old*vxcent_combined(iGroup) + C12_old*vycent_combined(iGroup) + C13_old*vzcent_combined(iGroup)
    non_vycent(iGroup) = C21_old*vxcent_combined(iGroup) + C22_old*vycent_combined(iGroup) + C23_old*vzcent_combined(iGroup)
    non_vzcent(iGroup) = C31_old*vxcent_combined(iGroup) + C32_old*vycent_combined(iGroup) + C33_old*vzcent_combined(iGroup)
    angv_roll(iGroup)  = C11_old*temp_angvx + C12_old*temp_angvy + C13_old*temp_angvz
    angv_yaw(iGroup)   = C21_old*temp_angvx + C22_old*temp_angvy + C23_old*temp_angvz   
    angv_pitch(iGroup) = C31_old*temp_angvx + C32_old*temp_angvy + C33_old*temp_angvz
    
     
    C11(iGroup) = R11(iGroup)*C11_old + R12(iGroup)*C21_old +R13(iGroup)*C31_old
    C12(iGroup) = R11(iGroup)*C12_old + R12(iGroup)*C22_old +R13(iGroup)*C32_old
    C13(iGroup) = R11(iGroup)*C13_old + R12(iGroup)*C23_old +R13(iGroup)*C33_old
    C21(iGroup) = R21(iGroup)*C11_old + R22(iGroup)*C21_old +R23(iGroup)*C31_old
    C22(iGroup) = R21(iGroup)*C12_old + R22(iGroup)*C22_old +R23(iGroup)*C32_old
    C23(iGroup) = R21(iGroup)*C13_old + R22(iGroup)*C23_old +R23(iGroup)*C33_old 
    C31(iGroup) = R31(iGroup)*C11_old + R32(iGroup)*C21_old +R33(iGroup)*C31_old
    C32(iGroup) = R31(iGroup)*C12_old + R32(iGroup)*C22_old +R33(iGroup)*C32_old  
    C33(iGroup) = R31(iGroup)*C13_old + R32(iGroup)*C23_old +R33(iGroup)*C33_old  
    
    ELSE
    
    C11(iGroup) = C11(iGroup)
    C12(iGroup) = C12(iGroup)
    C13(iGroup) = C13(iGroup)
    C21(iGroup) = C21(iGroup)
    C22(iGroup) = C22(iGroup)
    C23(iGroup) = C23(iGroup)
    C31(iGroup) = C31(iGroup)
    C32(iGroup) = C32(iGroup)
    C33(iGroup) = C33(iGroup)
    
    
    ENDIF 
    
    det = C11(iGroup)*(C22(iGroup)*C33(iGroup) - C23(iGroup)*C32(iGroup))  +           &
          C12(iGroup)*(C23(iGroup)*C31(iGroup) - C21(iGroup)*C33(iGroup))  +           &
          C13(iGroup)*(C21(iGroup)*C32(iGroup) - C22(iGroup)*C31(iGroup))  
          
    invC11(iGroup) = ( C22(iGroup)*C33(iGroup) - C23(iGroup)*C32(iGroup) ) / det
    invC12(iGroup) = ( C13(iGroup)*C32(iGroup) - C12(iGroup)*C33(iGroup) ) / det
    invC13(iGroup) = ( C12(iGroup)*C23(iGroup) - C13(iGroup)*C22(iGroup) ) / det
    invC21(iGroup) = ( C23(iGroup)*C31(iGroup) - C21(iGroup)*C33(iGroup) ) / det
    invC22(iGroup) = ( C11(iGroup)*C33(iGroup) - C13(iGroup)*C31(iGroup) ) / det
    invC23(iGroup) = ( C13(iGroup)*C21(iGroup) - C11(iGroup)*C23(iGroup) ) / det
    invC31(iGroup) = ( C21(iGroup)*C32(iGroup) - C22(iGroup)*C31(iGroup) ) / det
    invC32(iGroup) = ( C12(iGroup)*C31(iGroup) - C11(iGroup)*C32(iGroup) ) / det
    invC33(iGroup) = ( C11(iGroup)*C22(iGroup) - C12(iGroup)*C21(iGroup) ) / det
                           

!dddddddddddddddddddddddddddddddddddddddddd  
!IF (MOD(ntime,ndump) == 0 .AND. kimplicit_FIM == 1) THEN
IF(IMtheBOSS)THEN
 IF(derivative_flag == 0 .OR. k_derivative == 4)THEN
   WRITE(fname1,"('rotation_matrix.',i7.7,'.dat')") ntime  
   OPEN(UNIT=55555,FILE=fname1,STATUS='UNKNOWN') 
   WRITE(55555,*) '---------------------------------'    
   WRITE(55555,*) C11(iGroup),C12(iGroup),C13(iGroup)
   WRITE(55555,*) C21(iGroup),C22(iGroup),C23(iGroup)
   WRITE(55555,*) C31(iGroup),C32(iGroup),C33(iGroup)
   WRITE(55555,*) '=================================' 
   WRITE(55555,*) invC11(iGroup),invC12(iGroup),invC13(iGroup)
   WRITE(55555,*) invC21(iGroup),invC22(iGroup),invC23(iGroup)
   WRITE(55555,*) invC31(iGroup),invC32(iGroup),invC33(iGroup)
   WRITE(55555,*) '---------------------------------' 
   WRITE(55555,*) 'temp_angvx,temp_angvy,temp_angvz'    
   WRITE(55555,*) temp_angvx,temp_angvy,temp_angvz 
   WRITE(55555,*) 'x,y,z,theta'        
   WRITE(55555,*) x,y,z,theta_group(iGroup)  
   WRITE(55555,*) 'angv_roll,angv_yaw,angv_pitch'
   WRITE(55555,*) angv_roll(iGroup),angv_yaw(iGroup),angv_pitch(iGroup)
   WRITE(55555,*) 'non_vxcent,non_vycent,non_vzcent'     
   WRITE(55555,*) non_vxcent(iGroup), non_vycent(iGroup), non_vzcent(iGroup)
 ENDIF
ENDIF
!ENDIF     
!dddddddddddddddddddddddddddddddddddddddddd 
END SUBROUTINE Vector_Rotation_Matrix



SUBROUTINE calc_fixed_point_vel(u, v, w, fixed_body, fixed_point, mother_body)
!   ------------------------------------------------------------------
!    Only for ***Forced Motion***
!    Calculate the velocity of joint point fixed on the mother body (Same as the newly added 2nd order method to preserves length)  
!   ------------------------------------------------------------------
USE global_parameters 
USE flow_parameters
USE grid_arrays
USE boundary_arrays
USE unstructured_surface_arrays

IMPLICIT NONE
INTEGER         , INTENT(IN)    :: fixed_body, fixed_point, mother_body
REAL(KIND=CGREAL), INTENT(OUT)   :: u, v, w
REAL(KIND=CGREAL)   :: theta
REAL(KIND=CGREAL)   :: temp_uB, temp_vB, temp_wB, temp_angvx, temp_angvy, temp_angvz
REAL(KIND=CGREAL)   :: temp_uB_pre, temp_vB_pre, temp_wB_pre
REAL(KIND=CGREAL)   :: temp_uB_new, temp_vB_new, temp_wB_new
REAL(KIND=CGREAL)   :: average_uB,average_vB,average_wB
REAL(KIND=CGREAL)   :: temp_vel_magni_pre,temp_vel_magni_new
REAL(KIND=CGREAL)   :: xx,yy,zz    
REAL(KIND=CGREAL)   :: R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33 

temp_angvx  = 0.5_CGREAL*(angvx_old(mother_body)+angvx(mother_body))
temp_angvy  = 0.5_CGREAL*(angvy_old(mother_body)+angvy(mother_body))
temp_angvz  = 0.5_CGREAL*(angvz_old(mother_body)+angvz(mother_body))

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

temp_uB_pre =  ( temp_angvy*(zBodyMarker(fixed_body,fixed_point)-zcent(mother_body)) &
              - temp_angvz*(yBodyMarker(fixed_body,fixed_point)-ycent(mother_body)) )
 
temp_vB_pre = -( temp_angvx*(zBodyMarker(fixed_body,fixed_point)-zcent(mother_body)) &
              - temp_angvz*(xBodyMarker(fixed_body,fixed_point)-xcent(mother_body)) )
 
temp_wB_pre =  ( temp_angvx*(yBodyMarker(fixed_body,fixed_point)-ycent(mother_body)) &
              - temp_angvy*(xBodyMarker(fixed_body,fixed_point)-xcent(mother_body)) )
            

temp_vel_magni_pre = sqrt(temp_uB_pre**2 + temp_vB_pre**2 + temp_wB_pre**2)


! -----1. The First Step 

        temp_uB_new = temp_uB_pre*R_11 + temp_vB_pre*R_12 + temp_wB_pre*R_13

        temp_vB_new = temp_uB_pre*R_21 + temp_vB_pre*R_22 + temp_wB_pre*R_23
		
        temp_wB_new = temp_uB_pre*R_31 + temp_vB_pre*R_32 + temp_wB_pre*R_33
 

! -----2. The Second Step 

        average_uB = temp_uB_pre + temp_uB_new

        average_vB = temp_vB_pre + temp_vB_new
		
        average_wB = temp_wB_pre + temp_wB_new

		temp_vel_magni_new = sqrt(average_uB**2 + average_vB**2 + average_wB**2)
	

        temp_uB = average_uB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new

        temp_vB = average_vB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new
		
        temp_wB = average_wB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new


! -----3. The Third Step

        u =  vxcent(mother_body) + temp_uB       

        v =  vycent(mother_body) + temp_vB
 
        w =  vzcent(mother_body) + temp_wB        
     
END SUBROUTINE calc_fixed_point_vel
 



SUBROUTINE Calc_Single_Marker_Vel(iBody, theta,temp_angvx,temp_angvy,temp_angvz,    &
                                    xxx,yyy,zzz,u,v,w,                                &
                                    R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33)


    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)           :: iBody
    REAL(KIND=CGREAL),INTENT(IN)  :: theta
    REAL(KIND=CGREAL),INTENT(IN)  :: temp_angvx, temp_angvy, temp_angvz
    REAL(KIND=CGREAL),INTENT(IN)  :: xxx,yyy,zzz    
    REAL(KIND=CGREAL),INTENT(OUT) :: u,v,w      
    REAL(KIND=CGREAL),INTENT(IN)  :: R_11,R_12,R_13,R_21,R_22,R_23,R_31,R_32,R_33

  
    REAL(KIND=CGREAL)   :: uB,vB,wB
    REAL(KIND=CGREAL)   :: temp_uB, temp_vB, temp_wB
    REAL(KIND=CGREAL)   :: temp_uB_pre, temp_vB_pre, temp_wB_pre
    REAL(KIND=CGREAL)   :: temp_uB_new, temp_vB_new, temp_wB_new
    REAL(KIND=CGREAL)   :: average_uB,average_vB,average_wB
    REAL(KIND=CGREAL)   :: temp_vel_magni_pre,temp_vel_magni_new  
    
!  Newly added 2nd order method to preserves length during rotation by Zhang Chao
! -------------------------------------------------------------------
! ================================================================================================= 
! 1. The First Step
!        n+1                   n
!      v       =   R(n)  *   v   
!       - P                   - P
! 2. The Second Step                       
!                     n          n+1
!                   v       +  v           
!        n+1/2       - P        - P         |    n  |     2*sin(theta/2)
!      v       = ----------------------  *  |  v    |  * ---------------
!       - P      |    n          n+1  |     |   - P |         theta 
!                |  v       +  v      |
!                |   - P        - P   |

! 3. The Third Step
!        n+1        n+1     n+1/2
!      v       =  v     +  v
!      - i        - c      - P 
! ================================================================================================= 

	
  IF(boundary_motion_type(iBody) == FORCED .AND. nGroup_Combined == 0)THEN
  
        temp_uB_pre =  ( temp_angvy*(zzz - z_rot_cent(iBody)) &
                      - temp_angvz*(yyy - y_rot_cent(iBody)) )
 
        temp_vB_pre = -( temp_angvx*(zzz - z_rot_cent(iBody)) &
                      - temp_angvz*(xxx - x_rot_cent(iBody)) )
 
        temp_wB_pre =  ( temp_angvx*(yyy - y_rot_cent(iBody)) &
                      - temp_angvy*(xxx - x_rot_cent(iBody)) )  
                      
  ELSE
  
        temp_uB_pre =  ( temp_angvy*(zzz - zcent(iBody)) &
                      - temp_angvz*(yyy - ycent(iBody)) )
 
        temp_vB_pre = -( temp_angvx*(zzz - zcent(iBody)) &
                      - temp_angvz*(xxx - xcent(iBody)) )
 
        temp_wB_pre =  ( temp_angvx*(yyy - ycent(iBody)) &
                      - temp_angvy*(xxx - xcent(iBody)) )
  ENDIF                    

		temp_vel_magni_pre = sqrt(temp_uB_pre**2 + temp_vB_pre**2 + temp_wB_pre**2)
		
		
! -----1. The First Step 

        temp_uB_new = temp_uB_pre*R_11 + temp_vB_pre*R_12 + temp_wB_pre*R_13

        temp_vB_new = temp_uB_pre*R_21 + temp_vB_pre*R_22 + temp_wB_pre*R_23
		
        temp_wB_new = temp_uB_pre*R_31 + temp_vB_pre*R_32 + temp_wB_pre*R_33


! -----2. The Second Step 

        average_uB = temp_uB_pre + temp_uB_new

        average_vB = temp_vB_pre + temp_vB_new
		
        average_wB = temp_wB_pre + temp_wB_new

		temp_vel_magni_new = sqrt(average_uB**2 + average_vB**2 + average_wB**2)
                 
  		
		
		IF(temp_vel_magni_new > 0.0_CGREAL .AND. temp_vel_magni_pre > 0.0_CGREAL)THEN
	

        temp_uB = average_uB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new

        temp_vB = average_vB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new
		
        temp_wB = average_wB * temp_vel_magni_pre * sin(0.5_CGREAL*theta)/(0.5_CGREAL*theta) / temp_vel_magni_new
        
	        
        
        ELSE
        
        temp_uB = 0.0_CGREAL

        temp_vB = 0.0_CGREAL
		
        temp_wB = 0.0_CGREAL       
        
        ENDIF
     

! -----3. The Third Step

        uB =  vxcent(iBody) + temp_uB       

        vB =  vycent(iBody) + temp_vB
 
        wB =  vzcent(iBody) + temp_wB 
        
               
        u  =  uB
        v  =  vB
        w  =  wB        
  

END SUBROUTINE Calc_Single_Marker_Vel    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    