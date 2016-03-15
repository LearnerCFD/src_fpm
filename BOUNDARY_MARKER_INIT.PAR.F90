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
!  Filename: BOUNDARY_MARKER_INIT.PAR.F90
!  Latest Modification: December 30, 2007 (ver. 0.3.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     initialize_marker()
!     extend_cylinder_3D(iBody)
!     extend_cylinder_vel_3D(iBody)
!     get_edge_node(iBody,nv,tot_num_edge,num_of_edge_node)
! --------------------------------------------------------------------



SUBROUTINE initialize_marker()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!   Loop variables
!   --------------
    INTEGER :: i, j, k, iBody, nnt, iErr

!   Local variables
!   ---------------
    INTEGER :: m, nBodyMarkerIn, nPtsBodyMarkerIn, totNumTriElemIn, num_edge
    INTEGER, ALLOCATABLE, DIMENSION(:) :: edge_node_number

    INTEGER, DIMENSION(nBody) :: body_type_orig

    REAL(KIND=CGREAL) :: dMin, dxMin, dyMin,                        &
                         bodyResolution, bodyResolutionNormalized,  &
                         theta, phi, xTemp, yTemp

    LOGICAL           :: readMarkerFlag, openMemb 

    num_edge = 0

    IF (nRead /= 1) THEN
    
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A)') 'Setting up CANONICAL BODIES'
        WRITE(STDOUT,'(5X,A)') '==========================='
        WRITE(STDOUT,*)
      END IF

!     Save a copy of canonical body type
!     ----------------------------------
      body_type_orig(1:nBody) = canonical_body_type(1:nBody)

!     Initialize for  no restart
!     --------------------------
      uBodyMarker = 0.0_CGREAL
      vBodyMarker = 0.0_CGREAL
      wBodyMarker = 0.0_CGREAL

      xBodyMarker = 0.0_CGREAL
      yBodyMarker = 0.0_CGREAL
      zBodyMarker = 0.0_CGREAL
    
!     Select appropriate body type and determine marker locations
!     -----------------------------------------------------------
      IF (ImtheBOSS) OPEN(UNIT=ifuMarkerTRI, file='marker_unstruc_tri.dat', STATUS='UNKNOWN')

      DO iBody = 1, nBody

        SELECT CASE (canonical_body_type(iBody))

        CASE(ELLIPTIC_CYLINDER)
          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I3.3,A)') 'Body number',iBody,': ELLIPTICAL CYLINDER'
            WRITE(STDOUT,'(5X,A)') '------------------------------------'
          END IF

!         Test resolution for GCM
!         -----------------------
          IF ( boundary_formulation == GCM_METHOD ) THEN
            dxMin          = MINVAL(dx(1:nxc_GLBL))
            dyMin          = MINVAL(dy(1:nyc_GLBL))
            dMin           = MIN(dxMin,dyMin)
            bodyResolution = PI*( radiusx(iBody)+radiusy(iBody) ) /   &
                             REAL(nPtsBodyMarker(iBody),KIND=CGREAL)
            bodyResolutionNormalized = dMin/bodyResolution

            IF (monitorON) THEN
              WRITE(STDOUT,'(5X,A,1X,1PE15.7)') 'dx Min =',dxMin
              WRITE(STDOUT,'(5X,A,1X,1PE15.7)') 'dy Min =',dyMin
              WRITE(STDOUT,'(5X,A,1X,1PE15.7)') 'd  Min =',dMin
              WRITE(STDOUT,'(5X,A,1X,1PE15.7)') 'Body Resolution =',bodyResolution
              WRITE(STDOUT,'(5X,A,1X,1PE15.7)') 'Current Normalized Resolution for Body =',bodyResolutionNormalized
            END IF

            IF ( bodyResolutionNormalized < 2.0_CGREAL ) THEN
              WRITE(STDOUT,'(5X,A)') '!!! WARINING !!! Ideal Normalized Resolution Should be at LEAST 2 .. Aborting!'
              CALL flow_stop()
              STOP
            ENDIF ! bodyResolutionNormalized

          ENDIF ! boundary_formulation

          DO m = 1, nPtsBodyMarkerOrig(iBody)
            theta  = REAL(m-1,KIND=CGREAL)*2.0_CGREAL*PI/ &
                     REAL(nPtsBodyMarkerOrig(iBody),KIND=CGREAL)        
            xTemp                = radiusx(iBody)*COS(theta)
            yTemp                = radiusy(iBody)*SIN(theta)
            xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody) 
            yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody) 
            zBodyMarker(iBody,m) = z(1)
          ENDDO ! m

          CALL extend_cylinder_3D(iBody)
          
          IF (monitorON) WRITE(STDOUT,*)

        CASE(GENERAL_CYLINDER)
          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I3.3,A)') 'Body number',iBody,': GENERAL CYLINDER'
            WRITE(STDOUT,'(5X,A)') '---------------------------------'
          END IF
          
          READ(ifuMarkerIn,*) nBodyMarkerIn
          IF ( nBodyMarkerIn /= nPtsBodyMarkerOrig(iBody) ) THEN
            WRITE(STDOUT,'(5X,A,1X,I5)') 'Init_Marker: Inconsistent body_in.dat and marker_in.dat files for body =', iBody
            WRITE(STDOUT,'(5X,A,1X,I5)') 'Defined in body_in.dat, nPtsBodyMarker =', nPtsBodyMarkerOrig(iBody)
            WRITE(STDOUT,'(5X,A,1X,I5)') 'Defined in marker_in.dat, nBodyMarkerIn =', nBodyMarkerIn
            CALL flow_stop()
            STOP
          ENDIF ! nBodyMarkerIn
            
          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I5)') 'nBodyMarkerIn  =', nBodyMarkerIn
            WRITE(STDOUT,'(5X,A,1X,I5)') 'nPtsBodyMarker =', nPtsBodyMarker(iBody)
            WRITE(STDOUT,*)
          END IF  ! monitorON

          DO m = 1, nPtsBodyMarkerOrig(iBody)
            READ(ifuMarkerIn,*) xBodyMarker(iBody,m), yBodyMarker(iBody,m)
            zBodyMarker(iBody,m) = z(1)
          ENDDO ! m

          CALL extend_cylinder_3D(iBody)
           
        CASE(ELLIPSOID)
          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I3.3,A)') 'Body number',iBody,': ELLIPSOID'
            WRITE(STDOUT,'(5X,A)') '--------------------------'
          END IF

          m = 0
 
          DO i = 1, n_phi(iBody)
            phi  = REAL(i-1           ,KIND=CGREAL)*PI/ &
                   REAL(n_phi(iBody)-1,KIND=CGREAL)

            IF (i==1 .OR. i==n_phi(iBody)) THEN
              theta = 0.0_CGREAL
              m = m + 1
 
              xTemp                = radiusx(iBody)*COS(theta)*SIN(phi)
              yTemp                = radiusy(iBody)*SIN(theta)*SIN(phi)
              xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody) 
              yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody)
              zBodyMarker(iBody,m) = zcent(iBody) + radiusz(iBody)*COS(phi)
            ELSE
              DO j = 1,n_theta(iBody)
                theta = REAL(j-1,KIND=CGREAL)*2.0_CGREAL*PI/ &
                        REAL(n_theta(iBody),KIND=CGREAL)
                m = m + 1
                xTemp                = radiusx(iBody)*COS(theta)*SIN(phi)
                yTemp                = radiusy(iBody)*SIN(theta)*SIN(phi)
                xBodyMarker(iBody,m) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody)
                yBodyMarker(iBody,m) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody)
                zBodyMarker(iBody,m) = zcent(iBody) + radiusz(iBody)*COS(phi)
              ENDDO ! j
            END IF
          ENDDO ! i

          nPtsBodyMarker(iBody) = m
          IF(monitorON) THEN
           WRITE(STDOUT,'(9X,A,1X,I8)') 'iBody                         = ', iBody
           WRITE(STDOUT,'(9X,A,1X,I8)') 'Total number of marker points = ', m
          ENDIF
    
          i = 1
          j = 0

          DO k = 1, n_theta(iBody)
            j = j+1
            triElemNeig(iBody,1,j) = i
            triElemNeig(iBody,2,j) = i + k
            IF (i+k==n_theta(ibody) + 1) THEN
              triElemNeig(iBody,3,j) = i + 1
            ELSE 
              triElemNeig(iBody,3,j) = i + k + 1
            END IF
          END DO ! end k 

          nnt = 1

          DO i = 2, m-n_theta(iBody)-1
            j = j + 1
            triElemNeig(iBody,1,j) = i
            triElemNeig(iBody,2,j) = i + n_theta(iBody)

            IF ( i .EQ. nnt*n_theta(ibody)+1 ) THEN
              triElemNeig(iBody,3,j) = i + 1                 
            ELSE
              triElemNeig(iBody,3,j) = i + n_theta(iBody) + 1
            END IF

            j = j + 1 
            triElemNeig(iBody,1,j) = i
            IF ( i .EQ. nnt*n_theta(ibody) + 1 ) THEN
              triElemNeig(iBody,2,j) = i + 1
              triElemNeig(iBody,3,j) = i + 1 - n_theta(iBody)
              nnt = nnt + 1                   
            ELSE
              triElemNeig(iBody,2,j) = i + n_theta(iBody) + 1
              triElemNeig(iBody,3,j) = i + 1
            END IF
          ENDDO ! i

          DO i = m-n_theta(iBody), m-1
            j = j + 1
            triElemNeig(iBody,1,j) = i
            triElemNeig(iBody,2,j) = m

            IF (i+1 .EQ. m) THEN
              triElemNeig(iBody,3,j) = m-n_theta(iBody) 
            ELSE
              triElemNeig(iBody,3,j) = i + 1
            END IF
          END DO ! end i

          totNumTriElem(iBody) = j

          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I8)') 'Total number of marker points =', m
            WRITE(STDOUT,'(5X,A,1X,I8)') 'Total number of elements      =', totNumTriElem(iBody)
            WRITE(STDOUT,*)
          END IF

        CASE(UNSTRUCTURED_SURFACE)
          IF (monitorON) THEN
            WRITE(STDOUT,'(5X,A,1X,I3.3,A)') 'Body number',iBody,': UNSTRUCTURED SURFACE'
            WRITE(STDOUT,*)
          END IF

          IF ( boundary_motion /= MOVING_BOUNDARY .OR. nread /= 1) THEN
            READ(ifuUnstrucSurfIn,*)
            READ(ifuUnstrucSurfIn,*) nPtsBodyMarkerIn, totNumTriElemIn
            READ(ifuUnstrucSurfIn,*)

            IF ( nPtsBodyMarkerIn /= nPtsBodyMarker(iBody) ) THEN
              WRITE(STDOUT,'(9X,A,I5)') &
                'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
              WRITE(STDOUT,'(9X,A,I5)') &
                'Defined in canonical_body_in.dat,  nPtsBodyMarker   = ',nPtsBodyMarker(iBody)
              WRITE(STDOUT,'(9X,A,I5)') &
                'Defined in unstruc_surface_in.dat, nPtsBodyMarkerIn = ',nPtsBodyMarkerIn
              CALL flow_stop()
              STOP
            ENDIF ! nPtsBodyMarkerIn
             
            IF ( totNumTriElemIn /= totNumTriElem(iBody) ) THEN
              WRITE(STDOUT,'(9X,A,I5)') & 
                'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
              WRITE(STDOUT,'(9X,A,I5)') &
                'Defined in canonical_body_in.dat,   totNumTriElem   = ', totNumTriElem(iBody)
              WRITE(STDOUT,'(9X,A,I5)') &
                'Defined in unstruc_surface_in.dat,  totNumTriElemIn = ', totNumTriElemIn
              CALL flow_stop()
              STOP
            ENDIF !  totNumTriElemIn
          
            DO m=1,nPtsBodyMarker(iBody)
              READ(ifuUnstrucSurfIn,*) i, xBodyMarker(iBody,m), yBodyMarker(iBody,m), zBodyMarker(iBody,m)
            ENDDO
            READ(ifuUnstrucSurfIn,*)

            DO  j=1,totNumTriElem(iBody)
              READ(ifuUnstrucSurfIn,*) i, triElemNeig(iBody,1,j), triElemNeig(iBody,2,j), triElemNeig(iBody,3,j)
            ENDDO
            READ(ifuUnstrucSurfIn,*)

            READ(ifuUnstrucSurfIn,*) pointOutsideBodyX(iBody), pointOutsideBodyY(iBody), pointOutsideBodyZ(iBody)
          ENDIF

!         Write surface mesh data in Tecplot Format to check
!         --------------------------------------------------
          IF (ImtheBOSS) THEN
            WRITE(ifuMarkerTRI,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
            WRITE(ifuMarkerTRI,*) 'VARIABLES="X","Y","Z"'
            WRITE(ifuMarkerTRI,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'F=FEPOINT, ET=TRIANGLE'

            DO m=1,nPtsBodyMarker(iBody)
              WRITE(ifuMarkerTRI,*) xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
            ENDDO
            DO j=1,totNumTriElem(iBody)
              WRITE(ifuMarkerTRI,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
            ENDDO
          END IF

        END SELECT ! canonical_body_type

      ENDDO ! iBody

      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A)') 'End of CANONICAL BODY set up'
        WRITE(STDOUT,'(5X,A)') '============================'
        WRITE(STDOUT,*)
      END IF

!     set up initial marker velocities
!     --------------------------------
      DO iBody = 1, nBody

        SELECT CASE (canonical_body_type(iBody))
        CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)

          SELECT CASE (boundary_motion_type(iBody))
          CASE (STATIONARY)

            SELECT CASE (wall_type(iBody))
            CASE(POROUS_OR_SLIP)
              CALL wall_velocity(iBody)

            CASE(NONPOROUS_AND_NONSLIP)
              DO m = 1, nPtsBodyMarkerOrig(iBody)
                uBodyMarker(iBody,m) =  0.0_CGREAL
                vBodyMarker(iBody,m) =  0.0_CGREAL
                wBodyMarker(iBody,m) =  0.0_CGREAL
              ENDDO
            END SELECT ! wall_type

          CASE (FORCED)
            DO m = 1, nPtsBodyMarkerOrig(iBody)
              uBodyMarker(iBody,m) = vxcent(iBody)                                       &
                + ( -angvz(iBody)*(yBodyMarker(iBody,m)-ycent(iBody)) )
              vBodyMarker(iBody,m) = vycent(iBody)                                       &
                - ( -angvz(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
              wBodyMarker(iBody,m) = vzcent(iBody) 
            ENDDO

          CASE (FLOW_INDUCED:PRESCRIBED)
            DO m = 1, nPtsBodyMarkerOrig(iBody)
              uBodyMarker(iBody,m) =  0.0_CGREAL
              vBodyMarker(iBody,m) =  0.0_CGREAL
              wBodyMarker(iBody,m) =  0.0_CGREAL
            ENDDO

          END SELECT ! boundary_motion_type

          CALL extend_cylinder_vel_3d(iBody)
   
        CASE(ELLIPSOID:UNSTRUCTURED_SURFACE)

          IF (monitorON) WRITE(STDOUT,'(5X,A)') 'Setting up MOVING ELLIPSOID OR UNSTRUCTURED SURFACE : SSM'

!         set marker velocity
!         -------------------
          SELECT CASE (boundary_motion_type(iBody))
          CASE (STATIONARY)

            SELECT CASE (wall_type(iBody))
            CASE(POROUS_OR_SLIP)
              CALL wall_velocity(iBody)

            CASE(NONPOROUS_AND_NONSLIP)
              DO m = 1, nPtsBodyMarker(iBody)
                uBodyMarker(iBody,m) =  0.0_CGREAL
                vBodyMarker(iBody,m) =  0.0_CGREAL
                wBodyMarker(iBody,m) =  0.0_CGREAL
              ENDDO
            END SELECT ! wall_type

          CASE (FORCED)
            DO m = 1, nPtsBodyMarker(iBody)
              uBodyMarker(iBody,m)  = vxcent(iBody)   &
                + ( angvy(iBody)*(zBodyMarker(iBody,m)-zcent(iBody))  &
                - angvz(iBody)*(yBodyMarker(iBody,m)-ycent(iBody)) )
              vBodyMarker(iBody,m)  = vycent(iBody)   &
                - ( angvx(iBody)*(zBodyMarker(iBody,m)-zcent(iBody))  &
                - angvz(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
              wBodyMarker(iBody,m)  = vzcent(iBody)   &
                + ( angvx(iBody)*(yBodyMarker(iBody,m)-ycent(iBody))   &
                - angvy(iBody)*(xBodyMarker(iBody,m)-xcent(iBody)) )
            ENDDO

          CASE (FLOW_INDUCED:PRESCRIBED)
            DO m = 1, nPtsBodyMarker(iBody)
              uBodyMarker(iBody,m) =  0.0_CGREAL
              vBodyMarker(iBody,m) =  0.0_CGREAL
              wBodyMarker(iBody,m) =  0.0_CGREAL
            ENDDO
          END SELECT ! boundary_motion_type

        END SELECT ! canonical_body_type
      ENDDO ! iBody

!     write body markers into Output file
!     -----------------------------------
      IF (ImtheBOSS) WRITE(ifuBodyOut,'(3I13,2x,A)') nBody, nBody_Solid, nBody_Membrane, &
                                                  'nBody, nBody_Solid, nBody_Membrane'

      DO iBody = 1, nBody

        IF (ImtheBOSS) THEN
!         writing out unstructured surface
!         --------------------------------
          WRITE(ifuUnstrucSurfOut,*)
          WRITE(ifuUnstrucSurfOut,*) nPtsBodyMarker(iBody), totNumTriElem(iBody)
          WRITE(ifuUnstrucSurfOut,*)
          DO m=1,nPtsBodyMarker(iBody)
            WRITE(ifuUnstrucSurfOut,*) i,xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
          ENDDO
          WRITE(ifuUnstrucSurfOut,*)
          DO  j=1,totNumTriElem(iBody)
            WRITE(ifuUnstrucSurfOut,*) i,triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
          ENDDO
          WRITE(ifuUnstrucSurfOut,*)
          WRITE(ifuUnstrucSurfOut,*)pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)
        END IF

!       writing out body parameter file
!       -------------------------------
        IF(canonical_body_type(iBody) <= GENERAL_CYLINDER) THEN
          canonical_body_type(iBody) = UNSTRUCTURED_SURFACE
          body_dim(iBody)            = BODY_DIM2
          radiusz(iBody)             = 0.0_CGREAL
          zcent(iBody)               = 0.0_CGREAL
        ELSE 
          IF (canonical_body_type(iBody) == ELLIPSOID) canonical_body_type(iBody) = UNSTRUCTURED_SURFACE
        ENDIF

        IF (ImtheBOSS) THEN
          WRITE(ifuBodyOut,'(4I13        ,A)') canonical_body_type(iBody) , body_dim(iBody),                &
                                               boundary_motion_type(iBody), membrane_type(iBody),           &
                                               '  body_type, body_dim, motion_type, membrane_type'

          WRITE(ifuBodyOut,'(1I13        ,A)') wall_type(iBody),                                            &
                                               '  wall_type'

          WRITE(ifuBodyOut,'(2I13        ,A)') nPtsBodyMarker(iBody), totNumTriElem(iBody),                 &
                                               '  nPtsGCMBodyMarker, nTriElement'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') radiusx(iBody), radiusy(iBody), radiusz(iBody),              &
                                               '  radiusx, radiusy, radiusz'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') xcent(iBody), ycent(iBody), zcent(iBody),                    &
                                               '  xcent, ycent, zcent'

          WRITE(ifuBodyOut,'(  F21.15,5X ,A)') alpha(iBody),                                                &
                                               ' alpha'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') vxcentTrans(iBody), vycentTrans(iBody), vzcentTrans(iBody),  &
                                               ' vxcentTrans, vycentTrans, vzcentTrans'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') ampx(iBody), ampy(iBody), ampz(iBody),                       &
                                               ' ampx, ampy, ampz'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') freqx(iBody), freqy(iBody), freqz(iBody),                    &
                                               'freqx, freqy, freqz'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') angvx(iBody), angvy(iBody), angvz(iBody),                    &
                                               'angvx, angvy, angvz'

          WRITE(ifuBodyOut,'(  F21.15,5X ,A)') phase(iBody),                                                &
                                               'phase'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') ampangx(iBody), ampangy(iBody), ampangz(iBody),              &
                                               'ampangx, ampangy, ampangz'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') freqangx(iBody), freqangy(iBody), freqangz(iBody),           &
                                               ' freqangx, freqangy, freqangz'

          WRITE(ifuBodyOut,'(2I13        ,A)') mMinWallVel(iBody), mMaxWallVel(iBody),                      &
                                               ' mMinWallVel, mMaxWallVel'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') ampVelx(iBody), ampVelY(iBody), ampVelZ(iBody),              &
                                               ' ampVelx, ampVelY, ampVelZ'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') freqVelx(iBody), freqVelY(iBody), freqVelZ(iBody),           &
                                               ' freqVelx, freqVelY, freqVelZ'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') phaseVelx(iBody), phaseVelY(iBody), phaseVelZ(iBody),        &
                                               ' phaseVelx, phaseVelY, phaseVelz'

          WRITE(ifuBodyOut,'(2(F21.15,5X),A)') density_fluid, density_solid(iBody),                         &
                                               ' Density_Fluid, Density_Solid'

          WRITE(ifuBodyOut,'(3(F21.15,5X),A)') xcentConstr(iBody), ycentConstr(iBody), zcentConstr(iBody),  &
                                               ' Centroid Translation in X,Y,Z direction'
          IF(ibody == nBody) THEN
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,'(A)') 'body_type'
           WRITE(ifuBodyOut,'(A)') '        1: ELLIPTIC CYLINDER'
           WRITE(ifuBodyOut,'(A)') '        2: GENERAL CYLINDER'
           WRITE(ifuBodyOut,'(A)') '        3: ELLIPSOID'
           WRITE(ifuBodyOut,'(A)') '        4: UNSTRUCTURED SURFACE'
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,'(A)') 'body_dim'
           WRITE(ifuBodyOut,'(A)') '        2: 2D'
           WRITE(ifuBodyOut,'(A)') '        3: 3D'
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,'(A)') 'motion_type'
           WRITE(ifuBodyOut,'(A)') '        0: STATIONARY'
           WRITE(ifuBodyOut,'(A)') '        1: FORCED'
           WRITE(ifuBodyOut,'(A)') '        2: FLOW_INDUCED'
           WRITE(ifuBodyOut,'(A)') '        3: PRESCRIBED'
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,'(A)') 'membrane_type'
           WRITE(ifuBodyOut,'(A)') '        1: OPEN_MEMBRANE'
           WRITE(ifuBodyOut,'(A)') '        2: CLOSED_MEMBRANE'
           WRITE(ifuBodyOut,*)
           WRITE(ifuBodyOut,'(A)') 'wall_type'
           WRITE(ifuBodyOut,'(A)') '        0: NONPOROUS_AND_NONSLIP'
           WRITE(ifuBodyOut,'(A)') '        1: POROUS_OR_SLIP'
		  ENDIF !ibody == nbody
        END IF
      ENDDO

      DO iBody = 1, nBody
        IF(body_type_orig(iBody) <= ELLIPSOID ) THEN
          IF (ImtheBOSS) THEN
            PRINT*,'######################CODE NEEDS TO BE RERUN ##########################################'
            PRINT*,'Body parameters have been written out in      :        canonical_body_out.dat'
            PRINT*,'Cylinder Surface(s) have been written out in  :  unstructured_surface_out.dat'
            PRINT*,'execute the following statements'
            PRINT*,'cp unstruc_surface_out.dat  unstruc_surface_in.dat'
            PRINT*,'cp canonical_body_out.dat   canonical_body_in.dat'
            PRINT*,'Now run code again'
          END IF

          CALL flow_stop()
		  call mpi_finalize(iErr)
          STOP
        ENDIF
      ENDDO

    END IF
   
!   ------------------------------------------------------------------------------
!    SAMK:
!    From this point on the only body type the code knows is USTRUCTURED SURFACE.
!    This means canonical_body_type()=(always)=UNSTRUCTURED_SURFACE
!    In the future a pre-processor will be developed and then the code will only
!    know USTRUCTURED SURFACE.
!   ------------------------------------------------------------------------------


    IF (ImtheBOSS) THEN
      OPEN(ifuMarkerOut, FILE='marker_out.dat', STATUS='UNKNOWN')
      DO iBody = 1, nBody
        WRITE(ifuMarkerOut,*) 'Body:',iBody
        DO m = 1, nPtsBodyMarkerOrig(iBody)
          WRITE(ifuMarkerOut,*) xBodyMarker(iBody,m), yBodyMarker(iBody,m), zBodyMarker(iBody,m) 
        ENDDO
      ENDDO
    END IF

! Read edge node information for open membranes
!    openMemb = .FALSE.
!    DO iBody = 1, nBody
!      IF (nBody_membrane > 0 .AND. iBody > nbody_solid .AND. membrane_type(iBody) == OPEN_MEMBRANE) openMemb = .TRUE.
!    ENDDO
!    IF(openMemb == .TRUE.) &

!        OPEN(ifuOpenMembraneEdgeIn, FILE='open_membrane_edge_in.dat', STATUS='UNKNOWN')

!    DO iBody = 1, nBody
!      IF (nBody_membrane > 0 .AND. iBody > nbody_solid .AND. membrane_type(iBody) == OPEN_MEMBRANE) THEN
!        WRITE(STDOUT,'(5X,A)') 'Reading edge nodes'
!        READ(ifuOpenMembraneEdgeIn,*)num_edge
!        IF (num_edge > 0) THEN
!          DO i = 1,num_edge
!            READ(ifuOpenMembraneEdgeIn,*)edge_node_number
!            edge_node(iBody,edge_node_number) = 1
!          ENDDO

!        ENDIF
!      ENDIF
!    ENDDO

    IF (monitorON) THEN
      WRITE(STDOUT,'(5X,A)') 'Get Edge Node...'
    END IF

    DO iBody = 1, nBody
      IF (nBody_membrane > 0 .AND. iBody > nbody_solid .AND. membrane_type(iBody) == OPEN_MEMBRANE) THEN
        ALLOCATE(edge_node_number(nPtsBodyMarker(iBody)))
        edge_node_number = 0
        CALL get_edge_node(iBody,nPtsBodyMarker(iBody),num_edge,edge_node_number)
		
        IF (monitorON) THEN
         WRITE(*,*) iBody, 'OK', num_edge
        END IF		

        IF (num_edge > 0) THEN
            DO i = 1,num_edge
              edge_node(iBody,edge_node_number(i)) = 1
            ENDDO
        ENDIF

       DEALLOCATE(edge_node_number)

     ENDIF
    ENDDO

    IF (monitorON) THEN
      WRITE(STDOUT,'(5X,A)') 'Calculating solid surface parameters'
    END IF
	
    CALL calculate_arclength_norm_ds()

    IF (ImtheBOSS) CLOSE(ifuMarkerTRI)

END SUBROUTINE initialize_marker 
!---------------------------------------------------------------------

SUBROUTINE get_edge_node(iBody,nv,tot_num_edge,num_of_edge_node)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

IMPLICIT NONE
INTEGER, INTENT(IN) :: iBody, nv
INTEGER, INTENT(OUT) :: tot_num_edge
INTEGER, DIMENSION(nv), INTENT(INOUT) :: num_of_edge_node
INTEGER, ALLOCATABLE, DIMENSION(:) :: is_index, v2nv, v2he
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tris
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: opphes
INTEGER :: nvpe, nepe,  ntris, ne, v, vn, found, opp, ifound
INTEGER ::i, j, ii, jj, heid2fid, heid2leid, index, junk

nepe = 3
ntris = totNumTriElem(iBody)
ne = ntris*nepe
tot_num_edge = 0

ALLOCATE(is_index(nv+1))
ALLOCATE(tris(ntris,3))
ALLOCATE(v2nv(ne))
ALLOCATE(v2he(ne))
ALLOCATE(opphes(ntris,3))


DO i = 1, totNumTriElem(iBody)
  DO j = 1, 3
  tris(i,j) = triElemNeig(iBody, j, i)
  ENDDO
ENDDO
!_________________________________________________________
!Build is_index to store starting position for each vertex

is_index = 0
DO ii = 1, ntris
   DO i = 1, 3
   is_index(tris(ii,i)+1) = is_index(tris(ii,i)+1) + 1
   ENDDO
ENDDO
is_index(1) = 1
DO ii = 1, nv
   is_index(ii+1) = is_index(ii) + is_index(ii+1)
ENDDO

v2nv = 0    !Vertex to next vertex in each halfedge.
v2he = 0    !Vertex to half-edge
DO ii = 1, ntris
   DO i = 1, 3

     IF ( i == 3) THEN
     j = 1
     ELSE
     j = i+1
     ENDIF

   v2nv(is_index(tris(ii, i))) = tris(ii, j)
   v2he(is_index(tris(ii, i))) = 4*ii - 1 + i
   is_index(tris(ii,i)) = is_index(tris(ii,i)) + 1

   ENDDO
ENDDO

DO i = nv, 2, -1
is_index(i) = is_index(i-1)
ENDDO
is_index(1) = 1
!________________________________________________________
!Set opphes
opphes = 0
DO ii = 1, ntris

   DO jj = 1, 3

     IF ( jj == 3) THEN
     j = 1
     ELSE
     j = jj+1
     ENDIF

   v = tris(ii, jj)
   vn = tris(ii,j)

! LOCATE : Locate index col in v2nv(first to last)
    DO index = is_index(vn), is_index(vn+1)-1

        IF(v2nv(index) == v)THEN

        opp = v2he(index)
        opphes(ii,jj) = opp

        heid2fid = aint(real(opp)/4)
        heid2leid = mod(opp,4)+1

        opphes(heid2fid, heid2leid) = ii*4 + jj - 1

        ENDIF
    ENDDO

   ENDDO

ENDDO


DO ii = 1, ntris
   DO jj = 1, 3

   if(opphes(ii,jj) == 0)THEN
   tot_num_edge = tot_num_edge + 1
   num_of_edge_node(tot_num_edge) = tris(ii, jj)
   ENDIF

   ENDDO
ENDDO

DEALLOCATE(is_index)
DEALLOCATE(tris)
DEALLOCATE(v2nv)
DEALLOCATE(v2he)
DEALLOCATE(opphes)

END SUBROUTINE get_edge_node

!--------------------------------------------------------------------

SUBROUTINE extend_cylinder_3D(iBody)

! --------------------------------------------------------------------
!  This subroutine extends a cicular cylinder to make pseudo-3D body.
! --------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iBody

!   Loop variables
!   --------------
    INTEGER :: j, k, m, midx, mp
    INTEGER :: nodesMax

    REAL(KIND=CGREAL) :: dist_beg_end_marker, ds_marker



!   Determine if body if open or closed
!   -----------------------------------
    dist_beg_end_marker = SQRT(  (xBodyMarker(iBody,1)-xBodyMarker(iBody,nPtsBodymarkerOrig(iBody)))**2 &
                                +(yBodyMarker(iBody,1)-yBodyMarker(iBody,nPtsBodymarkerOrig(iBody)))**2 )
    ds_marker = SQRT(  (xBodyMarker(iBody,2)-xBodyMarker(iBody,1))**2 &
                      +(yBodyMarker(iBody,2)-yBodyMarker(iBody,1))**2   ) 
    IF (dist_beg_end_marker > 2.0_CGREAL*ds_marker) THEN
      IF (monitorON) WRITE(STDOUT,'(5X,A)'),'*** WARINING: Body seems to be an "open" body'
      nodesMax             = nPtsBodymarkerOrig(iBody)-1
      totNumTriElem(iBody) = 2*(nPtsBodyMarkerOrig(iBody)-1)*(nz-1)
    ELSE
      nodesMax = nPtsBodymarkerOrig(iBody)
    ENDIF

!   assign nodal cooridinates
!   -------------------------
    DO k=2,nz
      DO m=1,nPtsBodymarkerOrig(iBody)
        midx = (k-1)*nPtsBodymarkerOrig(iBody)+m
        xBodyMarker(iBody,midx) = xBodyMarker(iBody,m)
        yBodyMarker(iBody,midx) = yBodyMarker(iBody,m)
        zBodyMarker(iBody,midx) = z(k)
      ENDDO
    ENDDO

!   create neighbor list for elements
!   ---------------------------------
    j = 0
    DO k = 1, nz-1
      DO m = 1, nodesMax
        j = j+1

        mp = m+1
        IF (m == nPtsBodyMarkerOrig(iBody)) mp = 1

        triElemNeig(iBody,1,j) = (k-1)*nPtsBodymarkerOrig(iBody) + m
        triElemNeig(iBody,2,j) =  k*nPtsBodymarkerOrig(iBody)    + m
        triElemNeig(iBody,3,j) = (k-1)*nPtsBodymarkerOrig(iBody) + mp

        j = j+1
        triElemNeig(iBody,1,j) = (k-1)*nPtsBodymarkerOrig(iBody) + mp
        triElemNeig(iBody,2,j) =  k*nPtsBodymarkerOrig(iBody)    + m
        triElemNeig(iBody,3,j) =  k*nPtsBodymarkerOrig(iBody)    + mp

      ENDDO
    ENDDO

    IF (ImtheBOSS) THEN
      WRITE(ifuMarkerTRI,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(ifuMarkerTRI,*) 'VARIABLES="X","Y","Z"'
      WRITE(ifuMarkerTRI,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'F=FEPOINT, ET=TRIANGLE'
      DO m = 1,nPtsBodyMarker(iBody)
        WRITE(145,*) xBodyMarker(iBody,m),yBodyMarker(iBody,m),zBodyMarker(iBody,m)
      ENDDO
      DO j=1,totNumTriElem(iBody)
        WRITE(145,*) triElemNeig(iBody,1,j),triElemNeig(iBody,2,j),triElemNeig(iBody,3,j)
      ENDDO
    END IF  ! ImtheBOSS

    pointOutsideBodyX(iBody) =-10.0_CGREAL
    pointOutsideBodyy(iBody) =-10.0_CGREAL
    pointOutsideBodyz(iBody) =-10.0_CGREAL

END SUBROUTINE extend_cylinder_3D
!---------------------------------------------------------------------



SUBROUTINE extend_cylinder_vel_3D(iBody)

! -----------------------------------------------------------------
!  This subroutine prescribes surface velocity for pseudo-3D body.
! -----------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iBody

!   Loop variables
!   --------------
    INTEGER :: k, m, midx



    DO k= 2,nz
      DO m= 1,nPtsBodymarkerOrig(iBody)
        midx = (k-1)*nPtsBodymarkerOrig(iBody)+m
        uBodyMarker(iBody,midx) = uBodyMarker(iBody,m)
        vBodyMarker(iBody,midx) = vBodyMarker(iBody,m)
        wBodyMarker(iBody,midx) = wBodyMarker(iBody,m) 
      ENDDO
    ENDDO

END SUBROUTINE extend_cylinder_vel_3D
!---------------------------------------------------------------------
