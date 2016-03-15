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
!  Filename: BOUNDARY_SET_ARCLENGTH_NORM.PAR.F90
!  Latest Modification: December 07, 2007 (ver. 0.3.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     calculate_arclength_norm_ds()
!     determine_norm_dir_unstruc(iBody)
! --------------------------------------------------------------------



SUBROUTINE calculate_arclength_norm_ds()

! --------------------------------------------------------------------
!  This subroutine calculates the solid surface normal vectors.
!  The surface normal vectors point outwards the flow field.
!  For membranes, normal vectors point from iblank=0 toward iblank=1.
! --------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!   Loop variables
!   --------------
    INTEGER :: iBody, iElem

    REAL(KIND=CGREAL) :: normMag, tangMag
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    REAL(KIND=CGREAL) :: dxMin, dyMin, dzMin, cellAreaMin, triElemAreaAv, bodyResolutionNormalized



!   Check direction of surface normal and adjust ordering of triangle if
!   normal does not point into body this is done only during fresh start.
!   ---------------------------------------------------------------------
    IF (ntime == 0 .AND. nread == 0) THEN
      IF (monitorON) THEN
        WRITE(STDOUT,'(7X,A)') 'Sorting body database with regard to the point ouside'
        WRITE(STDOUT,'(7X,A)') '====================================================='
        WRITE(STDOUT,*)
      END IF

      DO iBody = 1, nBody
        IF (monitorON) THEN
          WRITE(STDOUT,'(7X,A,1X,I3.3)') 'Body number',iBody
          WRITE(STDOUT,'(7X,A)') '---------------'
        END IF

        CALL determine_norm_dir_unstruc(iBody)
      ENDDO ! iBody

      IF (monitorON) WRITE(STDOUT,*)
    ENDIF
 
    sBodyMarker = 0.0_CGREAL
    surfArea    = 0.0_CGREAL
    
!   Surface Normal Vector
!   ---------------------
    IF (monitorON) THEN
      WRITE(STDOUT,'(7X,A)') 'Calculating normal surface vector'
      WRITE(STDOUT,'(7X,A)') '================================='     
      WRITE(STDOUT,*)
    END IF

    DO iBody = 1, nBody

      IF (monitorON) THEN
        WRITE(STDOUT,'(7X,A,1X,I3.3)') 'Body number',iBody
        WRITE(STDOUT,'(7X,A)') '---------------'
      END IF

!     Sweep surface of each triangular element
!     ----------------------------------------
      surfArea(iBody)  = 0.0_CGREAL

      DO iElem=1,totNumTriElem(iBody)
                                    
!       First vector of each element
!       ----------------------------
        triElemVectorx(1)= xBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                          -xBodyMarker(iBody,triElemNeig(iBody,1,iElem))
        triElemVectory(1)= yBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                          -yBodyMarker(iBody,triElemNeig(iBody,1,iElem))
        triElemVectorz(1)= zBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                          -zBodyMarker(iBody,triElemNeig(iBody,1,iElem))

!       Second vector of each element
!       -----------------------------
        triElemVectorx(2)= xBodyMarker(iBody,triElemNeig(iBody,3,iElem)) &
                          -xBodyMarker(iBody,triElemNeig(iBody,2,iElem))
        triElemVectory(2)= yBodyMarker(iBody,triElemNeig(iBody,3,iElem)) &
                          -yBodyMarker(iBody,triElemNeig(iBody,2,iElem))
        triElemVectorz(2)= zBodyMarker(iBody,triElemNeig(iBody,3,iElem)) &
                          -zBodyMarker(iBody,triElemNeig(iBody,2,iElem))

!       Normal vector of the element
!       ----------------------------
        triElemNormx(iBody,iElem)= (triElemVectory(1)*triElemVectorz(2)  &
                                   -triElemVectorz(1)*triElemVectory(2))
        triElemNormy(iBody,iElem)= (triElemVectorz(1)*triElemVectorx(2)  &
                                   -triElemVectorx(1)*triElemVectorz(2))
        triElemNormz(iBody,iElem)= (triElemVectorx(1)*triElemVectory(2)  &
                                   -triElemVectory(1)*triElemVectorx(2))


        IF (body_dim(iBody) == BODY_DIM2) THEN
          IF( ABS(triElemNormz(iBody,iElem)) > 1.0E-6_CGREAL ) THEN
            IF( ABS(triElemNormz(iBody,iElem)) > 1.0E-4_CGREAL ) THEN
              WRITE(STDOUT,'(9X,A,2(1X,I6))') 'Body:Element = ',ibody,iElem
              WRITE(STDOUT,'(9X,A,3(1X,1PE15.7))') 'triElemNormx,y,z = ',&
                triElemNormx(iBody,iElem),triElemNormy(iBody,iElem),triElemNormz(iBody,iElem),&
                triElemNeig(iBody,1,iElem),triElemNeig(iBody,2,iElem),triElemNeig(iBody,3,iElem)
              WRITE(STDOUT,'(9X,A)') '2D Body Element has normal in the z-direction'
              WRITE(STDOUT,'(9X,A)') 'Aborting'

              CALL flow_stop()
              STOP
            ELSE
              IF (monitorON) THEN
                WRITE(STDOUT,'(9X,A,2(1X,I6))') 'Body:Element = ',ibody,iElem
                WRITE(STDOUT,'(9X,A,3(1X,1PE15.7))') 'triElemNormx,y,z = ',&
                  triElemNormx(iBody,iElem),triElemNormy(iBody,iElem),triElemNormz(iBody,iElem)
                WRITE(STDOUT,'(9X,A)') 'Zeroing out the small z-normal'
              END IF
              triElemNormz(iBody,iElem) = 0.0_CGREAL
            ENDIF
          ENDIF
        ENDIF
                 
        normMag  = SQRT(triElemNormx(iBody,iElem)**2   &
                       +triElemNormy(iBody,iElem)**2   &
                       +triElemNormz(iBody,iElem)**2)

!       Area of element
!       ---------------
        triElemArea(iBody,iElem)  = 0.5_CGREAL*normMag
        surfArea(iBody)           = surfArea(iBody) + triElemArea(iBody,iElem)

!       Unit Normal vector
!       ------------------
        triElemNormx(iBody,iElem) = triElemNormx(iBody,iElem)/normMag
        triElemNormy(iBody,iElem) = triElemNormy(iBody,iElem)/normMag
        triElemNormz(iBody,iElem) = triElemNormz(iBody,iElem)/normMag

!       Unit Tangents
!       Tangent-2 defined parallel to vector from vertex-1 to vertex-2
!       --------------------------------------------------------------
        triElemTang2x(iBody,iElem)  = xBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                     -xBodyMarker(iBody,triElemNeig(iBody,1,iElem))
        triElemTang2y(iBody,iElem)  = yBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                     -yBodyMarker(iBody,triElemNeig(iBody,1,iElem))
        triElemTang2z(iBody,iElem)  = zBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                     -zBodyMarker(iBody,triElemNeig(iBody,1,iElem))

        tangMag     = SQRT(  triElemTang2x(iBody,iElem)**2   &
                           + triElemTang2y(iBody,iElem)**2   &
                           + triElemTang2z(iBody,iElem)**2 )

        triElemTang2x(iBody,iElem)  = triElemTang2x(iBody,iElem)/tangMag
        triElemTang2y(iBody,iElem)  = triElemTang2y(iBody,iElem)/tangMag
        triElemTang2z(iBody,iElem)  = triElemTang2z(iBody,iElem)/tangMag

!       t1 = t2 x n
!       -----------
        triElemTang1x(iBody,iElem)  =  triElemTang2y(iBody,iElem)*triElemNormz(iBody,iElem)  &
                                     - triElemTang2z(iBody,iElem)*triElemNormy(iBody,iElem)
        triElemTang1y(iBody,iElem)  =- triElemTang2x(iBody,iElem)*triElemNormz(iBody,iElem)  &
                                     + triElemTang2z(iBody,iElem)*triElemNormx(iBody,iElem)
        triElemTang1z(iBody,iElem)  =  triElemTang2x(iBody,iElem)*triElemNormy(iBody,iElem)  &
                                     - triElemTang2y(iBody,iElem)*triElemNormx(iBody,iElem)

!       Centroid of the each element
!       ----------------------------
        triElemCentx(iBody,iElem)=(xBodyMarker(iBody,triElemNeig(iBody,1,iElem)) &
                                  +xBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                  +xBodyMarker(iBody,triElemNeig(iBody,3,iElem)))/3.0_CGREAL
        triElemCenty(iBody,iElem)=(yBodyMarker(iBody,triElemNeig(iBody,1,iElem)) &
                                  +yBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                  +yBodyMarker(iBody,triElemNeig(iBody,3,iElem)))/3.0_CGREAL
        triElemCentz(iBody,iElem)=(zBodyMarker(iBody,triElemNeig(iBody,1,iElem)) &
                                  +zBodyMarker(iBody,triElemNeig(iBody,2,iElem)) &
                                  +zBodyMarker(iBody,triElemNeig(iBody,3,iElem)))/3.0_CGREAL                                 

      ENDDO ! iElem

      dxMin                    = MINVAL(dx(1:nxc))
      dyMin                    = MINVAL(dy(1:nyc))
      dzMin                    = MINVAL(dz(1:nzc))
      cellAreaMin              = (dxMin*dyMIn*dzMin)**(2.0_CGREAL/3.0_CGREAL)
      triElemAreaAv            = surfArea(iBody)/totNumTriElem(iBody)
      bodyResolutionNormalized = cellAreaMin/triElemAreaAv

      IF (monitorON) THEN
        WRITE(STDOUT,'(7X,A,1PE15.7)') 'Min Cell Area             = ',cellAreaMin
        WRITE(STDOUT,'(7X,A,1PE15.7)') 'Surface area of body      = ',surfArea(iBody)
        WRITE(STDOUT,'(7X,A,1PE15.7)') 'Average Element Area      = ',triElemAreaAv
        WRITE(STDOUT,'(7X,A,1PE15.7)') 'Norm. Surface resolution  = ',bodyResolutionNormalized
        WRITE(STDOUT,*)
      END IF
    ENDDO ! iBody


END SUBROUTINE calculate_arclength_norm_ds
!---------------------------------------------------------------------



SUBROUTINE determine_norm_dir_unstruc(iBody)

! -----------------------------------------------------------------------------------------
!  This subroutine determines if for a given unstructured surface mesh, the normal
!  vectors point inward or outward the flow field and reorders vertex numbers accordingly.
! -----------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)::iBody

!   Loop variables
!   --------------
    INTEGER :: iNode, iElem
    INTEGER :: cNode, cElem

!   Local variables
!   ---------------
    INTEGER           :: tmpNode

    REAL(KIND=CGREAL) :: Dnode, Dmin
    REAL(KIND=CGREAL) :: vectorx, vectory, vectorz, vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2), triElemVectorY(2), triElemVectorZ(2)
    REAL(KIND=CGREAL) :: cTriElemNormx, cTriElemNormy, cTriElemNormz
    REAL(KIND=CGREAL) :: cTriElemCentx, cTriElemCenty, cTriElemCentz



    cNode = 0

!   Find node closest to outside point
!   ----------------------------------
    Dmin = 1.0E8_CGREAL
    DO iNode= 1, nPtsBodyMarker(iBody)
      Dnode= SQRT( (xBodyMarker(iBody,iNode)-pointOutsideBodyX(iBody))**2  &
                  +(yBodyMarker(iBody,iNode)-pointOutsideBodyY(iBody))**2  &
                  +(zBodyMarker(iBody,iNode)-pointOutsideBodyZ(iBody))**2   )
      IF(Dnode <= Dmin) THEN
        Dmin = Dnode
        cNode= iNode
      ENDIF
    ENDDO

!   Find the first element having cNode as a vertex
!   -----------------------------------------------
    DO iElem= 1, totNumTriElem(iBody)
      IF ( triElemNeig(iBody,1,iElem) == cNode .OR. &
           triElemNeig(iBody,2,iElem) == cNode .OR. &
           triElemNeig(iBody,3,iElem) == cNode       ) THEN

           cElem= iElem
           EXIT   ! The first element will be good, no need to search more (SAMK).
      END IF
    ENDDO

    IF (monitorON) WRITE(STDOUT,'(7X,A,I8)') 'Closest Node to Outside Point = ', cNode
    IF (monitorON) WRITE(STDOUT,'(7X,A,I8)') 'Element with Closest Node     = ', cElem


!     1 *-------------* 3
!        \           /
!         \         /
!          \       /
!           \     /
!            \   /
!             \ /
!            2 *
!
!   Sweep surface of CLOSEST element to determine normal direction of triangular elements
!   -------------------------------------------------------------------------------------
    triElemVectorx(1)= xBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                      -xBodyMarker(iBody,triElemNeig(iBody,1,cElem))
    triElemVectory(1)= yBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                      -yBodyMarker(iBody,triElemNeig(iBody,1,cElem))       ! V(1-2)
    triElemVectorz(1)= zBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                      -zBodyMarker(iBody,triElemNeig(iBody,1,cElem))

    triElemVectorx(2)= xBodyMarker(iBody,triElemNeig(iBody,3,cElem)) &
                      -xBodyMarker(iBody,triElemNeig(iBody,2,cElem))
    triElemVectory(2)= yBodyMarker(iBody,triElemNeig(iBody,3,cElem)) &
                      -yBodyMarker(iBody,triElemNeig(iBody,2,cElem))       ! V(2-3)
    triElemVectorz(2)= zBodyMarker(iBody,triElemNeig(iBody,3,cElem)) &
                      -zBodyMarker(iBody,triElemNeig(iBody,2,cElem))

!   Normal of closest element
!   -------------------------
    cTriElemNormx = ( triElemVectory(1)*triElemVectorz(2)  &
                     -triElemVectorz(1)*triElemVectory(2) )
    cTriElemNormy = ( triElemVectorz(1)*triElemVectorx(2)  &
                     -triElemVectorx(1)*triElemVectorz(2) )                 ! V(1-2)xV(2-3)
    cTriElemNormz = ( triElemVectorx(1)*triElemVectory(2)  &
                     -triElemVectory(1)*triElemVectorx(2) )

!   Centroid of the closest element
!   -------------------------------
    cTriElemCentx = (xBodyMarker(iBody,triElemNeig(iBody,1,cElem)) &
                    +xBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                    +xBodyMarker(iBody,triElemNeig(iBody,3,cElem)))/3.0_CGREAL
    cTriElemCenty = (yBodyMarker(iBody,triElemNeig(iBody,1,cElem)) &
                    +yBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                    +yBodyMarker(iBody,triElemNeig(iBody,3,cElem)))/3.0_CGREAL
    cTriElemCentz = (zBodyMarker(iBody,triElemNeig(iBody,1,cElem)) &
                    +zBodyMarker(iBody,triElemNeig(iBody,2,cElem)) &
                    +zBodyMarker(iBody,triElemNeig(iBody,3,cElem)))/3.0_CGREAL

    vectorx = cTriElemCentX - pointOutsideBodyX(iBody)
    vectory = cTriElemCentY - pointOutsideBodyY(iBody)
    vectorz = cTriElemCentZ - pointOutsideBodyZ(iBody)

    vectProduct = vectorx*cTriElemNormx  &
                + vectory*cTriElemNormy  &
                + vectorz*cTriElemNormz

    IF (monitorON) WRITE(STDOUT,'(7X,A,1PE15.7)') 'vectProduct                   = ',vectproduct

    normDirFlag = 1.0_CGREAL

    IF (vectProduct < 0.0_CGREAL)  THEN
      normDirFlag =-1.0_CGREAL

      IF (monitorON) WRITE(STDOUT,'(7X,A)') 'Normal vectors point towards flow field: Reordering triangle vertices...'

      DO iElem= 1,totNumTriElem(iBody)         ! changing vertex ordering
        tmpNode = triElemNeig(iBody,2,iElem)
        triElemNeig(iBody,2,iElem) = triElemNeig(iBody,3,iElem)
        triElemNeig(iBody,3,iElem) = tmpNode
      ENDDO
    ENDIF

END SUBROUTINE determine_norm_dir_unstruc
!---------------------------------------------------------------------
