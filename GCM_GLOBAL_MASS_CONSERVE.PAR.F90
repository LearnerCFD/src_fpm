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
!  Filename: GCM_GLOBAL_MASS_CONSERVE.PAR.F90
!  Latest Modification: October 22, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian

!  Latest Modification: Oct, 20 2010 (ver. PAT 1.2.0)
!  by JHSeo
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_enforce_p_compatibility(pres)
!     GCM_enforce_global_mass_consv()
!     remove_rhs_adjust_bc()
!     rhs_readjust_bc()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1



SUBROUTINE GCM_enforce_p_compatibility(pres)

! ------------------------------------------------------------------------------------------------------
!  compute mass flux at all boundaries and adjust outflow BC so as to satisfy global mass conservation.
! ------------------------------------------------------------------------------------------------------

  USE global_parameters
  USE flow_parameters
  USE flow_arrays
  USE boundary_arrays
  USE grid_arrays
  USE cutcell_arrays

  IMPLICIT NONE

  REAL(KIND=CGREAL), INTENT(IN) :: pres(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

  INTEGER :: i,j,k
  INTEGER :: iG, jG


  REAL(KIND=CGREAL) :: pgradflux,correction_pgrad
  REAL(KIND=CGREAL) :: pgxw,pgxe,pgys,pgyn,pgzb,pgzf



  pgradflux = 0.0_CGREAL

  IF ( pbcx1 == PBC_DIRICHLET .OR. &                                 ! SER_TO_PAR. QX. CH23
       pbcx2 == PBC_DIRICHLET .OR. &
       pbcy1 == PBC_DIRICHLET .OR. &
       pbcy2 == PBC_DIRICHLET .OR. &
       pbcz1 == PBC_DIRICHLET .OR. &
       pbcz2 == PBC_DIRICHLET       )  RETURN
    
  IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
       bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
       bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
       bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
       bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
       bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

!   add pressure gradient to face velocities
!   ----------------------------------------
  
    DO k = 1,nzc
    DO j = 1,nyc
      jG = L2GJ(j)
    DO i = 1,nxc
      iG = L2GI(i)

      pgxw     = ( pres(i  ,j,k)                                                                    &
                 - pres(i-1,j,k)*(1.0_CGREAL - ium(i,j,k)) - pGhost(i-1,j,k)*ium(i,j,k))*dxcinv(iG)

      pgxe     = ( pres(i+1,j,k)*(1.0_CGREAL - iup(i,j,k)) + pGhost(i+1,j,k)*iup(i,j,k)             &
                 - pres(i  ,j,k)                                                       )*dxcinv(iG+1)

      pgys     = ( pres(i,j  ,k)                                                                    &
                 - pres(i,j-1,k)*(1.0_CGREAL - jum(i,j,k)) - pGhost(i,j-1,k)*jum(i,j,k) )*dycinv(jG)

      pgyn     = ( pres(i,j+1,k)*(1.0_CGREAL - jup(i,j,k)) + pGhost(i,j+1,k)*jup(i,j,k)             &
                 - pres(i,j  ,k)                                                      )*dycinv(jG+1)

      pgzb     = ( pres(i,j,k  )                                                                    &
                  -pres(i,j,k-1)*(1.0_CGREAL - kum(i,j,k)) - pGhost(i,j,k-1)*kum(i,j,k) )*dzcinv(k)

      pgzf     = ( pres(i,j,k+1)*(1.0_CGREAL - kup(i,j,k)) + pGhost(i,j,k+1)*kup(i,j,k)             &
                  -pres(i,j,k  )                                                      )*dzcinv(k+1)


      pgradflux = pgradflux +        &
                (-pgxw*dy(jG)*dz(k)*AREA_W(i,j,k)*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)   &
                 +pgxe*dy(jG)*dz(k)*AREA_E(i,j,k)*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)   &
                 -pgys*dx(iG)*dz(k)*AREA_S(i,j,k)*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)   &
                 +pgyn*dx(iG)*dz(k)*AREA_N(i,j,k)*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)   &
                 -pgzb*dx(iG)*dy(jG)*AREA_B(i,j,k)*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)   &
                 +pgzf*dx(iG)*dy(jG)*AREA_F(i,j,k)*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC))  &
                 *REAL(1-iblank_solid(i,j,k),KIND=CGREAL)


    ENDDO
    ENDDO
    ENDDO
   
 
#   ifdef MPI
      CALL par_getSumReal(pgradflux)
#   endif

!   adjust BC at outflow to satisfy global mass conservation
!   --------------------------------------------------------
    correction_pgrad =-pgradflux/outflow_area  
        

    IF (myCoords(1)==0       .AND. bcx1 == BC_TYPE_ZERO_GRADIENT) &
      pgradx1(myJLL:myJUL,0    :nzc+1) = pgradx1(myJLL:myJUL,0    :nzc+1) - correction_pgrad

    IF (myCoords(1)==Np(1)-1 .AND. bcx2 == BC_TYPE_ZERO_GRADIENT) &
      pgradx2(myJLL:myJUL,0    :nzc+1) = pgradx2(myJLL:myJUL,0    :nzc+1) + correction_pgrad

    IF (myCoords(2)==0       .AND. bcy1 == BC_TYPE_ZERO_GRADIENT) & 
      pgrady1(myILL:myIUL,0    :nzc+1) = pgrady1(myILL:myIUL,0    :nzc+1) - correction_pgrad      

    IF (myCoords(2)==Np(2)-1 .AND. bcy2 == BC_TYPE_ZERO_GRADIENT) &
      pgrady2(myILL:myIUL,0    :nzc+1) = pgrady2(myILL:myIUL,0    :nzc+1) + correction_pgrad 

    IF (                          bcz1 == BC_TYPE_ZERO_GRADIENT) &
      pgradz1(myILL:myIUL,myJLL:myJUL) = pgradz1(myILL:myIUL,myJLL:myJUL) - correction_pgrad 
    
    IF (                          bcz2 == BC_TYPE_ZERO_GRADIENT) &
      pgradz2(myILL:myIUL,myJLL:myJUL) = pgradz2(myILL:myIUL,myJLL:myJUL) + correction_pgrad 

  END IF


! Compute boundary source term for pressure Poisson equation
! ----------------------------------------------------------
  boundPresSource = 0.0_CGREAL
  
  IF(myCoords(1) == 0 ) THEN
   i = 1
   iG = L2GI(i)
   DO k = 1,nzc
   DO j = 1,nyc
    boundPresSource(i-1,j,k) = - dxinv(iG)*(-pgradx1(j,k))
   ENDDO
   ENDDO
 ENDIF

 IF(myCoords(1) == Np(1)-1) THEN
  i = nxc 
  iG = L2GI(i)
  DO k = 1,nzc
  DO j = 1,nyc
    boundPresSource(i+1,j,k) = - dxinv(iG)*( pgradx2(j,k))
  ENDDO
  ENDDO
 ENDIF

 IF(myCoords(2) == 0) THEN
  j = 1
  jG = L2GJ(j)
  DO k = 1,nzc
  DO i = 1,nxc
    boundPresSource(i,j-1,k) = - dyinv(jG)*(-pgrady1(i,k))
  ENDDO
  ENDDO
 ENDIF

 IF(myCoords(2)==Np(2)-1) THEN
  j = nyc 
  jG = L2GJ(j)
  DO k = 1,nzc
  DO i = 1,nxc
    boundPresSource(i,j+1,k) = - dyinv(jG)*( pgrady2(i,k))
  ENDDO
  ENDDO
 ENDIF

  DO j = 1,nyc
  DO i = 1,nxc
    k = 1
    boundPresSource(i,j,k-1) = - dzinv(k)*(-pgradz1(i,j))
    k = nzc
    boundPresSource(i,j,k+1) = - dzinv(k)*( pgradz2(i,j))
 ENDDO
 ENDDO
 
   

END SUBROUTINE GCM_enforce_p_compatibility
!-------------------------------------------------------------------------



SUBROUTINE GCM_enforce_global_mass_consv()

! ------------------------------------------------------------------------------------------------------
!  compute mass flux at all boundaries and adjust outflow BC so as to satisfy global mass conservation.
! ------------------------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
	USE unstructured_surface_arrays
	USE cutcell_arrays	

    IMPLICIT NONE

    INTEGER :: i, j, k, n
    INTEGER :: iG, jG

    REAL(KIND=CGREAL) :: massflux, correction_vel, fluxBody
	
    INTEGER :: cElem, iBody, iKill
	REAL(KIND=CGREAL) xBI,yBI,zBI,uBI,vBI,wBI,xBIN,yBIN,zBIN,USIGMA,xp,yp,zp



    massflux = 0.0_CGREAL

    IF ( pbcx1 == PBC_DIRICHLET .OR. &                                 ! SER_TO_PAR. QX. CH23
         pbcx2 == PBC_DIRICHLET .OR. &
         pbcy1 == PBC_DIRICHLET .OR. &
         pbcy2 == PBC_DIRICHLET .OR. &
         pbcz1 == PBC_DIRICHLET .OR. &
         pbcz2 == PBC_DIRICHLET       )  RETURN    

    IF ( bcx1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcx2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcy1 == BC_TYPE_ZERO_GRADIENT .OR. & 
         bcy2 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz1 == BC_TYPE_ZERO_GRADIENT .OR. &
         bcz2 == BC_TYPE_ZERO_GRADIENT       ) THEN

      CALL face_vel()

      DO k=1,nzc
      DO j=1,nyc
        jG=L2GJ(j)
        
        DO i=1,nxc
          iG=L2GI(i)
          
		  IF(iCC) call Cutcell_norm_vel(i,j,k,usigma)
		  
		  iKill = isolid(i,j,k)*iCC + iblank(i,j,k)*(1-iCC)
          
          massflux = massflux +                      &
                   ( -face_uw(i,j,k)*dy(jG)*dz(k )*AREA_W(i,j,k)   &
                     +face_ue(i,j,k)*dy(jG)*dz(k )*AREA_E(i,j,k)   &
                     -face_vs(i,j,k)*dx(iG)*dz(k )*AREA_S(i,j,k)   &
                     +face_vn(i,j,k)*dx(iG)*dz(k )*AREA_N(i,j,k)   &
                     -face_wb(i,j,k)*dx(iG)*dy(jG)*AREA_B(i,j,k)   &
                     +face_wf(i,j,k)*dx(iG)*dy(jG)*AREA_F(i,j,k)   &
					 +CUT_AREA(i,j,k)*usigma  ) &
                  *REAL(1-iKill,KIND=CGREAL)
				  
        ENDDO ! i
      ENDDO ! j
      ENDDO ! k
      
#ifdef MPI
      CALL par_GetSumReal(massflux)
#endif

!     adjust BC at outflow to satisfy global mass conservation
!     --------------------------------------------------------
      correction_vel =-massflux/outflow_area  

      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A,D15.5)') ' SET_BC:massflux       = ',massflux
        WRITE(STDOUT,'(5X,A,D15.5)') ' SET_BC:outflow_area   = ',outflow_area
        WRITE(STDOUT,'(5X,A,D15.5)') ' SET_BC:correction_vel = ',correction_vel
        WRITE(STDOUT,*)
      END IF ! ntime

      IF (bcx1 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(1) == 0) THEN
        DO k=1,nzc
        DO j=1,nyc
          bcxu(1,j,k) = bcxu(1,j,k) - correction_vel
        ENDDO ! j
        ENDDO ! k
      END IF

      IF (bcx2 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(1) == Np(1)-1) THEN
        DO k=1,nzc
        DO j=1,nyc
          bcxu(nxc,j,k) = bcxu(nxc,j,k) + correction_vel
        ENDDO ! j
        ENDDO ! k
      END IF

      IF (bcy1 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(2) == 0) THEN
        DO k=1,nzc
        DO i=1,nxc
          bcyv(i,1,k) = bcyv(i,1,k) - correction_vel
        ENDDO ! i
        ENDDO ! k
      END IF

      IF (bcy2 == BC_TYPE_ZERO_GRADIENT .AND. myCoords(2) == Np(2)-1) THEN
        DO k=1,nzc
        DO i=1,nxc
          bcyv(i,nyc,k) = bcyv(i,nyc,k) + correction_vel
        ENDDO ! i
        ENDDO ! k
      END IF

      IF (bcz1 == BC_TYPE_ZERO_GRADIENT ) THEN
        DO j=1,nyc
        DO i=1,nxc
          bczw(i,j,1) = bczw(i,j,1) - correction_vel
        ENDDO ! i
        ENDDO ! j
      END IF

      IF (bcz2 == BC_TYPE_ZERO_GRADIENT ) THEN
        DO j=1,nyc
        DO i=1,nxc
          bczw(i,j,nzc) = bczw(i,j,nzc) + correction_vel
        ENDDO ! i
        ENDDO ! j
      END IF

      CALL set_outer_ghost_vel()

    ENDIF ! bcx1

END SUBROUTINE GCM_enforce_global_mass_consv
!-------------------------------------------------------------------------
