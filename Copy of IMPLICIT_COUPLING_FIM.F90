! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  The Johns Hopkins University
!
!  VICAR3Dp (ver. 1.5.9)
!
!  This is a continuously developing project.
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
!  Filename: IMPLICIT_COUPLING_FIM.F90
!  CREATED BY: Zhang, Chao, May 03, 2011
!
!  Latest Modification:  May 03, 2011
!  by Zhang, Chao
!
!  
! --------------------------------------------------------------------

SUBROUTINE FIM_implicit_coupling
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
IMPLICIT NONE
INTEGER :: iBody,iGroup
   
   CALL FIM_calculate_alphaimplicit
   
   DO iBody = 1,nbody
     CALL FIM_store_old_bodymarkers(iBody) !xold,uold = x,u broadcast
   END DO 
     
!============================Store original data at time n
   IF(kimplicit_FIM == 1)THEN
    
    CALL FIM_store_vicar_old_variables !uOld = u, Iblankold = iblank
!---Store all bodies' velocity       
    DO iBody = 1,nbody
     CALL FIM_store_body_old_variables(iBody)
     CALL FIM_store_location(iBody)
    END DO
!---Store all combined groups' velocity    
    DO iGroup = 1,nGroup_Combined
     CALL FIM_store_group_old_variables(iGroup)
    END DO
!---Store rotation matrix    
     DO iGroup = 1,nGroup_Combined     
      CALL FIM_store_rotation_matrix(iGroup)
     END DO       
!---Store all marker points for all bodies   
   DO iBody = 1,nbody
     CALL FIM_store_prev_bodymarkers(iBody) !xprev = x broadcast
   END DO    
    
   END IF


!============================Restore original data of time n for iteration

   IF(kimplicit_FIM > 1)THEN
   
     CALL FIM_restore_iblank_and_velocities !iblank = iblankold,u=uOld
!---Restore all bodies' velocity       
     DO iBody = 1,nbody
      CALL FIM_restore_body_old_variables(iBody)
      CALL FIM_restore_location(iBody)
     END DO
!---Restore all combined groups' velocity    
     DO iGroup = 1,nGroup_Combined
      CALL FIM_restore_group_old_variables(iGroup)
     END DO   
!---Restore rotation matrix    
     DO iGroup = 1,nGroup_Combined     
      CALL FIM_restore_rotation_matrix(iGroup)
     END DO        
!---Restore all marker points for all bodies   
     DO iBody = 1,nbody
      CALL FIM_restore_prev_bodymarkers(iBody) !x = xpre broadcast
     END DO  
              
     CALL FIM_calculate_new_iblank ! new values of iblank
	 
	 IF(iCC .eq. 1) THEN 
	  IF(monitorON) WRITE(STDOUT,'(3X,A)') 'Calc. CutCell fractions'	  
	  CALL init_CutCell	
      ENDIF ! iCC
	   
     CALL rhs_advec_diff() !Computing advection-diffusion RHS
     CALL set_solve_ad() ! advection eqs coefficents
     
  END IF !kimplicit_FIM > 1

END SUBROUTINE FIM_implicit_coupling



!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
!=====================Store information of previous time step
SUBROUTINE FIM_store_rotation_matrix(iGroup)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
IMPLICIT NONE
INTEGER :: iGroup

    C11_prev(iGroup) = C11(iGroup)
    C12_prev(iGroup) = C12(iGroup)
    C13_prev(iGroup) = C13(iGroup)
    C21_prev(iGroup) = C21(iGroup)
    C22_prev(iGroup) = C22(iGroup)
    C23_prev(iGroup) = C23(iGroup)
    C31_prev(iGroup) = C31(iGroup)
    C32_prev(iGroup) = C32(iGroup)
    C33_prev(iGroup) = C33(iGroup)

END SUBROUTINE FIM_store_rotation_matrix

SUBROUTINE FIM_store_force(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
USE usr_module
IMPLICIT NONE
INTEGER :: iBody

    sCx_old(iBody) = sCx(iBody)
    sCy_old(iBody) = sCy(iBody)
    sCz_old(iBody) = sCz(iBody)
    scmx_old(iBody) = scmx(iBody)
    scmy_old(iBody) = scmy(iBody)
    scmz_old(iBody) = scmz(iBody)

END SUBROUTINE FIM_store_force

SUBROUTINE FIM_store_location(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
USE usr_module
IMPLICIT NONE
INTEGER :: iBody

    xcent_old(iBody) = xcent(iBody)
    ycent_old(iBody) = ycent(iBody)
    zcent_old(iBody) = zcent(iBody)

END SUBROUTINE FIM_store_location


 
SUBROUTINE FIM_store_body_old_variables(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
IMPLICIT NONE
INTEGER :: iBody

   vxcentOld(iBody) = vxcent(iBody)
   vycentOld(iBody) = vycent(iBody)
   vzcentOld(iBody) = vzcent(iBody)
   angvxOld(iBody) = angvx(iBody)
   angvyOld(iBody) = angvy(iBody)
   angvzOld(iBody) = angvz(iBody)

END SUBROUTINE FIM_store_body_old_variables


SUBROUTINE FIM_store_group_old_variables(iGroup)
USE implicit_coupling_parameters 
USE flow_parameters 
IMPLICIT NONE
INTEGER :: iGroup

   vxcent_combinedOld(iGroup) = vxcent_combined(iGroup)
   vycent_combinedOld(iGroup) = vycent_combined(iGroup)
   vzcent_combinedOld(iGroup) = vzcent_combined(iGroup)
   angvx_combinedOld(iGroup) = angvx_combined(iGroup)
   angvy_combinedOld(iGroup) = angvy_combined(iGroup)
   angvz_combinedOld(iGroup) = angvz_combined(iGroup)

END SUBROUTINE FIM_store_group_old_variables


SUBROUTINE FIM_store_old_bodymarkers(iBody)
USE global_parameters
USE flow_parameters !nbody,nPtsBodyMarker
USE implicit_coupling_parameters !xBodyMarkerOld 
USE boundary_arrays !xbodymarker
USE unstructured_surface_arrays
USE derivative
IMPLICIT NONE

   INTEGER :: iBody,i
   REAL(KIND=CGREAL) :: dummy
   
  IF(ImtheBoss)THEN
    DO i = 1,nPtsBodyMarker(ibody)
     xBodyMarkerOld(ibody,i)  =  xBodyMarker(ibody,i) 
     yBodyMarkerOld(ibody,i)  =  yBodyMarker(ibody,i) 
     zBodyMarkerOld(ibody,i)  =  zBodyMarker(ibody,i)   
     uBodyMarkerOld(ibody,i)  =  uBodyMarker(ibody,i) 
     vBodyMarkerOld(ibody,i)  =  vBodyMarker(ibody,i) 
     wBodyMarkerOld(ibody,i)  =  wBodyMarker(ibody,i)            
    END DO
    
    
    IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN
    DO i = 1,totNumTriElem(ibody)
     triElemCentxOld(ibody,i) =  triElemCentx(ibody,i) 
     triElemCentyOld(ibody,i) =  triElemCenty(ibody,i) 
     triElemCentzOld(ibody,i) =  triElemCentz(ibody,i)  
     triElemAreaOld(ibody,i)  =  triElemArea(ibody,i)
    END DO   
    ENDIF !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0
    
  
  END IF !ImtheBoss
   
#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(uBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(vBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(wBodyMarkerOld, ibody, nbody, 0, nPtsBodyMarker(ibody))
 
 
 IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN
 CALL par_bcast_marker(triElemCentxOld, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemCentyOld, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemCentzOld, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemAreaOld, ibody, nbody, 0, totNumTriElem(ibody)) 
 ENDIF
#  endif
  
END SUBROUTINE FIM_store_old_bodymarkers
 


SUBROUTINE FIM_store_prev_bodymarkers(iBody)
USE global_parameters
USE flow_parameters !nbody,nPtsBodyMarker
USE implicit_coupling_parameters !xBodyMarkerOld 
USE boundary_arrays !xbodymarker
IMPLICIT NONE

   INTEGER :: iBody,i
   REAL(KIND=CGREAL) :: dummy
   
    IF(ImtheBoss)THEN
    DO i = 1,nPtsBodyMarker(ibody)
     xBodyMarkerPrev(ibody,i)  =  xBodyMarker(ibody,i) 
     yBodyMarkerPrev(ibody,i)  =  yBodyMarker(ibody,i) 
     zBodyMarkerPrev(ibody,i)  =  zBodyMarker(ibody,i)
    END DO
    END IF !ImtheBoss
   
#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarkerPrev, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarkerPrev, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarkerPrev, ibody, nbody, 0, nPtsBodyMarker(ibody))
#  endif
  
END SUBROUTINE FIM_store_prev_bodymarkers


!=====================Restore information of previous time step
SUBROUTINE FIM_restore_old_bodymarkers(iBody)
USE global_parameters
USE flow_parameters !nbody,nPtsBodyMarker
USE implicit_coupling_parameters !xBodyMarkerOld 
USE boundary_arrays !xbodymarker
USE unstructured_surface_arrays
USE derivative
IMPLICIT NONE

   INTEGER :: iBody,i
   REAL(KIND=CGREAL) :: dummy
   
  IF(ImtheBoss)THEN
    DO i = 1,nPtsBodyMarker(ibody)
     xBodyMarker(ibody,i)  =  xBodyMarkerOld(ibody,i) 
     yBodyMarker(ibody,i)  =  yBodyMarkerOld(ibody,i) 
     zBodyMarker(ibody,i)  =  zBodyMarkerOld(ibody,i)   
     uBodyMarker(ibody,i)  =  uBodyMarkerOld(ibody,i) 
     vBodyMarker(ibody,i)  =  vBodyMarkerOld(ibody,i) 
     wBodyMarker(ibody,i)  =  wBodyMarkerOld(ibody,i) 
    END DO

    IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN
    DO i = 1,totNumTriElem(ibody)
     triElemCentx(ibody,i) =  triElemCentxOld(ibody,i) 
     triElemCenty(ibody,i) =  triElemCentyOld(ibody,i) 
     triElemCentz(ibody,i) =  triElemCentzOld(ibody,i)  
     triElemArea(ibody,i)  =  triElemAreaOld(ibody,i)
    END DO    
    ENDIF 
  END IF !ImtheBoss
   
#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(uBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(vBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(wBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 
 
 IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN
 CALL par_bcast_marker(triElemCentx, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemCenty, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemCentz, ibody, nbody, 0, totNumTriElem(ibody))
 CALL par_bcast_marker(triElemArea, ibody, nbody, 0, totNumTriElem(ibody))  
 ENDIF
 
#  endif
  
END SUBROUTINE FIM_restore_old_bodymarkers


SUBROUTINE FIM_restore_rotation_matrix(iGroup)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
IMPLICIT NONE
INTEGER :: iGroup

    C11(iGroup) = C11_prev(iGroup)
    C12(iGroup) = C12_prev(iGroup)
    C13(iGroup) = C13_prev(iGroup)
    C21(iGroup) = C21_prev(iGroup)
    C22(iGroup) = C22_prev(iGroup)
    C23(iGroup) = C23_prev(iGroup)
    C31(iGroup) = C31_prev(iGroup)
    C32(iGroup) = C32_prev(iGroup)
    C33(iGroup) = C33_prev(iGroup)

END SUBROUTINE FIM_restore_rotation_matrix

SUBROUTINE FIM_restore_force(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
USE usr_module
IMPLICIT NONE
INTEGER :: iBody

    sCx(iBody) = sCx_old(iBody)
    sCy(iBody) = sCy_old(iBody)
    sCz(iBody) = sCz_old(iBody)
    scmx(iBody) = scmx_old(iBody)
    scmy(iBody) = scmy_old(iBody)
    scmz(iBody) = scmz_old(iBody)


END SUBROUTINE FIM_restore_force

SUBROUTINE FIM_restore_location(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
USE boundary_arrays
USE usr_module
IMPLICIT NONE
INTEGER :: iBody

    xcent(iBody) = xcent_old(iBody)
    ycent(iBody) = ycent_old(iBody)
    zcent(iBody) = zcent_old(iBody)

END SUBROUTINE FIM_restore_location


SUBROUTINE FIM_restore_body_old_variables(iBody)
USE implicit_coupling_parameters 
USE flow_parameters 
IMPLICIT NONE
INTEGER :: iBody

   vxcent(iBody) = vxcentOld(iBody)
   vycent(iBody) = vycentOld(iBody)
   vzcent(iBody) = vzcentOld(iBody)
   angvx(iBody) = angvxOld(iBody)
   angvy(iBody) = angvyOld(iBody)
   angvz(iBody) = angvzOld(iBody)

END SUBROUTINE FIM_restore_body_old_variables


SUBROUTINE FIM_restore_group_old_variables(iGroup)
USE implicit_coupling_parameters 
USE flow_parameters 
IMPLICIT NONE
INTEGER :: iGroup

   vxcent_combined(iGroup) = vxcent_combinedOld(iGroup)
   vycent_combined(iGroup) = vycent_combinedOld(iGroup)
   vzcent_combined(iGroup) = vzcent_combinedOld(iGroup)
   angvx_combined(iGroup) = angvx_combinedOld(iGroup)
   angvy_combined(iGroup) = angvy_combinedOld(iGroup)
   angvz_combined(iGroup) = angvz_combinedOld(iGroup)

END SUBROUTINE FIM_restore_group_old_variables

SUBROUTINE FIM_restore_prev_bodymarkers(iBody)
USE global_parameters
USE flow_parameters !nbody,nPtsBodyMarker
USE implicit_coupling_parameters !xBodyMarkerOld 
USE boundary_arrays !xbodymarker
IMPLICIT NONE
INTEGER :: iBody,i
REAL(KIND=CGREAL) :: dummy

    IF(ImtheBoss)THEN
    DO i = 1,nPtsBodyMarker(ibody)
     xBodyMarker(ibody,i)  =  xBodyMarkerPrev(ibody,i) 
     yBodyMarker(ibody,i)  =  yBodyMarkerPrev(ibody,i) 
     zBodyMarker(ibody,i)  =  zBodyMarkerPrev(ibody,i)    
    END DO
    END IF !ImtheBoss
   
#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))   
#  endif
  
END SUBROUTINE FIM_restore_prev_bodymarkers



!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

 SUBROUTINE FIM_store_vicar_old_variables
   USE global_parameters
   USE boundary_arrays
   USE flow_arrays   
   USE flow_parameters
   USE pressure_arrays
   USE implicit_coupling_parameters  
   USE tahoe_parameters
   USE multiuse_arrays
   USE nlold_arrays
   USE derivative
# ifdef MPI
    USE mpi
# endif
   IMPLICIT NONE
   INTEGER :: i,j,k,iMemb
   REAL(KIND=CGREAL) :: dummy
   
   nFresh_Old = nFresh
      
   PRINT*,'Store old values of u,v,w'
    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
     uOld(i,j,k) = u(i,j,k)
     vOld(i,j,k) = v(i,j,k)
     wOld(i,j,k) = w(i,j,k)
     uGhostOld(i,j,k) = uGhost(i,j,k) 
     vGhostOld(i,j,k) = vGhost(i,j,k) 
     wGhostOld(i,j,k) = wGhost(i,j,k)
     IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN 
     pGhostOld(i,j,k) = pGhost(i,j,k)         
     ENDIF  !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0    
    END DO
    END DO
    END DO
!   Store old values of Iblank

    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
     iblankOld(i,j,k)       = iblank(i,j,k)
	 iblank_solidOld(i,j,k) = iblank_solid(i,j,k)  	 	  
     pPrime_old(i,j,k) =pPrime(i,j,k) 
     
     !IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN     
     p_old(i,j,k) =p(i,j,k)  

       
     uTilde_Old(i,j,k) = uTilde(i,j,k) 
     vTilde_Old(i,j,k) = vTilde(i,j,k)
     wTilde_Old(i,j,k) = wTilde(i,j,k) 
     
      
     fresh_cell_old(i,j,k) = fresh_cell(i,j,k) 
     
     iup_old(i,j,k) = iup(i,j,k)  
     ium_old(i,j,k) = ium(i,j,k) 
     jup_old(i,j,k) = jup(i,j,k) 
     jum_old(i,j,k) = jum(i,j,k) 
     kup_old(i,j,k) = kup(i,j,k) 
     kum_old(i,j,k) = kum(i,j,k) 
       
     ghostCellMark_old(i,j,k) = ghostCellMark(i,j,k) 
     ghostCellMemb_old(i,j,k) = ghostCellMemb(i,j,k) 
     ghostCellSolid_old(i,j,k) = ghostCellSolid(i,j,k)  
    
	 bodyNum_Old(i,j,k) = bodyNum(i,j,k)
	 !ENDIF !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0
	                    
    END DO
    END DO
    END DO  
    
    !IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN   
    DO k = 1, nzc
    DO j = 1, nyc
    DO i = 1, nxc
     face_ue_Old(i,j,k) = face_ue(i,j,k)
     face_vn_Old(i,j,k) = face_vn(i,j,k) 	 	 
     face_wf_Old(i,j,k) = face_wf(i,j,k)
     face_uw_Old(i,j,k) = face_uw(i,j,k)
     face_vs_Old(i,j,k) = face_vs(i,j,k)
     face_wb_Old(i,j,k) = face_wb(i,j,k)    
     
     nlu_Old(i,j,k) = nlu(i,j,k)  
     nlv_Old(i,j,k) = nlv(i,j,k)  
     nlw_Old(i,j,k) = nlw(i,j,k)   
     
     nluold_Old(i,j,k) = nluold(i,j,k)   
     nlvold_Old(i,j,k) = nlvold(i,j,k) 
     nlwold_Old(i,j,k) = nlwold(i,j,k)        
     
     bcxu_Old(i,j,k) = bcxu(i,j,k)  
     bcxv_Old(i,j,k) = bcxv(i,j,k)  
     bcxw_Old(i,j,k) = bcxw(i,j,k)  
     bcyu_Old(i,j,k) = bcyu(i,j,k)  
     bcyv_Old(i,j,k) = bcyv(i,j,k)  
     bcyw_Old(i,j,k) = bcyw(i,j,k)  
     bczu_Old(i,j,k) = bczu(i,j,k) 
     bczv_Old(i,j,k) = bczv(i,j,k) 
     bczw_Old(i,j,k) = bczw(i,j,k)  
     
     iupp_Old(i,j,k) = iupp(i,j,k)  
     iumm_Old(i,j,k) = iumm(i,j,k)             
     jupp_Old(i,j,k) = jupp(i,j,k)  
     jumm_Old(i,j,k) = jumm(i,j,k)
     kupp_Old(i,j,k) = kupp(i,j,k)  
     kumm_Old(i,j,k) = kumm(i,j,k)                
    END DO
    END DO
    END DO    
    
    DO k = 0, nzc+1
    DO j = 0, nyc+1
    DO i = 0, nxc+1
     boundPresSource_old(i,j,k) = boundPresSource(i,j,k)                        
    END DO
    END DO
    END DO 
           
      
    
    
    pgradx1_old = pgradx1
    pgradx2_old = pgradx2
    pgrady1_old = pgrady1
    pgrady2_old = pgrady2
    pgradz1_old = pgradz1
    pgradz2_old = pgradz2
    
    !ENDIF !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0
    
    DO iMemb = 1,nBody_membrane
    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl 	 
	 iblank_membOld(i,j,k,iMemb) = iblank_memb(i,j,k,iMemb)  	
 		           
    END DO
    END DO
    END DO  
    END DO      

END SUBROUTINE FIM_store_vicar_old_variables


  SUBROUTINE FIM_restore_iblank_and_velocities
   USE global_parameters
   USE flow_arrays   
   USE boundary_arrays
   USE flow_parameters
   USE pressure_arrays
   USE implicit_coupling_parameters  
   USE multiuse_arrays   
   USE nlold_arrays
   USE derivative
# ifdef MPI
    use mpi
# endif      
   IMPLICIT NONE
   INTEGER :: i,j,k,iMemb
   REAL(KIND=CGREAL) :: dummy
   
   nFresh = nFresh_Old
!   Update values from old values of u,v,w
    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
     u(i,j,k) = uOld(i,j,k)
     v(i,j,k) = vOld(i,j,k)
     w(i,j,k) = wOld(i,j,k)
     uGhost(i,j,k) = uGhostOld(i,j,k)
     vGhost(i,j,k) = vGhostOld(i,j,k)
     wGhost(i,j,k) = wGhostOld(i,j,k)
     !IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN         
     pGhost(i,j,k) = pGhostOld(i,j,k)                     
     !ENDIF  !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0     
    END DO
    END DO
    END DO
!   Update values from old values of iblank
    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl
     iblank(i,j,k) = iblankOld(i,j,k)
	 iblank_solid(i,j,k) = iblank_solidOld(i,j,k)  	 	 
     pPrime(i,j,k) =pPrime_old(i,j,k)
     
     !IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN          
     p(i,j,k) =p_old(i,j,k)      
     
     uTilde(i,j,k) = uTilde_Old(i,j,k) 
     vTilde(i,j,k) = vTilde_Old(i,j,k)
     wTilde(i,j,k) = wTilde_Old(i,j,k)  
     
     fresh_cell(i,j,k) = fresh_cell_old(i,j,k) 
     iup(i,j,k) = iup_old(i,j,k)  
     ium(i,j,k) = ium_old(i,j,k) 
     jup(i,j,k) = jup_old(i,j,k) 
     jum(i,j,k) = jum_old(i,j,k) 
     kup(i,j,k) = kup_old(i,j,k) 
     kum(i,j,k) = kum_old(i,j,k) 
     ghostCellMark(i,j,k) = ghostCellMark_old(i,j,k) 
     ghostCellMemb(i,j,k) = ghostCellMemb_old(i,j,k) 
     ghostCellSolid(i,j,k) = ghostCellSolid_old(i,j,k)  
     
	 bodyNum(i,j,k) = bodyNum_Old(i,j,k)  
	 
	 !ENDIF !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0     
	                         
    END DO
    END DO
    END DO
    
    
    !IF(derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0)THEN     
    DO k = 1, nzc
    DO j = 1, nyc
    DO i = 1, nxc
     face_ue(i,j,k) = face_ue_Old(i,j,k)
     face_vn(i,j,k) = face_vn_Old(i,j,k) 	 	 
     face_wf(i,j,k) = face_wf_Old(i,j,k)
     face_uw(i,j,k) = face_uw_Old(i,j,k)
     face_vs(i,j,k) = face_vs_Old(i,j,k)
     face_wb(i,j,k) = face_wb_Old(i,j,k)  
     
     nlu(i,j,k) = nlu_Old(i,j,k)  
     nlv(i,j,k) = nlv_Old(i,j,k)  
     nlw(i,j,k) = nlw_Old(i,j,k)
     
     nluold(i,j,k) = nluold_Old(i,j,k)   
     nlvold(i,j,k) = nlvold_Old(i,j,k) 
     nlwold(i,j,k) = nlwold_Old(i,j,k) 
     
     bcxu(i,j,k) = bcxu_Old(i,j,k)  
     bcxv(i,j,k) = bcxv_Old(i,j,k)  
     bcxw(i,j,k) = bcxw_Old(i,j,k)  
     bcyu(i,j,k) = bcyu_Old(i,j,k)  
     bcyv(i,j,k) = bcyv_Old(i,j,k)  
     bcyw(i,j,k) = bcyw_Old(i,j,k)  
     bczu(i,j,k) = bczu_Old(i,j,k) 
     bczv(i,j,k) = bczv_Old(i,j,k) 
     bczw(i,j,k) = bczw_Old(i,j,k) 
     
     iupp(i,j,k) = iupp_Old(i,j,k)  
     iumm(i,j,k) = iumm_Old(i,j,k)             
     jupp(i,j,k) = jupp_Old(i,j,k)  
     jumm(i,j,k) = jumm_Old(i,j,k)
     kupp(i,j,k) = kupp_Old(i,j,k)  
     kumm(i,j,k) = kumm_Old(i,j,k)                                 
    END DO
    END DO
    END DO    
    
    DO k = 0, nzc+1
    DO j = 0, nyc+1
    DO i = 0, nxc+1
     boundPresSource(i,j,k) = boundPresSource_old(i,j,k)                        
    END DO
    END DO
    END DO    
    
 
      
    pgradx1 = pgradx1_old
    pgradx2 = pgradx2_old
    pgrady1 = pgrady1_old
    pgrady2 = pgrady2_old
    pgradz1 = pgradz1_old
    pgradz2 = pgradz2_old   
    
    !ENDIF !derivative_flag == 1 .AND. implicit_coupling_flag_combined_motion == 0     
    
    DO iMemb = 1,nBody_membrane
    DO k = 0, nzc+1
    DO j = 1-Ngl, nyc+Ngl
    DO i = 1-Ngl, nxc+Ngl 	 
	 iblank_memb(i,j,k,iMemb) = iblank_membOld(i,j,k,iMemb)  	          
    END DO
    END DO
    END DO  
    END DO        
      
  END SUBROUTINE FIM_restore_iblank_and_velocities
  
  
  SUBROUTINE FIM_calculate_new_iblank 
   USE global_parameters
   REAL(KIND=CGREAL) :: myCommTime
   CALL calculate_arclength_norm_ds()
   CALL set_boundary(myCommTime)  ! Calculates new IBLANK
  END SUBROUTINE FIM_calculate_new_iblank 
  
  
  
  SUBROUTINE FIM_WRITE_CONVERGENCE_DATA
 
    USE global_parameters
    USE flow_parameters
	USE boundary_arrays
    USE tahoe_parameters
    USE implicit_coupling_parameters
    USE usr_module

   IMPLICIT NONE
     
    INTEGER           :: i,ibody,k
    REAL(KIND=CGREAL) :: dummy,xmar,ymar,zmar,pmar

    IF(ntime == nstart_FIM + 1 .AND. kimplicit_FIM == 1)THEN
     OPEN (imconverg, file = 'implicit_convergence_FIM.dat',STATUS='REPLACE')
    ELSE
     OPEN (imconverg, file = 'implicit_convergence_FIM.dat',POSITION = 'APPEND')  
    ENDIF 
    WRITE (imconverg, 10) kimplicit_FIM,ImplicitResidual_FIM,ntime  !,num_fresh,num_dead  
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    WRITE (imconverg, *) '( density_fluid-density_solid(1) )*volume(bodyNumber)*Fr'
    WRITE (imconverg, *) ( density_fluid-density_solid(1) )*volume(1)*Fr    
    WRITE (imconverg, *) 'scx(1),scy(1),scz(1)'            
    WRITE (imconverg, *) scx(1),scy(1),scz(1)
    WRITE (imconverg, *) 'force_x_combined(1),force_y_combined(1),force_z_combined(1)'
    WRITE (imconverg, *) 0.5*force_x_combined(1),0.5*force_y_combined(1),0.5*force_z_combined(1)
    WRITE (imconverg, *) 'moment_x_combined(1),moment_y_combined(1),moment_z_combined(1)'
    WRITE (imconverg, *) 0.5*moment_x_combined(1),0.5*moment_y_combined(1),0.5*moment_z_combined(1)       
    WRITE (imconverg, *) 'xcent_combined(1),ycent_combined(1),zcent_combined(1)'    
    WRITE (imconverg, *) xcent_combined(1),ycent_combined(1),zcent_combined(1)
    WRITE (imconverg, *) 'vxcent_combined(1),vycent_combined(1),vzcent_combined(1)'        
    WRITE (imconverg, *) vxcent_combined(1),vycent_combined(1),vzcent_combined(1)  !,num_fresh,num_dead  
    WRITE (imconverg, *) 'vxcent_combinedOld(1),vycent_combinedOld(1),vzcent_combinedOld(1)'      
    WRITE (imconverg, *) vxcent_combinedOld(1),vycent_combinedOld(1),vzcent_combinedOld(1)
    WRITE (imconverg, *) 'angvx_combined(1),angvy_combined(1),angvz_combined(1)'    
    WRITE (imconverg, *) angvx_combined(1),angvy_combined(1),angvz_combined(1)
    WRITE (imconverg, *) 'angvx_combinedOld(1),angvy_combinedOld(1),angvz_combinedOld(1)'      
    WRITE (imconverg, *) angvx_combinedOld(1),angvy_combinedOld(1),angvz_combinedOld(1)  
    WRITE (imconverg, *) 'xBodyMarker(1,1),yBodyMarker(1,1),zBodyMarker(1,1)'
    WRITE (imconverg, *) xBodyMarker(1,1),yBodyMarker(1,1),zBodyMarker(1,1) 
    WRITE (imconverg, *) 'xBodyMarkerOld(1,1),yBodyMarkerOld(1,1),zBodyMarkerOld(1,1)'
    WRITE (imconverg, *) xBodyMarkerOld(1,1),yBodyMarkerOld(1,1),zBodyMarkerOld(1,1)     
    WRITE (imconverg, *) '---------------------------------'    
    WRITE (imconverg, *) C11(1),C12(1),C13(1)
    WRITE (imconverg, *) C21(1),C22(1),C23(1)
    WRITE (imconverg, *) C31(1),C32(1),C33(1)   
    WRITE (imconverg, *) '================================='    
    WRITE (imconverg, *) '================================='      
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddd      
    CLOSE (imconverg)

!10  FORMAT (i6,1x,1PE12.5,1x,i6,1x,i6,1x,i6)  
10  FORMAT (i6,1x,1PE12.5,1x,i6)  
 
ENDSUBROUTINE FIM_WRITE_CONVERGENCE_DATA


SUBROUTINE FIM_convergence_implicit
   USE implicit_coupling_parameters 
   USE tahoe_parameters
   USE flow_parameters 
   IMPLICIT NONE
   INTEGER :: iBody   
   REAL(KIND=CGREAL) :: myCommTime   
   
!---Update all the body's marker points (Include FIM body)   
  DO iBody=1,nBody 
	  CALL FIM_underrelax_xb_ub(iBody) !called by IMtheBoss and broadcast
  END DO !iBody

   CALL FIM_check_convergence_implicit_scheme !called by IMtheBoss and broadcast
 	
!---Update varibles about body mesh    
   CALL calculate_arclength_norm_ds()
   CALL set_boundary(myCommTime)   


END SUBROUTINE FIM_convergence_implicit


  SUBROUTINE FIM_check_convergence_implicit_scheme
   
   USE implicit_coupling_parameters
   USE flow_parameters !nbody,nPtsBodyMarker
   USE boundary_arrays !xbodymarker
   USE tahoe_parameters
   USE mpi
   IMPLICIT NONE
   INTEGER :: ibody,i,ierr
   REAL(KIND=CGREAL) :: resx,resxmax,resy,resymax,resz,reszmax
   
   IF(ImtheBOSS) THEN
   resx = 0.0
   resxmax = 0.0
   DO ibody = 1,nbody
   DO i = 1, nPtsBodyMarker(ibody)
      resx = sqrt((xBodyMarker(ibody,i) - xBodyMarkerOld(ibody,i))**2)
!	  IF(verbosetahoe==1) PRINT*,'x',xBodyMarker(ibody,i),xBodyMarkerOld(ibody,i)
      resxmax = max(resx,resxmax )
!	  IF(verbosetahoe==1) PRINT*,'x',resx,resxmax
   END DO
   END DO
   
   resy = 0.0
   resymax = 0.0
   DO ibody = 1,nbody
   DO i = 1, nPtsBodyMarker(ibody)
      resy = sqrt((yBodyMarker(ibody,i) - yBodyMarkerOld(ibody,i))**2)
!	  IF(verbosetahoe==1) PRINT*,'y',yBodyMarker(ibody,i),yBodyMarkerOld(ibody,i)
      resymax = max(resy,resymax )
!	  IF(verbosetahoe==1) PRINT*,'y',resy,resymax
   END DO
   END DO
   
   resz = 0.0
   reszmax = 0.0
   DO ibody = 1,nbody
      DO i = 1, nPtsBodyMarker(ibody)
      resz = sqrt((zBodyMarker(ibody,i) - zBodyMarkerOld(ibody,i))**2)
!	  IF(verbosetahoe==1) PRINT*,'z',zBodyMarker(ibody,i),zBodyMarkerOld(ibody,i)
      reszmax = max(resz,reszmax )
!	  IF(verbosetahoe==1) PRINT*,'z',resz,reszmax
   END DO
   END DO
   
   PRINT*,'MaxRes,xb,yb,zb',resxmax,resymax,reszmax
      
   ImplicitResidual_FIM = max(resxmax,resymax,reszmax)
   
   PRINT*,'MaxResTot',ImplicitResidual_FIM
   
   END IF

#  ifdef MPI
  CALL MPI_BCAST(ImplicitResidual_FIM, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
#  endif

  PRINT*,'Broadcasting done for ImplicitResidual_FIM'

  END SUBROUTINE FIM_check_convergence_implicit_scheme
  
  
  SUBROUTINE FIM_calculate_alphaimplicit
   USE global_parameters
   USE flow_parameters !nbody,nPtsBodyMarker
   USE implicit_coupling_parameters !xBodyMarkerOld 
   USE boundary_arrays !xbodymarker
   USE tahoe_parameters
   IMPLICIT NONE
   INTEGER :: iBody,i
   REAL(KIND=CGREAL) :: dummy,alp, alp1,b,termtime1,termtime2
   
   alp = alpha_underrelax_FIM
   
   IF(varyalphawithtime_FIM == 1)THEN
   
	   IF(ntime<timeimpl1_FIM) alp = alpha_underrelax_FIM
	   
	   IF(ntime>=timeimpl1_FIM .AND. ntime<=timeimpl2_FIM)THEN
	   termtime1 = timeimpl2_FIM - timeimpl1_FIM
	   termtime2 =  ntime - timeimpl1_FIM 
	   alp = alp + (alpha_underrelax2_FIM - alp)*termtime2/termtime1 
	   END IF
	   
	   IF(ntime>timeimpl2_FIM) alp = alpha_underrelax2
   
   END IF
   
   IF(flag_underrelax_FIM == 1 .AND.kimplicit_FIM > kimplicit_start_FIM)THEN
   alp = alpha_underrelax_FIM + slope_underrelax_FIM*(kimplicit_FIM-kimplicit_start_FIM)
   END IF
  
   IF(alp>1.0) alp = 1.0  
   
   alphaimplicit_FIM = alp
   
END SUBROUTINE FIM_calculate_alphaimplicit


SUBROUTINE FIM_underrelax_xb_ub(iBody)
   USE global_parameters
   USE flow_parameters !nbody,nPtsBodyMarker
   USE implicit_coupling_parameters !xBodyMarkerOld 
   USE boundary_arrays !xbodymarker
   USE tahoe_parameters
   IMPLICIT NONE
   INTEGER :: iBody,i
   REAL(KIND=CGREAL) :: dummy,alp, alp1,b,termtime1,termtime2
     
   alp = alphaimplicit_FIM
   alp1 = 1.0 - alp 
   
IF(ImtheBoss)THEN  
DO i = 1,nPtsBodyMarker(ibody)
xBodyMarker(ibody,i) = alp*xBodyMarker(ibody,i)+alp1*xBodyMarkerOld(ibody,i)
yBodyMarker(ibody,i) = alp*yBodyMarker(ibody,i)+alp1*yBodyMarkerOld(ibody,i)
zBodyMarker(ibody,i) = alp*zBodyMarker(ibody,i)+alp1*zBodyMarkerOld(ibody,i)
uBodyMarker(ibody,i) = alp*uBodyMarker(ibody,i)+alp1*uBodyMarkerOld(ibody,i)
vBodyMarker(ibody,i) = alp*vBodyMarker(ibody,i)+alp1*vBodyMarkerOld(ibody,i)
wBodyMarker(ibody,i) = alp*wBodyMarker(ibody,i)+alp1*wBodyMarkerOld(ibody,i)

END DO
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    WRITE (66666, *) '2222222222xBodyMarker(1,1),yBodyMarker(1,1),zBodyMarker(1,1)'
    WRITE (66666, *) xBodyMarker(1,1),yBodyMarker(1,1),zBodyMarker(1,1) 
close(66666)    
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd  

END IF !ImtheBoss
 
#   ifdef MPI 
 CALL par_bcast_marker(xBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(yBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(zBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(uBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(vBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
 CALL par_bcast_marker(wBodyMarker, ibody, nbody, 0, nPtsBodyMarker(ibody))
#  endif  

 

END SUBROUTINE FIM_underrelax_xb_ub


