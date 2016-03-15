! -------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  The George Washington University
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
!  Filename: BOUNDARY_MEMORY_ALLOCATE.PAR.F90
!  Latest Modification: Dec, 29 2010 (PAT 2.1.0)
!  by Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     BODY_allocate_memory()
!     BOUNDARY_allocate_memory()
!     MARKER_allocate_memory()
!     UNSTRUC_allocate_memory()
! --------------------------------------------------------------------
SUBROUTINE BODY_allocate_memory()

! -----------------------------------------------------
!  This subroutine allocates arrays for immersed body.
! -----------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE unstructured_surface_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    USE boundary_arrays
    USE usr_module 
    USE implicit_coupling_parameters
    USE derivative
    USE scalar
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
    
    IMPLICIT NONE
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER :: iBody,iGroup
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)     


    ALLOCATE(canonical_body_type(nBody))
    ALLOCATE(body_dim(nBody))
    ALLOCATE(boundary_motion_type(nBody))

!   new array for mixed bodies
!   --------------------------
    ALLOCATE(membrane_type(nBody))

    ALLOCATE(wall_type(nBody))
    ALLOCATE(nPtsBodyMarkerOrig(nBody))
    ALLOCATE(nPtsBodyMarker(nBody))
    ALLOCATE(totNumTriElem(nBody))
    ALLOCATE(unstruc_surface_type(nBody))

    ALLOCATE(n_phi(nBody))
    ALLOCATE(n_theta(nbody))
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ALLOCATE(combined_type(2))
    ALLOCATE(combined_Group_index(nBody))
    ALLOCATE(Surface_Integral(nBody)) 
    ALLOCATE(unsteady1(nBody))  
    ALLOCATE(unsteady1_prev(nBody)) 
    ALLOCATE(unsteady_term1(nBody)) 
    ALLOCATE(unsteady_term2(nBody))  
    ALLOCATE(unsteady_sum1(nBody))  
    ALLOCATE(unsteady_sum2(nBody))  
    ALLOCATE(kinematics_term2(nBody))
    ALLOCATE(unsteady_sum2p(nBody)) 
    ALLOCATE(unsteady_sum2v(nBody))         
    ALLOCATE(advection_sum(nBody))   
    ALLOCATE(Lamb_sum(nBody))  
    ALLOCATE(tau_Xgrad(nBody))  
    ALLOCATE(Grad_u_square_sum(nBody))   
    ALLOCATE(Grad_potential_u_square_sum(nBody)) 
    ALLOCATE(Grad_vorticity_u_square_sum(nBody))     
!dddddddddddddddddddddddddddddddddddddddddddddddddddd
    ALLOCATE(unsteady_sum_check(nBody)) 
    ALLOCATE(angle(nBody)) 
!dddddddddddddddddddddddddddddddddddddddddddddddddddd     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    ALLOCATE(radiusx(nBody))
    ALLOCATE(radiusy(nBody))
    ALLOCATE(radiusz(nBody))
    ALLOCATE(alpha(nBody))
    ALLOCATE(cosalpha(nBody))
    ALLOCATE(sinalpha(nBody))
    ALLOCATE(xcent(nBody))
    ALLOCATE(ycent(nBody))
    ALLOCATE(zcent(nBody))
    ALLOCATE(vxcent(nBody))
    ALLOCATE(vycent(nBody))
    ALLOCATE(vzcent(nBody))
    ALLOCATE(angvx(nBody))
    ALLOCATE(angvy(nBody))
    ALLOCATE(angvz(nBody))
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    ALLOCATE(x_rot_cent(nBody))
    ALLOCATE(y_rot_cent(nBody))
    ALLOCATE(z_rot_cent(nBody))
    ALLOCATE(force_x(nBody))
    ALLOCATE(force_y(nBody))
    ALLOCATE(force_z(nBody))
    ALLOCATE(moment_x(nBody))
    ALLOCATE(moment_y(nBody))
    ALLOCATE(moment_z(nBody))  
    
    ALLOCATE(vxcentOld(nBody))
    ALLOCATE(vycentOld(nBody))
    ALLOCATE(vzcentOld(nBody))
    ALLOCATE(angvxOld(nBody))
    ALLOCATE(angvyOld(nBody))
    ALLOCATE(angvzOld(nBody))
    
    ALLOCATE(vxcent_combined(nGroup_Combined))
    ALLOCATE(vycent_combined(nGroup_Combined))
    ALLOCATE(vzcent_combined(nGroup_Combined))
    ALLOCATE(angvx_combined(nGroup_Combined))
    ALLOCATE(angvy_combined(nGroup_Combined))
    ALLOCATE(angvz_combined(nGroup_Combined))          
    ALLOCATE(angvx_combined_old(nGroup_Combined))
    ALLOCATE(angvy_combined_old(nGroup_Combined))
    ALLOCATE(angvz_combined_old(nGroup_Combined)) 
    
    ALLOCATE(vxcent_combinedOld(nGroup_Combined))
    ALLOCATE(vycent_combinedOld(nGroup_Combined))
    ALLOCATE(vzcent_combinedOld(nGroup_Combined))
    ALLOCATE(angvx_combinedOld(nGroup_Combined))
    ALLOCATE(angvy_combinedOld(nGroup_Combined))
    ALLOCATE(angvz_combinedOld(nGroup_Combined)) 
    
    ALLOCATE(R11(nGroup_Combined))                           
    ALLOCATE(R12(nGroup_Combined))
    ALLOCATE(R13(nGroup_Combined))
    ALLOCATE(R21(nGroup_Combined))
    ALLOCATE(R22(nGroup_Combined))
    ALLOCATE(R23(nGroup_Combined))
    ALLOCATE(R31(nGroup_Combined))
    ALLOCATE(R32(nGroup_Combined))
    ALLOCATE(R33(nGroup_Combined))     
    
    ALLOCATE(C11(nGroup_Combined))                           
    ALLOCATE(C12(nGroup_Combined))
    ALLOCATE(C13(nGroup_Combined))
    ALLOCATE(C21(nGroup_Combined))
    ALLOCATE(C22(nGroup_Combined))
    ALLOCATE(C23(nGroup_Combined))
    ALLOCATE(C31(nGroup_Combined))
    ALLOCATE(C32(nGroup_Combined))
    ALLOCATE(C33(nGroup_Combined)) 
    
    ALLOCATE(C_11(nBody))                           
    ALLOCATE(C_12(nBody))
    ALLOCATE(C_13(nBody))
    ALLOCATE(C_21(nBody))
    ALLOCATE(C_22(nBody))
    ALLOCATE(C_23(nBody))
    ALLOCATE(C_31(nBody))
    ALLOCATE(C_32(nBody))
    ALLOCATE(C_33(nBody)) 
       

    ALLOCATE(non_vxcent(nGroup_Combined)) 
    ALLOCATE(non_vycent(nGroup_Combined)) 
    ALLOCATE(non_vzcent(nGroup_Combined)) 
        
    ALLOCATE(angv_roll(nGroup_Combined)) 
    ALLOCATE(angv_yaw(nGroup_Combined)) 
    ALLOCATE(angv_pitch(nGroup_Combined))    
    
    ALLOCATE(C11_prev(nGroup_Combined))                           
    ALLOCATE(C12_prev(nGroup_Combined))
    ALLOCATE(C13_prev(nGroup_Combined))
    ALLOCATE(C21_prev(nGroup_Combined))
    ALLOCATE(C22_prev(nGroup_Combined))
    ALLOCATE(C23_prev(nGroup_Combined))
    ALLOCATE(C31_prev(nGroup_Combined))
    ALLOCATE(C32_prev(nGroup_Combined))
    ALLOCATE(C33_prev(nGroup_Combined)) 
        
    ALLOCATE(theta_group(nGroup_Combined)) 
    
    ALLOCATE(vxcent_combined_init(nGroup_Combined))  
    ALLOCATE(vycent_combined_init(nGroup_Combined))  
    ALLOCATE(vzcent_combined_init(nGroup_Combined))  
    ALLOCATE(angvx_combined_init(nGroup_Combined))  
    ALLOCATE(angvy_combined_init(nGroup_Combined))  
    ALLOCATE(angvz_combined_init(nGroup_Combined))
    
    ALLOCATE(inertia_force_x_ref(nGroup_Combined)) 
    ALLOCATE(inertia_force_y_ref(nGroup_Combined)) 
    ALLOCATE(inertia_moment_z_ref(nGroup_Combined)) 
    ALLOCATE(inertia_force_x_deltaU(nGroup_Combined)) 
    ALLOCATE(inertia_force_y_deltaU(nGroup_Combined)) 
    ALLOCATE(inertia_moment_z_deltaU(nGroup_Combined)) 
    ALLOCATE(inertia_force_x_deltaV(nGroup_Combined)) 
    ALLOCATE(inertia_force_y_deltaV(nGroup_Combined)) 
    ALLOCATE(inertia_moment_z_deltaV(nGroup_Combined)) 
    ALLOCATE(inertia_force_x_deltaQ(nGroup_Combined)) 
    ALLOCATE(inertia_force_y_deltaQ(nGroup_Combined)) 
    ALLOCATE(inertia_moment_z_deltaQ(nGroup_Combined))     
    
    
    ALLOCATE(non_inertia_force_x_ref(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_y_ref(nGroup_Combined)) 
    ALLOCATE(non_inertia_moment_z_ref(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_x_deltaU(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_y_deltaU(nGroup_Combined)) 
    ALLOCATE(non_inertia_moment_z_deltaU(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_x_deltaV(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_y_deltaV(nGroup_Combined)) 
    ALLOCATE(non_inertia_moment_z_deltaV(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_x_deltaQ(nGroup_Combined)) 
    ALLOCATE(non_inertia_force_y_deltaQ(nGroup_Combined)) 
    ALLOCATE(non_inertia_moment_z_deltaQ(nGroup_Combined)) 
    ALLOCATE(invC11(nGroup_Combined)) 
    ALLOCATE(invC12(nGroup_Combined)) 
    ALLOCATE(invC13(nGroup_Combined)) 
    ALLOCATE(invC21(nGroup_Combined)) 
    ALLOCATE(invC22(nGroup_Combined)) 
    ALLOCATE(invC23(nGroup_Combined)) 
    ALLOCATE(invC31(nGroup_Combined)) 
    ALLOCATE(invC32(nGroup_Combined)) 
    ALLOCATE(invC33(nGroup_Combined)) 
    ALLOCATE(non_inertia_U(nGroup_Combined)) 
    ALLOCATE(non_inertia_V(nGroup_Combined)) 
    ALLOCATE(non_inertia_W(nGroup_Combined)) 
    ALLOCATE(non_inertia_angvx(nGroup_Combined)) 
    ALLOCATE(non_inertia_angvy(nGroup_Combined)) 
    ALLOCATE(non_inertia_angvz(nGroup_Combined)) 
    
    ALLOCATE(xcent_combined(nGroup_Combined))
    ALLOCATE(ycent_combined(nGroup_Combined))
    ALLOCATE(zcent_combined(nGroup_Combined))
    ALLOCATE(I_XX_COMBINED(nGroup_Combined))
    ALLOCATE(I_YY_COMBINED(nGroup_Combined))
    ALLOCATE(I_ZZ_COMBINED(nGroup_Combined))
    ALLOCATE(I_XY_COMBINED(nGroup_Combined))
    ALLOCATE(I_YZ_COMBINED(nGroup_Combined))
    ALLOCATE(I_XZ_COMBINED(nGroup_Combined))
    ALLOCATE(force_x_combined(nGroup_Combined))
    ALLOCATE(force_y_combined(nGroup_Combined))
    ALLOCATE(force_z_combined(nGroup_Combined))
    ALLOCATE(moment_x_combined(nGroup_Combined))
    ALLOCATE(moment_y_combined(nGroup_Combined))
    ALLOCATE(moment_z_combined(nGroup_Combined))
    ALLOCATE(non_dim_mass_combined(nGroup_Combined))
    ALLOCATE(invMI(3,3,nBody))   
    ALLOCATE(invMI_combined(3,3,nGroup_Combined))
    ALLOCATE(nonDimM_I_combined(3,3,nGroup_Combined))  
           
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
!   new arrays for rotation
!   -----------------------
    ALLOCATE(angvx_old(nBody))
    ALLOCATE(angvy_old(nBody))
    ALLOCATE(angvz_old(nBody))
    ALLOCATE(xcentinit(nBody))
    ALLOCATE(ycentinit(nBody))
    ALLOCATE(zcentinit(nBody))
    ALLOCATE(vxcentTrans(nBody))
    ALLOCATE(vycentTrans(nBody))
    ALLOCATE(vzcentTrans(nBody))
    ALLOCATE(ampx(nBody))
    ALLOCATE(ampy(nBody))
    ALLOCATE(ampz(nBody))
    ALLOCATE(freqx(nBody))
    ALLOCATE(freqy(nBody))
    ALLOCATE(freqz(nBody))

    ALLOCATE(phase(nBody))
    ALLOCATE(cosphase(nBody))
    ALLOCATE(sinphase(nBody))
    ALLOCATE(angvxinit(nBody))
    ALLOCATE(angvyinit(nBody))
    ALLOCATE(angvzinit(nBody))
    ALLOCATE(ampangx(nBody))
    ALLOCATE(ampangy(nBody))
    ALLOCATE(ampangz(nBody))
    ALLOCATE(freqangx(nBody))
    ALLOCATE(freqangy(nBody))
    ALLOCATE(freqangz(nBody))
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ALLOCATE(i_fixed(nBody)) 
    ALLOCATE(fixed_mother_body(nBody)) 
    ALLOCATE(pendulum_mother_body(nBody)) 
    ALLOCATE(hinge_x(nBody), hinge_y(nBody), hinge_z(nBody))  
    ALLOCATE(hinge_vx(nBody), hinge_vy(nBody), hinge_vz(nBody))    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    ALLOCATE(mMinWallVel(nBody))
    ALLOCATE(mMaxWallVel(nBody))
    ALLOCATE(ampVelx(nBody))
    ALLOCATE(ampVely(nBody))
    ALLOCATE(ampVelz(nBody))
    ALLOCATE(freqVelx(nBody))
    ALLOCATE(freqVely(nBody))
    ALLOCATE(freqVelz(nBody))
    ALLOCATE(phaseVelx(nBody))
    ALLOCATE(phaseVely(nBody))
    ALLOCATE(phaseVelz(nBody))

!   VEERA: Flow Induced Motion
!   --------------------------
    ALLOCATE(xcentConstr(nBody))
    ALLOCATE(ycentConstr(nBody))
    ALLOCATE(zcentConstr(nBody))
    ALLOCATE(density_solid(nBody))
    ALLOCATE(ntimePerCycle(nbody))
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    ALLOCATE(non_dim_mass(nBody))
    ALLOCATE(nonDimM_I(3,3,nBody))
    ALLOCATE(non_dim_volume(nBody))
    ALLOCATE(volume(nBody))
    ALLOCATE(I_XX(nBody))
    ALLOCATE(I_YY(nBody))  
    ALLOCATE(I_ZZ(nBody))
    ALLOCATE(I_XY(nBody)) 
    ALLOCATE(I_YZ(nBody)) 
    ALLOCATE(I_XZ(nBody)) 	
    
    ALLOCATE(scx_old(nBody))
    ALLOCATE(scy_old(nBody))  
    ALLOCATE(scz_old(nBody))
    ALLOCATE(scmx_old(nBody)) 
    ALLOCATE(scmy_old(nBody)) 
    ALLOCATE(scmz_old(nBody))     
    
    ALLOCATE(xcent_old(nBody))
    ALLOCATE(ycent_old(nBody))  
    ALLOCATE(zcent_old(nBody))    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    
   
    ALLOCATE(scx(nBody))
    ALLOCATE(scy(nBody))  
    ALLOCATE(scz(nBody))
    ALLOCATE(scmx(nBody)) 
    ALLOCATE(scmy(nBody)) 
    ALLOCATE(scmz(nBody)) 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)     
    ALLOCATE(scpw(nBody)) 
    ALLOCATE(ssspw(nBody))  
    ALLOCATE(adpw(nBody))   
    ALLOCATE(shear_x(nBody)) 
    ALLOCATE(shear_y(nBody)) 
    ALLOCATE(shear_z(nBody)) 
    ALLOCATE(pres_x(nBody)) 
    ALLOCATE(pres_y(nBody)) 
    ALLOCATE(pres_z(nBody))     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
  
!   Initialize arrays
!   -----------------
    canonical_body_type   = 0
    body_dim              = 0
    boundary_motion_type  = 0
    membrane_type         = 0

    wall_type             = 0
    nPtsBodyMarkerOrig    = 0
    nPtsBodyMarker        = 0
    totNumTriElem         = 0
    unstruc_surface_type  = 0

    n_phi                 = 0
    n_theta               = 0

    radiusx      = 0.0_CGREAL
    radiusy      = 0.0_CGREAL
    radiusz      = 0.0_CGREAL

    xcent        = 0.0_CGREAL 
    ycent        = 0.0_CGREAL 
    zcent        = 0.0_CGREAL 

    vxcent       = 0.0_CGREAL 
    vycent       = 0.0_CGREAL 
    vzcent       = 0.0_CGREAL 

    angvx        = 0.0_CGREAL 
    angvy        = 0.0_CGREAL 
    angvz        = 0.0_CGREAL 

    angvx_old    = 0.0_CGREAL
    angvy_old    = 0.0_CGREAL
    angvz_old    = 0.0_CGREAL

    xcentinit    = 0.0_CGREAL 
    ycentinit    = 0.0_CGREAL 
    zcentinit    = 0.0_CGREAL 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    angvx_combined_old  = 0.0_CGREAL 
    angvy_combined_old  = 0.0_CGREAL 
    angvz_combined_old  = 0.0_CGREAL   
    
    C11          = 1.0_CGREAL  
    C12          = 0.0_CGREAL  
    C13          = 0.0_CGREAL  
    C21          = 0.0_CGREAL  
    C22          = 1.0_CGREAL  
    C23          = 0.0_CGREAL  
    C31          = 0.0_CGREAL  
    C32          = 0.0_CGREAL  
    C33          = 1.0_CGREAL   
    
    C_11          = 1.0_CGREAL  
    C_12          = 0.0_CGREAL  
    C_13          = 0.0_CGREAL  
    C_21          = 0.0_CGREAL  
    C_22          = 1.0_CGREAL  
    C_23          = 0.0_CGREAL  
    C_31          = 0.0_CGREAL  
    C_32          = 0.0_CGREAL  
    C_33          = 1.0_CGREAL     
    
    R11          = 1.0_CGREAL  
    R12          = 0.0_CGREAL  
    R13          = 0.0_CGREAL  
    R21          = 0.0_CGREAL  
    R22          = 1.0_CGREAL  
    R23          = 0.0_CGREAL  
    R31          = 0.0_CGREAL  
    R32          = 0.0_CGREAL  
    R33          = 1.0_CGREAL  
                        
    force_x      = 0.0_CGREAL 
    force_y      = 0.0_CGREAL
    force_z      = 0.0_CGREAL
    
    moment_x     = 0.0_CGREAL    
    moment_y     = 0.0_CGREAL
    moment_z     = 0.0_CGREAL                

    hinge_vx     = 0.0_CGREAL   
    hinge_vy     = 0.0_CGREAL   
    hinge_vz     = 0.0_CGREAL   
      
    angle        = 0.0_CGREAL  

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)       
    vxcentTrans  = 0.0_CGREAL 
    vycentTrans  = 0.0_CGREAL 
    vzcentTrans  = 0.0_CGREAL 

    ampx         = 0.0_CGREAL 
    ampy         = 0.0_CGREAL
    ampz         = 0.0_CGREAL

    freqx        = 0.0_CGREAL 
    freqy        = 0.0_CGREAL
    freqz        = 0.0_CGREAL

    phase        = 0.0_CGREAL
    cosphase     = 0.0_CGREAL
    sinphase     = 0.0_CGREAL

    angvxinit    = 0.0_CGREAL 
    angvyinit    = 0.0_CGREAL 
    angvzinit    = 0.0_CGREAL 

    ampangx      = 0.0_CGREAL 
    ampangy      = 0.0_CGREAL
    ampangz      = 0.0_CGREAL

    freqangx    = 0.0_CGREAL 
    freqangy    = 0.0_CGREAL
    freqangz    = 0.0_CGREAL  
 
    mMinWallVel = 0
    mMaxWallVel = 0
    ampvelx     = 0.0_CGREAL 
    ampvely     = 0.0_CGREAL
    ampvelz     = 0.0_CGREAL        

    freqvelx    = 0.0_CGREAL 
    freqvely    = 0.0_CGREAL
    freqvelz    = 0.0_CGREAL        

    phasevelx   = 0.0_CGREAL 
    phasevely   = 0.0_CGREAL
    phasevelz   = 0.0_CGREAL        

    xcentConstr = 0.0_CGREAL
    ycentConstr = 0.0_CGREAL
    zcentConstr = 0.0_CGREAL       

    ntimePerCycle = 0  ! added for no-stop feature

END SUBROUTINE BODY_allocate_memory
!---------------------------------------------------------------------



SUBROUTINE BOUNDARY_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE unstructured_surface_arrays
	USE implicit_coupling_parameters  ! Added by Rajneesh
    
    IMPLICIT NONE

    INTEGER ::iBody, iErr

!=======================================================================
! arrays required by all internal boundaries
!=======================================================================

    ALLOCATE(iblank(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iblank'
      STOP
    ENDIF ! ierr 
	
	! Added by Rajneesh
	
	ALLOCATE(iblankOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iblank'
      STOP
    ENDIF ! ierr 
	
	ALLOCATE(iblank_solidOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iblank'
      STOP
    ENDIF ! ierr 
    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ALLOCATE(fresh_cell_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for fresh_cell_old'
      STOP
    ENDIF ! ierr   
    

    ALLOCATE(iup_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iup_old'
      STOP
    ENDIF ! ierr


    ALLOCATE(ium_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ium_old'
      STOP
    ENDIF ! ierr


    ALLOCATE(jup_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jup_old'
      STOP
    ENDIF ! ierr


    ALLOCATE(jum_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jum_old'
      STOP
    ENDIF ! ierr


    ALLOCATE(kup_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kup_old'
      STOP
    ENDIF ! ierr


    ALLOCATE(kum_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kum_old'
      STOP
    ENDIF ! ierr  
   
    ALLOCATE(ghostCellMark_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark_old'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(ghostCellSolid_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark_solid_old'
      STOP
    ENDIF ! ierr   
    
    ALLOCATE(bodyNum_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for bodyNum_Old'
      STOP
    ENDIF ! ierr   
    
    ALLOCATE(iupp_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iupp_Old'
      STOP
    ENDIF ! ierr  
    
    ALLOCATE(iumm_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iumm_Old'
      STOP
    ENDIF ! ierr  
 
    ALLOCATE(jupp_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jupp_Old'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(jumm_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jumm_Old'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(kupp_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kupp_Old'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(kumm_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kumm_Old'
      STOP
    ENDIF ! ierr    
    
      ALLOCATE(iblank_membOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,nBody_membrane),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
         'BOUNDARY_allocate_memory: Memory Allocation Error for iblank_membOld'
        STOP
      ENDIF ! ierr
      
   
    
    ALLOCATE(ghostCellMemb_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark_memb_old'
      STOP
    ENDIF ! ierr       
     
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)   
	
    
    IF (nBody_membrane > 0 ) THEN
      ALLOCATE(iblank_memb(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,nBody_membrane),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
         'BOUNDARY_allocate_memory: Memory Allocation Error for iblank_memb'
        STOP
      ENDIF ! ierr
             
    ENDIF
    
 
    
    ALLOCATE(iblank_solid(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iblank_solid'
      STOP
    ENDIF ! ierr
 
    ALLOCATE(ghostCellMark(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(ghostCellMemb(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark_memb'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(ghostCellSolid(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ghostCellMark_solid'
      STOP
    ENDIF ! ierr
 
    ALLOCATE(fresh_cell(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for fresh_cell'
      STOP
    ENDIF ! ierr
 
    ALLOCATE(exp_weight(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for exp_weight'
      STOP
    ENDIF ! ierr

    ALLOCATE(bodyNum(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for bodyNum'
      STOP
    ENDIF ! ierr

    ALLOCATE(bodyNumOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for bodyNumOld'
      STOP
    ENDIF ! ierr

    ALLOCATE(boundCell(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for boundCell'
      STOP
    ENDIF ! ierr

!    ALLOCATE(iup(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(iup(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iup'
      STOP
    ENDIF ! ierr

!    ALLOCATE(ium(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(ium(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for ium'
      STOP
    ENDIF ! ierr

!    ALLOCATE(jup(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(jup(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jup'
      STOP
    ENDIF ! ierr

!    ALLOCATE(jum(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(jum(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jum'
      STOP
    ENDIF ! ierr

!    ALLOCATE(kup(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(kup(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kup'
      STOP
    ENDIF ! ierr

!    ALLOCATE(kum(nxc,nyc,nzc),STAT=iErr)
    ALLOCATE(kum(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kum'
      STOP
    ENDIF ! ierr

    ALLOCATE(iBound(0:nx+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iBound'
      STOP
    ENDIF ! ierr

    ALLOCATE(jBound(0:ny+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jBound'
      STOP
    ENDIF ! ierr

    ALLOCATE(kBound(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kBound'
      STOP
    ENDIF ! ierr

!   Arrays relevant for 2nd Upwinding (Added by Rupesh)
!   ---------------------------------------------------
    ALLOCATE(iupp(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iupp'
      STOP
    ENDIF ! ierr  
    
    ALLOCATE(iumm(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for iumm'
      STOP
    ENDIF ! ierr  
 
    ALLOCATE(jupp(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jupp'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(jumm(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for jumm'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(kupp(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kupp'
      STOP
    ENDIF ! ierr  
  
    ALLOCATE(kumm(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'BOUNDARY_allocate_memory: Memory Allocation Error for kumm'
      STOP
    ENDIF ! ierr  

!   Initialize Arrays
!   -----------------
    iblank=0
    IF (nBody_membrane > 0 ) iblank_memb=0
    iblank_solid=0

    ghostCellMark=0
    ghostCellMemb=0
    ghostCellSolid=0

    fresh_cell=0
    bodyNum=0
    bodyNumOld=0
    boundCell=0

    exp_weight=0.0_CGREAL
    
    iup=0.0_CGREAL
    ium=0.0_CGREAL
    jup=0.0_CGREAL
    jum=0.0_CGREAL
    kup=0.0_CGREAL
    kum=0.0_CGREAL
    
    iBound=0.0_CGREAL
    jBound=0.0_CGREAL
    kBound=0.0_CGREAL
    
    IF (myCoords(1)==0)       iBound(1)   = 1.0_CGREAL
    IF (myCoords(1)==Np(1)-1) iBound(nxc) = 1.0_CGREAL

    IF (myCoords(2)==0)       jBound(1)   = 1.0_CGREAL
    IF (myCoords(2)==Np(2)-1) jBound(nyc) = 1.0_CGREAL

    kBound(1)   = 1.0_CGREAL
    kBound(nzc) = 1.0_CGREAL


    iupp=0.0_CGREAL
    iumm=0.0_CGREAL
    jupp=0.0_CGREAL
    jumm=0.0_CGREAL
    kupp=0.0_CGREAL
    kumm=0.0_CGREAL

!=======================================================================
! since elliptic & general cylinder will be converted into
! 3D unstruc surfaces, we need to determine memory requirement for these
!=======================================================================
    DO iBody = 1, nBody

      SELECT CASE (canonical_body_type(iBody))

        CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
           nPtsBodyMarker(iBody)    = nPtsBodyMarkerOrig(iBody)*nz
           totNumTriElem(iBody)     = 2*nPtsBodyMarkerOrig(iBody)*nzc

        CASE(ELLIPSOID)
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
           totNumTriElem(iBody)     = 2*nPtsBodyMarker(iBody) + 5

        CASE(UNSTRUCTURED_SURFACE)
           nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)

      END SELECT ! canonical_body_type

    ENDDO ! iBody
!=======================================================================
!   Allocate memory for markers and unstructured surfaces
!=======================================================================
    IF (monitorON) WRITE(STDOUT,'(5X,A)') 'Allocating Memory for Markers'
    CALL MARKER_allocate_memory()

    IF (monitorON) WRITE(STDOUT,'(5X,A)') 'Allocating Memory for unstructured Surfaces'
    CALL UNSTRUC_allocate_memory()

    IF ( NINT(gcmFLAG) == 1) THEN
      IF (monitorON) THEN
        WRITE(STDOUT,'(5X,A)') 'Allocating Memory for GCM static arrays'
        WRITE(STDOUT,*)
      END IF
      CALL GCM_allocate_static_arrays
    END IF ! gcmFLAG


END SUBROUTINE BOUNDARY_allocate_memory
!------------------------------------------------------------------------



SUBROUTINE MARKER_allocate_memory()

! ---------------------------------------------------------
!  Allocate Memory for various arrays pertinent to markers
! ---------------------------------------------------------

!!    USE global_parameters
    USE flow_parameters
!!    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays
!!    USE GCM_arrays
    USE implicit_coupling_parameters ! Added by Rajneesh
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    USE scalar
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    

    IMPLICIT NONE

!!    INTEGER :: iBody, iErr, nBodyPtsMax
    INTEGER :: iErr



!   Marker point arrays
!   -------------------
    nPtsMax = 1
    IF ( internal_boundary_present == INTR_BOUND_PRESENT .AND. &
         body_type == CANONICAL                                ) THEN
      nPtsMax = MAXVAL(nPtsBodyMarker(:))
    ELSE
      nPtsMax = 1
    ENDIF

    IF (monitorON) THEN
      WRITE(STDOUT,'(7X,A,1X,I8)') 'Maximum number of Points:  ', nPtsMax
      WRITE(STDOUT,*)
    END IF

!   Allocate Arrays
!   ---------------
    ALLOCATE(xBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for xBodyMarker'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for yBodyMarker'
      STOP
    ENDIF ! ierr

    ALLOCATE(zBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for zBodyMarker'
      STOP
    ENDIF ! ierr

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      
   ALLOCATE(xBodyMarkerPrev(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for xBodyMarkerPrev'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yBodyMarkerPrev(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for yBodyMarkerPrev'
      STOP
    ENDIF ! ierr

    ALLOCATE(zBodyMarkerPrev(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for zBodyMarkerPrev'
      STOP
    ENDIF ! ierr	
    
    ALLOCATE(X_Marker(nBody,nPtsMax,2),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for X_Marker'
      STOP
    ENDIF ! ierr   
    
    ALLOCATE(undot(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for undot'
      STOP
    ENDIF ! ierr     
    
    ALLOCATE(udotn(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for udotn'
      STOP
    ENDIF ! ierr    
    
    ALLOCATE(adpw_l(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for adpw_l'
      STOP
    ENDIF ! ierr       
    
    ALLOCATE(sspw_l(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for sspw_l'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(ppw_l(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for ppw_l'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(totalpw_l(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for totalpw_l'
      STOP
    ENDIF ! ierr            
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)	
	
	!   Added by Rajneesh

   ALLOCATE(xBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for xBodyMarkerOld'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for yBodyMarkerOld'
      STOP
    ENDIF ! ierr

    ALLOCATE(zBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for zBodyMarkerOld'
      STOP
    ENDIF ! ierr
	
	

   ALLOCATE(edge_node(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for edge_node'
      STOP
    ENDIF ! ierr


    ALLOCATE(sBodyMarker(nBody,nPtsMax+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for sBodyMarker'
      STOP
    ENDIF ! ierr
    
!!    ALLOCATE(dsBodyMarker(nBody,nPtsMax+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for dsBodyMarker'
!!      STOP
!!    ENDIF ! ierr
    
!!    ALLOCATE(xNormBodyMarker(nBody,nPtsMax+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for xNormBodyMarker'
!!      STOP
!!    ENDIF ! ierr
    
!!    ALLOCATE(yNormBodyMarker(nBody,nPtsMax+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for yNormBodyMarker'
!!      STOP
!!    ENDIF ! ierr
    
    ALLOCATE(uBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for uBodyMarker'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for vBodyMarker'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(wBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for wBodyMarker'
      STOP
    ENDIF ! ierr
	
	    ALLOCATE(pBodyMarker(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for pBodyMarker'
      STOP
    ENDIF ! ierr
	
	
	
	
!   Added by Rajneesh for implicit coupling

!____________	
	  ALLOCATE(uBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for uBodyMarker'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for vBodyMarker'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(wBodyMarkerOld(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for wBodyMarker'
      STOP
    ENDIF ! ierr 
!_____________
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
	  ALLOCATE(uBodyMarker_rel(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for uBodyMarker_relr'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vBodyMarker_rel(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for vBodyMarker_rel'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(wBodyMarker_rel(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for wBodyMarker_rel'
      STOP
    ENDIF ! ierr  
    
    
      ALLOCATE(uBodyMarker_rel_old(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for uBodyMarker_rel_old'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vBodyMarker_rel_old(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for vBodyMarker_rel_old'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(wBodyMarker_rel_old(nBody,nPtsMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'MARKER_allocate_memory: Memory Allocation Error for wBodyMarker_rel_old'
      STOP
    ENDIF ! ierr 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

!!    ALLOCATE(axBodyMarker(nBody,nPtsMax),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for axBodyMarker'
!!      STOP
!!    ENDIF ! ierr
    
!!    ALLOCATE(ayBodyMarker(nBody,nPtsMax),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for ayBodyMarker'
!!      STOP
!!    ENDIF ! ierr
     
!!    ALLOCATE(azBodyMarker(nBody,nPtsMax),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!       'MARKER_allocate_memory: Memory Allocation Error for azBodyMarker'
!!      STOP
!!    ENDIF ! ierr 

!   Initialize Arrays
!   -----------------
    xBodyMarker     = 0.0_CGREAL
    yBodyMarker     = 0.0_CGREAL
    zBodyMarker     = 0.0_CGREAL

    edge_node       = 0


!!    sBodyMarker     = 0.0_CGREAL
!!    dsBodyMarker    = 0.0_CGREAL
!!    xNormBodyMarker = 0.0_CGREAL
!!    yNormBodyMarker = 0.0_CGREAL

    uBodyMarker     = 0.0_CGREAL
    vBodyMarker     = 0.0_CGREAL
    wBodyMarker     = 0.0_CGREAL
    pBodyMarker     = 0.0_CGREAL

!!    axBodyMarker    = 0.0_CGREAL
!!    ayBodyMarker    = 0.0_CGREAL
!!    azBodyMarker    = 0.0_CGREAL

END SUBROUTINE MARKER_allocate_memory
!-----------------------------------------------------------------------



SUBROUTINE UNSTRUC_allocate_memory()
    
! -----------------------------------------------------------------------
!  Allocate Memory for various arrays pertinent to unstructured surface.
! -----------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE unstructured_surface_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)    
    USE implicit_coupling_parameters
    USE scalar
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    
    IMPLICIT NONE

    INTEGER :: nTriElemMax, iBody, iErr
    LOGICAL :: unstruc



!   Set initial values
!   ------------------
    nTriElemMax = 0
    unstruc     = .FALSE.

!   Set logical flag if unstructured surface is present
!   ---------------------------------------------------
    DO iBody = 1, nBody
      IF ( canonical_body_type(iBody) == ELLIPTIC_CYLINDER  .OR.   &
           canonical_body_type(iBody) == GENERAL_CYLINDER   .OR.   &
           canonical_body_type(iBody) == ELLIPSOID          .OR.   &
           canonical_body_type(iBody) == UNSTRUCTURED_SURFACE ) unstruc  = .TRUE. 
    ENDDO ! iBody

!   Allocate arrays and initialize
!   ------------------------------
    IF (unstruc) THEN
      nTriElemMax = MAXVAL(totNumTriElem(:))
      IF (monitorON) THEN
        WRITE(STDOUT,'(7X,A,1X,I8)') 'Maximum number of Elements:',nTriElemMax
        WRITE(STDOUT,*)
      END IF

      ALLOCATE(triElemNeig(nBody,3,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemNeig'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemtang1X(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang1X'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemtang1Y(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang1Y'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemtang1Z(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang1Z'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemtang2X(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang2X'
        STOP
      ENDIF ! ierr

      ALLOCATE(triElemtang2Y(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang2Y'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemtang2Z(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemtang2Z'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemNormX(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemNormX'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemNormY(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemNormY'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemNormZ(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemNormZ'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemCentX(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentX'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemCentY(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentY'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemCentZ(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentZ'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemArea(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemArea'
        STOP
      ENDIF ! ierr
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
!dddddddddddddddddddddddddddddddddddddddddddddddddd
      ALLOCATE(vorTriElemCent0(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vorTriElem0'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(vorTriElemCent1(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vorTriElem1'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(vorNormCent0(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vorNorm0'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(vorNormCent1(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vorNorm1'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(ssBodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for ssBodyMarker'
        STOP
      ENDIF ! ierr 

      ALLOCATE(dssBodyMarker(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for dssBodyMarker'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(pTriElemCent0(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for pTriElemCent0'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(pTriElemCent1(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for pTriElemCent1'
        STOP
      ENDIF ! ierr 
      
      !ALLOCATE(p_force0(nBody,nPtsMax),STAT=iErr)
      !IF ( iErr /= ERR_NONE ) THEN
      !  WRITE(STDOUT,*) &
      !  'UNSTRUC_allocate_memory: Memory Allocation Error for p_force0'
      !  STOP
      !ENDIF ! ierr 
      
      !ALLOCATE(p_force1(nBody,nPtsMax),STAT=iErr)
      !IF ( iErr /= ERR_NONE ) THEN
      !  WRITE(STDOUT,*) &
      !  'UNSTRUC_allocate_memory: Memory Allocation Error for p_force1'
      !  STOP
      !ENDIF ! ierr 
      
      ALLOCATE(x_force(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for x_force'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(y_force(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for y_force'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(z_force(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for z_force'
        STOP
      ENDIF ! ierr 
!dddddddddddddddddddddddddddddddddddddddddddddddddd

      ALLOCATE(XgradxBodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for XgradxBodyMarker'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(XgradyBodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for XgradyBodyMarker'
        STOP
      ENDIF ! ierr    

      ALLOCATE(XgradzBodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for XgradzBodyMarker'
        STOP
      ENDIF ! ierr       
        
      ALLOCATE(X_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for X_BodyMarker'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(X_BodyMarker_prev(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for X_BodyMarker_prev'
        STOP
      ENDIF ! ierr    

      ALLOCATE(advection_x_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for advection_x_BodyMarker'
        STOP
      ENDIF ! ierr 

      ALLOCATE(advection_y_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for advection_y_BodyMarker'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(advection_z_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for advection_z_BodyMarker'
        STOP
      ENDIF ! ierr   
      
      ALLOCATE(Lamb_x_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Lamb_x_BodyMarker'
        STOP
      ENDIF ! ierr 

      ALLOCATE(Lamb_y_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Lamb_y_BodyMarker'
        STOP
      ENDIF ! ierr 
      
      ALLOCATE(Lamb_z_BodyMarker(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Lamb_z_BodyMarker'
        STOP
      ENDIF ! ierr        
      
 
      ALLOCATE(Grad_u_square_x(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_u_square_x'
        STOP
      ENDIF ! ierr           

      ALLOCATE(Grad_u_square_y(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_u_square_y'
        STOP
      ENDIF ! ierr  

      ALLOCATE(Grad_u_square_z(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_u_square_z'
        STOP
      ENDIF ! ierr  
      
      
      ALLOCATE(Grad_potential_u_square_x(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_potential_u_square_x'
        STOP
      ENDIF ! ierr           

      ALLOCATE(Grad_potential_u_square_y(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_potential_u_square_y'
        STOP
      ENDIF ! ierr  

      ALLOCATE(Grad_potential_u_square_z(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_potential_u_square_z'
        STOP
      ENDIF ! ierr
      
      
      ALLOCATE(Grad_vorticity_u_square_x(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_vorticity_u_square_x'
        STOP
      ENDIF ! ierr           

      ALLOCATE(Grad_vorticity_u_square_y(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_vorticity_u_square_y'
        STOP
      ENDIF ! ierr  

      ALLOCATE(Grad_vorticity_u_square_z(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_vorticity_u_square_z'
        STOP
      ENDIF ! ierr      

      
      
      ALLOCATE(u_n(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_n'
        STOP
      ENDIF ! ierr

      
      ALLOCATE(u_n_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_n_prev'
        STOP
      ENDIF ! ierr
           
            
      ALLOCATE(u_n_dot(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_n_dot'
        STOP
      ENDIF ! ierr            
      
      ALLOCATE(u_n1_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_n1_prev'
        STOP
      ENDIF ! ierr   
      
      ALLOCATE(u_n2_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_n2_prev'
        STOP
      ENDIF ! ierr          
      
              
      ALLOCATE(triElemCentXOld(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentXOld'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemCentYOld(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentYOld'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemCentZOld(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemCentZOld'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(triElemAreaOld(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for triElemAreaOld'
        STOP
      ENDIF ! ierr
      

      ALLOCATE(u_dot_n(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for u_dot_n'
        STOP
      ENDIF ! ierr
      
      ALLOCATE(up_dot_n(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for up_dot_n'
        STOP
      ENDIF ! ierr 
      
      
      ALLOCATE(uv_dot_n(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for uv_dot_n'
        STOP
      ENDIF ! ierr           

      ALLOCATE(advection_n(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for advection_n'
        STOP
      ENDIF ! ierr
      
      ALLOCATE(Lamb_normal(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Lamb_normal'
        STOP
      ENDIF ! ierr      
    
      ALLOCATE(Grad_u_square_n(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_u_square_n'
        STOP
      ENDIF ! ierr     
      
      ALLOCATE(Grad_potential_u_square_n(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_potential_u_square_n'
        STOP
      ENDIF ! ierr  
      
      ALLOCATE(Grad_vorticity_u_square_n(nBody,nTriElemMax,2),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for Grad_vorticity_u_square_n'
        STOP
      ENDIF ! ierr              
       
      
      ALLOCATE(uBodyCent_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for uBodyCent_prev'
        STOP
      ENDIF ! ierr
      
      
      ALLOCATE(vBodyCent_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vBodyCent_prev'
        STOP
      ENDIF ! ierr
      
 
      ALLOCATE(wBodyCent_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for wBodyCent_prev'
        STOP
      ENDIF ! ierr      
      

      ALLOCATE(uBodyCent_p_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for uBodyCent_p_prev'
        STOP
      ENDIF ! ierr
      
      
      ALLOCATE(vBodyCent_p_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vBodyCent_p_prev'
        STOP
      ENDIF ! ierr
      
 
      ALLOCATE(wBodyCent_p_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for wBodyCent_p_prev'
        STOP
      ENDIF ! ierr          
      
      ALLOCATE(uBodyCent_v_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for uBodyCent_v_prev'
        STOP
      ENDIF ! ierr
      
      
      ALLOCATE(vBodyCent_v_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for vBodyCent_v_prev'
        STOP
      ENDIF ! ierr
      
 
      ALLOCATE(wBodyCent_v_prev(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for wBodyCent_v_prev'
        STOP
      ENDIF ! ierr         
      
      ALLOCATE(adpwl(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for adpwl'
        STOP
      ENDIF ! ierr      
      
      ALLOCATE(sspwl(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for sspwl'
        STOP
      ENDIF ! ierr        
      
      ALLOCATE(ppwl(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for ppwl'
        STOP
      ENDIF ! ierr   
      
      ALLOCATE(totalpwl(nBody,nTriElemMax),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for totalpwl'
        STOP
      ENDIF ! ierr             

!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)   
       
      ALLOCATE(pointOutsideBodyX(nBody),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for pointOutsideBodyX'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(pointOutsideBodyY(nBody),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for pointOutsideBodyY'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(pointOutsideBodyZ(nBody),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for pointOutsideBodyZ'
        STOP
      ENDIF ! ierr
       
      ALLOCATE(surfArea(nBody),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
        WRITE(STDOUT,*) &
        'UNSTRUC_allocate_memory: Memory Allocation Error for surfArea'
        STOP
      ENDIF ! ierr
       
!     Initialize Arrays
!     -----------------
      triElemNeig = 0

      triElemtang1X = 0.0_CGREAL
      triElemtang1Y = 0.0_CGREAL
      triElemtang1Z = 0.0_CGREAL

      triElemtang2X = 0.0_CGREAL
      triElemtang2Y = 0.0_CGREAL
      triElemtang2Z = 0.0_CGREAL

      triElemNormX  = 0.0_CGREAL
      triElemNormY  = 0.0_CGREAL
      triElemNormZ  = 0.0_CGREAL

      triElemCentX  = 0.0_CGREAL
      triElemCentY  = 0.0_CGREAL
      triElemCentZ  = 0.0_CGREAL

      triElemArea       = 0.0_CGREAL
      pointOutsideBodyX = 0.0_CGREAL
      pointOutsideBodyY = 0.0_CGREAL
      pointOutsideBodyZ = 0.0_CGREAL
      surfArea          = 0.0_CGREAL
      
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
      IF(nread == 0)THEN
         X_BodyMarker_prev      = 0.0_CGREAL  
         advection_x_BodyMarker = 0.0_CGREAL
         advection_y_BodyMarker = 0.0_CGREAL
         advection_z_BodyMarker = 0.0_CGREAL
         Lamb_x_BodyMarker      = 0.0_CGREAL
         Lamb_y_BodyMarker      = 0.0_CGREAL
         Lamb_z_BodyMarker      = 0.0_CGREAL         
         unsteady1_prev         = 0.0_CGREAL 
         uBodyCent_prev         = 0.0_CGREAL 
         vBodyCent_prev         = 0.0_CGREAL
         wBodyCent_prev         = 0.0_CGREAL     
         uBodyCent_p_prev       = 0.0_CGREAL 
         vBodyCent_p_prev       = 0.0_CGREAL
         wBodyCent_p_prev       = 0.0_CGREAL       
         uBodyCent_v_prev       = 0.0_CGREAL 
         vBodyCent_v_prev       = 0.0_CGREAL
         wBodyCent_v_prev       = 0.0_CGREAL          
      ENDIF
      
      X_BodyMarker = 0.0_CGREAL  
      adpwl        = 0.0_CGREAL        
      sspwl        = 0.0_CGREAL
      ppwl         = 0.0_CGREAL    
      totalpwl     = 0.0_CGREAL                
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)      

    END IF ! unstruc 

END SUBROUTINE UNSTRUC_allocate_memory
