! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  Johns Hopkins University
!
!  ParVICAR3D, a parallelized version of VICAR3D.
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
!  Filename: CORE_MEMORY_ALLOCATE.PAR.F90
!  Latest Modification: Dec, 29 2010 (PAT 2.1.0)
!  by Rajneesh Bhardwaj
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     allocate_memory()
!     deallocate_memory()
! --------------------------------------------------------------------



SUBROUTINE allocate_memory()

!   ---------------------------------------------------------
!    This subroutine allocate memory for the main variables.
!   ---------------------------------------------------------
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE boundary_arrays
    USE nlold_arrays
    USE multiuse_arrays
    USE grid_arrays
    USE solver_arrays
    USE solver_ad_arrays
    USE stat_arrays
    USE mg_arrays
    USE implicit_coupling_parameters  ! Added by Rajneesh
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
	USE scalar	
	USE mpi
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)


    IMPLICIT NONE
    
    INTEGER :: iErr, nmax



!   Global grid variables
!   ---------------------
    ALLOCATE(x(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for x'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(y(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for y'
      CALL flow_stop
      STOP
    ENDIF ! iErr   
      
    ALLOCATE(z(0:nz+1) ,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for z'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
       
    ALLOCATE(xc(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for xc'
      CALL flow_stop
      STOP
    ENDIF ! iErr    

    ALLOCATE(yc(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for yc'
      CALL flow_stop
      STOP
    ENDIF ! iErr    

    ALLOCATE(zc(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for zc'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
       
    ALLOCATE(dx(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dx'
      CALL flow_stop
      STOP
    ENDIF ! iErr    

    ALLOCATE(dy(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dy'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
       
    ALLOCATE(dz(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dz'
      CALL flow_stop
      STOP
    ENDIF ! iErr
        
    ALLOCATE(dxc(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dxc'
      CALL flow_stop
      STOP
    ENDIF ! iErr   

    ALLOCATE(dyc(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dyc'
      CALL flow_stop
      STOP
    ENDIF ! iErr
       
    ALLOCATE(dzc(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dzc'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
      
    ALLOCATE(dxinv(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dxinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(dyinv(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dyinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
     
    ALLOCATE(dzinv(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dzinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
     
    ALLOCATE(dxcinv(0:nx_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dxcinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(dycinv(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dycinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(dzcinv(0:nz+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dzcinv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(fx(0:nxc_GLBL+1) ,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for fx'
      CALL flow_stop
      STOP
    ENDIF ! iErr   

    ALLOCATE(fy(0:ny_GLBL+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for fy'
      CALL flow_stop
      STOP
    ENDIF ! iErr    

    ALLOCATE(fz(0:nz+1) ,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for fz'
      CALL flow_stop
      STOP
    ENDIF ! iErr   

    ALLOCATE(damper(0:nxc_GLBL+1) ,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for fx'
      CALL flow_stop
      STOP
    ENDIF ! iErr   

!   Flow field variables - Velocity
!   -------------------------------
    ALLOCATE(u(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for u'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(v(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for v'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(w(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for w'
      CALL flow_stop
      STOP
    ENDIF ! iErr
	
	
!   Added by Rajneesh for implicit coupling
!________________________________	
	ALLOCATE(uOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for uOld'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(vOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for vOld'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(wOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for wOld'
      CALL flow_stop
      STOP
    ENDIF ! iErr
	
!__________________________	
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    ALLOCATE(uGhostOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pGhostOld'
      STOP
    ENDIF ! iErr  
    
    
    ALLOCATE(vGhostOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pGhostOld'
      STOP
    ENDIF ! iErr  

    ALLOCATE(wGhostOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pGhostOld'
      STOP
    ENDIF ! iErr            
    
    
    ALLOCATE(pPrime_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pPrime_old'
      CALL flow_stop
      STOP
    ENDIF ! iErr   
    
    ALLOCATE(p_old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for p_old'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(pGhostOld(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pGhostOld'
      STOP
    ENDIF ! iErr 
    
    
    ALLOCATE(face_ue_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_u_Old'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_vn_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_v_Old'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_wf_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_w_Old'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_uw_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_u_Old'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_vs_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_v_Old'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_wb_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_w_Old'
      STOP
    ENDIF ! iErr   
    
    ALLOCATE(nlu_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlu_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nlv_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlv_Old'
      STOP
    ENDIF ! iErr
   
    ALLOCATE(nlw_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlw_Old'
      STOP
    ENDIF ! iErr     
    
    ALLOCATE(uTilde_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for uTilde_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(vTilde_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for vTilde_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(wTilde_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for wTilde_Old'
      STOP
    ENDIF ! iErr       
    
    ALLOCATE(boundPresSource_old(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for boundPresSource_old'
      STOP
    ENDIF ! iErr            
    
    ALLOCATE(bcxu_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxu_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcxv_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxv_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcxw_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxw_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyu_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyu_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyv_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyv_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyw_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyw_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczu_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczu_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczv_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczv_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczw_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczw_Old'
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(pgradx1_Old(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradx1_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradx2_Old(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradx2_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgrady1_Old(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgrady1_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgrady2_Old(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgrady2_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradz1_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradz1_Old'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradz2_Old(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradz2_Old'
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(nluold_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nluold_Old'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nlvold_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlvold_Old'
      STOP
    ENDIF ! iErr
   
    ALLOCATE(nlwold_Old(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlwold_Old'
      STOP
    ENDIF ! iErr    
    

    ALLOCATE(omega_cross_u(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for omega_cross_u'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(Div_Lamb(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Div_Lamb'
      CALL flow_stop
      STOP
    ENDIF ! iErr     
  
    ALLOCATE(u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr     
    
    ALLOCATE(u_squareGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for u_squareGhost'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(potential_u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for potential_u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(vorticity_u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for vorticity_u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr                 
  
  
    ALLOCATE(Grad_u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Grad_u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    
    ALLOCATE(Grad_potential_u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Grad_potential_u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(Grad_vorticity_u_square(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Grad_vorticity_u_square'
      CALL flow_stop
      STOP
    ENDIF ! iErr  
    
    
    ALLOCATE(AB2_advection_Xx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for AB2_advection_Xx'
      CALL flow_stop
      STOP
    ENDIF ! iErr  
   
   
    ALLOCATE(Gradu2_Xx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Gradu2_Xx'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
             
    
    ALLOCATE(Div_AB2_Xx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Div_AB2_Xx'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(Div_AB2_Xx_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Div_AB2_Xx_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr   
    
    ALLOCATE(Div_Gradu2_Xx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Div_Gradu2_Xx'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    
    ALLOCATE(Div_Gradu2_Xx_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Div_Gradu2_Xx_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr         
    
    
    ALLOCATE(Grad_u_square_dot_Grad_Dx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Grad_u_square_dot_Grad_Dx'
      CALL flow_stop
      STOP
    ENDIF ! iErr              
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
    ALLOCATE(Advection_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Advection_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    
    ALLOCATE(Unsteady_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Unsteady_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(Shear_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Shear_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
           

    ALLOCATE(F_i(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for F_i'
      CALL flow_stop
      STOP
    ENDIF ! iErr


    ALLOCATE(Vortex_force(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Vortex_force'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(diff_prev(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for diff_prev'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(diff_star(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for diff_star'
      CALL flow_stop
      STOP
    ENDIF ! iErr   
       
    ALLOCATE(gradP(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for gradP'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    
    !CALL mpi_barrier(MPI_COMM_WORLD, ierr)    
    ALLOCATE(gradu2(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for gradu2'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    !CALL mpi_barrier(MPI_COMM_WORLD, ierr) 
    
    ALLOCATE(PuPt(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for PuPt'
      CALL flow_stop
      STOP
    ENDIF ! iErr            
    
    ALLOCATE(PvPt(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for PvPt'
      CALL flow_stop
      STOP
    ENDIF ! iErr         
 
    ALLOCATE(PwPt(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for PwPt'
      CALL flow_stop
      STOP
    ENDIF ! iErr        
!ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd  
    ALLOCATE(AB2_advection(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for AB2_advection'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
      
    ALLOCATE(advection(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for advection'
      CALL flow_stop
      STOP
    ENDIF ! iErr        
        
    ALLOCATE(advection_prev(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1,3),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for advection_prev'
      CALL flow_stop
      STOP
    ENDIF ! iErr         
        
        
    ALLOCATE(scalarX(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for scalarX'
      CALL flow_stop
      STOP
    ENDIF ! iErr   
    
    
    ALLOCATE(potential(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for potential'
      CALL flow_stop
      STOP
    ENDIF ! iErr     
    
    ALLOCATE(Xgradx(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradx'
      CALL flow_stop
      STOP
    ENDIF ! iErr      
    ALLOCATE(Xgrady(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradx'
      CALL flow_stop
      STOP
    ENDIF ! iErr    
    ALLOCATE(Xgradz(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradx'
      CALL flow_stop
      STOP
    ENDIF ! iErr  
    
    ALLOCATE(Potential_u(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Potential_u'
      CALL flow_stop
      STOP
    ENDIF ! iErr  
    ALLOCATE(Potential_v(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Potential_v'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    ALLOCATE(Potential_w(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Potential_w'
      CALL flow_stop
      STOP
    ENDIF ! iErr           
               
               
    ALLOCATE(Vorticity_u(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Vorticity_u'
      CALL flow_stop
      STOP
    ENDIF ! iErr  
    ALLOCATE(Vorticity_v(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Vorticity_v'
      CALL flow_stop
      STOP
    ENDIF ! iErr 
    ALLOCATE(Vorticity_w(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Vorticity_w'
      CALL flow_stop
      STOP
    ENDIF ! iErr                 
    
    ALLOCATE(Xgradx1(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradx1'
      STOP
    ENDIF ! iErr 
    ALLOCATE(Xgradx2(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradx2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Xgrady1(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgrady1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Xgrady2(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgrady2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Xgradz1(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradz1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Xgradz2(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Xgradz2'
      STOP
    ENDIF ! iErr
    
    
    ALLOCATE(Phigradx1(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigradx1'
      STOP
    ENDIF ! iErr 
    ALLOCATE(Phigradx2(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigradx2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Phigrady1(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigrady1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Phigrady2(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigrady2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Phigradz1(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigradz1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(Phigradz2(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for Phigradz2'
      STOP
    ENDIF ! iErr    

    ALLOCATE(boundScalarSource(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for boundScalarSource'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(boundPhiSource(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for boundPhiSource'
      STOP
    ENDIF ! iErr    

!   Flow field variables- AD RHS
!   ----------------------------
    ALLOCATE(zero(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for zero'
      CALL flow_stop
      STOP
    ENDIF ! iErr       
    
    ALLOCATE(XGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for XGhost'
      STOP
    ENDIF ! iErr 
    
    ALLOCATE(PhiGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for PhiGhost'
      STOP
    ENDIF ! iErr         
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)	

    ALLOCATE(face_ue(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_u'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_vn(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_v'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_wf(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_w'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_uw(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_u'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_vs(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_v'
      STOP
    ENDIF ! iErr

    ALLOCATE(face_wb(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face_w'
      STOP
    ENDIF ! iErr

!   Flow field variables - Pressure
!   -------------------------------
    ALLOCATE(p(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for p'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pPrime(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pPrime'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradx1(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradx1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradx2(1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradx2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgrady1(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgrady1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgrady2(1-Ngl:nxc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgrady2'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradz1(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradz1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pgradz2(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradz2'
      STOP
    ENDIF ! iErr

    ALLOCATE(boundPresSource(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pgradz2'
      STOP
    ENDIF ! iErr

!   Flow field variables- AD RHS
!   ----------------------------
    ALLOCATE(nlu(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlu'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nlv(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlv'
      STOP
    ENDIF ! iErr
   
    ALLOCATE(nlw(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlw'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nluold(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nluold'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nlvold(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlvold'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(nlwold(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for nlwold'
      STOP
    ENDIF ! iErr

!   Flow field variables- Van-Kan velocity
!   --------------------------------------
    ALLOCATE(uTilde(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for uTilde'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(vTilde(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for vTilde'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(wTilde(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for wTilde'
      STOP
    ENDIF ! iErr

!   Boundary variables
!   ------------------
    ALLOCATE(bcxu(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxu'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcxv(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcxw(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxw'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyu(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyu'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyv(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyv'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bcyw(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyw'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczu(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczu'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczv(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczv'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(bczw(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczw'
      STOP
    ENDIF ! iErr

!   AD Solver coefficients
!   ----------------------
    ALLOCATE(amx(nxc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amx'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(acx(nxc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for acx'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(apx(nxc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for apx'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(amy(nyc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amy'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(acy(nyc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for acy'
      STOP
    ENDIF ! iErr

    ALLOCATE(apy(nyc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for apy'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(amz(nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amz'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(acz(nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for acz'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(apz(nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for acz'
      STOP
    ENDIF ! iErr

    nmax = max(nxc,nyc,nzc)
    ALLOCATE(rhs(nmax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for rhs'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(dummy(nmax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for dummy'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(face1(nmax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face1'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(face2(nmax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for face2'
      STOP
    ENDIF ! iErr

    ALLOCATE(amx_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amx_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(apx_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for apx_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(amy_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amy_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(apy_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for apy_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr

    ALLOCATE(amz_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for amz_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr
    
    ALLOCATE(apz_ad(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for apz_ad'
      CALL flow_stop
      STOP
    ENDIF ! iErr

!   Statistics variables
!   --------------------
    IF (nStat > STATS_NONE ) THEN
      ALLOCATE(uav(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for uav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(vav(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for vav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(wav(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for wav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(pav(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for pav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(uvAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for uvav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(uwAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for uwav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(vwAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for vwav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(uuAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for uuav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(vvAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for vvav'
         STOP
      ENDIF ! iErr
      
      ALLOCATE(wwAv(nxc,nyc,nzc),STAT=iErr)
      IF ( iErr /= ERR_NONE ) THEN
         WRITE(STDOUT,*) &
         'Allocate_memory: Memory Allocation Error for wwav'
         STOP
      ENDIF ! iErr
    ENDIF ! nStat

!   New unified formulation of SSM and GCM variables
!   ------------------------------------------------
!!    ALLOCATE(imv(0:nx+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for imv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(ipv(0:nx+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for ipv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(immv(0:nx+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for immv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(ippv(0:nx+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for ippv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(jmv(0:ny+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for jmv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(jpv(0:ny+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for jpv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(jmmv(0:ny+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for jmmv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(jppv(0:ny+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for jppv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(kmv(0:nz+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for kmv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(kpv(0:nz+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for kpv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(kmmv(0:nz+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for kmmv'
!!      STOP
!!    ENDIF ! iErr
    
!!    ALLOCATE(kppv(0:nz+1),STAT=iErr)
!!    IF ( iErr /= ERR_NONE ) THEN
!!      WRITE(STDOUT,*) &
!!      'Allocate_memory: Memory Allocation Error for kppv'
!!      STOP
!!    ENDIF ! iErr

    ALLOCATE(uGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for uGhost'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(vGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for vGhost'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(wGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for wGhost'
      STOP
    ENDIF ! iErr
    
    ALLOCATE(pGhost(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for pGhost'
      STOP
    ENDIF ! iErr
  
!   Viscosity variables
!   -------------------
    ALLOCATE(viscTot(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(bcxvisc(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcxvisc'
      STOP
    ENDIF ! iErr
      
    ALLOCATE(bcyvisc(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bcyvisc'
      STOP
    ENDIF ! iErr

    ALLOCATE(bczvisc(0:nxc+1,0:nyc+1,0:nzc+1),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for bczvisc'
      STOP
    ENDIF ! iErr

!   Temporary variables
!   -------------------
    ALLOCATE(div(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(temp1(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(temp2(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(temp3(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

    ALLOCATE(vor(nxc,nyc,nzc),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'Allocate_memory: Memory Allocation Error for viscTot'
      STOP
    ENDIF ! iErr 

        
END SUBROUTINE allocate_memory
!------------------------------------------------------------------------------
