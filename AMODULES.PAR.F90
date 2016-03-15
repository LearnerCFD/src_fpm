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
!     Qian Xue
!     Rupesh Babu K. A.
!     Xudong Zheng 
!     Reza Ghias 
!     S. A. Mohsen Karimian
!
!  Filename: AMODULE.PAR.F90
!  Latest Modification: Jan, 25 2011 (PAT 2.1.1)
!  by xudong zheng 
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following modules:
!     global_parameters
!     flow_parameters
!     grid_arrays
!     flow_arrays
!     boundary_arrays
!     multiuse_arrays
!     pressure_arrays
!     nlold_arrays
!     solver_arrays
!     solver_ad_arrays
!     GCM_arrays
!     unstructured_surface_arrays
!     probe_parameters
!     stat_arrays
!     blasius_profile
!     mg_parameters
!     mg_arrays
!     usr_module
!     stat_vort_arrays
!     tahoe_parameters
! --------------------------------------------------------------------




MODULE global_parameters
   
    IMPLICIT NONE
 
    INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=14,R=30) 
   
    INTEGER, PARAMETER :: VISCOUS_FLOW             = 1, &
                          POTENTIAL_FLOW           = 2
 
    INTEGER, PARAMETER :: UNIFORM_GRID             = 1, &
                          NONUNIFORM_GRID          = 2
    
    INTEGER, PARAMETER :: BC_TYPE_DIRICHLET        = 1, & 
                          BC_TYPE_ZERO_GRADIENT    = 2, &
                          BC_TYPE_PULSATILE_INFLOW = 3, &  
                          BC_TYPE_SYMMETRY         = 4, &
                          BC_TYPE_PERIODIC         = 5, &
                          BC_TYPE_USER_SPECIFIED   = 6, &
                          BC_TYPE_SHEAR            = 7
      
    INTEGER, PARAMETER :: IT_SOLVER_TYPE_LSOR      = 1, & 
                          IT_SOLVER_TYPE_MG        = 2

    INTEGER, PARAMETER :: FIXED_BOUNDARY           = 1, & 
                          MOVING_BOUNDARY          = 2 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER, PARAMETER :: UNCOMBINED               = 0, & 
                          COMBINED                 = 1 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)                           

    INTEGER, PARAMETER :: STATIONARY               = 0, &
                          FORCED                   = 1, & 
                          FLOW_INDUCED             = 2, &
                          PRESCRIBED               = 3, &
                          FEA_FLOW_STRUC_INTERACTION = 4, & !SER_TO_PAR, QX, CH1     
                          EIGEN_MOTION             = 5, &
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)                          
			              FEA_TAHOE                = 6, &      ! Added by Rajneesh
			              PENDULUM                 = 7
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)			              
 

    INTEGER, PARAMETER :: PBC_DIRICHLET            = 1, &!added by H. Luo
                          PBC_NEUMANN              = 2

    INTEGER, PARAMETER :: ADAMS_BASHFORTH2         = 1, &!added by H. Luo
                          CRANK_NICOLSON1          = 2, &
                          CRANK_NICOLSON2          = 3   

    INTEGER, PARAMETER :: NONPOROUS_AND_NONSLIP    = 0, &
                          POROUS_OR_SLIP           = 1

    INTEGER, PARAMETER :: NONE                     = 0, &
                          GENERAL                  = 1, &
                          CANONICAL                = 2

    INTEGER, PARAMETER :: ELLIPTIC_CYLINDER        = 1, &
                          GENERAL_CYLINDER         = 2, &
                          ELLIPSOID                = 3, &
                          UNSTRUCTURED_SURFACE     = 4

    INTEGER, PARAMETER :: OPEN_MEMBRANE            = 1, &
                          CLOSED_MEMBRANE          = 2

    INTEGER, PARAMETER :: SOLID_BODY               = 1, &
                          MEMBRANE                 = 2     
    
    INTEGER, PARAMETER :: BODY_DIM2                = 2, &
                          BODY_DIM3                = 3
    
    INTEGER, PARAMETER :: INTR_BOUND_NONE          = 0, &
                          INTR_BOUND_PRESENT       = 1

    INTEGER, PARAMETER :: NO_VAN_KAN               = 0, &
                          VAN_KAN                  = 1

!    INTEGER, PARAMETER :: TECPLOT                  = 1, & 
!                          FIELDVIEW                = 2 

!    INTEGER, PARAMETER :: IBLANK_READ              = 1

    INTEGER, PARAMETER :: DIM_2D                   = 2, &
                          DIM_3D                   = 3

!    INTEGER, PARAMETER :: IBLANK_USED              = 1

    INTEGER, PARAMETER :: IBLANK_SLOW              = 0, &
                          IBLANK_FAST              = 1
    
    INTEGER, PARAMETER :: NO_INTERNAL_BOUNDARY     = 0, &
                          SSM_METHOD               = 1, &
                          GCM_METHOD               = 2
    
!    INTEGER, PARAMETER :: INVISCID                 = 1                  

    INTEGER, PARAMETER :: ICOORD                   = 1, &
                          JCOORD                   = 2, &
                          KCOORD                   = 3

    INTEGER, PARAMETER :: STATS_NONE               = 0

    INTEGER, PARAMETER :: STDOUT                   = 6

    INTEGER, PARAMETER :: INACTIVE                 = 0, &
                          ACTIVE                   = 1, &
                          ERR_NONE                 = 0

!   File Unit Numbers
!   -----------------
    INTEGER, PARAMETER :: ifuMarkerTRI             = 145
                          
!--------Added For FEA---------------------------------  ! SER_TO_PAR, QX, CH2
    INTEGER, PARAMETER :: PLANE_STRESS             = 1,&
                          PLANE_STRAIN             = 2,&
                          GENERAL_3DBODY           = 3

    INTEGER, PARAMETER :: STATIC_ANALYSIS          = 1,&
                          DYNAMIC_RESPONSE_CD      = 2,&
                          DYNAMIC_RESPONSE_NEWMARK = 3,&
                          DYNAMIC_CHARACTER_INVERSE = 4,&
                          DYNAMIC_CHARACTER_SURFACE = 5
    INTEGER, PARAMETER :: PENALTY_CONTACT_MODEL     = 1,&
                          GEO_CONTACT_MODEL         = 2
                           
    INTEGER, PARAMETER :: SCALAR_ON                   = 1, &
                          SCALAR_OFF                  = 0
    INTEGER, PARAMETER :: ADIABATIC                = 1, &
                          ISOTHERMAL               = 2
    INTEGER, PARAMETER :: UNIFORM_INIT_T           = 0, &
                          NONUNIFRORM_INIT_T       = 1

   END MODULE global_parameters
!------------------------------------------------------

   MODULE flow_parameters 

    USE global_parameters
    
    IMPLICIT NONE
 
    LOGICAL :: monitorON, monitorIT, dumpIT
    LOGICAL :: ImtheBOSS

    REAL(KIND=CGREAL), PARAMETER :: sidw  = 2.0_CGREAL ! parameter used in IDW interpolation

    REAL(KIND=CGREAL) :: pi
    REAL(KIND=CGREAL) :: KillFor2D
    
#   ifdef MPI
      REAL(KIND=CGREAL) :: MainTimeStart, MainTimeEnd
#   else
      INTEGER :: MainTimeStart, MainTimeEnd
#   endif
    REAL(KIND=CGREAL) :: system_wall_time, exit_safe_time
    INTEGER :: last_stop_ntime, iSAFE_EXIT
    LOGICAL :: exit_safe

!   Common domain and MPI variables for parallel and serial versions (SAMK)
!   -----------------------------------------------------------------------
    INTEGER :: myRank
    INTEGER :: Np(2)
    INTEGER :: myCoords(2)
    INTEGER :: myIs, myIe, myJs, myJe
    INTEGER :: myILL, myIUL, myJLL, myJUL
    INTEGER :: myImin, myImax, myJmin, myJmax
    INTEGER :: Ngl
!   -----------------------------------------------------------------------

    INTEGER  :: iacou, Nsa, iacou_rest   ! additional input parameters
    INTEGER  :: iCC, iSSMP, iDragLift    ! added for CutCell
    INTEGER  :: nread
    INTEGER  :: ndim
    INTEGER  :: flow_type
    INTEGER  :: nx_GLBL,ny_GLBL
    INTEGER  :: nxc_GLBL,nyc_GLBL
    INTEGER  :: nx,ny,nz
    INTEGER  :: nxc,nyc,nzc
    INTEGER  :: xgrid_unif,ygrid_unif,zgrid_unif
    INTEGER  :: FDprof
    INTEGER  :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
    INTEGER  :: no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag,ninit, qdump
!   INTEGER  :: format_dump
    INTEGER  :: ntime, ntime_start, ntime_skip
    INTEGER  :: it_solver_type, itermax
!   INTEGER  :: mlev,iwmg,mlw
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER  :: restart_FIM
    REAL(KIND=CGREAL) :: perturb_switch
    REAL(KIND=CGREAL) :: Fr
    INTEGER  :: nstart_FIM
    INTEGER  :: nFIM_dump
    INTEGER  :: boundary_motion, nBody, nBody_Solid, nBody_Membrane, nGroup_Combined  !nGroup_Combined
    INTEGER, DIMENSION(:), ALLOCATABLE  :: combined_Group_index
    INTEGER, DIMENSION(:), ALLOCATABLE  :: Surface_Integral    
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
!   INTEGER  :: ssm_Run

    INTEGER  :: boundary_formulation, extended_outflow
    INTEGER  :: iterMaxPoisson
    INTEGER  :: body_type

    character(11) :: ParLogName

    INTEGER  :: ifuInput, ifuParLog, ifuIblankIn, &
                ifuRstrtFlowIn, ifuRstrtFlowOut,  &
                ifuBodyIn, ifuBodyOut,            &
                ifuRstrtBodyIn, ifuRstrtBodyOut,  &
                ifuMarkerIn, ifuMarkerOut,        &
                ifuProbeIn, ifuProbeOut,          &
                ifuUnstrucSurfIn, ifuUnstrucSurfOut, &
                ifuOpenMembraneEdgeIn, ifurunningstatus, &
                ifuDyeIn, ifuScalarIn
                

    INTEGER  :: ifuStatOut
    INTEGER  :: ifuDragOut, ifuFreshCellOut
    
    INTEGER  :: ifuFlux ! used for FEA 
                                                                                      ! SER_to_PAR. QX. CH3

    INTEGER  :: internal_boundary_present, nPtsMax, iblankFast
    INTEGER  :: idxRstrt, idxStat!!, indexStatVort
    INTEGER  :: frac_step_type
   
    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarkerOrig,canonical_body_type, &
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
                                           body_dim,boundary_motion_type,combined_type,wall_type, & !combined_type
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
                                           membrane_type,unstruc_surface_type
 
    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarker

    INTEGER, DIMENSION(:), ALLOCATABLE  :: n_theta,n_phi

    INTEGER, DIMENSION(:), ALLOCATABLE  :: ntimePerCycle ! added by Haibo

!------------- new arrays -- added by Haibo
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: imv,ipv,immv,ippv
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: jmv,jpv,jmmv,jppv
!    INTEGER, DIMENSION(:),        ALLOCATABLE :: kmv,kpv,kmmv,kppv
!--------------------------------

!    LOGICAL            :: readIblankFlag

    REAL(KIND=CGREAL)  :: xout,yout,zout,Xext,DampFact
    REAL(KIND=CGREAL)  :: uinit,vinit,winit
    REAL(KIND=CGREAL)  :: ux1,ux2,vx1,vx2,wx1,wx2
    REAL(KIND=CGREAL)  :: uy1,uy2,vy1,vy2,wy1,wy2
    REAL(KIND=CGREAL)  :: uz1,uz2,vz1,vz2,wz1,wz2
    REAL(KIND=CGREAL)  :: freq_ux1,freq_vx1,freq_wx1
    REAL(KIND=CGREAL)  :: freq_ux2,freq_vx2,freq_wx2
    REAL(KIND=CGREAL)  :: freq_uy1,freq_vy1,freq_wy1
    REAL(KIND=CGREAL)  :: freq_uy2,freq_vy2,freq_wy2
    REAL(KIND=CGREAL)  :: freq_uz1,freq_vz1,freq_wz1
    REAL(KIND=CGREAL)  :: freq_uz2,freq_vz2,freq_wz2
  
    REAL(KIND=CGREAL)  :: re,dt,reinv,dtinv
    REAL(KIND=CGREAL)  :: restol, restolPoisson, omega, omega_adv
    REAL(KIND=CGREAL)  :: time
    REAL(KIND=CGREAL)  :: bodyInterceptWeight, imagePointWeight, probeLengthNormalized
    REAL(KIND=CGREAL)  :: gcmFlag

    REAL(KIND=CGREAL)  :: areax1,areax2,  &
                          areay1,areay2,  &
                          areaz1,areaz2

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: radiusx,radiusy,radiusz
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: alpha,cosalpha,sinalpha
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: phase,cosphase,sinphase
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: vxcent,vycent,vzcent,  &
                                                   angvx,angvy,angvz,     &
                                                   xcent,ycent,zcent
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: x_rot_cent, y_rot_cent, z_rot_cent
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcentOld,vycentOld,vzcentOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvxOld,angvyOld,angvzOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcent_combinedOld,vycent_combinedOld,vzcent_combinedOld
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combinedOld,angvy_combinedOld,angvz_combinedOld 
       
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: xcent_combined, ycent_combined, zcent_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcent_combined, vycent_combined, vzcent_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combined, angvy_combined, angvz_combined
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: vxcent_combined_init, vycent_combined_init, vzcent_combined_init
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combined_init, angvy_combined_init, angvz_combined_init
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE  :: angvx_combined_old, angvy_combined_old, angvz_combined_old             
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)                                                    
! for new motion -- Added by Haibo
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvx_old,angvy_old,angvz_old
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: xcentinit,ycentinit,zcentinit,        &
                                                   vxcentTrans,vycentTrans,vzcentTrans,  &
                                                   ampx,ampy,ampz,freqx,freqy,freqz
                                                        
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvxinit,angvyinit,angvzinit,  &
                                                   ampangx,ampangy,ampangz,        &
                                                   freqangx,freqangy,freqangz
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
    INTEGER,DIMENSION(:), ALLOCATABLE          :: i_fixed, fixed_mother_body   
    INTEGER,DIMENSION(:), ALLOCATABLE          :: pendulum_mother_body 
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: hinge_x, hinge_y, hinge_z       
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: hinge_vx, hinge_vy, hinge_vz                        
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
! wall velocity valiables
    INTEGER,           DIMENSION(:),     ALLOCATABLE  :: mMinWallVel, mMaxWallVel
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: ampVelX, ampVelY, ampVelZ
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: freqVelX, freqVelY, freqVelZ
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE  :: phaseVelX, phaseVelY, phaseVelZ

    REAL(KIND=CGREAL)  :: area_left,area_right,     &
                          area_bot,area_top,        &
                          area_back,area_front,     &
                          outflow_area, inflow_area, prim_left, outflow_area_FP

    REAL(KIND=CGREAL)  :: alfa       ! Weighting factor for hybrid scheme - added by Rupesh

    REAL(KIND=CGREAL)  :: vPert2Dto3D ! for 2D-3D initial random perturbations
  
! For Flow-Induced Motion. --- Added by veera

    REAL(KIND=CGREAL) ,DIMENSION(:),ALLOCATABLE  :: xcentConstr, ycentConstr, zcentConstr ! Centroid Constraint Flag
    
    REAL(KIND=CGREAL) :: density_fluid 
    
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: density_solid
   
    INTEGER            :: pbcx1,pbcx2, pbcy1,pbcy2, pbcz1,pbcz2  ! Added by H. Luo
    REAL(KIND=CGREAL)  :: pppx1,pppx2, pppy1,pppy2, pppz1,pppz2  !
    INTEGER            :: advec_scheme                           ! added for implicit scheme

    INTEGER            :: ifluxplane, jfluxplane, kfluxplane     ! used for FEA
                                                                 ! SER_TO_PAR. QX. CH4

    INTEGER            :: is_scalar_on !xzheng flag for adding dyes into flow

   END MODULE flow_parameters
!------------------------------------------------------

   MODULE grid_arrays

    USE global_parameters
   
    IMPLICIT NONE

!   Global Grid Variables
!   ---------------------
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: x,y,z
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: xc,yc,zc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dx,dy,dz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxinv,dyinv,dzinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxc,dyc,dzc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxcinv,dycinv,dzcinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: fx,fy,fz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: damper
   
   END MODULE grid_arrays
!------------------------------------------------------

   MODULE flow_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u,v,w
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_ue,face_vn,face_wf
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_uw,face_vs,face_wb
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxu,bcxv,bcxw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcyu,bcyv,bcyw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bczu,bczv,bczw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: viscTot
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxvisc,bcyvisc,bczvisc

    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: uGhost, vGhost, wGhost
    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: pGhost

    REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: div, temp1, temp2, temp3,vor

   END MODULE flow_arrays
!------------------------------------------------------

   MODULE boundary_arrays

    USE global_parameters
    
    IMPLICIT NONE

    INTEGER  :: nFresh, num_fresh, num_dead

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblank, fresh_cell, ghostCellMark, boundCell
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ghostCellSolid, ghostCellMemb
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: bodyNum, bodyNumOld

    INTEGER, DIMENSION(:,:,:),   ALLOCATABLE :: iblank_solid
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iblank_memb

!    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: pot_flag
    REAL(KIND=CGREAL), DIMENSION(:),     ALLOCATABLE :: iBound, jBound, kBound
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iup , ium , jup , jum , kup,  kum, exp_weight
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iupp, iumm, jupp, jumm, kupp, kumm              ! For 2nd Upwinding - Added by Rupesh
!    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarker,yBodyMarker,zBodyMarker, &
!                                                        uBodyMarker,vBodyMarker,wBodyMarker, &
!                                                        axBodyMarker,ayBodyMarker,azBodyMarker, &
!                                                        sBodyMarker,dsBodyMarker,xNormBodyMarker, &
!                                                        yNormBodyMarker
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarker,yBodyMarker,zBodyMarker, &
                                                        uBodyMarker,vBodyMarker,wBodyMarker, &
                                                        sBodyMarker,pBodyMarker!xudong added June 3rd 
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)  
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: theta_group
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_vxcent, non_vycent, non_vzcent  
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: angv_roll,angv_yaw,angv_pitch
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: R11,R12,R13,R21,R22,R23,R31,R32,R33 
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: C11,C12,C13,C21,C22,C23,C31,C32,C33  
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: C_11,C_12,C_13,C_21,C_22,C_23,C_31,C_32,C_33                                                           
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyMarker_rel,vBodyMarker_rel,wBodyMarker_rel
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyMarker_rel_old,vBodyMarker_rel_old,wBodyMarker_rel_old
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)                                                          

    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradx1,pgradx2
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgrady1,pgrady2 
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradz1,pgradz2 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundPresSource

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblankUndecided, iblankTemp

   END MODULE boundary_arrays
!------------------------------------------------------

   MODULE multiuse_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu,nlv,nlw 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uTilde,vTilde,wTilde
   
   END MODULE multiuse_arrays
!------------------------------------------------------

   MODULE pressure_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: p, pPrime
   
   END MODULE pressure_arrays
!------------------------------------------------------

   MODULE nlold_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nluold,nlvold,nlwold
   
   END MODULE nlold_arrays
!------------------------------------------------------

   MODULE solver_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amx,apx,acx
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amy,apy,acy
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amz,apz,acz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: rhs,dummy
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: face1, face2
   
   END MODULE solver_arrays
!------------------------------------------------------

   MODULE solver_ad_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad,apx_ad
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad,apy_ad
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad,apz_ad
   
   END MODULE solver_ad_arrays
!------------------------------------------------------

   MODULE GCM_arrays
   
    USE global_parameters
   
    IMPLICIT NONE
   
    INTEGER                                :: iRowMax, nGhost
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: incI, incJ, incK, iPvt
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarker,            &
                                              iGhost,jGhost,kGhost,     &
                                              iCellIndex,jCellIndex,kCellIndex
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iFresh,jFresh,kFresh      &
                                             ,iFreshCellIndex,jFreshCellIndex,kFreshCellIndex
!    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarkerFresh
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestElementFresh
!    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iBodyRank
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iCellIndexS, jCellIndexS, kCellIndexS

!    REAL(KIND=CGREAL), DIMENSION(2)                :: det  
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: work
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: closestMarkerRatio, &
                                                      xBodyInterceptTang, &
                                                      yBodyInterceptTang, &
                                                      zBodyInterceptTang, &
                                                      xBodyInterceptNorm, &
                                                      yBodyInterceptNorm, &
                                                      zBodyInterceptNorm, &
                                                      xBodyIntercept,     &
                                                      yBodyIntercept,     &
                                                      zBodyIntercept,     &
                                                      xImagePoint,        &
                                                      yImagePoint,        &
                                                      zImagePoint,        &
                                                      probeLength,        &
                                                      uBodyIntercept,     &
                                                      vBodyIntercept,     &
                                                      wBodyIntercept,     &
                                                      pBodyIntercept,     &
                                                      dpdnBodyIntercept,  &
                                                      dpdtBodyIntercept

!    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: xBIG,yBIG,ZBIG

    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: &!!closestMarkerRatioFresh, &
                                                      xBodyInterceptFresh,     &
                                                      yBodyInterceptFresh,     &
                                                      zBodyInterceptFresh,     &
                                                      uBodyInterceptFresh,     &
                                                      vBodyInterceptFresh,     &
                                                      wBodyInterceptFresh

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMD, coeffGCMN,         &
                                                      vanMatrixD, vanMatrixN,       &
                                                      coeffGCMFreshD!,dSFaceProject, &   
!                                                      xBodyCentroid,yBodyCentroid,  &
!                                                      sBodyCentroid,                &
!                                                      xCentroidTang,yCentroidTang, &
!                                                      xCentroidNorm,yCentroidNorm
 
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: probeLengthS,probeLengthNormalizedS, &
                                                      imagePointWeightS,                   &
                                                      xImagePointS,yImagePointS,zImagePointS

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMDS, coeffGCMNS!, &
!                                                      vanMatrixDS, vanMatrixNS

   END MODULE GCM_arrays
!------------------------------------------------------

   MODULE unstructured_surface_arrays

    USE global_parameters
    
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemNormX,      triElemNormY,      triElemNormZ,      triElemArea
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemCentX,      triElemCentY,      triElemCentZ
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang1X,     triElemTang1Y,     triElemTang1Z
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang2X,     triElemTang2Y,     triElemTang2Z
    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: pointOutsideBodyX, pointOutsideBodyY, pointOutsideBodyZ, surfArea

    INTEGER, DIMENSION(:,:,:),         ALLOCATABLE :: triElemNeig
    INTEGER, DIMENSION(:),             ALLOCATABLE :: totNumTriElem
    INTEGER, DIMENSION(:),             ALLOCATABLE :: closestElementGC!,cElementG
    INTEGER, DIMENSION(:,:),           ALLOCATABLE :: edge_node

    REAL (KIND=CGREAL)                             :: normDirFlag

   END MODULE unstructured_surface_arrays
!------------------------------------------------------

   MODULE probe_parameters

    IMPLICIT NONE
    
    INTEGER                           :: nProbe  
    INTEGER, DIMENSION(:),ALLOCATABLE :: iProbe, jProbe, kProbe
    
    LOGICAL, DIMENSION(:),ALLOCATABLE :: myProbe

   END MODULE probe_parameters
!------------------------------------------------------

   MODULE stat_arrays

    USE global_parameters

        IMPLICIT NONE

         INTEGER                                        :: statCtr
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uAv,vAv,wAv,pAv
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uvAv,vwAv,uwAv
         REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uuAv,vvAv,wwAv

   END MODULE stat_arrays
!------------------------------------------------------

    MODULE blasius_profile

        USE global_parameters

        IMPLICIT NONE

!        REAL(KIND=CGREAL)                          :: cavity_H,slot_H,ddratio,d,delta,uinf
!        REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:) :: eta,u_blasius
!        INTEGER                                    :: l,junk,i_start

    END MODULE blasius_profile
!------------------------------------------------------

    MODULE mg_parameters
 
    USE global_parameters
    USE flow_parameters
 
    IMPLICIT NONE
 
    INTEGER :: mgLevels_X, mgLevels_Y, mgLevels_Z
!    INTEGER :: mgcyclex, mgcycley, mgcyclez, infoconv, incrlev
    INTEGER :: mgcyclex, mgcycley, mgcyclez, infoconv
    INTEGER :: iterFinest, iterInter, iterCoarsest
!    INTEGER :: ittt1, nCount
 
    INTEGER :: iRedBlack, TNcolorX, TNcolorY, TNcolorZ, iStep, jStep, kStep    !new for Redblack LSOR
    
    integer :: nCoarseGrids

 
   END MODULE mg_parameters
 
!------------------------------------------------------
   MODULE mg_arrays

    USE global_parameters
 
    IMPLICIT NONE
 
    TYPE :: MGXtype
!   ===============
      INTEGER :: nx_GLBL , nx
      INTEGER :: nxc_GLBL, nxc

      INTEGER :: myIs, myIe

      INTEGER, DIMENSION(:,:,:), POINTER :: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER :: x, xc, dxinv, dxcinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef MPI
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
    END TYPE MGXtype
    
    TYPE :: MGYtype
!   ===============
      INTEGER :: ny_GLBL , ny
      INTEGER :: nyc_GLBL, nyc

      INTEGER :: myJs, myJe

      INTEGER, DIMENSION(:,:,:),POINTER  :: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER :: y, yc, dyinv, dycinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef MPI
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
    END TYPE MGYtype
    
    TYPE :: MGZtype
!   ===============
      INTEGER :: nz, nzc

      INTEGER, DIMENSION(:,:,:), POINTER:: iblank

      REAL(KIND=CGREAL), DIMENSION(:),     POINTER:: z, zc, dzinv, dzcinv
      REAL(KIND=CGREAL), DIMENSION(:,:,:), POINTER :: rhs, res, phi

#ifdef MPI
      INTEGER :: parVecWE , ParVecSN
      INTEGER :: parIVecWE, ParIVecSN
#endif
    END TYPE MGZtype
    
    TYPE(MGXtype), ALLOCATABLE :: MGX(:)
    TYPE(MGYtype), ALLOCATABLE :: MGY(:)
    TYPE(MGZtype), ALLOCATABLE :: MGZ(:)
    
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE              :: iblank_MG
!!    INTEGER, DIMENSION(:,:,:), ALLOCATABLE              :: ghostcellMark_MG, iblank_MG
 
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE    :: ium_mg, iup_mg, &
                                                           jum_mg, jup_mg, &
                                                           kum_mg, kup_mg

!------- new arrays ----------
!   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: iup_mg_Outer, jup_mg_Outer, kup_mg_Outer
!   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: ium_mg_Outer, jum_mg_Outer, kum_mg_Outer
 
   END MODULE mg_arrays
     
!------------------------------------------------------
MODULE usr_module

    USE global_parameters
     
   REAL(KIND=CGREAL) :: density_ratio
   REAL(KIND=CGREAL) :: lScale,vScale
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: I_XX_COMBINED,I_YY_COMBINED,I_ZZ_COMBINED,I_XY_COMBINED,I_YZ_COMBINED,I_XZ_COMBINED
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_volume,volume
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: I_XX,I_YY,I_ZZ,I_XY,I_YZ,I_XZ
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_mass
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: moment_x_combined,moment_y_combined,moment_z_combined
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: force_x_combined,force_y_combined,force_z_combined
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: non_dim_mass_combined
   REAL(KIND=CGREAL),DIMENSION(:,:,:),ALLOCATABLE :: nonDimM_I_combined, invMI_combined
   REAL(KIND=CGREAL) :: vxcent_combined_prev,vycent_combined_prev,vzcent_combined_prev
   REAL(KIND=CGREAL) :: angvx_combined_prev, angvy_combined_prev, angvz_combined_prev   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
   REAL(KIND=CGREAL) :: vxcent_prev,vycent_prev,vzcent_prev
   REAL(KIND=CGREAL) :: angvx_prev, angvy_prev, angvz_prev
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)   
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: force_x,force_y,force_z
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: moment_x,moment_y,moment_z   
   REAL(KIND=CGREAL),DIMENSION(:,:,:),ALLOCATABLE :: invMI
   REAL(KIND=CGREAL),DIMENSION(:,:,:),ALLOCATABLE :: nonDimM_I
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: scx,scy,scz
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: scmx,scmy,scmz
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)   
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: scpw, ssspw, adpw
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: shear_x, shear_y, shear_z
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: pres_x, pres_y, pres_z
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)

END MODULE usr_module 

!------------------------------------------------------
MODULE stat_vort_arrays

    USE global_parameters

    IMPLICIT NONE

    INTEGER                                        :: statCtrv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxAv,oyAv,ozAv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxoxAv,oyoyAv,ozozAv

END MODULE stat_vort_arrays
!------------------------------------------------------

!=====Added for FEA====================================

   MODULE finite_element_parameters  ! SER_TO_PAR. QX. CH5

    USE global_parameters

    IMPLICIT NONE

    INTEGER :: fea_nbody, fea_nnode, fea_nelementmax, fea_nptmax, fea_mbandmax, fea_nmat
    INTEGER :: fea_nfixptmax, fea_nloadmax, fea_neignvmax
    INTEGER :: fea_nndof, fea_nnstress, fea_solvertype, fea_probtype, fea_outputtype, fea_nva
    INTEGER :: fea_nedof, fea_ngdof
    INTEGER :: fea_ntime, fea_itime
    INTEGER :: fea_inputfile, fea_meshfile, fea_outputfile, fea_conloadfile, fea_solutionfile,&
               fea_staticoutput, fea_initdisvel, fea_correstablefile, fea_restartfilein, fea_restartfileout,&
               fea_stressfile, fea_probein, fea_ifuProbeOut, fea_Mfile, fea_keffectfile, fea_afile
    INTEGER :: fea_boundary_flag, fea_contactprobflag
    INTEGER :: fea_nd
    INTEGER :: fea_nfilterstart, fea_nfilterend, fea_dtratio
    INTEGER :: fea_readM, fea_readkeffect, fea_readinita
    INTEGER :: fea_contact_model,fea_trans_dir

    REAL(KIND=CGREAL)  :: fea_gravkey, fea_grav, fea_time
    REAL(KIND=CGREAL)  :: fea_freq, fea_cc1, fea_cc2, fea_dt, fea_beta, fea_gamma
    REAL(KIND=CGREAL)  :: fea_nmc0, fea_nmc1, fea_nmc2, fea_nmc3, fea_nmc4, fea_nmc5, fea_nmc6, fea_nmc7
    REAL(KIND=CGREAL)  :: fea_omegaload
    REAL(KIND=CGREAL)  :: fea_coeff_penalty

   END MODULE  finite_element_parameters

!------------------------------------------------------
   MODULE finite_element_arrays  ! SER_TO_PAR. QX. CH6

    USE global_parameters

    IMPLICIT NONE

    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_nelement, fea_npt, fea_mband,fea_nfixpt, fea_nload, fea_muv, fea_neignv
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: fea_ielement, fea_ifixed, fea_iload
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: markertofeapointtable, featomarkerpointtable


    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_cood, fea_vfixed, fea_vload
    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_gmm, fea_gkm, fea_gcm, fea_keffect
    REAL(KIND=CGREAL),DIMENSION(:, :, :), ALLOCATABLE :: fea_veignv
    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_vmati, fea_ieignv
    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_d, fea_v, fea_a, fea_gp, fea_gp_old, fea_gu

    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_icontactdir, fea_icontactplane, fea_icontactdirnorm
    INTEGER, DIMENSION(:), ALLOCATABLE :: fea_contactflag, fea_ncontactpoint, fea_ncontactsurf, fea_ncontactside
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: fea_icontactpoint, fea_icontactsurf, fea_icontactside

    REAL(KIND=CGREAL),DIMENSION(:, :),    ALLOCATABLE :: fea_original_d, fea_original_v, fea_original_a
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: fea_penaltycoeff
   
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: tempt,keffecttempt
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: tempta
 
   END MODULE finite_element_arrays

!------------------------------------------------------
   MODULE fea_probe_parameters  ! SER_TO_PAR. QX. CH7

    IMPLICIT NONE

    INTEGER                           :: fea_nProbe
    INTEGER, DIMENSION(:),ALLOCATABLE :: fea_iprobebody, fea_iprobenode

   END MODULE fea_probe_parameters

!------------------------------------------------------
   MODULE eigen_motion_parameters  !SER_TO_PAR. QX.CH8

    USE global_parameters

    IMPLICIT NONE

    INTEGER ::  neig_body, maxnmodes, fea_ifeigenin

   ENDMODULE eigen_motion_parameters

!-----------------------------------------------------
   MODULE eigen_motion_arrays  !SER_TO_PAR. QX. CH9

    USE global_parameters

    IMPLICIT NONE

    INTEGER, ALLOCATABLE, DIMENSION(:) :: ieig_body_table, nmodes_eig, ncontact_eig
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: icon_dir_eig, icon_plane_eig, icon_surnor_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:, :, :) :: x_eig, y_eig, z_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:, :) :: fre_eig, coe_eig, phase_eig
    REAL(KIND=CGREAL), ALLOCATABLE, DIMENSION(:,:) :: x_coor_org, y_coor_org, z_coor_org

  ENDMODULE eigen_motion_arrays
 
!========End FEA====================================

   MODULE tahoe_parameters  ! Added by Rajneesh
   
   USE global_parameters

   IMPLICIT NONE

   INTEGER :: ndim_tahoe,elementType,ConstitutiveModel, TahoeSolver, nstart_tahoe,ndump_tahoe,nprobe_tahoe,&
   dtratio_tahoe,nmax,emax,FBCmax,KBCmax,Markmax
   
!  Flags
   INTEGER :: fea_boundary_flag_tahoe,verbosetahoe,tahoe_dataprobe,itime_tahoe,KBC_flag,FBC_flag,&
   iTahoeOutput,itahoedraglift,iTahoeRestart,iTahoeProbe,iflag_restart_tahoe

!  Files strings
   INTEGER :: tahoe_inputfile, tahoe_fbcnodes,tahoe_kbcnodes,tahoe_markers2fem,tahoe_u,tahoe_Du,tahoe_xml,&
   tahoe_bodydat,tahoe_geom,implicit_inputfile,imconverg,time_data1,time_data2,tahoe_DDu,tahoe_elem0,tahoe_rs,&
   restart_tahoebodyin,restart_tahoebodyout

   REAL(KIND=CGREAL)  :: SolidDen,YmodulusEQ,YmodulusNEQ,PoissonR,etaDamping,muNEQ,muEQ,kappaNEQ,kappaEQ,tauBulk,tauShear,dtTahoe
   
   INTEGER, ALLOCATABLE, DIMENSION(:) :: numNodes,numElts,numMarkers,FBCnodes,KBCnodes
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: FBCnode2marker,Markers2femGridpt,nodesFBC,nodesKBC,BodyProbe,marker2elem
   REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:) ::utahoe,Dutahoe,DDutahoe
   REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:) :: x0i,y0i,z0i

   
   
  ENDMODULE tahoe_parameters

!-----------------------------------------------------
MODULE cutcell_arrays
 USE global_parameters
 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: AREA_W, AREA_E, AREA_S, AREA_N, AREA_B, AREA_F
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: VOL_CELL, CUT_AREA, cent_x,cent_y,cent_z
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: isolid

END MODULE  
!-----------------------------------------------------
MODULE implicit_coupling_parameters  ! Added by Rajneesh
 
  USE global_parameters
   IMPLICIT NONE
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)       
   INTEGER :: implicit_coupling_flag_combined_motion
   INTEGER :: fim_inputfile
   INTEGER :: max_ntime_implicit_FIM
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarkerPrev,yBodyMarkerPrev,zBodyMarkerPrev 
   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: C11_prev,C12_prev,C13_prev,C21_prev,C22_prev,C23_prev,C31_prev,C32_prev,C33_prev  
   INTEGER :: kimplicit_FIM,kimplicitMax_FIM,flag_underrelax_FIM,kimplicit_start_FIM,varyalphawithtime_FIM
   REAL(KIND=CGREAL) :: alpha_underrelax_FIM,slope_underrelax_FIM,ImplicitResidual_FIM,ImplicitResidualMax_FIM
   REAL(KIND=CGREAL) :: timeimpl1_FIM,timeimpl2_FIM,alpha_underrelax2_FIM, alphaimplicit_FIM   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)         
   INTEGER :: implicit_coupling_flag,kimplicit,kimplicitMax,flag_underrelax,kimplicit_start,varyalphawithtime
   REAL(KIND=CGREAL) :: alpha_underrelax,slope_underrelax,ImplicitResidual,ImplicitResidualMax
   REAL(KIND=CGREAL) :: timeimpl1,timeimpl2,alpha_underrelax2, alphaimplicit
 
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iblankOld,iblank_solidOld
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)     
   INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: iblank_membOld
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uOld,vOld,wOld,uGhostOld,vGhostOld,wGhostOld,pPrime_old, p_old, pGhostOld
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_ue_Old, face_vn_Old, face_wf_Old, face_uw_Old, face_vs_Old, face_wb_Old 
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu_Old, nlv_Old, nlw_Old, uTilde_Old, vTilde_Old, wTilde_Old 
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: fresh_cell_old 
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iup_old, ium_old, jup_old, jum_old, kup_old, kum_old  
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: ghostCellMark_old, ghostCellMemb_old, ghostCellSolid_old
   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE     :: sCx_old, sCy_old, sCz_old, scmx_old, scmy_old, scmz_old
   REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE     :: xcent_old, ycent_old, zcent_old
   REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE   :: triElemCentxOld, triElemCentyOld, triElemCentzOld, triElemAreaOld      
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundPresSource_old 
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxu_Old, bcxv_Old, bcxw_Old, bcyu_Old, bcyv_Old, bcyw_Old, bczu_Old, bczv_Old, bczw_Old 
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bodyNum_Old
   INTEGER :: nFresh_Old
   REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradx1_old, pgradx2_old, pgrady1_old, pgrady2_old, pgradz1_old, pgradz2_old
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nluold_Old, nlvold_Old, nlwold_Old
   REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: iupp_Old,iumm_Old,jupp_Old,jumm_Old,kupp_Old,kumm_Old
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)    
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarkerOld,yBodyMarkerOld,zBodyMarkerOld
   REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: uBodyMarkerOld,vBodyMarkerOld,wBodyMarkerOld

ENDMODULE implicit_coupling_parameters
   
!-----------------------------------------------------
MODULE diffusive_material_transport
 
 USE global_parameters
 IMPLICIT NONE
 
 INTEGER :: nscalars
 
 INTEGER, DIMENSION(:), ALLOCATABLE :: T_inner_BCtype
 INTEGER, DIMENSION(:), ALLOCATABLE :: bcx1_scalar, bcy1_scalar, bcz1_scalar   !BOUNDARY CONDITION TYPE
 INTEGER, DIMENSION(:), ALLOCATABLE :: bcx2_scalar, bcy2_scalar, bcz2_scalar
 INTEGER, DIMENSION(:), ALLOCATABLE :: T_init_type

 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Tbcinner
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: sc
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: T_init_val
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: Tbcx1, Tbcx2, Tbcy1, Tbcy2, Tbcz1, Tbcz2

 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: T_scalar, nlT, nlTold, bcxT, bcyT, bczT
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: TGhost

 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amx_ad_T,apx_ad_T
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amy_ad_T,apy_ad_T
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: amz_ad_T,apz_ad_T
 

ENDMODULE diffusive_material_transport
   
   
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START) 
MODULE derivative

 USE global_parameters
 IMPLICIT NONE
 
INTEGER :: derivative_flag
INTEGER :: non_inertia_coord
INTEGER :: k_derivative
REAL(KIND=CGREAL) :: deltaU
REAL(KIND=CGREAL) :: deltaV
REAL(KIND=CGREAL) :: deltaQ
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_x_ref     
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_y_ref   
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_moment_z_ref
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_x_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_y_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_moment_z_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_x_deltaV 
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_y_deltaV
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_moment_z_deltaV
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_x_deltaQ
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_force_y_deltaQ
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_moment_z_deltaQ

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_x_ref     
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_y_ref   
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_moment_z_ref
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_x_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_y_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_moment_z_deltaU
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_x_deltaV 
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_y_deltaV
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_moment_z_deltaV
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_x_deltaQ
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_force_y_deltaQ
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: inertia_moment_z_deltaQ

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC11
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC12
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC13
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC21
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC22
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC23
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC31
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC32
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: invC33

REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_U
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_V
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_W
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_angvx
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_angvy
REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: non_inertia_angvz

ENDMODULE derivative

 
 
MODULE scalar

 USE global_parameters
 IMPLICIT NONE
 
 
 INTEGER, PARAMETER :: XBC_DIRICHLET            = 1, &
                          XBC_NEUMANN              = 2 
                          
 INTEGER, PARAMETER :: PhiBC_DIRICHLET            = 1, &
                          PhiBC_NEUMANN              = 2                           
                          
                          
 INTEGER             :: VDVactive, PFactive
 INTEGER             :: n_visual, n_potential
 INTEGER             :: idirection
 INTEGER             :: iterMaxLaplace
 INTEGER             :: iterMaxPotential
 INTEGER             :: POSSION_INDEX
 INTEGER,PARAMETER :: index_pressure   = 1
 INTEGER,PARAMETER :: index_scalar     = 2
 INTEGER,PARAMETER :: index_potential  = 3
 REAL(KIND=CGREAL)   :: restolLaplace
 REAL(KIND=CGREAL)   :: restolPotential
 INTEGER            :: Xbcx1, Xbcx2, Xbcy1, Xbcy2, Xbcz1, Xbcz2 
 INTEGER            :: Phibcx1, Phibcx2, Phibcy1, Phibcy2, Phibcz1, Phibcz2 
 REAL(KIND=CGREAL)   :: xxxx1,xxxx2, xxxy1,xxxy2, xxxz1,xxxz2
 REAL(KIND=CGREAL)   :: phix1,phix2, phiy1,phiy2, phiz1,phiz2
 REAL(KIND=CGREAL)   :: fx1, fx2, fy1, fy2, fz1, fz2 
 REAL(KIND=CGREAL)   :: botsurf_kinetic_energy, topsurf_kinetic_energy
 REAL(KIND=CGREAL)   :: kinetic_energy_difference
 REAL(KIND=CGREAL)   :: botsurf_pressure, topsurf_pressure
 REAL(KIND=CGREAL)   :: pressure_difference
 REAL(KIND=CGREAL)   :: outer_shearforce
 REAL(KIND=CGREAL)   :: Grad_u_square_dot_Grad_Dx_total
 REAL(KIND=CGREAL)   :: outer_Dx_dot_advection, outer_Dx_dot_gradu2
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: scalarX
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: potential
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Xgradx, Xgrady, Xgradz
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Potential_u, Potential_v, Potential_w
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Vorticity_u, Vorticity_v, Vorticity_w
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgradx1,Xgradx2
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgrady1,Xgrady2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Xgradz1,Xgradz2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigradx1,Phigradx2
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigrady1,Phigrady2 
 REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: Phigradz1,Phigradz2 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundScalarSource
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: boundPhiSource
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: zero
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: XGhost
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: PhiGhost 
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: omega_cross_u
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: XgradxBodyMarker, XgradyBodyMarker, XgradzBodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: X_BodyMarker, X_BodyMarker_prev
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: advection_x_BodyMarker, advection_y_BodyMarker, advection_z_BodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Lamb_x_BodyMarker, Lamb_y_BodyMarker, Lamb_z_BodyMarker
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_term1, unsteady_term2, tau_Xgrad 
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady1, unsteady1_prev
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_sum1, unsteady_sum2, unsteady_sum2p, unsteady_sum2v, advection_sum, Lamb_sum
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: kinematics_term2
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: X_Marker 
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: ssBodyMarker
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: dssBodyMarker 
 !REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: xforce, yforce, zforce
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: x_force, y_force, z_force
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: pTriElemCent0,pTriElemCent1
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: vorTriElemCent0,vorTriElemCent1
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: vorNormCent0,vorNormCent1
 !REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: p_force0,p_force1
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: u_n, u_n_prev, u_n1_prev, u_n2_prev, u_n_dot, undot, udotn
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: advection, advection_prev, AB2_advection
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: F_i
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Advection_force, Unsteady_force, Shear_force, Vortex_force
 REAL(KIND=CGREAL) :: F_i_total, Advection_force_total, Unsteady_force_total, Shear_force_total, Vortex_force_total, Div_AB2_Xx_force_total, Div_Gradu2_Xx_force_total
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: diff_prev, diff_star, gradP, gradu2
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: PuPt, PvPt, PwPt
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_prev, vBodyCent_prev, wBodyCent_prev  
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_p_prev, vBodyCent_p_prev, wBodyCent_p_prev
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: uBodyCent_v_prev, vBodyCent_v_prev, wBodyCent_v_prev  
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: adpwl, sspwl, ppwl, totalpwl
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: adpw_l, sspw_l, ppw_l, totalpw_l
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: u_dot_n, up_dot_n, uv_dot_n
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: advection_n, Lamb_normal
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Div_Lamb
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u_square, u_squareGhost
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: potential_u_square, vorticity_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: Grad_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: Grad_potential_u_square, Grad_vorticity_u_square
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_u_square_x, Grad_u_square_y, Grad_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_potential_u_square_x, Grad_potential_u_square_y, Grad_potential_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_vorticity_u_square_x, Grad_vorticity_u_square_y, Grad_vorticity_u_square_z
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_u_square_n
 REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: Grad_potential_u_square_n, Grad_vorticity_u_square_n
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: Grad_u_square_sum, Grad_potential_u_square_sum, Grad_vorticity_u_square_sum 
 REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: AB2_advection_Xx, Gradu2_Xx
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Div_AB2_Xx, Div_Gradu2_Xx
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Div_AB2_Xx_force, Div_Gradu2_Xx_force
 REAL(KIND=CGREAL), DIMENSION(:,:,:),   ALLOCATABLE :: Grad_u_square_dot_Grad_Dx
!dddddddddddddddddddddddddddddddddddddddd
 REAL(KIND=CGREAL) :: outer_unsteady
 REAL(KIND=CGREAL) :: outer_shearforce_down
 REAL(KIND=CGREAL) :: outer_shearforce_up
 REAL(KIND=CGREAL) :: outer_shearforce_front
 REAL(KIND=CGREAL) :: outer_shearforce_back
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: unsteady_sum_check 
 REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: angle
!dddddddddddddddddddddddddddddddddddddddd

ENDMODULE scalar
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 
