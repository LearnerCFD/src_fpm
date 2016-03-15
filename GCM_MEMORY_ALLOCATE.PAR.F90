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
!  Filename: GCM_MEMORY_ALLOCATE.PAR.F90
!  Latest Modification: October 21, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_Allocate_static_arrays()
!     GCM_AllocateGhostCellArrays()
!     GCM_AllocateFreshCellArrays()
! --------------------------------------------------------------------



SUBROUTINE GCM_Allocate_static_arrays()

! ----------------------------------------------------------
!  static arrays for GCM that need to be declared only once
! ----------------------------------------------------------

    USE global_parameters
    USE GCM_arrays

    IMPLICIT NONE

    INTEGER :: iErr



    iRowMax = 8

    ALLOCATE(incI(iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for incI'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(incJ(iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for incJ'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(incK(iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for incK'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(iPvt(iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for iPvt'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(work(iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for work'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vanMatrixD(iRowMax,iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for vanMatrixD'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vanMatrixN(iRowMax,iRowMax),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_Allocate_static_arrays: Memory Allocation Error for vanMatrixN'
      STOP
    ENDIF ! ierr

!=======================================================================
! These allow us to define stencil image point
! Assumed clockwise from lower left corner.
!=======================================================================

    incI(:) = (/0, 0, 1, 1, 0, 0, 1, 1/)
    incJ(:) = (/0, 1, 1, 0, 0, 1, 1, 0/)
    incK(:) = (/0, 0, 0, 0, 1, 1, 1, 1/)

END SUBROUTINE GCM_Allocate_static_arrays
!---------------------------------------------------------------------



SUBROUTINE GCM_AllocateGhostCellArrays()

!   ----------------------------------------------------------------------------------------
!    This subroutine allocates memory for various arrays pertinent to GCM ghost cell arrays
!   ----------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE GCM_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

    INTEGER :: iErr, nFaceMax



    iRowMax  = 8
    nFaceMax = 6


   ALLOCATE(iGhost(1:nGhost),STAT=iErr)
   IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &   
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for iGhost'
       STOP
   ENDIF ! ierr 

    ALLOCATE(jGhost(1:nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for jGhost'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(kGhost(1:nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for kGhost'
      STOP
    ENDIF ! ierr 

    ALLOCATE(coeffGCMD(iRowMax,nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for coeffGCMD'
      STOP
    ENDIF ! ierr 

    ALLOCATE(coeffGCMN(iRowMax,nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for coeffGCMN'
      STOP
    ENDIF ! ierr 

    
    ALLOCATE(xBodyInterceptTang(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for xBodyInterceptTang'
      STOP
    ENDIF ! ierr 

    ALLOCATE(yBodyInterceptTang(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for yBodyInterceptTang'
      STOP
    ENDIF ! ierr 

    ALLOCATE(zBodyInterceptTang(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for zBodyInterceptTang'
      STOP
    ENDIF ! ierr 

    ALLOCATE(xBodyInterceptNorm(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for xBodyInterceptNorm'
      STOP
    ENDIF ! ierr 

    ALLOCATE(yBodyInterceptNorm(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for yBodyInterceptNorm'
      STOP
    ENDIF ! ierr 

    ALLOCATE(zBodyInterceptNorm(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for zBodyInterceptNorm'
      STOP
    ENDIF ! ierr 
  
    ALLOCATE(xBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for xBodyIntercept'
      STOP
    ENDIF ! ierr 

    ALLOCATE(yBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for yBodyIntercept'
      STOP
    ENDIF ! ierr 

    ALLOCATE(zBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
      'GCM_AllocateGhostCellArrays: Memory Allocation Error for xBodyIntercept'
      STOP
    ENDIF ! ierr 

    ALLOCATE(uBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for uBodyIntercept'
      STOP
    ENDIF ! ierr 

    ALLOCATE(vBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for vBodyIntercept'
      STOP
    ENDIF ! ierr 
 
    ALLOCATE(wBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for wBodyIntercept'
      STOP
    ENDIF ! ierr 

    ALLOCATE(pBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for pBodyIntercept'
      STOP
    ENDIF ! ierr 
 
    ALLOCATE(dpdnBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for dpdnBodyIntercept'
      STOP
    ENDIF ! ierr 
 
    ALLOCATE(dpdtBodyIntercept(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for dpdtBodyIntercept'
      STOP
    ENDIF ! ierr 
                        
    ALLOCATE(closestMarker(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for closestMarker'
      STOP
    ENDIF ! ierr 

    ALLOCATE(closestMarkerRatio(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for closestMarkerRatio'
      STOP
    ENDIF ! ierr 
 
    ALLOCATE(closestElementGC(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for closestElementGC'
      STOP
    ENDIF ! ierr 

    ALLOCATE(iCellIndex(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for iCellIndex'
      STOP
    ENDIF ! ierr 

    ALLOCATE(jCellIndex(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for jCellIndex'
      STOP
    ENDIF ! ierr 

    ALLOCATE(kCellIndex(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for kCellIndex'
      STOP
    ENDIF ! ierr

    ALLOCATE(xImagePoint(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for xImagePoint'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yImagePoint(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for yImagePoint'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(zImagePoint(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for zImagePoint'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(probeLength(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for probeLength'
      STOP
    ENDIF ! ierr 

!   Allocate infrastructure for shear stress
!   ----------------------------------------
    ALLOCATE(coeffGCMDS(iRowMax,nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for coeffGCMDS'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(coeffGCMNS(iRowMax,nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for coeffGCMNS'
      STOP
    ENDIF ! ierr

    ALLOCATE(iCellIndexS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for iCellIndexS'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(jCellIndexS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for jCellIndexS'
      STOP
    ENDIF ! ierr 
       
    ALLOCATE(kCellIndexS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for kCellIndexS'
      STOP
    ENDIF ! ierr    

    ALLOCATE(xImagePointS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for xImagePointS'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yImagePointS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for yImagePointS'
      STOP
    ENDIF ! ierr
     
    ALLOCATE(zImagePointS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for zImagePointS'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(probeLengthS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for probeLengthS'
      STOP
    ENDIF ! ierr 
    
    ALLOCATE(probeLengthNormalizedS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for probeLengthNormalizedS'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(imagePointWeightS(nGhost),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateGhostCellArrays: Memory Allocation Error for imagePointWeightS'
      STOP
    ENDIF ! ierr

!   Initialize arrays
!   -----------------
    iGhost             = 0
    jGhost             = 0 
    kGhost             = 0

    coeffGCMD          = 0.0_CGREAL
    coeffGCMN          = 0.0_CGREAL

    xBodyInterceptTang = 0.0_CGREAL
    yBodyInterceptTang = 0.0_CGREAL
    zBodyInterceptTang = 0.0_CGREAL

    xBodyInterceptNorm = 0.0_CGREAL
    yBodyInterceptNorm = 0.0_CGREAL
    zBodyInterceptNorm = 0.0_CGREAL
      
    xBodyIntercept     = 0.0_CGREAL
    yBodyIntercept     = 0.0_CGREAL 
    zBodyIntercept     = 0.0_CGREAL 
      
    uBodyIntercept     = 0.0_CGREAL
    vBodyIntercept     = 0.0_CGREAL
    wBodyIntercept     = 0.0_CGREAL 

    pBodyIntercept     = 0.0_CGREAL 
    dpdnBodyIntercept  = 0.0_CGREAL 
    dpdtBodyIntercept  = 0.0_CGREAL 
        
    xImagePoint        = 0.0_CGREAL
    yImagePoint        = 0.0_CGREAL
    zImagePoint        = 0.0_CGREAL
    probeLength        = 0.0_CGREAL

    closestMarker      = 0
    closestMarkerRatio = 0.0_CGREAL
    closestElementGC   = 0.0_CGREAL

    coeffGCMDS         = 0.0_CGREAL
    coeffGCMNS         = 0.0_CGREAL

    xImagePointS       = 0.0_CGREAL
    yImagePointS       = 0.0_CGREAL 
    zImagePointS       = 0.0_CGREAL

    probeLengthS           = 0.0_CGREAL
    probeLengthNormalizedS = 0.0_CGREAL
    imagePointWeightS      = 0.0_CGREAL
    
END SUBROUTINE GCM_AllocateGhostCellArrays
!---------------------------------------------------------------------



SUBROUTINE GCM_AllocateFreshCellArrays()

    USE global_parameters
!!    USE flow_parameters
!!    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE
    
    INTEGER :: iErr

!=======================================================================
! Allocate Memory for various arrays pertinent to GCM Fresh Cells
!=======================================================================

    ALLOCATE(iFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for iFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(jFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for jFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(kFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for kFresh'
      STOP
    ENDIF ! ierr
    
!&&&&&&&&&&&&&& HAS NEVER BEEN USED (SAMK) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    ALLOCATE(closestMarkerFresh(nFresh),STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!      WRITE(STDOUT,*) &
!       'GCM_AllocateFreshCellArrays: Memory Allocation Error for closestMarkerFresh'
!      STOP
!    ENDIF ! ierr

    ALLOCATE(closestElementFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for closestElementFresh'
      STOP
    ENDIF ! ierr
    
!&&&&&&&&&&&&&& HAS NEVER BEEN USED (SAMK) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    ALLOCATE(closestMarkerRatioFresh(nFresh),STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!      WRITE(STDOUT,*) &
!       'GCM_AllocateFreshCellArrays: Memory Allocation Error for closestMarkerRatioFresh'
!      STOP
!    ENDIF ! ierr
    
    ALLOCATE(xBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for xBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(yBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for yBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(zBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for zBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(iFreshCellIndex(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for iFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(jFreshCellIndex(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for jFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(kFreshCellIndex(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for kFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(coeffGCMFreshD(iRowMax,nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for coeffGCMFreshD'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(uBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for uBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(vBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for vBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    ALLOCATE(wBodyInterceptFresh(nFresh),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_AllocateFreshCellArrays: Memory Allocation Error for wBodyInterceptFresh'
      STOP
    ENDIF ! ierr

!=======================================================================
!   Initialize arrays
!=======================================================================

    iFresh                   = 0
    jFresh                   = 0
    kFresh                   = 0
!    closestMarkerFresh       = 0
    closestElementFresh      = 0
!    closestMarkerRatioFresh  = 0.0_CGREAL
    xBodyInterceptFresh      = 0.0_CGREAL
    yBodyInterceptFresh      = 0.0_CGREAL
    zBodyInterceptFresh      = 0.0_CGREAL
    iFreshCellIndex          = 0
    jFreshCellIndex          = 0
    kFreshCellIndex          = 0
    coeffGCMFreshD           = 0.0_CGREAL
    uBodyInterceptFresh      = 0.0_CGREAL
    vBodyInterceptFresh      = 0.0_CGREAL
    wBodyInterceptFresh      = 0.0_CGREAL

END SUBROUTINE GCM_AllocateFreshCellArrays
!------------------------------------------------------------------------
