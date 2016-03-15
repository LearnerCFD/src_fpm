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
!  Filename: GCM_MEMORY_DEALLOCATE.PAR.F90
!  Latest Modification: October 21, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_DeallocateGhostCellArrays()
!     GCM_DeallocateFreshCellArrays()
! --------------------------------------------------------------------



SUBROUTINE GCM_DeallocateGhostCellArrays()

! -------------------------------------------------------------------------------------
!  This subroutine deallocates memory for various arrays pertinent to GCM Body Markers
! -------------------------------------------------------------------------------------

    USE global_parameters
    USE GCM_arrays
    USE unstructured_surface_arrays
    
    IMPLICIT NONE

    INTEGER :: iErr



    DEALLOCATE(iGhost,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for iGhost'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(jGhost,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for jGhost'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(kGhost,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for kGhost'
      CALL flow_stop
      STOP
    ENDIF ! ierr 

    DEALLOCATE(coeffGCMD,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for coeffGCMD'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(coeffGCMN,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for coeffGCMN'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(xBodyInterceptTang,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for xBodyInterceptTang'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(yBodyInterceptTang,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for yBodyInterceptTang'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(zBodyInterceptTang,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for zBodyInterceptTang'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(xBodyInterceptNorm,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for xBodyInterceptNorm'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(yBodyInterceptNorm,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for yBodyInterceptNorm'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(zBodyInterceptNorm,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for zBodyInterceptNorm'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(xBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for xBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(yBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for yBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(zBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for zBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(uBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for uBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(vBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for vBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(wBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for wBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr
     
    DEALLOCATE(pBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for pBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr
     
    DEALLOCATE(dpdnBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for dpdnBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(dpdtBodyIntercept,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for dpdtBodyIntercept'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
                        
    DEALLOCATE(closestMarker,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for closestMarker'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(closestMarkerRatio,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for closestMarkerRatio'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(closestElementGC,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for closestElementGC'
      CALL flow_stop
      STOP
    ENDIF ! ierr
          
    DEALLOCATE(iCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for iCellIndex'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(jCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for jCellIndex'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(kCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for kCellIndex'
      CALL flow_stop
      STOP
    ENDIF ! ierr
 
    DEALLOCATE(xImagePoint,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for xImagePoint'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(yImagePoint,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for yImagePoint'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(zImagePoint,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for zImagePoint'
      CALL flow_stop
      STOP
    ENDIF ! ierr 
    
    DEALLOCATE(probeLength,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for probeLength'
      CALL flow_stop
      STOP
    ENDIF ! ierr

!   Deallocate infrastructure for shear stress
!   ------------------------------------------
    DEALLOCATE(coeffGCMDS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for coeffGCMDS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(coeffGCMNS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for coeffGCMNS'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(iCellIndexS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for iCellIndexS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(jCellIndexS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for jCellIndexS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
        
    DEALLOCATE(kCellIndexS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for kCellIndexS'
      CALL flow_stop
      STOP
    ENDIF ! ierr    

    DEALLOCATE(xImagePointS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for xImagePointS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
     
    DEALLOCATE(yImagePointS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for yImagePointS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
 
    DEALLOCATE(zImagePointS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for zImagePointS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
 
    DEALLOCATE(probeLengthS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for probeLengthS'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(probeLengthNormalizedS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for probeLengthNormalizedS'
      CALL flow_stop
      STOP
    ENDIF ! ierr

    DEALLOCATE(imagePointWeightS,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateGhostCellArrays: Memory Deallocation Error for imagePointWeightS'
      CALL flow_stop
      STOP
    ENDIF ! ierr
  
END SUBROUTINE GCM_DeallocateGhostCellArrays
!---------------------------------------------------------------------



SUBROUTINE GCM_DeallocateFreshCellArrays()

    USE global_parameters
!!    USE flow_parameters
!!    USE grid_arrays
!!    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE
    
    INTEGER :: iErr
    
!=======================================================================
! Deallocate Memory for various arrays pertinent to GCM Fresh Cells
!=======================================================================

    DEALLOCATE(iFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for iFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(jFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for jFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(kFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for kFresh'
      STOP
    ENDIF ! ierr
    
!    DEALLOCATE(closestMarkerFresh,STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!      WRITE(STDOUT,*) &
!       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for closestMarkerFresh'
!      STOP
!    ENDIF ! ierr
    
    DEALLOCATE(closestElementFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for closestElementFresh'
      STOP
    ENDIF ! ierr
    
!    DEALLOCATE(closestMarkerRatioFresh,STAT=iErr)
!    IF ( iErr /= ERR_NONE ) THEN
!      WRITE(STDOUT,*) &
!       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for closestMarkerRatioFresh'
!      STOP
!    ENDIF ! ierr
    
    DEALLOCATE(xBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for xBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(yBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for yBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(zBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for zBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(iFreshCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for iFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(jFreshCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for jFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(kFreshCellIndex,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for kFreshCellIndex'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(coeffGCMFreshD,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for coeffGCMFreshD'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(uBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for uBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(vBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for vBodyInterceptFresh'
      STOP
    ENDIF ! ierr
    
    DEALLOCATE(wBodyInterceptFresh,STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'GCM_DeallocateFreshCellArrays: Memory Deallocation Error for wBodyInterceptFresh'
      STOP
    ENDIF ! ierr

END SUBROUTINE GCM_DeallocateFreshCellArrays
!------------------------------------------------------------------------
