! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  Johns Hopkins University
!
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
!  Filename: GCM_INTERCEPT_KERNELS_UNSTRUC.F90
!  Last Modification: June ,17 2010 (ver. 1.5.8)
!  By X.Zheng
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     GCM_calc_bodyIntercept_Unstruc()
!     determine_BIValidityMembrane(iBody,n1,n2,valid)
!     GCM_calc_BIVelocity_Unstruc(iGBI,jGBI,kGBI,xGBI,yGBI,zGBI,closestElementGBI,uGBI,vGBI,wGBI)
! --------------------------------------------------------------------



SUBROUTINE GCM_calc_bodyIntercept_Unstruc(iGP, jGP, kGP, xGP, yGP, zGP, xBI, yBI, zBI , closestElementGP)

! ---------------------------------------------------------------------
!  This subroutine computes the coordinates of the intercept points for
!  any generic point onto the unstructured surface mesh.
! Description: none.
! ---------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE gcm_arrays
    USE grid_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER,           INTENT(IN)  :: iGP, jGP, kGP
    REAL(KIND=CGREAL), INTENT(IN)  :: xGP, yGP, zGP
    INTEGER,           INTENT(OUT) :: closestElementGP
    REAL(KIND=CGREAL), INTENT(OUT) :: xBI, yBI, zBI
    
!... loop variables

    INTEGER :: iEdge,m,n,nc

!... local variables
 
    INTEGER, PARAMETER       :: NSIZE = 1000
    INTEGER, PARAMETER       :: MSIZE = 20
    INTEGER, DIMENSION(:),ALLOCATABLE   :: NeighElemInd

    REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE   :: distMarker

    INTEGER                  :: nCheck
    INTEGER,DIMENSION(1)     :: iDummy(1)
    INTEGER                  :: iBody,numNeighElement
    INTEGER                  :: elemInd,node1,node2,node3,nMarker,iErr,nodeV
    INTEGER                  :: nEdges,nodeCV,nodeA,nodeB
    INTEGER                  :: cElementGP(1:MSIZE),closestNodeGP(1:MSIZE)
    INTEGER                  :: cEdgeNode1(1:MSIZE),cEdgeNode2(1:MSIZE)
    INTEGER                  :: shortestProbe
    INTEGER                  :: edgeNode1,edgeNode2
    INTEGER                  :: closestEdgeNode1,closestEdgeNode2

    REAL(KIND=CGREAL)        :: cNormal,dMin,dsIntercept,xM,yM,zM
    REAL(KIND=CGREAL)        :: area123,areaDiff,distBIElem,distBIElemMin
    REAL(KIND=CGREAL)        :: epsiArea,distGPEI,distGPEIMin
    REAL(KIND=CGREAL)        :: xBITemp, yBITemp, zBITemp
    REAL(KIND=CGREAL)        :: xCV, yCV, zCV
    REAL(KIND=CGREAL)        :: xEI, yEI, zEI
    REAL(KIND=CGREAL)        :: xEITemp, yEITemp, zEITemp
    REAL(KIND=CGREAL)        :: vec01x, vec01y, vec01z
    REAL(KIND=CGREAL)        :: vec12x, vec12y, vec12z
    REAL(KIND=CGREAL)        :: magnitude12,magnitude12Inv,projectedLength
    REAL(KIND=CGREAL)        :: xBIT(1:MSIZE),yBIT(1:MSIZE),zBIT(1:MSIZE),dist(1:MSIZE)

    LOGICAL                  :: validIntercept



    iBody   = bodyNum(iGP,jGP,kGP)

	IF(iBody == 0) THEN
	  iBody = max(bodyNum(iGP+1,jGP,kGP),bodyNum(iGP-1,jGP,kGP),bodyNum(iGP,jGP+1,kGP),bodyNum(iGP,jGP-1,kGP))
	  IF( ndim == 3 ) iBody = max(iBody,bodyNum(iGP,jGP,kGP+1),bodyNum(iGP,jGP,kGP-1))
	ENDIF	
	
    IF (iBody == 0 ) THEN
      PRINT*,'bodyNum = 0 for',iGP,jGP,kGP,'in GCM_INTERCEPT_KERNELS_UNSTRUC' 
      STOP
    ENDIF
    nMarker = nPtsBodyMarker(iBody)


!   NCheck:  Number of closesest nodes to check
!   high values of this variable increases robustness of procedure 
!   and also CPU time for finding body intercept.
!   --------------------------------------------------------------
    nCheck = 3 

    IF (nCheck > MSIZE) THEN
       PRINT*,'nCheck in GCM_calc_bodyIntercept_Unstruc is limited to', MSIZE
       PRINT*,'Increase array size'
       STOP
    ENDIF

!   Allocate local array
!   --------------------
    ALLOCATE(distMarker(nMarker),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) 'search_vertex_dotNorm: Memory Allocation Error for distMarker'
      STOP
    ENDIF ! ierr 
    ALLOCATE(NeighElemInd(NSIZE),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) 'search_vertex_dotNorm: Memory Allocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr 

!   Get closestMarker for generic point 
!   -----------------------------------
    dMin = 1.0E+5_CGREAL

    DO m = 1, nMarker
      xM = xBodyMarker(iBody,m)
      yM = yBodyMarker(iBody,m)
      zM = zBodyMarker(iBody,m)
          
      distMarker(m) = (xM-xGP)**2 + (yM-yGP)**2 + (zM-zGP)**2

    ENDDO ! m 

    DO nc = 1,NCheck
      iDummy                        = MINLOC(distmarker(1:nMarker))
      closestNodeGP(nc)             = iDummy(1)
      distmarker(closestNodeGP(nc)) = 1.0E20_CGREAL
    ENDDO 

!   Find elements that share closest node/marker
!   --------------------------------------------
    DO nc = 1,nCheck
      numNeighElement = 0
      DO m=1,totNumTriElem(iBody)
        IF ( triElemNeig(iBody,1,m) == closestNodeGP(nc) .OR. &
            triElemNeig(iBody,2,m) == closestNodeGP(nc) .OR. &
            triElemNeig(iBody,3,m) == closestNodeGP(nc) ) THEN
          numNeighElement               = numNeighElement + 1
          NeighElemInd(numNeighElement) = m
        ENDIF
      ENDDO

!     Trap error if array NeighElemenInd overflows
!     --------------------------------------------
      IF ( numNeighElement > NSIZE ) THEN
        WRITE(STDOUT,*) 'GCM_calc_bodyIntercept_Unstruc: Memory Overflow Error for NeighElemInd'
        WRITE(STDOUT,*) ' Allocated size = ',NSIZE
        WRITE(STDOUT,*) ' Current size   = ',numNeighElement
        WRITE(STDOUT,*) ' Aborting Run'

        STOP
       ENDIF ! NeighElemInd

!     Determine which element contains normal intercept
!     -------------------------------------------------
      closestElementGP = 0

      epsiArea = 1.0E-04_CGREAL
      distBIElemMin = 1.0E16_CGREAL

      DO n = 1,numNeighElement

        elemInd = NeighElemInd(n)

        node1   = triElemNeig(iBody,1,elemInd)
        node2   = triElemNeig(iBody,2,elemInd)
        node3   = triElemNeig(iBody,3,elemInd)
      
!       Check if BITemp is located inside triangle of surface element
!         through area difference
!       -------------------------------------------------------------
        CALL check_BIInsideTriangle(iBody,elemInd,node1,node2,node3,xGP,yGP,zGP,&
                                       xBITemp,yBITemp,zBITemp,area123,areaDiff)
    
!       Select closest Elem and BI coordinates:
!       If BI falls inside the element use that
!       Else Base the selection on the minimum distance
!       between BI and either the norm to closest side or vertices of side
!       ------------------------------------------------------------------
        IF ( ABS(areaDiff) < epsiArea*area123) THEN
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp
          closestElementGP = elemInd
          closestEdgeNode1 =  0
          closestEdgeNode2 =  0
          GOTO 999
        ELSE 
          CALL calc_BIOutsideTriangle(iBody,elemInd,node1,node2,node3,  &
                                   xBITemp,yBITemp,zBITemp,distBIElem)
       
          IF (distBIElem <= distBIElemMin) THEN
            distBIElemMin    = distBIElem
            closestElementGP = elemInd
            xBI              = xBITemp
            yBI              = yBITemp
            zBI              = zBITemp
          ENDIF ! distBIElem
        ENDIF ! areaDiff

!DEBUG
!     IF (iGP == 102 .and. jGP == 93 .and. kGP == 1) then
!        WRITE(355,*)'ZONE'
!        WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!        WRITE(355,*)xBodyMarker(iBody,node2),yBodyMarker(iBody,node2),zBodyMarker(iBody,node2)
!        WRITE(355,*)xBodyMarker(iBody,node3),yBodyMarker(iBody,node3),zBodyMarker(iBody,node3)
!        WRITE(355,*)xBodyMarker(iBody,node1),yBodyMarker(iBody,node1),zBodyMarker(iBody,node1)
!        WRITE(355,*)'ZONE'
!        WRITE(355,*)xGP,yGP,zGP
!        WRITE(355,*)xBItemp,yBItemp,zBItemp
!        WRITE(356,*)areadiff,closestElementGP,epsiArea*area123
!     ENDIF
!DEBUG

      ENDDO ! n

!     Compute coordinates of Body Intercept in a robust manner
!     for the case where the temporary BI is located outside 
!     all the surface elements
!     --------------------------------------------------------

!     1. Load coordinates of closest vertex (CV)
!     ------------------------------------------
      xCV = xBodyMarker(iBody,closestNodeGP(nc))
      yCV = yBodyMarker(iBody,closestNodeGP(nc))
      zCV = zBodyMarker(iBody,closestNodeGP(nc))

      distGPEIMin = 1.0E+16_CGREAL

!     2. Detemine the indices of the 2 vertices connected to CV
!     ---------------------------------------------------------
      node1   = triElemNeig(iBody,1,closestElementGP)
      node2   = triElemNeig(iBody,2,closestElementGP)
      node3   = triElemNeig(iBody,3,closestElementGP)

      IF ( node1 == closestNodeGP(nc) ) THEN
        nodeCV = node1
        nodeA  = node2
        nodeB  = node3
      ELSEIF ( node2 == closestNodeGP(nc) ) THEN
        nodeCV = node2
        nodeA  = node3
        nodeB  = node1
      ELSEIF ( node3 == closestNodeGP(nc) ) THEN
        nodeCV = node3
        nodeA  = node1
        nodeB  = node2
      END IF ! node1 

!     3. Compute edge01 (CV-->GP), edge12 (CV-->A), edge13 (CV-->B) vectors
!        Project vector GP-CV onto CV-A or CV-B to find temporary edge intercept
!        If the projectedLength is < 0 or > edgeLength, EI is outside
!        Else EI is inside then compute its Location and distance to GP
!     --------------------------------------------------------------------------
      vec01x = xGP - xCV
      vec01y = yGP - yCV
      vec01z = zGP - zCV

      nEdges = 2

      DO iEdge = 1, nEdges

        SELECT CASE(iEdge)

        CASE(1)
          nodeV = nodeA
        CASE(2)
          nodeV = nodeB
        END SELECT ! iEdge

        vec12x = xBodyMarker(iBody,nodeV) - xCV
        vec12y = yBodyMarker(iBody,nodeV) - yCV
        vec12z = zBodyMarker(iBody,nodeV) - zCV

        magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)

        magnitude12Inv = 1.0_CGREAL/magnitude12
     
        vec12x = vec12x*magnitude12Inv
        vec12y = vec12y*magnitude12Inv
        vec12z = vec12z*magnitude12Inv

        projectedLength = vec01x*vec12x +vec01y*vec12y +vec01z*vec12z 
 
!       -----------------------------------------------------------------------
!        Edge-Intercept (EI) point is outside Edge if pL < 0 or pL>magnitude12
!        else EI is inside and compute coordinates and distance
!        No need to take SQRT for distGPEI to save computations
!        Load EI into temporary value
!       -----------------------------------------------------------------------
        IF ( projectedLength < 0.0_CGREAL .OR. &
             projectedLength > magnitude12     ) THEN
          xEITemp = xCV
          yEITemp = yCV
          zEITemp = zCV
          edgeNode1 = nodeCV
          edgeNode2 = nodeCV
        ELSE
          xEITemp = xCV + projectedLength*vec12x
          yEITemp = yCV + projectedLength*vec12y
          zEITemp = zCV + projectedLength*vec12z
          edgeNode1 = nodeCV
          edgeNode2 = nodeV
        END IF ! projectedLength

!       Find mininum value of |GP-EI| and corresponding EI
!       --------------------------------------------------
        distGPEI = (xEITemp-xGP)**2.0_CGREAL &
                 + (yEITemp-yGP)**2.0_CGREAL &
                 + (zEITemp-zGP)**2.0_CGREAL 
  
        IF ( distGPEI < distGPEIMin ) THEN
          distGPEIMin = distGPEI
          xEI = xEITemp
          yEI = yEITemp
          zEI = zEITemp
          closestEdgeNode1 = edgeNode1 
          closestEdgeNode2 = edgeNode2 
        END IF ! distGPEI
      END DO ! iEdge

      xBI = xEI 
      yBI = yEI
      zBI = zEI


!DEBUG
!    IF (iGP == 126 .and. jGP == 23 .and. kGP == 12) then
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xGP,yGP,zGP
!       WRITE(355,*)xBI,yBI,zBI
!    ENDIF
!DEBUG

! TEMPORARY
!    WRITE(365,*)xBI,yBI,zBI
! END TEMPORARY

999   CONTINUE

      xBIT(nc) = xBI
      yBIT(nc) = yBI
      zBIT(nc) = zBI
      cElementGP(nc) = closestElementGP
      cEdgeNode1(nc) = closestEdgeNode1
      cEdgeNode2(nc) = closestEdgeNode2

    ENDDO  ! nc

    DO nc = 1,nCheck
      dist(nc) =  (xBIT(nc)-xGP)**2.0_CGREAL &
                + (yBIT(nc)-yGP)**2.0_CGREAL &
                + (zBIT(nc)-zGP)**2.0_CGREAL 
	  IF(ndim == DIM_2D) THEN
       IF(abs(zBIT(nc)-zGP) > dz(1)/2.0_CGREAL) THEN
          dist(nc) = 1.0e26_CGREAL
       ENDIF
      ENDIF    			 
    ENDDO


    iDummy           = MINLOC(dist(1:nCheck))
    shortestProbe    = iDummy(1)
    xBI              = xBIT(shortestProbe)
    yBI              = yBIT(shortestProbe)
    zBI              = zBIT(shortestProbe)
    closestElementGP = cElementGP(shortestProbe)

!DEBUG
!  IF (igp == 102 .and. jgp == 93 .and. kgp == 1) then
!     print*,'n1,n2'
!     print*,cEdgeNode1(shortestProbe),cEdgeNode2(shortestProbe)
!     print*,iBody,membrane_type(iBody),OPEN_MEMBRANE
!  ENDIF
!DEBUG

    IF (membrane_type(iBody) == OPEN_MEMBRANE) THEN
      CALL determine_BIValidityMembrane(iBody,cEdgeNode1(shortestProbe),cEdgeNode2(shortestProbe),validIntercept)     
!DEBUG
!    print*,'valid',iGP,jGP,kGP,validIntercept,cEdgeNode1(shortestProbe),cEdgeNode2(shortestProbe)
!DEBUG
      IF (validIntercept .EQV. .FALSE.) THEN
        xBI = -1.0_CGREAL
        yBI = -1.0_CGREAL
        zBI = -1.0_CGREAL
        closestElementGP = -1
      ENDIF
    ENDIF

    DEALLOCATE(NeighElemInd)
    DEALLOCATE(distMarker)

END SUBROUTINE GCM_calc_bodyIntercept_Unstruc
!---------------------------------------------------------------------



SUBROUTINE determine_BIValidityMembrane(iBody,n1,n2,valid)          

! ----------------------------------------------------------------------------------------
!  Purpose: Subroutine  is designed to determine if a grid cell qualifies as a true ghost
!  cell for an open membrane
! 
!  Input: iBody : body corresponding to closest element to cell
!         n1,n2 : node number of two nodes that define the edge on which
!                 lies the body interceprt
!         n1=n2=0 : there is a true normal intercept
!         n1=n2   : no normal intercpet. Closest intercept point is on a element node
!         n1/=n2  : no normal intercpet. Closest intercept point is edge defined by nodes n1 and n2
!
!  Output:         valid : .TRUE. => grid cell is a true ghost point.
! ----------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN)  :: iBody,n1,n2
    LOGICAL, INTENT(OUT) :: valid
!   
!... loop variables

    INTEGER :: j,jj

!... local variables

    INTEGER :: nn1,nn2,nn3,edge_shared


    valid = .TRUE.

    IF (n1 == 0 .AND. n2 == 0) THEN
       valid = .TRUE.
       RETURN
    ENDIF

    IF (edge_node(iBody,n1) == 1 .AND. edge_node(iBody,n2) == 1) THEN
!     PRINT*,'Closest node or edge is on edge of membrane '
!     PRINT*,'Therefore this is not a true ghost cell'
      valid = .FALSE.
      RETURN
    ENDIF


 END SUBROUTINE determine_BIValidityMembrane




SUBROUTINE GCM_calc_BIVelocity_Unstruc( iGBI, jGBI, kGBI, xGBI, yGBI, zGBI, closestElementGBI,          &
                                        uGBI, vGBI, wGBI )

! ----------------------------------------------------------------------------
!  Purpose: generalized kernel to compute the velocity components of
!           the intercept points for any generic point onto the body markers.
!
!   Input: iGP, jGP, kGP      = indices of Generic Point,
!             closestMarkerGP = closest Marker for Generic point,
!
!   Output: uGP, vGP, wGP = velocity components of Generic Point.
! ----------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE unstructured_surface_arrays

    IMPLICIT NONE

!... parameters

    INTEGER,           INTENT(IN)  :: iGBI, jGBI, kGBI, closestElementGBI 
    REAL(KIND=CGREAL), INTENT(IN)  :: xGBI, yGBI,zGBI
    REAL(KIND=CGREAL), INTENT(OUT) :: uGBI, vGBI, wGBI
!   
!... loop variables

    INTEGER :: i

!... local variables

    INTEGER                           :: iBody,node1,node2,node3,nMarker
    INTEGER                           :: info
    REAL(KIND=CGREAL)                 :: cX, cY, cZ, cC, rCond
    REAL(KIND=CGREAL), DIMENSION(4,4) :: vanTri
    REAL(KIND=CGREAL), DIMENSION(4)   :: rhsTri

! use the following approach
!
!  u = a x + b y + c z + d 
!  (a,b,c,d) determined by using four conditions
!   u = u(i) at ith node for i=1,3
!   GRAD(u) . n = 0  where n is normal to plane of triangle.
!  
!******************************************************************************

    iBody   = bodyNum(iGBI,jGBI,kGBI)
	
	IF(iBody == 0) THEN
	  iBody = max(bodyNum(iGBI+1,jGBI,kGBI),bodyNum(iGBI-1,jGBI,kGBI),bodyNum(iGBI,jGBI+1,kGBI),bodyNum(iGBI,jGBI-1,kGBI))
	  IF( ndim == 3 ) iBody = max(iBody,bodyNum(iGBI,jGBI,kGBI+1),bodyNum(iGBI,jGBI,kGBI-1))
	ENDIF

	
    IF (iBody == 0) THEN
       PRINT*,'bodyNum = 0 for',iGBI,jGBI,kGBI,'in GCM_calc_BIVelocity_Unstruc'
       STOP
    ENDIF
    nMarker = nPtsBodyMarker(iBody)
!  
!   assume linear variation of velocity across element and then compute value at intercept.

    node1   = triElemNeig(iBody,1,closestElementGBI)
    node2   = triElemNeig(iBody,2,closestElementGBI)
    node3   = triElemNeig(iBody,3,closestElementGBI)
!
    vanTri(1,1) = xBodyMarker(iBody,node1) 
    vanTri(1,2) = yBodyMarker(iBody,node1) 
    vanTri(1,3) = zBodyMarker(iBody,node1) 
    vanTri(1,4) = 1.0_CGREAL

    vanTri(2,1) = xBodyMarker(iBody,node2) 
    vanTri(2,2) = yBodyMarker(iBody,node2) 
    vanTri(2,3) = zBodyMarker(iBody,node2) 
    vanTri(2,4) = 1.0_CGREAL

    vanTri(3,1) = xBodyMarker(iBody,node3) 
    vanTri(3,2) = yBodyMarker(iBody,node3) 
    vanTri(3,3) = zBodyMarker(iBody,node3) 
    vanTri(3,4) = 1.0_CGREAL

    vanTri(4,1) = triElemNormx(iBody,closestElementGBI) 
    vanTri(4,2) = triElemNormy(iBody,closestElementGBI) 
    vanTri(4,3) = triElemNormz(iBody,closestElementGBI) 
    vanTri(4,4) = 0.0_CGREAL

    CALL DGETRF(4, 4, vanTri,4,iPvt, info)
    CALL DGETRI(4, vanTri,4,iPvt,work, 4, info)

! compute uGBI
    rhsTri(1) = uBodyMarker(iBody,node1)
    rhsTri(2) = uBodyMarker(iBody,node2)
    rhsTri(3) = uBodyMarker(iBody,node3)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    uGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute vGBI
    rhsTri(1) = vBodyMarker(iBody,node1)
    rhsTri(2) = vBodyMarker(iBody,node2)
    rhsTri(3) = vBodyMarker(iBody,node3)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    vGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute wGBI
    rhsTri(1) = wBodyMarker(iBody,node1)
    rhsTri(2) = wBodyMarker(iBody,node2)
    rhsTri(3) = wBodyMarker(iBody,node3)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO
  
    wGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

END SUBROUTINE GCM_calc_BIVelocity_Unstruc
!---------------------------------------------------------------------
