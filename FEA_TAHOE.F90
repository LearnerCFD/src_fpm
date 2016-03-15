! --------------------------------------------------------------------
!  Flow Simulations and Analysis Group
!  The Johns Hopkins University
!
!  ViCar3D_PAT (ver. 2.0.0)
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
!  Filename: FEA_TAHOE.F90
!  Created by Rajneesh Bhardwaj on 02 April 2010,


!  Latest Modification: Feb, 25 2011 
!  by Rajneesh Bhardwaj
!
!  This file contains following subroutines:
!  initialize_tahoe
!  create_tahoe_geom
!  tahoe_solver
!  update_markers_positions_velo_tahoe
!  write_data_files_tahoe
!  write_output_tahoe
! --------------------------------------------------------------------
  
! --------------------------------------------------------------------
! Define global to local variables
# define G2LI(i)      i-(myIs-1)                             
# define G2LJ(j)      j-(myJs-1)

 SUBROUTINE initialize_tahoe 
   USE flow_parameters
   USE tahoe_parameters
   USE global_parameters
   USE unstructured_surface_arrays
   USE implicit_coupling_parameters     
   IMPLICIT NONE
   INTEGER :: i,j,a,iBody,k,kk,n,nelem,nn1,nn2,nn3
   REAL(KIND=CGREAL) :: dummy
   INTEGER, DIMENSION(:),ALLOCATABLE :: dummyint
   CHARACTER (LEN=80) :: string1,string2,string3,string4,string5
   CHARACTER (LEN=120) :: st1
   
!  Assign Tahoe file numbers
   tahoe_DDu = 740
   tahoe_elem0 = 750
   tahoe_rs = 760 
   
   OPEN(tahoe_inputfile,FILE='input_tahoe.dat',STATUS='UNKNOWN')
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)nstart_tahoe,ndump_tahoe
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)ConstitutiveModel,TahoeSolver
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)dtratio_tahoe
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)KBC_flag,FBC_flag
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)SolidDen,YmodulusEQ,YmodulusNEQ,PoissonR,etaDamping
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)elementType,iTahoeOutput,itahoedraglift,iTahoeRestart,iTahoeProbe 
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)nprobe_tahoe,verbosetahoe
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)implicit_coupling_flag,kimplicitMax,ImplicitResidualMax 
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)alpha_underrelax
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)flag_underrelax
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)kimplicit_start
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)slope_underrelax
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)varyalphawithtime,alpha_underrelax2,timeimpl1,timeimpl2
   
   ALLOCATE(numNodes(nbody),numelts(nbody),FBCnodes(nbody),KBCnodes(nbody))
   ALLOCATE(numMarkers(nbody),BodyProbe(nbody,2))
   
   DO n = 1, nbody
   SELECT CASE (boundary_motion_type(n))
   CASE(FEA_TAHOE)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*) numNodes(n),numelts(n),FBCNodes(n),KBCNodes(n)
   READ(tahoe_inputfile,*)
   READ(tahoe_inputfile,*) (BodyProbe(n,i),i = 1,nprobe_tahoe)
   READ(tahoe_inputfile,*)
   END SELECT
   END DO 
    
   CLOSE(tahoe_inputfile)
   
   WRITE(STDOUT,'(3X,A)')'Reading of input_tahoe.dat is complete'
   
   IF( nDim == DIM_2D) ndim_tahoe = 2  !ndim_tahoe is needed to pass to maintahoe
   IF( nDim == DIM_3D) ndim_tahoe = 3
   
   
IF(ImtheBoss)THEN
     
   IF(FBC_flag ==1 .OR. FBC_flag == 0 .OR. KBC_flag ==1 .OR. KBC_flag == 0)THEN
   PRINT*,"KBC and FBC flags are:",KBC_flag,FBC_flag
   ELSE
   PRINT*,"WRONG CHOICE OF KBC or FBC flag"
   STOP
   END IF
   IF(ConstitutiveModel >=1 .OR. ConstitutiveModel <= 6 )THEN
   PRINT*,"ConstitutiveModel is:",ConstitutiveModel
   ELSE
   PRINT*,"WRONG CHOICE OF ConstitutiveModel"
   STOP
   END IF
   IF(TahoeSolver ==1 .OR. TahoeSolver == 2 )THEN
   PRINT*,"Tahoe Solver is:",TahoeSolver
   ELSE
   PRINT*,"Wrong choice option of TahoeSolver"
   STOP
   END IF
   
   DO iBody = 1,nBody
   SELECT CASE (boundary_motion_type(iBody))
   CASE(FEA_TAHOE)
   IF(ndim_tahoe==2) numMarkers(iBody) = nPtsBodyMarker(iBody)/3 ! because there are three layers
   IF(ndim_tahoe==3) numMarkers(iBody) = nPtsBodyMarker(iBody)
   PRINT*,'Numbers of markers on body',iBody,'is',numMarkers(iBody)
   END SELECT
   END DO 
   
   nmax = 0
   emax = 0
   FBCmax = 0
   KBCmax = 0
   Markmax = 0
   DO n = 1, nbody
   SELECT CASE (boundary_motion_type(n))
   CASE(FEA_TAHOE)
     WRITE(*, *) 'Number of nodes in Tahoe Body',n,'is',numNodes(n)
	 WRITE(*, *) 'Number of elements in Tahoe Body',n,'is',numelts(n)
	 WRITE(*, *) 'Number of FBC nodes in Tahoe Body',n,'is',FBCNodes(n)
	 WRITE(*, *) 'Number of KBC nodes in Tahoe Body',n,'is',KBCNodes(n)
	 WRITE(*, *) 'Number of markers in Tahoe Body',n,'is',numMarkers(n)
     IF(nmax   < numNodes(n))    nmax   =  numNodes(n)
     IF(emax   < numelts(n))     emax   =  numelts(n)
     IF(FBCmax < FBCNodes(n))    FBCmax =  FBCNodes(n)
     IF(KBCmax < KBCNodes(n))    KBCmax =  KBCNodes(n)
     IF(Markmax < numMarkers(n)) Markmax = numMarkers(n)
   END SELECT
   END DO
   
   WRITE(*, *) 'Total Number of bodies is', nbody
   WRITE(*, *) 'Maximum number of elements among all Tahoe bodies is',emax
   WRITE(*, *) 'Maximum number of nodes among all Tahoe bodies is',nmax
   WRITE(*, *) 'Maximum number of FBC nodes among Tahoe bodies is',FBCmax
   WRITE(*, *) 'Maximum number of KBC nodes among Tahoe bodies is',KBCmax
   WRITE(*, *) 'Maximum number of Markers among Tahoe bodies is',Markmax
   
   WRITE(STDOUT,'(3X,A)')'Allocating memory now...'

!  Needs to be allocated with Tahoebodies   
   ALLOCATE(x0i(nmax,nbody),y0i(nmax,nbody),z0i(nmax,nbody))
   ALLOCATE(FBCnode2marker(FBCmax,nbody),Markers2femGridpt(Markmax,nbody))
   ALLOCATE(nodesFBC(FBCmax,nbody),nodesKBC(KBCmax,nbody),dummyint(Markmax),marker2elem(Markmax,nbody))
   ALLOCATE(utahoe(2*nmax),Dutahoe(2*nmax),DDutahoe(2*nmax))

   
   WRITE(STDOUT,'(3X,A)')'Memory allocation complete'
    
   dummy = dtratio_tahoe
   dtTahoe = dt/dummy
!  Calculate viscoelastic parameters
   muNEQ = YmodulusNEQ/(2.0*(1+PoissonR))
   muEQ = YmodulusEQ/(2.0*(1+PoissonR))
   kappaNEQ = YmodulusNEQ/(3.0*(1-2.0*PoissonR))
   kappaEQ = YmodulusEQ/(3.0*(1-2.0*PoissonR))
   tauBulk =  etaDamping/kappaNEQ
   tauShear = etaDamping/muNEQ
   
! Read data for boundary conditions (FBC and KBC)
  
DO iBody = 1,nBody
   
 SELECT CASE (boundary_motion_type(iBody))
 CASE(FEA_TAHOE)
 
   WRITE(string1,'("TahoeFBCnodesBody",i1.1,".dat")'), iBody
   PRINT*,'OPENED FILE', string1
   WRITE(string2,'("TahoeMarkers2FEMgridpointsBody",i1.1,".dat")'), iBody
   PRINT*,'OPENED FILE', string2
   WRITE(string3,'("TahoeKBCnodesBody",i1.1,".dat")'), iBody
   PRINT*,'OPENED FILE', string3
   
   OPEN (tahoe_fbcnodes, FILE=string1)
  
   DO i = 1,FBCnodes(iBody)
   READ (tahoe_fbcnodes,*)nodesFBC(i,iBody)
   END DO
   CLOSE(tahoe_fbcnodes)
   
   PRINT*,'READING OF FILE', string1,'IS COMPLETE'
   
   OPEN (tahoe_markers2fem, FILE=string2)
   DO i = 1,numMarkers(iBody) 
   READ (tahoe_markers2fem,*)dummyint(i), Markers2femGridpt(i,iBody)
   END DO
   CLOSE(tahoe_markers2fem)
   
   PRINT*,'READING OF FILE', string2,'IS COMPLETE'
   
   DO i = 1,numMarkers(iBody)-1 
   IF(dummyint(i) .gt. dummyint(i+1))THEN
   PRINT*, 'ERROR IN FORMAT OF',string2,'for Body',iBody
   STOP
   END IF
   END DO
   
   OPEN (tahoe_kbcnodes, FILE=string3)
   DO i = 1,KBCNodes(iBody)
   READ (tahoe_kbcnodes,*)nodesKBC(i,iBody)
   END DO
   CLOSE(tahoe_kbcnodes)
   
   PRINT*,'READING OF FILE', string3,'IS COMPLETE'
   
   DO i = 1,FBCnodes(iBody)
   FBCnode2marker(i,iBody)=0
   END DO
   
  DO i = 1,FBCnodes(iBody)
    k = nodesFBC(i,iBody)
   DO j = 1,numMarkers(iBody)
    kk = Markers2femGridpt(j,iBody)
    IF(k==kk)THEN
    FBCnode2marker(i,iBody)=j
    END IF
   END DO
     IF(FBCnode2marker(i,iBody)==0)THEN
     PRINT*,"ERROR IN FINDING FBC NODE"
     STOP
     END IF
  END DO
 END SELECT
  
END DO

IF(nread /= 1)THEN

PRINT*,"INITIALIZING u AND Du FILES"

dummy = 0.0
DO iBody = 1,nBody
   
   SELECT CASE (boundary_motion_type(iBody))
   CASE(FEA_TAHOE)
   
 WRITE(string1,'("tahoe_body",i1.1,".rs1of1")'), iBody
 WRITE(string2,'("tahoe_body",i1.1,".rs1of1.displacement.u")'), iBody
 WRITE(string3,'("tahoe_body",i1.1,".rs1of1.displacement.Du")'), iBody
 WRITE(string4,'("tahoe_body",i1.1,".rs1of1.displacement.DDu")'), iBody
 WRITE(string5,'("tahoe_body",i1.1,".rs1of1.elem0")'), iBody
   
 OPEN (tahoe_rs,    FILE=string1)
 OPEN (tahoe_u,     FILE=string2)
 OPEN (tahoe_Du,    FILE=string3)
 OPEN (tahoe_DDu,   FILE=string4)
 OPEN (tahoe_elem0, FILE=string5)
 
 DO i = 1,numNodes(iBody)
 IF(ndim_tahoe==2)WRITE (tahoe_u,*)  dummy,dummy
 IF(ndim_tahoe==2)WRITE (tahoe_Du,*) dummy,dummy
 IF(ndim_tahoe==2)WRITE (tahoe_DDu,*)dummy,dummy
 IF(ndim_tahoe==3)WRITE (tahoe_u,*)  dummy,dummy,dummy
 IF(ndim_tahoe==3)WRITE (tahoe_Du,*) dummy,dummy,dummy
 IF(ndim_tahoe==3)WRITE (tahoe_DDu,*)dummy,dummy,dummy
 END DO
 
 WRITE (tahoe_rs,*)0.0
 WRITE (tahoe_rs,*)0.0
 WRITE (tahoe_rs,*)dtTahoe
 
 DO i = 1,numelts(iBody)
 WRITE (tahoe_elem0,*)1 
 END DO
 
 CLOSE(tahoe_rs)
 CLOSE(tahoe_u)
 CLOSE(tahoe_Du)
 CLOSE(tahoe_DDu)
 CLOSE(tahoe_elem0)
 
 
END SELECT
END DO

! Open probe files
IF(iTahoeProbe==1)THEN
PRINT*,"INITIALIZING PROBE FILES...."
DO iBody = 1,nBody
SELECT CASE (boundary_motion_type(iBody))
CASE(FEA_TAHOE)
 DO i = 1,nprobe_tahoe
   WRITE(string1,'("data_body",i1.1,"_probe",i1.1,".dat")'), iBody,i
   OPEN (tahoe_dataprobe, file = string1)
   WRITE (tahoe_dataprobe, *) 'ntime time X Y Z P num_fresh num_dead'
   CLOSE (tahoe_dataprobe)
 END DO 
END SELECT
END DO !nBody
END IF

ELSE 

DO iBody = 1,nBody
SELECT CASE (boundary_motion_type(iBody))
CASE(FEA_TAHOE)
WRITE(st1,'("cp restart_tahoe_body",i1.1,".rs1of1                   tahoe_body",i1.1,".rs1of1")'), iBody,iBody
CALL system (st1)  
WRITE(st1,'("cp restart_tahoe_body",i1.1,".rs1of1.displacement.u    tahoe_body",i1.1,".rs1of1.displacement.u")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp restart_tahoe_body",i1.1,".rs1of1.displacement.Du   tahoe_body",i1.1,".rs1of1.displacement.Du")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp restart_tahoe_body",i1.1,".rs1of1.displacement.DDu  tahoe_body",i1.1,".rs1of1.displacement.DDu")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp restart_tahoe_body",i1.1,".rs1of1.elem0             tahoe_body",i1.1,".rs1of1.elem0")'), iBody,iBody
CALL system (st1)
END SELECT
END DO


END IF!nread


WRITE(STDOUT,'(3X,A)')'CREATING GEOM FILES....'
CALL create_tahoe_geom
  
IF(ndim_tahoe==2)THEN

WRITE(STDOUT,'(3X,A)')'Calculating correspondence table for 2D case...'

!markers2elem needed for all bodies for calculating left and drag
DO iBody = 1,nBody
    DO i=1,numMarkers(iBody)-1 
    nelem = 0 
    DO j = 1,totNumTriElem(iBody)
    nn1 = triElemNeig(iBody,1,j) !node1 of the element
    nn2 = triElemNeig(iBody,2,j) !node2 of the element
    nn3 = triElemNeig(iBody,3,j) !node3 of the element 
     IF((nn1.eq.i .AND. nn2.eq.i+1).OR.(nn2.eq.i .AND. nn1.eq.i+1).OR. &
     (nn2.eq.i .AND. nn3.eq.i+1).OR.(nn3.eq.i .AND. nn2.eq.i+1).OR. &
     (nn1.eq.i .AND. nn3.eq.i+1).OR.(nn3.eq.i .AND. nn1.eq.i+1))THEN
     nelem = j
!	PRINT*,'element #',i,nelem
     END IF
     END DO
     IF(nelem.eq.0)THEN
     PRINT*,"ELEMENT NOT FOUND....STOP...."
     STOP
   END IF
    marker2elem(i,iBody) = nelem 
   END DO !numMarkers
 END DO !iBody
END IF !ndim_tahoe

! Open data files if needed 

IF(implicit_coupling_flag .AND. nread == 0)THEN
OPEN (imconverg, file = 'implicit_convergence.dat')
WRITE (imconverg, *) 'ntime, ImplicitResidual kimplicit'
CLOSE (imconverg)
END IF !implicit_coupling_flag

WRITE(STDOUT,'(3X,A)')"ALL TAHOE RELATED DATASTRUCTURES ARE SUCCESSFULLY CREATED"
  
END IF !ImtheBoss
  
END SUBROUTINE initialize_tahoe
! --------------------------------------------------------------------

SUBROUTINE tahoe_solver (iBody)
    USE global_parameters
    USE flow_parameters 
    USE grid_arrays
    USE boundary_arrays
    USE unstructured_surface_arrays 
    USE tahoe_parameters
!  This subroutine calls tahoe to calculate the deformation of a solid body. 
!  First it creates a xml file and then call tahoe main solver. The formulation 
!  of calculating the nodal forces is given in Xudong thesis at page 49 and 
!  Appendix II. 

   IMPLICIT NONE
   REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: forcex,forcey,forcez,xx,yy,zz,psn,ps
   INTEGER           :: iBody,i1,nn1,nn2,nn3
   REAL(KIND=CGREAL) :: areatri,a,peq1,peq2,peq3,L,dummy,normal(3),xx3,yy3,zz3
   CHARACTER (LEN=80) :: command1
   INTEGER ::i,j,k,n,n2,nelem
   CHARACTER text
   CHARACTER (LEN=30) :: string1,string2,string3
       
   ALLOCATE(forcex(numMarkers(iBody)),forcey(numMarkers(iBody)),forcez(numMarkers(iBody)),&
   xx(numMarkers(iBody)),yy(numMarkers(iBody)),zz(numMarkers(iBody)),ps(numMarkers(iBody)))
   
IF(verbosetahoe==1) PRINT*,'Initialize forces on the body',iBody

DO i=1,numMarkers(iBody)
     xx(i) = xBodyMarker(iBody,i)  
     yy(i) = yBodyMarker(iBody,i) 
     zz(i) = zBodyMarker(iBody,i) 
     ps(i) = pBodyMarker(iBody,i)  !   ps is pressure at grid point
     forcex(i) = 0.0_CGREAL
     forcey(i) = 0.0_CGREAL 
     forcez(i) = 0.0_CGREAL       
END DO !i


IF(ndim_tahoe==2)THEN  
    
  DO i=1,numMarkers(iBody)-1 
    nelem = marker2elem(i,iBody)
    L = sqrt((xx(i) - xx(i+1))**2.0_CGREAL + (yy(i) - yy(i+1))**2.0_CGREAL) 
    peq1 = L * (ps(i)/3.0_CGREAL + ps(i+1)/6.0_CGREAL)
    peq2 = L * (ps(i+1)/3.0_CGREAL + ps(i)/6.0_CGREAL)
    forcex(i) = (peq1+peq2)*triElemNormx(iBody,nelem)
    forcey(i) = (peq1+peq2)*triElemNormy(iBody,nelem)
  ENDDO ! i
  
  !for last point
  i = numMarkers(iBody)
  nelem = marker2elem(i,iBody)
  L = sqrt((xx(i) - xx(1))**2.0_CGREAL + (yy(i) - yy(1))**2.0_CGREAL)
  peq1 = L * (ps(i)/3.0_CGREAL + ps(1)/6.0_CGREAL)
  peq2 = L * (ps(1)/3.0_CGREAL + ps(i)/6.0_CGREAL)
  forcex(i) =  (peq1+peq2)*triElemNormx(iBody,nelem)
  forcey(i) =  (peq1+peq2)*triElemNormy(iBody,nelem)
  END IF

IF(ndim_tahoe==3)THEN
DO i = 1,totNumTriElem(iBody)
    nn1 = triElemNeig(iBody,1,i) !node1 of the element
    nn2 = triElemNeig(iBody,2,i) !node2 of the element
    nn3 = triElemNeig(iBody,3,i) !node3 of the element

    CALL fea_compute_tri_area(xx(nn1), yy(nn1), zz(nn1), & 
                             xx(nn2), yy(nn2),zz(nn2), &
                             xx(nn3), yy(nn3), zz(nn3), areatri)
   
      peq1 = areatri *(ps(nn1) / 6.0_CGREAL + ps(nn2)/ 12.0_CGREAL &
                + ps(nn3)/ 12.0_CGREAL)
      peq2 = areatri *(ps(nn2) / 6.0_CGREAL + ps(nn3)/ 12.0_CGREAL &
                + ps(nn1)/ 12.0_CGREAL)
      peq3 = areatri *(ps(nn3) / 6.0_CGREAL + ps(nn1)/ 12.0_CGREAL &
                + ps(nn2)/ 12.0_CGREAL)

!Calculate nodal forces
     forcex(nn1) =  forcex(nn1) + peq1*triElemNormx(iBody,i)             
     forcey(nn1) =  forcey(nn1) + peq1*triElemNormy(iBody,i) 
     forcez(nn1) =  forcez(nn1) + peq1*triElemNormz(iBody,i)
     forcex(nn2) =  forcex(nn2) + peq2*triElemNormx(iBody,i)
     forcey(nn2) =  forcey(nn2) + peq2*triElemNormy(iBody,i)
     forcez(nn2) =  forcez(nn2) + peq2*triElemNormz(iBody,i)
     forcex(nn3) =  forcex(nn3) + peq3*triElemNormx(iBody,i)
     forcey(nn3) =  forcey(nn3) + peq3*triElemNormy(iBody,i)
     forcez(nn3) =  forcez(nn3) + peq3*triElemNormz(iBody,i)  

     END DO ! i = totNumTriElem(iBody)

    END IF ! ndim = DIM_3D


IF(verbosetahoe==1) PRINT*,'Writing XML file'


IF(ntime<nstart_tahoe+2)THEN ! FOR COMPILING TOGETHER

   WRITE(string1,'("tahoe_body",i1.1,".xml")'),iBody
   WRITE(string2,'("tahoe_body",i1.1,".geom")'),iBody
   WRITE(string3,'("tahoe_body",i1.1,".rs1of1")'),iBody
   OPEN (tahoe_xml, FILE=string1)
   WRITE (tahoe_xml,5)'<?xml version="1.0" encoding="UTF-8"?>'
   WRITE (tahoe_xml,*)'<tahoe echo_input="false" geometry_file="',string2,'"' 
   WRITE (tahoe_xml,*)'output_format="TecPlot"' 
   
   
   IF(ntime>nstart_tahoe .AND. dtratio_tahoe==1)THEN
   WRITE(tahoe_xml,*)'restart_file="',string3,'"'
   END IF 
    
   IF(dtratio_tahoe>1 .AND. iflag_restart_tahoe==1)THEN
   WRITE(tahoe_xml,*)'restart_file="',string3,'"'
   END IF 
   
   
    WRITE (tahoe_xml,*)'restart_output_inc="1"' 
    WRITE (tahoe_xml,*)'xmlns:x0="http://www.w3.org/2001/XMLSchema">'
    WRITE (tahoe_xml,*)'<time max_step_cuts = "0" num_steps="1" output_inc="1"' 
    WRITE (tahoe_xml,*)'time_step="',dtTahoe,'">'
    WRITE (tahoe_xml,*)'</time>'
    WRITE (tahoe_xml,*)'<nodes>'
    WRITE (tahoe_xml,*)'<field field_name="displacement"' 

    IF(TahoeSolver==1)THEN 
    WRITE (tahoe_xml,*)'integrator="linear_HHT">'
    END IF
    IF(TahoeSolver==2)THEN 
    WRITE (tahoe_xml,*)'integrator="nonlinear_HHT">'
    END IF

    WRITE (tahoe_xml,*)'<dof_labels>'
    WRITE (tahoe_xml,*)'<String value="D_X"/>'
    WRITE (tahoe_xml,*)'<String value="D_Y"/>'
    IF(ndim==DIM_3D)THEN
    WRITE (tahoe_xml,*)'<String value="D_Z"/>'
    END IF
    WRITE (tahoe_xml,*)'</dof_labels>'

    IF(KBC_flag==1)THEN
    i = 1
    WRITE (tahoe_xml,*)'<kinematic_BC dof="1" node_ID="',i,'"/>'
    WRITE (tahoe_xml,*)'<kinematic_BC dof="2" node_ID="',i,'"/>'
    IF(ndim==DIM_3D)THEN
    WRITE (tahoe_xml,*)'<kinematic_BC dof="3" node_ID="',i,'"/>'
    END IF
    END IF !KBC_flag

    IF(FBC_flag==1)THEN
    !Node ID i+1 is used because the first nodeset is with kinematic BC
    DO i = 1 , FBCNodes(iBody) 
    k = FBCnode2marker(i,iBody)
    WRITE (tahoe_xml,*)'<force_BC dof="1" node_ID="',i+1,'" value="',forcex(k),'"/>'
    END DO
    DO i = 1 , FBCNodes(iBody) 
    k = FBCnode2marker(i,iBody)
    WRITE (tahoe_xml,*)'<force_BC dof="2" node_ID="',i+1,'" value="',forcey(k),'"/>'
    END DO
    IF(ndim==DIM_3D)THEN
    DO i = 1 , FBCNodes(iBody) 
    k = FBCnode2marker(i,iBody)
    WRITE (tahoe_xml,*)'<force_BC dof="3" node_ID="',i+1,'" value="',forcez(k),'"/>'
    END DO
    END IF !DIM_3D
    END IF !FBC_flag

    WRITE (tahoe_xml,*)'</field>'
    WRITE (tahoe_xml,*)'</nodes>'
    WRITE (tahoe_xml,*)'<element_list>'
   
    IF(ConstitutiveModel==1)THEN !   Small strain St Venant model
    WRITE (tahoe_xml,*)'<small_strain field_name="displacement">'
	IF(elementType==1)WRITE (tahoe_xml,*)'<triangle/>'
    IF(elementType==2)WRITE (tahoe_xml,*)'<quadrilateral/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<small_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<small_strain_material_2D>'
    WRITE (tahoe_xml,*)'<small_strain_StVenant_2D'
    WRITE (tahoe_xml,*)'constraint_2D="plane_strain" density="',SolidDen,'">'
    WRITE (tahoe_xml,*)'<E_and_nu Poisson_ratio="',PoissonR,'"' 
    WRITE (tahoe_xml,*)'Young_modulus="',YmodulusEQ,'"/>'
    WRITE (tahoe_xml,*)'</small_strain_StVenant_2D>'
    WRITE (tahoe_xml,*)'</small_strain_material_2D>'
    WRITE (tahoe_xml,*)'</small_strain_element_block>'
    WRITE (tahoe_xml,*)'</small_strain>'
    END IF

    IF(ConstitutiveModel==2)THEN !   Small strain Linear Viscoelatic model
    WRITE (tahoe_xml,*)'<small_strain field_name="displacement">'
    IF(elementType==1)WRITE (tahoe_xml,*)'<triangle/>'
    IF(elementType==2)WRITE (tahoe_xml,*)'<quadrilateral/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<small_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<small_strain_material_2D>'
    WRITE (tahoe_xml,*)'<linear_viscoelastic_2D constraint_2D="plane_strain"'
    WRITE (tahoe_xml,*)'density="',SolidDen,'" kappa_EQ="',kappaEQ,'"' 
    WRITE (tahoe_xml,*)'kappa_NEQ="',kappaNEQ,'"'
    WRITE (tahoe_xml,*)'mu_EQ="',muEQ ,'" mu_NEQ="',muNEQ ,'"' 
    WRITE (tahoe_xml,*)'tau_bulk="',tauBulk,'" tau_shear="',tauShear,'"/>'
    WRITE (tahoe_xml,*)'</small_strain_material_2D>'
    WRITE (tahoe_xml,*)'</small_strain_element_block>'
    WRITE (tahoe_xml,*)'</small_strain>'
    END IF

    IF(ConstitutiveModel==3)THEN !   Large strain model
    WRITE (tahoe_xml,*)'<updated_lagrangian field_name="displacement">'
    IF(elementType==1)WRITE (tahoe_xml,*)'<triangle/>'
    IF(elementType==2)WRITE (tahoe_xml,*)'<quadrilateral/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<large_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<large_strain_material_2D>'
    WRITE (tahoe_xml,*)'<large_strain_StVenant_2D'
    WRITE (tahoe_xml,*)'constraint_2D="plane_strain" density="',SolidDen,'">'
    WRITE (tahoe_xml,*)'<E_and_nu Poisson_ratio="',PoissonR,'" Young_modulus="',YmodulusEQ,'"/>'
    WRITE (tahoe_xml,*)'</large_strain_StVenant_2D>'
    WRITE (tahoe_xml,*)'</large_strain_material_2D>'
    WRITE (tahoe_xml,*)'</large_strain_element_block>'
    WRITE (tahoe_xml,*)'</updated_lagrangian>'
    END IF

    IF(ConstitutiveModel==4)THEN !   Small strain St Venant model
    WRITE (tahoe_xml,*)'<small_strain field_name="displacement">'
	IF(elementType==3)WRITE (tahoe_xml,*)'<tetrahedron/>'
    IF(elementType==4)WRITE (tahoe_xml,*)'<hexahedron/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<small_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<small_strain_material_3D>'
    WRITE (tahoe_xml,*)'<small_strain_StVenant'
    WRITE (tahoe_xml,*)'density="',SolidDen,'">'
    WRITE (tahoe_xml,*)'<E_and_nu Poisson_ratio="',PoissonR,'"' 
    WRITE (tahoe_xml,*)'Young_modulus="',YmodulusEQ,'"/>'
    WRITE (tahoe_xml,*)'</small_strain_StVenant>'
    WRITE (tahoe_xml,*)'</small_strain_material_3D>'
    WRITE (tahoe_xml,*)'</small_strain_element_block>'
    WRITE (tahoe_xml,*)'</small_strain>'
    END IF

    IF(ConstitutiveModel==5)THEN !   Small strain Linear Viscoelatic model
    WRITE (tahoe_xml,*)'<small_strain field_name="displacement">'
    IF(elementType==3)WRITE (tahoe_xml,*)'<tetrahedron/>'
    IF(elementType==4)WRITE (tahoe_xml,*)'<hexahedron/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<small_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<small_strain_material_3D>'
    WRITE (tahoe_xml,*)'<linear_viscoelastic'
    WRITE (tahoe_xml,*)'density="',SolidDen,'" kappa_EQ="',kappaEQ,'" kappa_NEQ="',kappaNEQ,'"'
    WRITE (tahoe_xml,*)'mu_EQ="',muEQ ,'" mu_NEQ="',muNEQ ,'" tau_bulk="',tauBulk,'" tau_shear="',tauShear,'"/>'
    WRITE (tahoe_xml,*)'</small_strain_material_3D>'
    WRITE (tahoe_xml,*)'</small_strain_element_block>'
    WRITE (tahoe_xml,*)'</small_strain>'
    END IF

    IF(ConstitutiveModel==6)THEN !   Large strain model
    WRITE (tahoe_xml,*)'<updated_lagrangian field_name="displacement">'
    IF(elementType==3)WRITE (tahoe_xml,*)'<tetrahedron/>'
    IF(elementType==4)WRITE (tahoe_xml,*)'<hexahedron/>'
    WRITE (tahoe_xml,*)'<solid_element_nodal_output displacements="1" stress="1"/>'
    WRITE (tahoe_xml,*)'<solid_element_element_output linear_momentum="1"/>'
    WRITE (tahoe_xml,*)'<large_strain_element_block>'
    WRITE (tahoe_xml,*)'<block_ID_list>'
    WRITE (tahoe_xml,*)'<String value="1"/>'
    WRITE (tahoe_xml,*)'</block_ID_list>'
    WRITE (tahoe_xml,*)'<large_strain_material_3D>'
    WRITE (tahoe_xml,*)'<large_strain_StVenant'
    WRITE (tahoe_xml,*)'density="',SolidDen,'">'
    WRITE (tahoe_xml,*)'<E_and_nu Poisson_ratio="',PoissonR,'" Young_modulus="',YmodulusEQ,'"/>'
    WRITE (tahoe_xml,*)'</large_strain_StVenant>'
    WRITE (tahoe_xml,*)'</large_strain_material_3D>'
    WRITE (tahoe_xml,*)'</large_strain_element_block>'
    WRITE (tahoe_xml,*)'</updated_lagrangian>'
    END IF

    WRITE (tahoe_xml,*)'</element_list>'

    IF(TahoeSolver==1)THEN
    WRITE (tahoe_xml,*)'<linear_solver>'
    WRITE (tahoe_xml,*)'<profile_matrix/>'
    WRITE (tahoe_xml,*)'</linear_solver>'
    END IF
    IF(TahoeSolver==2)THEN
    WRITE (tahoe_xml,*)'<nonlinear_solver abs_tolerance="1.0e-8"' 
    WRITE (tahoe_xml,*)'check_code = "small_pivots" divergence_tolerance="1.0e+04"' 
    WRITE (tahoe_xml,*)'max_iterations="12" output_inc = "0" rel_tolerance="1.0e-8">'
    WRITE (tahoe_xml,*)'<SPOOLES_matrix/>'
    WRITE (tahoe_xml,*)'</nonlinear_solver>'
    END IF
    
    WRITE (tahoe_xml,*)'</tahoe>'
 
    CLOSE(tahoe_xml)
5   FORMAT(a)

END IF ! for ntime<nstart_tahoe+2

     DO i = 1,FBCNodes(iBody)
     k = FBCnode2marker(i,iBody)
     forcex(i) = forcex(k)
     forcey(i) = forcey(k)
     forcez(i) = forcez(k)
     END DO
	 IF(iBody.gt.5)THEN
	 WRITE(*,*)'Need modifications in maintahoe, stopping now'
	 STOP
	 END IF
     CALL maintahoe(i1,text,ndim_tahoe,numNodes(iBody),FBCNodes(iBody),&
     forcex,forcey,forcez,iBody,time,dtTahoe,utahoe,Dutahoe,DDutahoe)
 
 
     IF(dtratio_tahoe>1 .AND. ntime==nstart_tahoe)iflag_restart_tahoe=1
 
    DEALLOCATE(forcex,forcey,forcez,xx,yy,zz,ps) 

END SUBROUTINE tahoe_solver
 ! --------------------------------------------------------------------  
   
  SUBROUTINE update_markers_positions_velo_tahoe(iBody)
    USE flow_parameters
    USE boundary_arrays
    USE tahoe_parameters
    IMPLICIT NONE
    INTEGER ::i,k,n,iBody,nk 
	REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: u_x,v_y,w_z,d_x,d_y,d_z
    CHARACTER (LEN=80) :: string1,string2

    IF(ntime>=nstart_tahoe .AND. ImtheBOSS) THEN
	
	IF(verbosetahoe==1) PRINT*,'Allocating with nmax =',nmax
	ALLOCATE(d_x(nmax),d_y(nmax),u_x(nmax),v_y(nmax),d_z(nmax),w_z(nmax))

!   Need equal to sign (>=) because markers are updated in 
!   each outer iterations during implicit coupling.

       WRITE(string1,'("tahoe_body",i1.1,".rs1of1.displacement.Du")'), iBody 
       WRITE(string2,'("tahoe_body",i1.1,".rs1of1.displacement.u")'),  iBody 
   
       IF(verbosetahoe==1) PRINT*,'OPENED FILE',string1
       IF(verbosetahoe==1) PRINT*,'OPENED FILE',string2
   
       OPEN (tahoe_Du, FILE = string1)
       DO i = 1,numNodes(iBody)
       IF(ndim_tahoe==2)READ (tahoe_Du,*)u_x(i),v_y(i)
       IF(ndim_tahoe==3)READ (tahoe_Du,*)u_x(i),v_y(i),w_z(i)
       END DO
       CLOSE(tahoe_Du)
   
       OPEN (tahoe_u, FILE = string2)
       DO i = 1,numNodes(iBody)
       IF(ndim_tahoe==2)READ (tahoe_u,*)d_x(i),d_y(i)
       IF(ndim_tahoe==3)READ (tahoe_u,*)d_x(i),d_y(i),d_z(i)
       END DO
       CLOSE(tahoe_u)
   
       IF(verbosetahoe==1) PRINT*,'Reading of files complete'

    IF(ndim_tahoe==2)THEN
       DO i=1,numMarkers(iBody)
       k = Markers2femGridpt(i,iBody) ! convert MARKER to GRID POINT
       ubodymarker(iBody,i) = u_x(k)
       vbodymarker(iBody,i) = v_y(k)
       ubodymarker(iBody,i+numMarkers(iBody)) = u_x(k)
       vbodymarker(iBody,i+numMarkers(iBody)) = v_y(k)
       ubodymarker(iBody,i+2*numMarkers(iBody)) = u_x(k)
       vbodymarker(iBody,i+2*numMarkers(iBody)) = v_y(k)
       xBodyMarker(iBody,i) =  x0i(k,iBody) + d_x(k)
       yBodyMarker(iBody,i) =  y0i(k,iBody) + d_y(k)
       xBodyMarker(iBody,i+numMarkers(iBody)) = x0i(k,iBody) + d_x(k)
       yBodyMarker(iBody,i+numMarkers(iBody)) = y0i(k,iBody) + d_y(k)
       xBodyMarker(iBody,i+2*numMarkers(iBody)) = x0i(k,iBody) + d_x(k)
       yBodyMarker(iBody,i+2*numMarkers(iBody)) = y0i(k,iBody) + d_y(k)
       END DO
    END IF

    IF(ndim_tahoe==3)THEN
      DO i=1,numMarkers(iBody)
       k = Markers2femGridpt(i,iBody)
       ubodymarker(iBody,i) = u_x(k)
       vbodymarker(iBody,i) = v_y(k)
       wbodymarker(iBody,i) = w_z(k)
       xBodyMarker(iBody,i) =  x0i(k,iBody) + d_x(k)
       yBodyMarker(iBody,i) =  y0i(k,iBody) + d_y(k)
       zBodyMarker(iBody,i) =  z0i(k,iBody) + d_z(k)
     END DO
    END IF
  
  DEALLOCATE(d_x,d_y,d_z,u_x,v_y,w_z)
  END IF!nstart_tahoe !Imtheboss 
  
  IF(verbosetahoe==1) PRINT*,'Updated bodymarkers....'

#   ifdef MPI 
 CALL  par_bcast_marker(xBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
 CALL  par_bcast_marker(yBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
 CALL  par_bcast_marker(zBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
 CALL  par_bcast_marker(uBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
 CALL  par_bcast_marker(vBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
 CALL  par_bcast_marker(wBodyMarker, iBody, nBody, 0, nPtsBodyMarker(iBody))
#  endif
  
  

END SUBROUTINE update_markers_positions_velo_tahoe
! --------------------------------------------------------------------
SUBROUTINE write_probe_files_tahoe
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE tahoe_parameters
    IMPLICIT NONE
    INTEGER           :: i,iBody,k,numprobe
    REAL(KIND=CGREAL) :: dummy,xmar,ymar,zmar,pmar
    CHARACTER (LEN=80) :: string1
	
IF(verbosetahoe==1)WRITE(*,*)'Writing probe data of Tahoe...'

DO iBody = 1,nBody
   SELECT CASE (boundary_motion_type(iBody))
   CASE(FEA_TAHOE)
   DO i = 1,nprobe_tahoe
    WRITE(string1,'("data_body",i1.1,"_probe",i1.1,".dat")'), iBody,i
    xmar = xBodyMarker(iBody,BodyProbe(iBody,i))
    ymar = yBodyMarker(iBody,BodyProbe(iBody,i))
    zmar = zBodyMarker(iBody,BodyProbe(iBody,i))
    pmar = pBodyMarker(iBody,BodyProbe(iBody,i))
    OPEN (tahoe_dataprobe, file = string1, POSITION = 'APPEND')
    WRITE (tahoe_dataprobe, 11) ntime, time,xmar,ymar,zmar,pmar,num_fresh,num_dead
    CLOSE (tahoe_dataprobe)
  END DO
END SELECT 
END DO 

IF(verbosetahoe==1)WRITE(*,*)'Writing of probe data of Tahoe is complete'

11 FORMAT (i6,1x,1PE12.5,1x,1PE14.7,1x,1PE14.7,1x,1PE14.7,1x,1PE12.5,1x,i6,1x,i6)

ENDSUBROUTINE write_probe_files_tahoe

 SUBROUTINE create_tahoe_geom
    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE tahoe_parameters
    IMPLICIT NONE
    INTEGER           :: i,iBody,k,nen,kcounter
    REAL(KIND=CGREAL) :: dummy
    INTEGER, ALLOCATABLE, DIMENSION(:) ::aelts,belts,celts,delts
    INTEGER, ALLOCATABLE, DIMENSION(:) ::eelts,felts,gelts,helts
    CHARACTER (LEN=80) :: string1

IF(elementType==1)nen = 3 !traingles
IF(elementType==2)nen = 4 !quads
IF(elementType==3)nen = 4 !tets
IF(elementType==4)nen = 8 !bricks

IF(elementType==1)ALLOCATE(aelts(emax),belts(emax),celts(emax))
IF(elementType==2)ALLOCATE(aelts(emax),belts(emax),celts(emax),delts(emax))
IF(elementType==3)ALLOCATE(aelts(emax),belts(emax),celts(emax),delts(emax))
IF(elementType==4)ALLOCATE(aelts(emax),belts(emax),celts(emax),delts(emax),&
                           eelts(emax),felts(emax),gelts(emax),helts(emax))

DO iBody = 1,nBody

 SELECT CASE (boundary_motion_type(iBody))
 CASE(FEA_TAHOE)

WRITE(string1,'("TahoeBody",i1.1,".dat")'), iBody 
OPEN (tahoe_bodydat, FILE=string1)
READ(tahoe_bodydat,*)
READ(tahoe_bodydat,*)
READ(tahoe_bodydat,*)
READ(tahoe_bodydat,*)
DO i = 1,numNodes(iBody)
 IF(ndim_tahoe==2)READ (tahoe_bodydat,*)x0i(i,iBody),y0i(i,iBody)
 IF(ndim_tahoe==3)READ (tahoe_bodydat,*)x0i(i,iBody),y0i(i,iBody),z0i(i,iBody)
END DO
DO i = 1,numElts(iBody) 
 IF(elementType==1) READ (tahoe_bodydat,*)aelts(i),belts(i),celts(i) 
 IF(elementType==2) READ (tahoe_bodydat,*)aelts(i),belts(i),celts(i),delts(i)
 IF(elementType==3) READ (tahoe_bodydat,*)aelts(i),belts(i),celts(i),delts(i)
 IF(elementType==4) READ (tahoe_bodydat,*)aelts(i),belts(i),celts(i),delts(i),&
                                            eelts(i),felts(i),gelts(i),helts(i)
END DO
CLOSE(tahoe_bodydat)
 
PRINT*,'Reading of tahoe_body',iBody,'.dat is complete'
    
WRITE(string1,'("tahoe_body",i1.1,".geom")'), iBody    
OPEN (tahoe_geom, FILE=string1) !WRITE TAHOE compatible geometry file
    WRITE (tahoe_geom,*) '*version'
    WRITE (tahoe_geom,*) '1.0'
    WRITE (tahoe_geom,*) '*title'
    WRITE (tahoe_geom,*) 'geom_file'
    WRITE (tahoe_geom,*) '*dimensions'
    WRITE (tahoe_geom,*) numNodes(iBody),'# number of nodes'
    WRITE (tahoe_geom,*) ndim_tahoe,'# number of spatial dimensions'
    WRITE (tahoe_geom,*) '1 # number of element sets'
    WRITE (tahoe_geom,*) '# [ID] [nel] [nen]'
    WRITE (tahoe_geom,*) '1', numElts(iBody), nen  !Number of element sets is 1
 !Number of node sets = sets of nodes where BCs are applied
 !Node sets: 1 node set for inner surface (KBC)
 !Several (numFBCnodes) node sets for outer surface (FBC)
 !total node sets = numFBCnodes + 1
    WRITE (tahoe_geom,*) FBCNodes(iBody)+1,' # number of node sets'
    WRITE (tahoe_geom,*) '# [ID] [nnd]'  ! nnd = number of nodes in a particular node set
    WRITE (tahoe_geom,*) '1',KBCNodes(iBody)
    DO i = 2,FBCNodes(iBody)+1
    WRITE (tahoe_geom,*) i,'1'
    END DO
    WRITE (tahoe_geom,*) '0 # number of side sets'
    WRITE (tahoe_geom,*) '# [ID] [element set ID] [ns]'
    WRITE (tahoe_geom,*) '# end dimensions'
    WRITE (tahoe_geom,*) '*nodesets'
!first set 
    WRITE (tahoe_geom,*) '*set'
    WRITE (tahoe_geom,*) KBCNodes(iBody), '# number of nodes'
    DO i = 1,KBCNodes(iBody)
    WRITE (tahoe_geom,*) nodesKBC(i,iBody)
    END DO
!second set.....end of set
    DO i = 1,FBCnodes(iBody)
    WRITE (tahoe_geom,*) '*set'
    WRITE (tahoe_geom,*) '1 # number of nodes'
    WRITE (tahoe_geom,*)nodesFBC(i,iBody)
    END DO
    WRITE (tahoe_geom,*) '*#end node sets'
    WRITE (tahoe_geom,*) '*sidesets'
!WRITE (tahoe_geom,*) '*set'
!WRITE (tahoe_geom,*) '1 # number of sides'
!WRITE (tahoe_geom,*) '# [element] [face]'
    WRITE (tahoe_geom,*) '*elements'
    WRITE (tahoe_geom,*) '*set'
    WRITE (tahoe_geom,*) numElts(iBody),'# number of elements'
!Number of element nodes = Number of nodes per element 
    WRITE (tahoe_geom,*) nen, '# number of element nodes'
    kcounter = 0
    DO i=1,numElts(iBody)
    kcounter = kcounter + 1
    IF(elementType==1)WRITE (tahoe_geom,*) kcounter, aelts(i),belts(i),celts(i)
    IF(elementType==2)WRITE (tahoe_geom,*) kcounter, aelts(i),belts(i),celts(i),delts(i)
    IF(elementType==3)WRITE (tahoe_geom,*) kcounter, aelts(i),belts(i),celts(i),delts(i)
    IF(elementType==4)WRITE (tahoe_geom,*) kcounter, aelts(i),belts(i),celts(i),delts(i),&
                                                 eelts(i),felts(i),gelts(i),helts(i)

    END DO
    WRITE (tahoe_geom,*) '# end elements'
    WRITE (tahoe_geom,*) '*nodes'
    WRITE (tahoe_geom,*) numNodes(iBody),'# number of nodes'
    WRITE (tahoe_geom,*) ndim_tahoe,'# number of spatial dimensions'
    kcounter = 0
    DO i = 1,numNodes(iBody)
    kcounter = kcounter + 1
    IF(ndim_tahoe==2)WRITE (tahoe_geom,*) kcounter, x0i(i,iBody),y0i(i,iBody)
    IF(ndim_tahoe==3)WRITE (tahoe_geom,*) kcounter, x0i(i,iBody),y0i(i,iBody),z0i(i,iBody)
    END DO
    CLOSE(tahoe_geom)

    PRINT*,'Writing of tahoe_body',iBody,'.geom is complete'
 
END SELECT
END DO !ibody

IF(elementType==1)DEALLOCATE(aelts,belts,celts)
IF(elementType==2)DEALLOCATE(aelts,belts,celts,delts)
IF(elementType==3)DEALLOCATE(aelts,belts,celts,delts)
IF(elementType==4)DEALLOCATE(aelts,belts,celts,delts,eelts,felts,gelts,helts)

ENDSUBROUTINE create_tahoe_geom
! --------------------------------------------------------------------
SUBROUTINE write_output_tahoe
   USE global_parameters
   USE flow_parameters
   USE tahoe_parameters
   USE implicit_coupling_parameters 
   IMPLICIT NONE
   INTEGER  :: iBody 
   CHARACTER (LEN=80) :: command1,string1
   
   DO iBody = 1,nBody
    SELECT CASE (boundary_motion_type(iBody))
    CASE(FEA_TAHOE)
     IF(mod(ntime,ndump_tahoe)==0 .AND. ntime>=nstart_tahoe)THEN
     WRITE(command1,'("cp tahoe_body",i1.1,".io0.ps0000.dat tahoe_out_body",i1.1,".",i7.7)'), iBody,iBody,ntime
     CALL system (command1)
     END IF
    END SELECT
   END DO
   
END SUBROUTINE write_output_tahoe
! --------------------------------------------------------------------
SUBROUTINE write_restart_files_tahoe
   USE flow_parameters
   USE tahoe_parameters
   IMPLICIT NONE
   INTEGER  :: iBody ,i ,j
   CHARACTER (LEN=120) :: st1
     
IF(mod(ntime,nrestart)==0 .AND. ntime>=nstart_tahoe)THEN
WRITE(STDOUT,'(3X,A,I8)') 'Writing out restart file of Tahoe at ntime ' ,ntime

DO iBody = 1,nBody
SELECT CASE (boundary_motion_type(iBody))
CASE(FEA_TAHOE)
SELECT CASE(idxRstrt)
CASE(2)

WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1                  restart_tahoe_body",i1.1,".rs1of1.1")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.u   restart_tahoe_body",i1.1,".rs1of1.displacement.u.1")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.Du  restart_tahoe_body",i1.1,".rs1of1.displacement.Du.1")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.DDu restart_tahoe_body",i1.1,".rs1of1.displacement.DDu.1")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.elem0            restart_tahoe_body",i1.1,".rs1of1.elem0.1")'), iBody,iBody
CALL system (st1)

CASE(1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1                  restart_tahoe_body",i1.1,".rs1of1.2")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.u   restart_tahoe_body",i1.1,".rs1of1.displacement.u.2")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.Du  restart_tahoe_body",i1.1,".rs1of1.displacement.Du.2")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.displacement.DDu restart_tahoe_body",i1.1,".rs1of1.displacement.DDu.2")'), iBody,iBody
CALL system (st1)
WRITE(st1,'("cp tahoe_body",i1.1,".rs1of1.elem0            restart_tahoe_body",i1.1,".rs1of1.elem0.2")'), iBody,iBody
CALL system (st1)

END SELECT ! indexRstrt
END SELECT ! boundary_motion
END DO !iBody
END IF !nrestart 

ENDSUBROUTINE write_restart_files_tahoe

