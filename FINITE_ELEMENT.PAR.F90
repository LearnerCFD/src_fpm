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
!
!  Filename: FINITE_ELEMENT.PAR.F90
!  Latest Modification: March 10, 2010 (ver. P1.5.7)
!  Made by X. Zheng
! --------------------------------------------------------------------
!--------------------------------------------------------------
!  This code simulated structure by using Fnite Element Method.
!  This file contains the following subroutines
!  fea_initial
!  fea_readinput
!  fea_allocatememory
!  fea_deallocatememory
!  fea_openfiles
!  fea_readmesh
!  fea_readconload
!  fea_getcorresponancetable
!  fea_readcorresponancetable
!  fea_assemble
!  fea_elementmatrix
!  fea_hammerintegeration3D
!  fea_shapetetrahedron 
!  fea_guassintegration
!  fea_shapequadra8
!  fea_shapequadra9
!  fea_hammerintegeration
!  fea_shapetriangle
!  fea_elementmatrixD
!  fea_elementmatrixA
!  fea_swap
!  fea_elementmatrixB
!  fea_elementjacobi
!  fea_decomposition
!  fea_backsubstitution
!  fea_dynamicsetup
!  fea_contactsetup
!  fea_dynamicsolver
!  fea_calload
!  fea_compute_tri_area
!  fea_compute_normal_3D
!  fea_compute_normal_2D
!  fea_getpressure_SSM
!  fea_getpressure_GCM
!  fea_collisondetction
!  fea_contact3D
!  fea_contact2D
!  fea_contactsurf
!  fea_contactforce3D
!  fea_contactside
!  fea_contactforce2D
!  fea_setnewmark
!  fea_dynamicsolutionoutput
!  fea_open_probe_files
!  fea_read_probe_inputs
!  fea_readrestart
!  fea_writerestart
!  fea_readMfile
!  fea_writeMfile
!  fea_readkeffectfile
!  fea_writekeffectfile
!  fea_body_marker_pressure
!  fea_compute_initial_acceleration
!  fea_writeafile
!  fea_readafile
!------------------------------------------------------------
# define G2LI(i)      i-(myIs-1)                                !SER_TO_PAR. QX. CH32
# define G2LJ(j)      j-(myJs-1)

 SUBROUTINE fea_initial	
   USE finite_element_parameters
   USE flow_parameters

   IMPLICIT NONE

   INTEGER             :: ibody
	
   CALL fea_openfiles
   CALL fea_readinput
   CALL fea_allocatememory
   CALL fea_readmesh
   CALL fea_readconload
     
   IF (nread == 0) THEN
     CALL fea_getcorresponancetable
   ELSE 
     CALL fea_readcorresponancetable 
   ENDIF

   DO ibody = 1, fea_nbody
     CALL fea_assemble(ibody)
     CALL fea_dynamicsetup(ibody)
     IF (fea_contactprobflag==1) THEN
        CALL fea_contactsetup(ibody)
     ENDIF
   END DO
 
   IF(nread == 1) CALL fea_readrestart()
   IF(monitorON)  THEN
     CALL fea_read_probe_inputs()
     CALL fea_open_probe_files()
   ENDIF
 ENDSUBROUTINE fea_initial 
!----------------------------------------------------------

 SUBROUTINE fea_readinput

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
  
   IMPLICIT NONE

   INTEGER :: n, i, j

   READ(fea_inputfile,*)
   READ(fea_inputfile,*)
   READ(fea_inputfile,*) fea_nbody, fea_nmat, fea_nd             
   READ(fea_inputfile,*) 
   READ(fea_inputfile,*) fea_nnode, fea_probtype, fea_solvertype,&
                         fea_grav, fea_gravkey
   READ(fea_inputfile,*) 
   READ(fea_inputfile,*) fea_omegaload, fea_nfilterstart, &
                         fea_nfilterend
   READ(fea_inputfile,*) 
   READ(fea_inputfile,*) fea_coeff_penalty, fea_dtratio

   READ(fea_inputfile,*)
   READ(fea_inputfile,*) fea_readM, fea_readkeffect, fea_readinita, fea_contact_model

   ALLOCATE(fea_nelement(fea_nbody), &
            fea_npt(fea_nbody), &
            fea_mband(fea_nbody))
   ALLOCATE(fea_nfixpt(fea_nbody), &
            fea_nload(fea_nbody), &
            fea_muv(fea_nbody), &
            fea_neignv(fea_nbody))
   ALLOCATE(fea_icontactplane(fea_nbody), &
            fea_icontactdir(fea_nbody),&
            fea_icontactdirnorm(fea_nbody))
  ALLOCATE(fea_vmati(fea_nmat, 7))
   DO n = 1, fea_nbody
     READ(fea_inputfile,*)
     READ(fea_inputfile,*)
     READ(fea_inputfile,*) fea_nelement(n), fea_npt(n), fea_mband(n)
     READ(fea_inputfile,*)
     READ(fea_inputfile,*) fea_nfixpt(n), fea_nload(n)
     READ(fea_inputfile,*)
     READ(fea_inputfile,*) fea_muv(n), fea_ntime, fea_dt, fea_freq,&
                           fea_cc1, fea_cc2, fea_beta, fea_gamma 
     READ(fea_inputfile,*)
     READ(fea_inputfile,*) fea_neignv(n)
     READ(fea_inputfile,*)
     READ(fea_inputfile,*) fea_contactprobflag, fea_icontactplane(n),&
                           fea_icontactdir(n), fea_icontactdirnorm(n)
   END DO 
  
   fea_nelementmax = 0
   fea_nptmax      = 0
   fea_mbandmax    = 0
   fea_nfixptmax   = 0
   fea_nloadmax    = 0
   fea_neignvmax   = 0
   DO n = 1, fea_nbody
     IF(fea_nelementmax < fea_nelement(n)) fea_nelementmax = fea_nelement(n)
     IF(fea_nptmax      < fea_npt(n))      fea_nptmax      = fea_npt(n)
     IF(fea_mbandmax    < fea_mband(n))    fea_mbandmax    = fea_mband(n)
     IF(fea_nfixptmax   < fea_nfixpt(n))   fea_nfixptmax   = fea_nfixpt(n)
     IF(fea_nloadmax    < fea_nload(n))    fea_nloadmax    = fea_nload(n)
     IF(fea_neignvmax   < fea_neignv(n))   fea_neignvmax   = fea_neignv(n)
   END DO
   
   WRITE(*, *) 'Number of finite element body is', fea_nbody
   WRITE(*, *) 'Maxium number of nodes in a single element is', fea_nnode
   WRITE(*, *) 'Maxium number of elements among all of feabodies is',fea_nelementmax
   WRITE(*, *) 'Maxium number of nodes among all of feabodies is',fea_nptmax
   WRITE(*, *) 'Maxium matrix band width among all of feabodies is',fea_mbandmax
   WRITE(*, *) 'Maxium contraint nodes among all of feabodies is', fea_nfixptmax
   WRITE(*, *) 'Maxium load number among all of feabodies is', fea_nloadmax 
 
   WRITE(*, *) 'Total number of materials is', fea_nmat

   READ(fea_inputfile,*)
   READ(fea_inputfile,*)

   DO n = 1, fea_nmat
      READ(fea_inputfile,*) j, (fea_vmati(n, i), i = 1, 7) 
   END DO 
   READ(fea_inputfile,*)
   READ(fea_inputfile,*) fea_trans_dir


 END SUBROUTINE fea_readinput
!----------------------------------------------------------

 SUBROUTINE fea_allocatememory
   
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE flow_parameters

   IMPLICIT NONE

   ALLOCATE(fea_ielement(fea_nelementmax, 2+fea_nd+fea_nnode, fea_nbody))
   ALLOCATE(fea_cood(fea_nptmax, fea_nd, fea_nbody))
   
   IF(fea_nfixptmax > 0) ALLOCATE&
                         (fea_vfixed(fea_nfixptmax, fea_nd, fea_nbody), &
                          fea_ifixed(fea_nfixptmax, fea_nd+1, fea_nbody))
   IF(fea_nloadmax > 0)  ALLOCATE&
                         (fea_vload(fea_nloadmax, fea_nd, fea_nbody), &
                          fea_iload(fea_nloadmax, fea_nd+1, fea_nbody))
   
   ALLOCATE(fea_gmm(fea_nptmax*fea_nd, fea_mbandmax, fea_nbody))
   ALLOCATE(fea_gkm(fea_nptmax*fea_nd, fea_mbandmax, fea_nbody))
   ALLOCATE(fea_gcm(fea_nptmax*fea_nd, fea_mbandmax, fea_nbody))
   ALLOCATE(fea_keffect(fea_nptmax*fea_nd, fea_mbandmax, fea_nbody))
  
   ALLOCATE(fea_gu(fea_nptmax*fea_nd, fea_nbody), &
            fea_gp(fea_nptmax*fea_nd, fea_nbody), &
            fea_gp_old(fea_nptmax*fea_nd, fea_nbody))
   ALLOCATE(fea_d(fea_nptmax*fea_nd, fea_nbody), &
            fea_v(fea_nptmax*fea_nd, fea_nbody), &
            fea_a(fea_nptmax*fea_nd, fea_nbody))
   ALLOCATE(markertofeapointtable(nPtsMax, 2, nbody), &
            featomarkerpointtable(fea_nptmax,2, fea_nbody))
   ALLOCATE(fea_contactflag(fea_nbody), fea_ncontactpoint(fea_nbody), &
            fea_ncontactsurf(fea_nbody), fea_ncontactside(fea_nbody))
   ALLOCATE(fea_original_d(fea_nptmax*fea_nd, fea_nbody), &
            fea_original_v(fea_nptmax*fea_nd, fea_nbody), &
            fea_original_a(fea_nptmax*fea_nd, fea_nbody))
   ALLOCATE(fea_icontactpoint(fea_nptmax, fea_nbody), &
            fea_icontactsurf(fea_nloadmax, fea_nbody),&
            fea_icontactside(fea_nloadmax, fea_nbody))
   ALLOCATE(fea_penaltycoeff(fea_nbody))
   IF(fea_neignvmax > 0) &
      ALLOCATE (fea_veignv(fea_nptmax*fea_nd, fea_neignvmax, fea_nbody), &
                fea_ieignv(fea_neignvmax, fea_nbody))                     
  ALLOCATE(keffecttempt(fea_nptmax*fea_nd,fea_mbandmax))
  !initial arrays
  
  fea_gmm  = 0.0_CGREAL
  fea_gkm  = 0.0_CGREAL
  fea_gcm  = 0.0_CGREAL
  fea_keffect = 0.0_CGREAL
  fea_gu  = 0.0_CGREAL
  fea_gp = 0.0_CGREAL
  fea_gp_old = 0.0_CGREAL
  fea_d = 0.0_CGREAL
  fea_v = 0.0_CGREAL
  fea_a = 0.0_CGREAL
  markertofeapointtable = 0
  featomarkerpointtable = 0
  fea_contactflag = 0 
  fea_ncontactpoint =0
  fea_ncontactsurf = 0
  fea_ncontactside = 0

  fea_original_d = 0.0_CGREAL
  fea_original_v = 0.0_CGREAL
  fea_original_a = 0.0_CGREAL
  fea_icontactpoint = 0
  fea_icontactsurf = 0
 fea_icontactside = 0
 fea_penaltycoeff = 0.0_CGREAL  
 END SUBROUTINE fea_allocatememory
 !----------------------------------------------------------

 SUBROUTINE fea_deallocatememory
   
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   DEALLOCATE(fea_nelement, fea_npt, fea_mband)
   DEALLOCATE(fea_nfixpt, fea_nload, fea_muv)
   DEALLOCATE(fea_ielement)
   DEALLOCATE(fea_cood)
   IF(fea_nfixptmax > 0) DEALLOCATE(fea_vfixed, fea_ifixed)
   IF(fea_nloadmax > 0)  DEALLOCATE(fea_vload, fea_iload)
   DEALLOCATE(fea_gmm, fea_gkm, fea_gcm)
   DEALLOCATE(fea_gu, fea_gp)

 END SUBROUTINE fea_deallocatememory
!----------------------------------------------------------

 SUBROUTINE fea_openfiles
  
   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters

   IMPLICIT NONE

 ! Added for fea, open files  for flow structure interaction
   OPEN(fea_inputfile,      FILE='fea_input.dat',           STATUS='UNKNOWN')
   OPEN(fea_meshfile,       FILE='fea_mesh.dat',            STATUS='UNKNOWN')
   OPEN(fea_conloadfile,    FILE='fea_conload.dat',         STATUS='UNKNOWN')
   OPEN(fea_outputfile,     FILE='fea_output.dat',          STATUS='UNKNOWN')
   OPEN(fea_staticoutput,   FILE='fea_staticoutput.dat',    STATUS='UNKNOWN')
   OPEN(fea_initdisvel  ,   FILE='fea_initdisvel.dat',      STATUS='UNKNOWN')
   OPEN(fea_correstablefile,FILE='fea_correspond_table.dat',STATUS='UNKNOWN')
   OPEN(fea_restartfilein,  FILE='fea_restart_in.dat'      ,FORM='UNFORMATTED')
   OPEN(fea_Mfile,          FILE='fea_Matrix.dat',          FORM='UNFORMATTED',STATUS='UNKNOWN')
   OPEN(fea_keffectfile,    FILE='fea_keffect.dat',         FORM='UNFORMATTED',STATUS='UNKNOWN')
   OPEN(fea_afile,          FILE='fea_a.dat',               FORM='UNFORMATTED',STATUS='UNKNOWN') 
END SUBROUTINE fea_openfiles
!------------------------------------------------------------------------------------

 SUBROUTINE fea_readmesh

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER :: n, i, j

   DO n = 1, fea_nbody
     READ(fea_meshfile,*) 
     READ(fea_meshfile,*)
     DO i = 1, fea_npt(n)
       READ(fea_meshfile,*) (fea_cood(i, j, n), j = 1, fea_nd)
     END DO
     DO i = 1, fea_nelement(n)
       READ(fea_meshfile,*) (fea_ielement(i, j, n), j = 1, 2+fea_nd+fea_nnode)
     END DO
   END DO !n
  
 END SUBROUTINE fea_readmesh
!----------------------------------------------------------

 SUBROUTINE fea_readconload

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER :: n, i, j, ii
   
   DO n = 1, fea_nbody	
     READ(fea_conloadfile,*)

     IF(fea_nfixpt(n) /= 0) THEN
       READ(fea_conloadfile,*)
       READ(fea_conloadfile,*)
       DO i = 1, fea_nfixpt(n)
         READ(fea_conloadfile,*) ii, (fea_ifixed(i, j, n), j = 1, fea_nd+1), &
                                     (fea_vfixed(i, j, n), j = 1, fea_nd)
       END DO !i
     END IF

     IF(fea_nload(n) /= 0) THEN
       READ(fea_conloadfile,*)
       READ(fea_conloadfile,*)
       DO i = 1, fea_nload(n)
         READ(fea_conloadfile,*) ii, (fea_iload(i, j, n), j = 1, fea_nd+1),&
                                     (fea_vload(i, j, n), j = 1, fea_nd)
       END DO !i
     END IF
   END DO !n
    	
END SUBROUTINE fea_readconload
!----------------------------------------------------------

 SUBROUTINE fea_getcorresponancetable
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE flow_parameters
   USE boundary_arrays

   IMPLICIT NONE

   INTEGER :: i,  k, ibody, ifeabody, FLAG
   REAL(KIND=CGREAL) x1, x2, y1, y2, z1, z2, DIS, eps

   FLAG = 0
   eps = 1.0e-5_CGREAL

   DO ibody = 1,nbody
   DO i =nPtsBodyMarker(ibody), 1, -1
     x1 = xBodyMarker(ibody, i)
     y1 = yBodyMarker(ibody, i)
     z1 = zBodyMarker(ibody, i)

     IF (fea_nd == DIM_3D) THEN

       DO ifeabody = 1, fea_nbody
       DO k = 1, fea_npt(ifeabody)
         x2 = fea_cood(k, 1, ifeabody)
         y2 = fea_cood(k, 2, ifeabody)
         z2 = fea_cood(k, 3, ifeabody)
         DIS = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
         IF(DIS < eps) THEN
           markertofeapointtable(i, 1, ibody)    = k
           markertofeapointtable(i, 2, ibody)    = ifeabody
           featomarkerpointtable(k, 1, ifeabody) = i
           featomarkerpointtable(k, 2, ifeabody) = ibody
           FLAG = 1
           GOTO 1000
         ENDIF
       ENDDO ! k
       ENDDO  ! ifeabody

     ELSE IF (fea_nd == DIM_2D)THEN

       DO ifeabody = 1, fea_nbody
       DO k = 1, fea_npt(ifeabody)
         x2 = fea_cood(k, 1, ifeabody)
         y2 = fea_cood(k, 2, ifeabody)
         DIS = sqrt((x1 - x2)**2 + (y1 - y2)**2 )
         IF(DIS < eps) THEN
           markertofeapointtable(i, 1, ibody)    = k
           markertofeapointtable(i, 2, ibody)    = ifeabody
           featomarkerpointtable(k, 1, ifeabody) = i
           featomarkerpointtable(k, 2, ifeabody) = ibody
           FLAG = 1
           GOTO 1000
         ENDIF
       ENDDO ! k
       ENDDO  ! ifeabody

     ENDIF
1000 CONTINUE

     IF(FLAG == 0) THEN
       markertofeapointtable(i, 1, ibody)    = -1
       markertofeapointtable(i, 2, ibody)    = -1
     ENDIF

     FLAG = 0

   ENDDO !i
   ENDDO !ibody

   DO ibody = 1, nbody
   DO i = 1, nPtsBodyMarker(ibody)
     write(fea_correstablefile, *) markertofeapointtable(i, 1, ibody),&
                                   markertofeapointtable(i, 2, ibody)
   ENDDO
   ENDDO
   
   write (fea_correstablefile, *)
   DO ifeabody = 1, fea_nbody
   DO k = 1, fea_npt(ifeabody)
     WRITE (fea_correstablefile, *) featomarkerpointtable(k, 1, ifeabody),&
                                    featomarkerpointtable(k, 2, ifeabody)
   ENDDO
   ENDDO
 
 END SUBROUTINE fea_getcorresponancetable
!-------------------------------------------------------------

 SUBROUTINE fea_readcorresponancetable
   USE finite_element_parameters
   USE finite_element_arrays
   USE flow_parameters

   INTEGER :: ibody, i, k, ifeabody

   DO ibody = 1, nbody
   DO i = 1, nPtsBodyMarker(ibody)
     READ(fea_correstablefile, *) markertofeapointtable(i, 1, ibody),&
                                  markertofeapointtable(i, 2, ibody)
   ENDDO
   ENDDO

   READ(fea_correstablefile, *)

   DO ifeabody = 1, fea_nbody
   DO k = 1, fea_npt(ifeabody)
     READ(fea_correstablefile, *) featomarkerpointtable(k, 1, ifeabody),&
                                  featomarkerpointtable(k, 2, ifeabody)
   ENDDO
   ENDDO

 ENDSUBROUTINE
!-------------------------------------------------------------

 SUBROUTINE fea_assemble(ibody)
   
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   
   IMPLICIT NONE
   
   INTEGER , INTENT (IN) :: ibody
   
   INTEGER :: ie, i, j, k, ii, jj, kk, ih, iv, ihh, ivv
   INTEGER :: fea_node, fea_intx, fea_inty, fea_imati, fea_intz
   INTEGER :: fea_iel(fea_nnode)
   
   REAL(KIND=CGREAL) :: fea_vmae(7)
   REAL(KIND=CGREAL),  DIMENSION(:, :), ALLOCATABLE :: fea_xy ,fea_lmm,&
                                                       fea_lkm, fea_lcm
  
  IF(fea_readM == 0) THEN
   
   DO i = 1 , fea_nnode * fea_nd
   DO j = 1, fea_mband(ibody)
     fea_gkm(i, j, ibody) = 0.0_CGREAL
     fea_gmm(i, j, ibody) = 0.0_CGREAL
     fea_gcm(i, j, ibody) = 0.0_CGREAL
   END DO !j
   END DO !i
 
   DO ie = 1, fea_nelement(ibody)
     fea_node  = fea_ielement(ie, 1, ibody)
     fea_imati = fea_ielement(ie, 2, ibody)
     fea_intx  = fea_ielement(ie, 3, ibody)
     fea_inty  = fea_ielement(ie, 4, ibody)

     IF (fea_nd == 3) fea_intz  = fea_ielement(ie, 5, ibody)

     DO j = 1,7 
       fea_vmae(j) = fea_vmati(fea_imati, j)
     END DO !j 
    
     ALLOCATE(fea_xy (fea_node, fea_nd),&
              fea_lmm(fea_node*fea_nd, fea_node*fea_nd), &
              fea_lkm(fea_node*fea_nd, fea_node*fea_nd), &
              fea_lcm(fea_node*fea_nd, fea_node*fea_nd))

     DO i = 1, fea_node
       fea_iel(i) = fea_ielement(ie, i+fea_nd+2, ibody)
       DO j = 1, fea_nd
         fea_xy(i, j) = 0.0_CGREAL
         IF(fea_iel(i) > 0) fea_xy(i,j) = fea_cood(fea_iel(i), j, ibody)
       END DO !j
     END DO !i

!Compute element matrix fea_lmm and stiffiness matrix fea_lkm
     CALL fea_elementmatrix(ie, ibody, fea_xy, fea_iel, fea_lmm, &
                            fea_lkm, fea_lcm, fea_vmae, fea_node,&
                            fea_intx, fea_inty, fea_imati)

!Aseemeble element matrixs to form the global mass matrix fea_gmm
!and global stiffness matrix fea_gkm
     DO i = 1, fea_node
     DO j = 1, fea_node
     DO ii = 1, fea_nd
     DO jj = 1, fea_nd
       ih  = fea_nd * (i - 1) + ii
       iv  = fea_nd * (j - 1) + jj
       ihh = fea_nd * (fea_iel(i) - 1) + ii
       ivv = fea_nd * (fea_iel(j) - 1) + jj
       ivv = ivv - ihh + 1
       IF(ihh > 0 .and. ivv > 0) THEN
         fea_gmm(ihh, ivv, ibody) = fea_gmm(ihh, ivv, ibody)&
                                  + fea_lmm(ih, iv)
         fea_gkm(ihh, ivv, ibody) = fea_gkm(ihh, ivv, ibody)&
                                  + fea_lkm(ih, iv) 
         fea_gcm(ihh, ivv, ibody) = fea_gcm(ihh, ivv, ibody)&
                                  + fea_lcm(ih, iv) 
       END IF 
     END DO !jj
     END DO !ii	 
     END DO !j
     END DO !i

     DEALLOCATE(fea_xy, fea_lmm, fea_lkm, fea_lcm)

   END DO !ie

! write out the matrix gmm, gcm, gkm 
  CALL fea_writeMfile(ibody)

ELSE !fea_readM

! read in the matrix gmm, gcm ,gkm 
  CALL fea_readMfile(ibody)

ENDIF !fea_readM   

! Draw nodal constraint
   IF(fea_solvertype == 1) THEN

     DO i = 1, fea_nfixpt(ibody)
     DO j = 1, fea_nd
       IF(fea_ifixed(i, j+1, ibody) /= 0) THEN
         ii = fea_nd * (fea_ifixed(i, 1, ibody) - 1) + j
         fea_gp(ii, ibody) = fea_vfixed(i, j, ibody) &
                           * fea_gkm(ii, 1, ibody) * 1.0e20_CGREAL
         fea_gkm(ii, 1, ibody) = fea_gkm(ii, 1, ibody) * 1.0e20_CGREAL
         IF(fea_gkm(ii, 1, ibody) <= 1.0e-20_CGREAL) &
            fea_gkm(ii, 1, ibody) = 1.0e-20_CGREAL 
       END IF
     END DO !j
     END DO !i

   END IF

   IF(fea_solvertype > 1) THEN

     DO i = 1, fea_nfixpt(ibody)
     DO j = 1, fea_nd

       IF(fea_ifixed(i, j+1, ibody) /= 0) THEN
         ii = fea_nd * (fea_ifixed(i, 1, ibody) - 1) + j
         fea_gu(ii, ibody) = 0.0_CGREAL
         fea_gkm(ii, 1, ibody) = 1.0_CGREAL
         fea_gmm(ii, 1, ibody) = 1.0_CGREAL
         fea_gcm(ii, 1, ibody) = 1.0_CGREAL 
        IF(fea_solvertype == DYNAMIC_CHARACTER_INVERSE &
            .or. fea_solvertype == DYNAMIC_CHARACTER_SURFACE) THEN
           fea_gmm(ii, 1, ibody) = 0.0_CGREAL
         END IF
         DO k = 2, fea_mband(ibody)
           fea_gkm(ii, k, ibody) = 0.0_CGREAL
           fea_gmm(ii, k, ibody) = 0.0_CGREAL
           fea_gcm(ii, k, ibody) = 0.0_CGREAL
           IF(ii >= k) THEN
             fea_gkm(ii-k+1, k, ibody) = 0.0_CGREAL
             fea_gmm(ii-k+1, k, ibody) = 0.0_CGREAL
             fea_gcm(ii-k+1, k, ibody) = 0.0_CGREAL
           END IF
         END DO !k
       END IF

     END DO !j
     END DO !i
  
   END IF

 END SUBROUTINE fea_assemble
!----------------------------------------------------------
 SUBROUTINE fea_elementmatrix(ie, ibody, fea_xy, fea_iel, fea_lmm,  &
                             fea_lkm, fea_lcm, fea_vmae, fea_node,  &
                             fea_intx, fea_inty, fea_imati)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   
   IMPLICIT NONE

   INTEGER            :: ie, ibody, fea_node, fea_intx, fea_inty, &
                         fea_imati, fea_iel(fea_nnode)
   REAL(KIND=CGREAL)  :: fea_xy(fea_node, fea_nd), fea_vmae(7)
   REAL(KIND=CGREAL)  :: fea_lmm(fea_node*fea_nd, fea_node*fea_nd), &
                         fea_lkm(fea_node*fea_nd, fea_node*fea_nd),  & 
                         fea_lcm(fea_node*fea_nd, fea_node*fea_nd)
   
   INTEGER            :: i, j, k, ii, jj, kk, ii1, jj1
   REAL(KIND=CGREAL)  :: fea_wt, fea_xt, fea_yt, fea_zt, fea_wxy, &
                         fea_sj, fea_sr, fea_density
   REAL(KIND=CGREAL)  :: fea_vN((fea_nd+1)*3), &
                         fea_vdN(fea_nd+1, (fea_nd+1)*3), &
                         fea_vd0(fea_nd+1, (fea_nd+1)*3), &
                         fea_vd(6, 6),&
                         fea_va(6, 6)
   REAL(KIND=CGREAL)  :: fea_vB(6, 27), fea_vsg(6, fea_nnode*fea_nd),&
                         fea_vcg(6, fea_nnode*fea_nd)
  REAL(KIND=CGREAL) :: mu_tmp
   
! Form elastic matrix [d]
   CALL fea_elementmatrixD(fea_probtype, fea_vd, fea_vmae, fea_nd,fea_trans_dir)   
   CALL fea_elementmatrixA(fea_probtype, fea_va, fea_vmae)
   
   DO i = 1, fea_node*fea_nd
   DO j = 1, fea_node*fea_nd
     fea_lmm(i, j) = 0.0_CGREAL
     fea_lkm(i, j) = 0.0_CGREAL
     fea_lcm(i, j) = 0.0_CGREAL
   END DO !j
   END DO !i

   IF(fea_nd == 3) THEN
      IF(fea_node== 4) fea_inty = 1
   ENDIF

   IF (fea_nd == 2) THEN
     IF(fea_node== 3 .or. fea_node == 6) fea_inty = 1
   ENDIF

   DO i = 1, fea_intx
   DO j = 1, fea_inty
     IF(fea_nd == 3 .and. fea_node == 4) THEN
       CALL fea_hammerintegeration3D(fea_intx, i, fea_wt, &
                                     fea_xt, fea_yt, fea_zt, fea_wxy)
       CALL fea_shapetetrahedron(fea_node, fea_wt, fea_xt, fea_yt,&
                                 fea_zt,fea_iel, fea_vN, fea_vdN)
     ENDIF
	    
     IF (fea_nd == 2) THEN    
       IF(fea_node == 4 .or. fea_node ==8) THEN
         CALL fea_guassintegration(fea_intx, fea_inty,&
                                   i, j, fea_xt, fea_yt, fea_wxy)
         CALL fea_shapequadra8(fea_node, fea_xt, fea_yt, fea_iel,&
                               fea_vN, fea_vdN)
       END IF
       IF(fea_node == 9) THEN
         CALL fea_guassintegration(fea_intx, fea_inty, i, j,&
                                   fea_xt, fea_yt, fea_wxy)
         CALL fea_shapequadra9(fea_node, fea_xt, fea_yt, &
                               fea_iel, fea_vN, fea_vdN)
       END IF
       IF(fea_node == 3 .or. fea_node ==6) THEN
         CALL fea_hammerintegeration(fea_intx, i, fea_xt,&
                                     fea_yt, fea_zt, fea_wxy)
         CALL fea_shapetriangle(fea_node, fea_xt, fea_yt, &
                                fea_zt,fea_iel, fea_vN, fea_vdN)
       END IF
     ENDIF  

! Form the jacobi matrix at integartion points
     CALL fea_elementjacobi(fea_nd, fea_node, fea_xy, fea_vdN,&
                            fea_sj, fea_vd0)

     IF(fea_sj <= 0.0_CGREAL) THEN
       WRITE(*, 999) ie, i, j, fea_sj
999    format(/3x,' fea-sj <= 0 in element =', i4, 3x,'INTX=', i2,3x,&
              'INTY=', i2, 3x, 'SJ=', e11.4) 
     END IF
	   
! Form B matrix at integeration points	[B] 
     CALL fea_elementmatrixB(fea_probtype, fea_nd, fea_node, fea_xy,&
                            fea_vN, fea_vd0, fea_vB, fea_sr) 

! Form the element stress matrix [S] = [D] * [B]
     fea_vsg = 0.0_CGREAL
     DO ii = 1, 6
     DO jj = 1, fea_node*fea_nd
     DO kk = 1, 6
       fea_vsg(ii, jj) = fea_vsg(ii, jj) &
                       + fea_vd(ii, kk) * fea_vb(kk, jj)
     END DO !k
     END DO !jj
     END DO !ii

! Form elemental stiffiness matrix [k] = [b] * [s]
     DO ii = 1, fea_node*fea_nd
     DO jj = 1, fea_node*fea_nd
     DO kk = 1, 6
       fea_lkm(ii, jj) = fea_lkm(ii, jj) &
                       + fea_vB(kk, ii) * fea_vsg(kk, jj) &
                       * fea_wxy * fea_sj * fea_sr
     END DO !kk
     END DO !jj
     END DO !ii
	   
! Form the element stress matrix [C] = [A] * [B]
     fea_vcg = 0.0_CGREAL
     DO ii = 1, 6
     DO jj = 1, fea_node*fea_nd
     DO kk = 1, 6
       fea_vcg(ii, jj) = fea_vcg(ii, jj) &
                       + fea_va(ii, kk) * fea_vb(kk, jj)
     END DO !k
     END DO !jj
     END DO !ii

! Form elemental damping matrix [C] = [b] * [c]
     DO ii = 1, fea_node*fea_nd
     DO jj = 1, fea_node*fea_nd
     DO kk = 1, 6
       fea_lcm(ii, jj) = fea_lcm(ii, jj) &
                       + fea_vB(kk, ii) * fea_vcg(kk, jj) &
                       * fea_wxy * fea_sj * fea_sr
     END DO !kk
     END DO !jj
     END DO !ii	   
!debug
   ! mu_tmp =  0.5_CGREAL  * fea_vmae(2) / (1.0_CGREAL + fea_vmae(3))  
   ! DO ii = 1, fea_node*fea_nd   
   ! DO jj = 1, fea_node*fea_nd
   !   fea_lcm(ii, jj) = fea_lkm(ii,jj)/ mu_tmp * fea_vmae(7)
   ! END DO !jj
   ! END DO !ii

   
! Form elemental mass matrix [m]
     fea_density = fea_vmae(1)
     DO ii = 1, fea_node
     DO jj = 1, fea_node
     DO kk = 1, fea_nd
       ii1 = fea_nd * (ii - 1) + kk
       jj1 = fea_nd * (jj - 1) + kk
       fea_lmm(ii1, jj1) = fea_lmm(ii1, jj1)& 
                         + fea_density * fea_vN(ii) * fea_vN(jj)&
                         * fea_wxy * fea_sj * fea_sr
     END DO !kk
     END DO !jj
     END DO !ii

   END DO !j
   END DO !i

 END SUBROUTINE fea_elementmatrix
!----------------------------------------------------------

 SUBROUTINE fea_hammerintegeration3D(intx, i, wt, xt, yt, zt, wx)

   USE global_parameters

   IMPLICIT NONE

   INTEGER           :: intx,  i, j
   REAL(KIND=CGREAL) :: wt, xt, yt, zt, wx
   REAL(KIND=CGREAL) ::a, b
  
   SELECT CASE(intx)
     CASE(1)
   ! Integeration constants of one point
       wt = 1.0_CGREAL / 4.0_CGREAL
       xt = 1.0_CGREAL / 4.0_CGREAL
       yt = 1.0_CGREAL / 4.0_CGREAL
       zt = 1.0_CGREAL / 4.0_CGREAL
       wx     = 1.0_CGREAL/6.0_CGREAL
     CASE(4)
       a = 0.58541020_CGREAL
       b = 0.13819660_CGREAL
       wx = 0.25_CGREAL
       wx = wx / 6.0_CGREAL 
       SELECT CASE(i)
         CASE(1)
           wt = a
	   xt = b
	   yt = b
	   zt = b
         CASE(2)
	   wt = b
	   xt = a
	   yt = b
	   zt = b
         CASE(3)
	   wt = b
	   xt = b
	   yt = a
	   zt = b
         CASE(4)
	   wt = b
	   xt = b
	   yt = b
	   zt = a
       END SELECT !i
   END SELECT !	intx

 END SUBROUTINE fea_hammerintegeration3D
!----------------------------------------------------------

 SUBROUTINE fea_shapetetrahedron(node, wt, xt, yt, zt, iel, N, DN)

   USE global_parameters
   
   IMPLICIT NONE
   
   INTEGER           :: node, iel(node)
   REAL(KIND=CGREAL) :: xt, yt, zt, wt, N(12), DN(4, 12)
   INTEGER           :: i

   N  = 0.0_CGREAL
   DN = 0.0_CGREAL

! Set function value for 4-node tetrahedron element
   N(1)    = wt
   N(2)    = xt
   N(3)    = yt
   N(4)    = zt
   DN(1,1) = 1.0_CGREAL
   DN(2,2) = 1.0_CGREAL
   DN(3,3) = 1.0_CGREAL
   DN(4,4) = 1.0_CGREAL

   DO i = 1, node
     DN(1, i) = DN(1, i) - DN(4, i)
     DN(2, i) = DN(2, i) - DN(4, i)
     DN(3, i) = DN(3, i) - DN(4, i)
   END DO !i

 ENDSUBROUTINE fea_shapetetrahedron
!----------------------------------------------------------------------------------

 SUBROUTINE fea_guassintegration(intx, inty, i, j, xt, yt, wxy)

   USE global_parameters

   IMPLICIT NONE

   INTEGER           :: intx, inty, i, j
   REAL(KIND=CGREAL) :: xt, yt, wxy

   REAL(KIND=CGREAL) :: gxy(3, 3), wg(3, 3)

! Guass integration constants for 1, 2 and 3 pints
   gxy(1, 1) =  0.0_CGREAL
   wg(1, 1)  =  2.0_CGREAL
   gxy(1, 2) = -0.577350269189626_CGREAL
   gxy(2, 2) =  0.577350269189626_CGREAL
   wg(1, 2)  =  1.0_CGREAL
   wg(2, 2)  =  1.0_CGREAL
   gxy(1, 3) = -0.774596669241483_CGREAL
   gxy(2, 3) =  0.0_CGREAL
   gxy(3, 3) =  0.774596669241483_CGREAL
   wg(1, 3)  =  0.555555555555556_CGREAL
   wg(2, 3)  =  0.888888888888889_CGREAL
   wg(3, 3)  =  0.555555555555556_CGREAL

! Get parameters of integration pint
   xt  = gxy(i, intx)
   yt  = gxy(j, inty)
   wxy = wg(i, intx) * wg(j, inty)
    
 ENDSUBROUTINE
!----------------------------------------------------------

 SUBROUTINE fea_shapequadra8(node, xt, yt, iel, N, DN)
  
   USE global_parameters

   IMPLICIT NONE
 
   INTEGER           :: node, iel(node)
   REAL(KIND=CGREAL) :: xt, yt, N(9), DN(3,9)
   INTEGER           :: i

   N  = 0.0_CGREAL
   DN = 0.0_CGREAL

! Set function value for quadradic element of 4 nodes
   N(1) = (1 + xt) * (1 + yt) / 4
   N(2) = (1 - xt) * (1 + yt) / 4
   N(3) = (1 - xt) * (1 - yt) / 4
   N(4) = (1 + xt) * (1 - yt) / 4
   
   DN(1, 1) =  (1 + yt) / 4
   DN(1, 2) = -(1 + yt) / 4
   DN(1, 3) = -(1 - yt) / 4
   DN(1, 4) =  (1 - yt) / 4
   DN(2, 1) =  (1 + xt) / 4
   DN(2, 2) =  (1 - xt) / 4
   DN(2, 3) = -(1 - xt) / 4
   DN(2, 4) = -(1 + xt) / 4
! Set function value for quadradic element of 8 nodes
   IF(node == 8) THEN
     N(5) = (1- xt * xt) * (1 + yt) / 2
     N(6) = (1- yt * yt) * (1 - xt) / 2
     N(7) = (1- xt * xt) * (1 - yt) / 2  
     N(8) = (1- yt * yt) * (1 + xt) / 2

     DN(1, 5) = (-2 * xt) * (1 + yt) / 2
     DN(1, 6) = (1 - yt * yt) * (-1) / 2
     DN(1, 7) = (-2 * xt) * (1 - yt) / 2
     DN(1, 8) = (1 - yt * yt) * (+1) / 2
     DN(2, 5) = (1 - xt * xt) * (+1) / 2
     DN(2, 6) = (-2 * yt) * (1 - xt) / 2
     DN(2, 7) = (1 - xt * xt) * (-1) / 2
     DN(2, 8) = (-2 * yt) * (1 + xt) / 2

     DO i = 1, 4
       IF(iel(4+i) == 0) THEN
         N(4+i)     = 0.0_CGREAL
         DN(1, 4+i) = 0.0_CGREAL
         DN(2, 4+i) = 0.0_CGREAL
       END IF
     END DO !i

     N(1) = N(1) - (N(5) + N(8)) / 2
     N(2) = N(2) - (N(5) + N(6)) / 2
     N(3) = N(3) - (N(6) + N(7)) / 2
     N(4) = N(4) - (N(7) + N(8)) / 2

     DO i = 1, 2
       DN(i, 1) = DN(i, 1) - (DN(i, 5) + DN(i, 8)) / 2
       DN(i, 2) = DN(i, 3) - (DN(i, 5) + DN(i, 6)) / 2
       DN(i, 3) = DN(i, 3) - (DN(i, 6) + DN(i, 7)) / 2
       DN(i, 4) = DN(i, 4) - (DN(i, 7) + DN(i, 8)) / 2
     END DO  !i  

   ENDIF

 END SUBROUTINE fea_shapequadra8
!----------------------------------------------------------

 SUBROUTINE fea_shapequadra9(node, xt, yt, iel, N, DN)

   USE global_parameters

   IMPLICIT NONE
 
   INTEGER           :: node, iel(node)
   REAL(KIND=CGREAL) :: xt, yt, N(9), DN(3,9)
   INTEGER           :: i, iix, iiy, ix(9), iy(9)
   REAL(KIND=CGREAL) :: vl0x(3), vl1x(3), vl0y(3), vl1y(3)

   ix(1) =  1
   ix(2) = -1
   ix(3) = -1
   ix(4) =  1
   ix(5) =  0
   ix(6) = -1
   ix(7) =  0
   ix(8) =  1
   ix(9) =  0


   ix(1) =  1
   ix(2) =  1
   ix(3) = -1
   ix(4) = -1
   ix(5) =  1
   ix(6) =  0
   ix(7) = -1
   ix(8) =  0
   ix(9) =  0

   N  = 0.0_CGREAL
   DN = 0.0_CGREAL
   
! Set the value of lagrange function and its derivatives
   vl0x(1) = 0.5_CGREAL * xt * (xt - 1.0_CGREAL)
   vl0x(2) = 1.0_CGREAL - xt * xt    
   vl0x(3) = 0.5_CGREAL * xt * (xt + 1.0_CGREAL)
   vl1x(1) = xt - 0.5_CGREAL
   vl1x(2) = -2.0_CGREAL * xt
   vl1x(3) = xt + 0.5_CGREAL
   vl0y(1) = 0.5_CGREAL * yt * (yt - 1.0_CGREAL)
   vl0y(2) = 1.0_CGREAL - yt * yt
   vl0y(3) = 0.5_CGREAL * yt * (yt + 1.0_CGREAL)
   vl1y(1) = yt - 0.5_CGREAL
   vl1y(2) = -2.0_CGREAL * yt
   vl1y(3) = yt + 0.5_CGREAL
! Set shape function for quadratic element of 9 nodes
   
   DO i = 1, 9
     iix      = ix(i) + 2
     iiy      = iy(i) + 2
     N(i)     = vl0x(iix) * vl0y(iiy)
     DN(1, i) = vl1x(iix) * vl0y(iiy)
     DN(2, i) = vl0x(iix) * vl1y(iiy)
   END DO !i
   
 END SUBROUTINE fea_shapequadra9
!----------------------------------------------------------

 SUBROUTINE fea_hammerintegeration(intx, i, xt, yt, zt, wx)

   USE global_parameters

   IMPLICIT NONE

   INTEGER           :: intx,  i, j
   REAL(KIND=CGREAL) :: xt, yt, zt, wx

   INTEGER           :: indextable(7)
   REAL(KIND=CGREAL) :: A1, B1, A2, B2
   REAL(KIND=CGREAL) :: hxy(3, 5), wh(5) 
! hxy - coordinate of hammer integration point
! wh  - weight of hammer integeration point 

   indextable(1) = 1
   indextable(3) = 2
   indextable(4) = 3
   indextable(7) = 4

! Integeration constants of one point
   hxy(1, 1) = 1.0_CGREAL / 3.0_CGREAL
   hxy(2, 1) = 1.0_CGREAL / 3.0_CGREAL
   hxy(3, 1) = 1.0_CGREAL / 3.0_CGREAL
   wh(1)     = 1.0_CGREAL

! Integeration constants of 3 points
   hxy(1, 2) = 2.0_CGREAL / 3.0_CGREAL
   hxy(2, 2) = 1.0_CGREAL / 6.0_CGREAL
   hxy(3, 2) = 1.0_CGREAL / 6.0_CGREAL
   wh(2)     = 1.0_CGREAL / 3.0_CGREAL
    
! Integeration constants of 4 points
   hxy(1, 3) = 0.6_CGREAL 
   hxy(2, 3) = 0.2_CGREAL
   hxy(3, 3) = 0.2_CGREAL
   wh(3)     = 25.0_CGREAL / 48.0_CGREAL

!Integration constants of 7 points 
   A1 = 0.0597158717_CGREAL
   B1 = 0.4701420641_CGREAL
   A2 = 0.7974269853_CGREAL
   B2 = 0.1012865073_CGREAL
   hxy(1,4) = A1
   hxy(2,4) = B1
   hxy(3,4) = B1
   wh(4) = 0.1323941527_CGREAL

   hxy(1,5) = A2
   hxy(2,5) = B2
   hxy(3,5) = B2
   wh(5) = 0.1259391805_CGREAL

! Get parameters of integration points
   xt = hxy(i+0-(i-1)/3*3, indextable(intx))
   yt = hxy(i+2-(i+1)/3*3, indextable(intx))
   zt = hxy(i+1-(i+0)/3*3, indextable(intx))
   wx = wh(indextable(intx)) / 2.0_CGREAL

   IF(intx == 7 .and. i >= 4) THEN
     j = i - 3
     xt = hxy(j+0-(j-1)/3*3, 5)
     yt = hxy(j+2-(j+1)/3*3, 5)
     zt = hxy(j+1-(j+0)/3*3, 5)
     wx = wh(5) / 2.0_CGREAL
   END IF

   IF(intx == 4 .and. i == 4) THEN
    xt = hxy(i+0-(i-1)/3*3, 1)
    yt = hxy(i+1-(i+0)/3*3, 1)
    zt = hxy(i+2-(i+1)/3*3, 1)
    wx = -27.0_CGREAL / 48.0_CGREAL / 2.0_CGREAL
   END IF

   IF(intx == 7 .and. i == 7) THEN
	 xt = hxy(1, 1)
     yt = hxy(2, 1)
     zt = hxy(3, 1)
     wx = 0.9_CGREAL / 0.8_CGREAL
   END IF

 END SUBROUTINE fea_hammerintegeration
!----------------------------------------------------------

 SUBROUTINE fea_shapetriangle(node, xt, yt, zt, iel, N, DN)

   USE global_parameters
   
   IMPLICIT NONE
   
   INTEGER           :: node, iel(node)
   REAL(KIND=CGREAL) :: xt, yt, zt, N(9), DN(3, 9)
   INTEGER           :: i

   N  = 0.0_CGREAL
   DN = 0.0_CGREAL

! Set function value for 3-node trianlge element
   N(1)    = xt
   N(2)    = yt
   N(3)    = zt
   DN(1,1) = 1.0_CGREAL
   DN(2,2) = 1.0_CGREAL
   DN(3,3) = 1.0_CGREAL

! Set function value for triangle element of 4-6 nodes
   IF(node == 6) THEN
     N(4)     = 4 * xt * yt
     N(5)     = 4 * yt * zt
     N(6)     = 4 * zt * xt
     DN(1, 4) = 4 * yt
     DN(1, 6) = 4 * zt
     DN(2, 4) = 4 * xt
     DN(2, 5) = 4 * zt
     DN(3, 5) = 4 * yt
     DN(3, 6) = 4 * xt

     DO i = 1, 3
       IF(iel(3+i) == 0) THEN
         N(3+i)    = 0.0_CGREAL
         DN(1,3+i) = 0.0_CGREAL
         DN(2,3+i) = 0.0_CGREAL
         DN(3,3+i) = 0.0_CGREAL 
       END IF
     END DO !i

     N(1) = N(1) - (N(4) - N(6)) / 2
     N(2) = N(2) - (N(4) + N(5)) / 2
     N(3) = N(3) - (N(5) + N(6)) / 2

     DO i = 1, 3
       DN(i, 1) = DN(i, 1) - (DN(i, 4) + DN(i, 6)) / 2
       DN(i, 2) = DN(i, 2) - (DN(i, 4) + DN(i, 5)) / 2
       DN(i, 3) = DN(i, 3) - (DN(i, 5) + DN(i, 6)) / 2
     END DO !i
   END IF

   DO i = 1, node
     DN(1, i) = DN(1, i) - DN(3, i)
     DN(2, i) = DN(2, i) - DN(3, i)
   END DO !i

 END SUBROUTINE fea_shapetriangle
!----------------------------------------------------------
 SUBROUTINE fea_elementmatrixD(typeofprob, vd, vmae, fea_nd, fea_trans_dir)
   
   USE global_parameters 

   IMPLICIT NONE
 
   INTEGER, INTENT(IN)  :: typeofprob, fea_nd, fea_trans_dir
   REAL(KIND=CGREAL) :: vd(6, 6), vmae(7)
   REAL(KIND=CGREAL) :: ep, mup, ez, mupz, muzp, gzp, delta, D0, tempt1, tempt2
   INTEGER :: i, j
   ep   = vmae(2)
   mup  = vmae(3)
   ez   = vmae(4)
   muzp = vmae(5)
   gzp  = vmae(6)

                            ! According to compliance condition muzp 
   mupz = muzp * ep / ez                       ! muzp / Ez = mupz / Ep
   delta = (1 + mup) * (1 - mup - 2 * mupz * muzp) / (ep**2 * ez)
   VD  = 0.0_CGREAL

   SELECT CASE(typeofprob)
     CASE(PLANE_STRESS)
       D0 = ep / (1 - mup **2)
       VD(1, 1) = D0
       VD(2, 2) = D0
       VD(3, 3) = D0 * (1 - mup) / 2
       VD(1, 2) = D0 * mup
       VD(2, 1) = D0 * mup
     CASE(PLANE_STRAIN)
       VD(1, 1) = (1 - mupz * muzp) / (ep * ez * delta)
       VD(2, 2) = VD(1, 1)
       VD(3, 3) = ep / (1 + mup) /2.0_CGREAL
       VD(1, 2) = (mup + mupz * muzp) / (ep * ez * delta) 
       VD(2, 1) = VD(1, 2)
     CASE(GENERAL_3DBODY)
       IF (1==1) THEN
         VD(1, 1) = (1-mupz*muzp)/(ep*ez*delta)
         VD(1, 2) = (mup+muzp*mupz) /(ep*ez*delta)
         VD(1, 3) = (muzp+mup*muzp)/(ep*ez*delta)
         VD(2, 1) = (mup+mupz*muzp)/(ep*ez*delta)
         VD(2, 2) = (1-mupz*muzp)/ (ep*ez*delta)
         VD(2, 3) = (muzp+muzp*mup)/(ep*ez*delta)
         VD(3, 1) = (mupz+mup*mupz)/(ep**2*delta)
         VD(3, 2) = (mupz+mup*mupz)/(ep**2*delta)
         VD(3, 3) = (1-mup**2)/(ep**2*delta)
         VD(4, 4) = gzp
         VD(5, 5) = gzp
         VD(6, 6) = ep/(2*(1+mup))
       ENDIF

       SELECT CASE(fea_trans_dir)
        CASE(ICOORD)
         DO i = 1, 6
          CALL fea_swap(VD(i, 1), VD(i, 3))
          CALL fea_swap(VD(i, 4), VD(i, 6))
         ENDDO
         DO j = 1, 6
          CALL fea_swap(VD(1, j), VD(3, j))
          CALL fea_swap(VD(4, j), VD(6, j))
         ENDDO
        CASE(JCOORD)
         DO i = 1, 6
          CALL fea_swap(VD(i, 2), VD(i, 3))
          CALL fea_swap(VD(i, 5), VD(i, 6))
         ENDDO
         DO j = 1, 6
          CALL fea_swap(VD(2, j), VD(3, j))
          CALL fea_swap(VD(5, j), VD(6, j))
         ENDDO
        CASE(KCOORD)
         CONTINUE
       ENDSELECT

       IF(1 == 0) THEN 
	 D0 = ep * (1 - mup) / ((1 + mup) *(1 - 2 * mup))
         tempt1 = mup / (1 - mup)
         tempt2 = 0.5_CGREAL * (1 - 2 *mup) / (1 - mup)
         VD(1, 1) = D0
         VD(2, 2) = D0
         VD(3, 3) = D0
         VD(4, 4) = D0 * tempt2
         VD(5, 5) = D0 * tempt2
         VD(6, 6) = D0 * tempt2
         VD(1, 2) = D0 * tempt1
         VD(2, 1) = D0 * tempt1
         VD(1, 3) = D0 * tempt1
         VD(3, 1) = D0 * tempt1
         VD(2, 3) = D0 * tempt1
         VD(3, 2) = D0 * tempt1
       ENDIF	
   END SELECT
   
 END SUBROUTINE fea_elementmatrixD
!----------------------------------------------------------

 SUBROUTINE fea_elementmatrixA(typeofprob, vA, vmae)
   
   USE global_parameters 
   USE flow_parameters
   IMPLICIT NONE
 
   INTEGER           :: typeofprob
   REAL(KIND=CGREAL) :: va(6, 6), vmae(7)
   REAL :: eta 
   
   eta   = vmae(7)
   va = 0.0_CGREAL
   IF(ndim == DIM_2D) THEN
    va(1, 1) = eta
    va(2, 2) = eta
    va(3, 3) = eta/2.0_CGREAL 
   ELSE
    Va(1,1) = eta
    va(2,2) = eta
    va(3,3) = eta
    va(4,4) = eta/2.0_CGREAL
    va(5,5) = eta/2.0_CGREAL
    va(6,6) = eta/2.0_CGREAL 
   ENDIF
 

 END SUBROUTINE fea_elementmatrixA
 !-------------------------------------------------------------------------
 
 SUBROUTINE fea_swap(var1, var2)
  USE global_parameters
  IMPLICIT NONE
  REAL(KIND=CGREAL) :: var1, var2
  REAL(KIND=CGREAL) :: temp
  temp = var1
  var1 = var2
  var2 = temp
 ENDSUBROUTINE fea_swap
!----------------------------------------------------------

 SUBROUTINE fea_elementmatrixB(typeofprob,NF, MND, XY, N, D0, B, sr)
   
   USE global_parameters 

   IMPLICIT NONE

 
   INTEGER           :: typeofprob, MND, NF
   REAL(KIND=CGREAL) :: XY(MND, NF), N((NF+1)*3), D0(NF+1, (NF+1)*3), B(6,27), sr
				        !XY-NODAL COORDINATES OF ELEMENT
				        !N - SHAPE FUNCTION
				        !D0 - GLOBAL DERIVATIVE OF N
				        !B -ELEMENT STRAIN MATRIX

   INTEGER           :: i, j

   B = 0.0_CGREAL

   SELECT CASE(typeofprob)
     CASE(PLANE_STRESS : PLANE_STRAIN)
       sr = 1.0_CGREAL
       DO i = 1, 9
         j = (i -1) * 2 
         B(1, j+1) = D0(1, i)
         B(2, j+2) = D0(2, i)
         B(3, j+1) = D0(2, i)
         B(3, j+2) = D0(1, i)
       END DO !i
     CASE(GENERAL_3DBODY)
       sr = 1.0_CGREAL 
       DO i = 1, 4 ! four node tetrahedron only !!!
         j = (i - 1) * 3
          j = (i - 1) * 3
		 B(1, j+1) = D0(1, i)
         B(2, j+2) = D0(2, i)
         B(3, j+3) = D0(3, i)
         B(4, j+2) = D0(3, i)
         B(4, j+3) = D0(2, i)
         B(5, j+1) = D0(3, i)
         B(5, J+3) = D0(1, i)
         B(6, j+1) = D0(2, i)
         B(6, j+2) = D0(1, i)
       ENDDO
    END SELECT

 END SUBROUTINE fea_elementmatrixB
!----------------------------------------------------------

 SUBROUTINE fea_elementjacobi(NF, MND, XY, DN, sj, VD0)
   
   USE global_parameters 

   IMPLICIT NONE

   INTEGER             :: MND, NF
   REAL(KIND = CGREAL) :: XY(MND, NF), DN(NF+1, (NF+1)*3), &
                          VD0(NF+1, (NF+1)*3), SJ
   INTEGER             :: ii, jj, kk
   REAL(KIND = CGREAL) :: VJJ(NF,NF), VJ1(NF,NF)

! Form jacobi matrix
   DO ii = 1, NF
   DO jj = 1, NF
     VJJ(ii,jj) = 0.0_CGREAL
     DO  kk = 1, MND
       VJJ(ii, jj) = VJJ(ii, jj) + DN(ii, kk) * XY(kk, jj)
     END DO !kk
   END DO !jj
   END DO !ii
! Calculate the determination of jacobi matrix
 
  IF(NF ==3) THEN
    sj =      VJJ(1, 1) * (VJJ(2,2) * VJJ(3,3) - VJJ(2,3) * VJJ(3,2))
    sj = sj - VJJ(1, 2) * (VJJ(2,1) * VJJ(3,3) - VJJ(2,3) * VJJ(3,1))
    sj = sj + VJJ(1, 3) * (VJJ(2,1) * VJJ(3,2) - VJJ(2,2) * VJJ(3,1))
    VJ1(1, 1) = (VJJ(2,2) * VJJ(3,3) - VJJ(2,3) * VJJ(3,2)) / sj
    VJ1(1, 2) = (VJJ(1,3) * VJJ(3,2) - VJJ(1,2) * VJJ(3,3)) / sj
    VJ1(1, 3) = (VJJ(1,2) * VJJ(2,3) - VJJ(1,3) * VJJ(2,2)) / sj
    VJ1(2, 1) = (VJJ(2,3) * VJJ(3,1) - VJJ(2,1) * VJJ(3,3)) / sj
    VJ1(2, 2) = (VJJ(1,1) * VJJ(3,3) - VJJ(1,3) * VJJ(3,1)) / sj
    VJ1(2, 3) = (VJJ(1,3) * VJJ(2,1) - VJJ(1,1) * VJJ(2,3)) / sj
    VJ1(3, 1) = (VJJ(2,1) * VJJ(3,2) - VJJ(2,2) * VJJ(3,1)) / sj
    VJ1(3, 2) = (VJJ(1,2) * VJJ(3,1) - VJJ(1,1) * VJJ(3,2)) / sj
    VJ1(3, 3) = (VJJ(1,1) * VJJ(2,2) - VJJ(1,2) * VJJ(2,1)) / sj

 !  Compute cartesian derivatives of shape functions
    DO ii = 1,NF
    DO jj = 1, 4 ! for 4 nodes only temperary
      VD0(ii, jj) = 0.0_CGREAL
      DO kk = 1, NF
        VD0(ii, jj) = VD0(ii, jj) + VJ1(ii, kk) * DN(kk, jj)
      END DO !kk
    END DO !jj
    END DO !ii
  END IF
  
  IF (NF == 2) THEN
    sj = VJJ(1,1) * VJJ(2,2) - VJJ(1,2) * VJJ(2,1)
    VJ1(1, 1) = +VJJ(2, 2) / sj
    VJ1(1, 2) = -VJJ(1, 2) / sj
    VJ1(2, 1) = -VJJ(2, 1) / sj
    VJ1(2, 2) = +VJJ(1, 1) / sj
    
! Compute cartesian derivatives of shape functions
    DO ii = 1,NF
    DO jj = 1,9
      VD0(ii, jj) = 0.0_CGREAL
      DO kk = 1, NF
        VD0(ii, jj) = VD0(ii, jj) + VJ1(ii, kk) * DN(kk, jj)
      END DO !kk
    END DO !jj
    END DO !ii
  ENDIF

 END SUBROUTINE
!----------------------------------------------------------

 SUBROUTINE fea_decomposition(numpt2, mband, array)

   USE global_parameters

   IMPLICIT NONE

   INTEGER           :: numpt2, mband
   REAL(KIND=CGREAL) :: array(numpt2, mband)

   INTEGER           :: i, j, mj, m
   
   DO j = 1, numpt2
     mj = j - mband + 1
     IF(mj < 1) mj = 1
     DO i = mj, j
     DO m = mj, i-1
       array(i, j-i+1) = array(i, j-i+1) &
                      - array(m, j-m+1) * array(m, i-m+1) / array(m, 1)
     END DO !m
     END DO !i
   END DO !j
   	  
 END SUBROUTINE  fea_decomposition
!----------------------------------------------------------

 SUBROUTINE fea_backsubstitution(numpt2, mband, array1, array2)
  
   USE global_parameters

   IMPLICIT NONE

   INTEGER           :: numpt2, mband
   REAL(KIND=CGREAL) :: array1(numpt2, mband), array2(numpt2) 

   INTEGER           :: i, j, mi, m
   
   DO i = 1, numpt2
     mi = i - mband + 1
     IF(mi < 1) mi = 1
     DO m = mi, i-1
       array2(i) = array2(i)&
                 - array1(m, i-m+1) * array2(m) / array1(m, 1)
     END DO !m
   END DO !i

   DO i = numpt2, 1, -1
     mi = i + mband - 1
     IF(mi > numpt2) mi = numpt2
     DO j = i+1, mi
       array2(i) = array2(i) - array1(i, j-i+1) * array2(j)
     END DO !j
     array2(i) = array2(i) / array1(i, 1)
   END DO !i
   
 END SUBROUTINE fea_backsubstitution
!----------------------------------------------------------

 SUBROUTINE fea_dynamicsetup(ibody)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER           :: ibody

   INTEGER           :: i, j, ii, ik, n
   !REAL(KIND=CGREAL) :: tempt(fea_npt(ibody)*fea_nd, fea_mband(ibody)),&
   !                     tempta(fea_npt(ibody)*fea_nd)

! Form Damping matrix fea_gcm
!  DO i = 1, fea_npt(ibody)*fea_nd
!   DO j = 1, fea_mband(ibody)
! fea_gcm(i, j, ibody) = fea_cc1 * fea_gmm(i, j, ibody) &
!                          + fea_cc2 * fea_gkm(i, j, ibody)
!  END DO !j
!  END DO !i

   DO i = 1, fea_npt(ibody)*fea_nd
     fea_d(i, ibody) = 0.0_CGREAL
     fea_v(i, ibody) = 0.0_CGREAL
     fea_a(i, ibody) = 0.0_CGREAL
   END DO !i
!assume at the begining thre is no acceleration occuring, this part has been commmented out to save computational time

IF(1==0) THEN
   ALLOCATE(tempt(fea_npt(ibody)*fea_nd, fea_mband(ibody)),tempta(fea_npt(ibody)*fea_nd))
   
   IF(fea_muv(ibody) == 1 .OR. fea_muv(ibody) == 3) THEN
     READ(fea_initdisvel, *)
     DO i = 1, fea_npt(ibody)*fea_nd
       READ(fea_initdisvel, *) ii, fea_d(i, ibody)
     END DO !i
   END IF
           
   IF(fea_muv(ibody) == 2 .OR. fea_muv(ibody) == 3) THEN
     READ(fea_initdisvel, *)
     DO i = 1, fea_npt(ibody)*fea_nd
       READ(fea_initdisvel, *) ii, fea_v(i, ibody) 
     END DO !i
   END IF
! Compute initial acceleration
   DO i = 1, fea_npt(ibody)*fea_nd
     tempta(i) = fea_gp(i, ibody)
   END DO !i
   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i - fea_mband(ibody) + 1
     IF(ik < 1) ik = 1
     DO j = ik, i
       tempta(i) = tempta(i)&
                 - fea_gkm(j, i-j+1, ibody) * fea_d(j, ibody)&
                 - fea_gcm(j, i-j+1, ibody) * fea_v(j, ibody)
     END DO !j
   END DO !i

   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i + fea_mband(ibody) - 1
     IF(ik > fea_npt(ibody)*fea_nd) ik = fea_npt(ibody) * fea_nd
     DO j = i+1, ik
       tempta(i) = tempta(i)&
                 - fea_gkm(i, j-i+1, ibody) * fea_d(j, ibody)&
                 - fea_gcm(i, j-i+1, ibody) * fea_v(j, ibody)
     END DO !j
   END DO !i

   DO i = 1, fea_npt(ibody)*fea_nd
   DO j = 1, fea_mband(ibody)
     tempt(i, j) = fea_gmm(i, j, ibody)
   END DO !j
   END DO !i

   CALL fea_decomposition(fea_npt(ibody)*fea_nd, fea_mband(ibody), tempt)
   CALL fea_backsubstitution(fea_npt(ibody)*fea_nd, &
                             fea_mband(ibody), tempt, tempta)

   DO i =1, fea_npt(ibody)*fea_nd
     fea_a(i, ibody) = tempta(i)
   END DO
 DEALLOCATE(tempt,tempta)
ENDIF !1==0

! Set keffect and other parameters for newmark scheme
   CALL fea_setnewmark(ibody)
  
 ENDSUBROUTINE fea_dynamicsetup
!----------------------------------------------------------

 SUBROUTINE fea_contactsetup(ibody)
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE grid_arrays


   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: i, j, k
   REAL(KIND=CGREAL) :: kr

   kr = 0.0_CGREAL
   DO i = 1, fea_npt(ibody)*fea_nd
   DO j =1, fea_mband(ibody)
     IF(kr < fea_gkm(i,j,ibody)) kr = fea_gkm(i,j,ibody)
   ENDDO
   ENDDO
   fea_penaltycoeff(ibody) = kr * fea_coeff_penalty

 ENDSUBROUTINE fea_contactsetup
!---------------------------------------------------------------------------------------------

 SUBROUTINE fea_dynamicsolver(ibody,n)
 
   USE global_parameters
   USE flow_parameters
   USE pressure_arrays
   USE finite_element_parameters
   USE finite_element_arrays
  
   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: i, j, ii, ik, n
!   REAL(KIND=CGREAL) :: tempt(fea_npt(ibody)*fea_nd, fea_mband(ibody)), &
!                         tempta(fea_npt(ibody)*fea_nd), & 
!                        rhsfea(fea_npt(ibody)*fea_nd), &
!                         kefffea(fea_npt(ibody)*fea_nd, fea_mband(ibody))
   REAL(KIND=CGREAL), POINTER, DIMENSION(:) :: rhsfea
   REAL(KIND=CGREAL),POINTER, DIMENSION(:,:) :: kefffea
   REAL(KIND=CGREAL) TEMP
  
   ALLOCATE(rhsfea(fea_npt(ibody)*fea_nd), &
            kefffea(fea_npt(ibody)*fea_nd, fea_mband(ibody)))


   rhsfea = 0.0_CGREAL
   kefffea =0.0_CGREAL

  IF(fea_contactprobflag==1 .and. fea_contact_model == PENALTY_CONTACT_MODEL) THEN 
 ! store old values for contact problem
   DO i = 1, fea_npt(ibody)*fea_nd
     fea_original_d(i, ibody) = fea_d(i, ibody)
     fea_original_v(i, ibody) = fea_v(i, ibody)
     fea_original_a(i, ibody) = fea_a(i, ibody)
   END DO !i
  ENDIF

   DO i = 1, fea_npt(ibody)*fea_nd
   DO j = 1, fea_mband(ibody)
     kefffea(i, j) = fea_keffect(i, j, ibody)
   END DO !j
   END DO !i
   !CALL fea_calload(ibody)
   CALL fea_compute_eqload(ibody)
   pi = 4.0_CGREAL * atan(1.0_CGREAL)

   DO i = 1, fea_nfixpt(ibody)
     DO j = 1, fea_nd

       IF(fea_ifixed(i, j+1, ibody) /= 0) THEN
         ii = fea_nd * (fea_ifixed(i, 1, ibody) - 1) + j
         fea_gp(ii, ibody) = 0.0_CGREAL
       ENDIF
     ENDDO
   ENDDO 
   DO i = 1, fea_npt(ibody)*fea_nd
     rhsfea(i) = fea_gp(i, ibody)
   END DO !i
  
   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i - fea_mband(ibody) + 1
     IF(ik < 1) ik = 1
     DO j = ik, i
       rhsfea(i) = rhsfea(i) &
                 + fea_gmm(j, i-j+1, ibody) * (fea_nmc0 * fea_d(j,ibody) &
                 + fea_nmc2 * fea_v(j,ibody) &
                 + fea_nmc3 * fea_a(j, ibody))&
                 + fea_gcm(j, i-j+1, ibody) * (fea_nmc1 * fea_d(j, ibody) &
                 + fea_nmc4 *fea_v(j, ibody)&
                 + fea_nmc5 * fea_a(j, ibody))

 
     END DO !j
   END DO !i
     
   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i + fea_mband(ibody) - 1
     IF(ik > fea_npt(ibody)*fea_nd) ik = fea_npt(ibody)*fea_nd
     DO j = i+1, ik
       rhsfea(i) = rhsfea(i) &
                 + fea_gmm(i, j-i+1, ibody) * (fea_nmc0 * fea_d(j,ibody)&
                 + fea_nmc2 * fea_v(j,ibody) &
                 + fea_nmc3 * fea_a(j, ibody))&
                 + fea_gcm(i, j-i+1, ibody) * (fea_nmc1 * fea_d(j, ibody) &
                                           + fea_nmc4 * fea_v(j, ibody)&
                                           + fea_nmc5 * fea_a(j, ibody))
 
 
     END DO !j
   END DO !i

   CALL fea_backsubstitution(fea_npt(ibody)*fea_nd, &
                             fea_mband(ibody), kefffea, rhsfea)

   DO i = 1, fea_npt(ibody)*fea_nd
     temp =  fea_nmc0 * (rhsfea(i) - fea_d(i, ibody)) - fea_nmc2 *          &
             fea_v(i, ibody) - fea_nmc3 * fea_a(i, ibody)

     fea_v(i, ibody) = fea_v(i, ibody) + fea_nmc6 * fea_a(i, ibody)&
                     + fea_nmc7 * temp
     fea_a(i, ibody) = temp
     fea_d(i, ibody) = rhsfea(i)
   END DO !i
 DO i = 1, fea_nfixpt(ibody)
     DO j = 1, fea_nd

       IF(fea_ifixed(i, j+1, ibody) /= 0) THEN
         ii = fea_nd * (fea_ifixed(i, 1, ibody) - 1) + j
         fea_d(ii, ibody) = 0.0_CGREAL
         fea_v(ii, ibody) = 0.0_CGREAL
         fea_a(ii, ibody) = 0.0_CGREAL
       ENDIF
     ENDDO
   ENDDO

!contact problem
   IF(fea_contactprobflag==1) then
     CALL fea_collisondetction(ibody)
     IF(fea_contactflag(ibody) == 1) THEN
       WRITE(*,*)  "BODY ", ibody, " is found collision"
      IF(fea_contact_model == PENALTY_CONTACT_MODEL) THEN
!      calcuate contact force
       IF (fea_nd == DIM_3D) THEN
         CALL fea_contact3D(ibody, rhsfea)
       ELSE IF (fea_nd ==DIM_2D) THEN
         CALL fea_contact2D(ibody, rhsfea)
       ENDIF !fea_nd
       DO i = 1, fea_npt(ibody)*fea_nd
         fea_d(i, ibody) = fea_original_d(i, ibody)
         fea_v(i, ibody) = fea_original_v(i, ibody)
         fea_a(i, ibody) = fea_original_a(i, ibody)
       END DO !i

       DO i = 1, fea_npt(ibody)*fea_nd
         ik = i - fea_mband(ibody) + 1
         IF(ik < 1) ik = 1
         DO j = ik, i
           rhsfea(i) = rhsfea(i) + fea_gmm(j, i-j+1, ibody) &
                     * (fea_nmc0 * fea_d(j,ibody) &
                     + fea_nmc2 * fea_v(j,ibody)  &
                     + fea_nmc3 * fea_a(j, ibody)) &
                     + fea_gcm(j, i-j+1, ibody) &
                     * (fea_nmc1 * fea_d(j, ibody) + fea_nmc4  &
                     * fea_v(j, ibody) + fea_nmc5 * fea_a(j, ibody))
         END DO !j
       END DO !i

       DO i = 1, fea_npt(ibody)*fea_nd
         ik = i + fea_mband(ibody) - 1
         IF(ik > fea_npt(ibody)*fea_nd) ik = fea_npt(ibody)*fea_nd
         DO j = i+1, ik
           rhsfea(i) = rhsfea(i) &
                     + fea_gmm(i, j-i+1, ibody) &
                     * (fea_nmc0 * fea_d(j,ibody) &
                     + fea_nmc2 * fea_v(j,ibody)  &
                     + fea_nmc3 * fea_a(j, ibody))&
                     + fea_gcm(i, j-i+1, ibody)&
                     * (fea_nmc1 * fea_d(j, ibody)&
                     + fea_nmc4 * fea_v(j, ibody) & 
                     + fea_nmc5 * fea_a(j, ibody))
         END DO !j
       END DO !i

       CALL fea_backsubstitution(fea_npt(ibody)*fea_nd,&
                                 fea_mband(ibody), kefffea, rhsfea)

       DO i = 1, fea_npt(ibody)*fea_nd
         temp =  fea_nmc0 * (rhsfea(i) - fea_d(i, ibody)) - fea_nmc2 *          &
                 fea_v(i, ibody) - fea_nmc3 * fea_a(i, ibody)
         fea_v(i, ibody) = fea_v(i, ibody) &
                         + fea_nmc6 * fea_a(i, ibody) + fea_nmc7 * temp
         fea_a(i, ibody) = temp
         fea_d(i, ibody) = rhsfea(i)
       END DO !i
      ENDIF !fea_contact_model
     ENDIF !contact detection
   ENDIF !contact flag
   deallocate(rhsfea, kefffea)

 ENDSUBROUTINE fea_dynamicsolver
!----------------------------------------------------------

 SUBROUTINE fea_calload(ibody)
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE boundary_arrays
   USE flow_parameters
   USE grid_arrays
   USE mpi

   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: ie(3)
   INTEGER           :: i, j, k, j1, j2, j3, j4, ii, ierr
   REAL(KIND=CGREAL) :: xx(3), yy(3), zz(3), xx3, yy3,&
                        xx4, yy4, zz4, A, ps(3), peq(3), normal(3), L

   DO i = 1, fea_nptmax*fea_nd
     fea_gp(i, ibody) = 0.0_CGREAL
   END DO !i

   DO i = 1, fea_nload(ibody)
     DO j = 1,fea_nd
       ie(j) = fea_iload(i, j, ibody)
     ENDDO

     DO j = 1, fea_nd  
       xx(j) = fea_cood(ie(j), 1, ibody)&
             + fea_d(fea_nd * (ie(j) - 1) + 1, ibody)
       yy(j) = fea_cood(ie(j), 2, ibody)&
             + fea_d(fea_nd * (ie(j) - 1) + 2, ibody)
       IF (fea_nd == 3) THEN
         zz(j) = fea_cood(ie(j), 3, ibody)&
               + fea_d(fea_nd * (ie(j) - 1) + 3, ibody)
       ELSE
         zz(j) = zc(1) 
       ENDIF
       IF(boundary_formulation == GCM_METHOD) THEN
         CALL fea_getpressure_GCM(xx(j), yy(j), zz(j), ps(j))
       ELSE
       !  CALL fea_getpressure_SSM(xx(j), yy(j), zz(j), ps(j))
       ENDIF
       
     ENDDO!j
    
     CALL par_getSumRealArray(ps, fea_nd)                                   !SER_TO_PAR. QX. CH31

     IF (fea_nd == DIM_3D) THEN

       CALL fea_compute_tri_area(xx(1), yy(1), zz(1), &
                             xx(2), yy(2),zz(2), &
                             xx(3), yy(3), zz(3), A)
       j = fea_iload(i, 4, ibody)
       xx4 = fea_cood(j, 1, ibody) + fea_d(fea_nd * (j - 1) + 1, ibody)
       yy4 = fea_cood(j, 2, ibody) + fea_d(fea_nd * (j - 1) + 2, ibody)
       zz4 = fea_cood(j, 3, ibody) + fea_d(fea_nd * (j - 1) + 3, ibody)
! Pei = A(pi/6+pj/12+pk/12)
       DO j = 1, 3
         j1=  j
         j2 = mod(j+1, 3)
         IF (j2==0) j2=3
         j3 = mod(j+2, 3)
         IF (j3==0) j3=3
         peq(j) = A *(ps(j1) / 6.0_CGREAL + ps(j2)/ 12.0_CGREAL &
                + ps(j3)/ 12.0_CGREAL)
       ENDDO
       CALL fea_compute_normal_3D(xx(1), yy(1), zz(1), xx(2), yy(2),&
                            zz(2), xx(3), yy(3), zz(3), &
                            XX4,YY4,ZZ4,NORMAL)
       DO j = 1, 3
       DO k = 1, fea_nd
         fea_gp(fea_nd*(ie(j)-1)+k, ibody) &
          = fea_gp(fea_nd*(ie(j)-1)+k, ibody) +  peq(j) * normal(k)
       ENDDO
       ENDDO

     ELSE IF (fea_nd == DIM_2D) THEN
 
       L = sqrt((xx(1) - xx(2))**2.0_CGREAL + (yy(1) - yy(2))**2.0_CGREAL)
       
       j = fea_iload(i, 3, ibody)
       xx3 = fea_cood(j, 1, ibody) + fea_d(fea_nd * (j - 1) + 1, ibody)
       yy3 = fea_cood(j, 2, ibody) + fea_d(fea_nd * (j - 1) + 2, ibody)

       peq(1) = L * (ps(1)/3.0_CGREAL + ps(2)/6.0_CGREAL)
       peq(2) = L * (ps(2)/3.0_CGREAL + ps(1)/6.0_CGREAL)
        
       WRITE(3333,*) i, xx(1), yy(1), ps(1)
       CALL fea_compute_normal_2D(xx(1), yy(1), xx(2), yy(2), xx3, yy3, NORMAL)
       DO j = 1,2 
       DO k = 1, fea_nd
         fea_gp(fea_nd*(ie(j)-1)+k, ibody) &
         = fea_gp(fea_nd*(ie(j)-1)+k, ibody)&
         +  peq(j) * normal(k)
       ENDDO
       ENDDO

     ENDIF

     DO j = 1, fea_nd
     DO k = 1, fea_nfixpt(ibody)
       IF ( ie(j) == fea_ifixed(k,1,ibody) ) THEN
         DO ii = 1, fea_nd
           fea_gp(fea_nd *(ie(j)-1)+ii, ibody)=0.0_CGREAL
         ENDDO
         EXIT
       ENDIF
     ENDDO
     ENDDO        

   END DO!i

   IF(ntime >= fea_nfilterstart .and. ntime <= fea_nfilterend) THEN
     DO i =1,  fea_nload(ibody)
       DO j = 1,fea_nd 
         ie(j) = fea_iload(i, j, ibody)
       ENDDO

       IF(ntime == fea_nfilterstart)THEN
         DO j = 1, fea_nd
         DO k = 1, fea_nd
           fea_gp_old(fea_nd*(ie(j)-1)+k, ibody) = 0.0_CGREAL
 !fea_gp(fea_nd*(ie(j)-1)+k, ibody)
         ENDDO
         ENDDO !j
       ENDIF

       DO j = 1, fea_nd
       DO k = 1, fea_nd
        fea_gp(fea_nd*(ie(j)-1)+k, ibody) &
        = fea_gp(fea_nd*(ie(j)-1)+k, ibody) * fea_omegaload &
        + (1 - fea_omegaload) * fea_gp_old(fea_nd*(ie(j)-1)+k, ibody)
       ENDDO
       ENDDO

       DO j = 1, fea_nd
       DO k = 1, fea_nd
         fea_gp_old(fea_nd*(ie(j)-1)+k, ibody) &
         = fea_gp(fea_nd*(ie(j)-1)+k, ibody)
       ENDDO
       ENDDO !j

     ENDDO !i
   ENDIF

 END SUBROUTINE fea_calload

!----------------------------------------------------------
 SUBROUTINE fea_compute_tri_area(X1, Y1, Z1, X2, Y2,Z2, X3, Y3, Z3, S)

   USE global_parameters

   IMPLICIT NONE

   REAL(KIND=CGREAL):: X1, Y1, Z1, X2, Y2, Z2, X3, Y3,Z3, S
   REAL(KIND=CGREAL) :: A, B, C

   A = SQRT((X1-X2)**2.0_CGREAL+(Y1-Y2)**2.0_CGREAL &
     + (Z1-Z2)**2.0_CGREAL)
   B = SQRT((X2-X3)**2.0_CGREAL+(Y2-Y3)**2.0_CGREAL &
     + (Z2-Z3)**2.0_CGREAL)
   C = SQRT((X3-X1)**2.0_CGREAL+(Y3-Y1)**2.0_CGREAL &
     + (Z3-Z1)**2.0_CGREAL)
   S = 0.5_CGREAL*(A+B+C)
   S = SQRT(S*(S-A)*(S-B)*(S-C))

 ENDSUBROUTINE fea_compute_tri_area
!----------------------------------------------------------

 SUBROUTINE fea_compute_normal_3D(X1,Y1,Z1, X2,Y2, Z2, &
                            X3,Y3, Z3, X4,Y4,Z4,NORMAL)

   USE global_parameters

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: X1,Y1,Z1,X2,Y2,Z2,&
                        X3,Y3,Z3,X4,Y4,Z4,NORMAL(3)
   REAL(KIND=CGREAL) :: XX1, YY1, ZZ1, XX2, YY2, ZZ2, S, NOR
   REAL(KIND=CGREAL) :: x_centroid, y_centroid, z_centroid
   REAL(KIND=CGREAL) :: XXC, YYC, ZZC
   INTEGER           :: I

   XX1 = X1 - X2
   XX2 = X2 - X3
   YY1 = Y1 - Y2
   YY2 = Y2 - Y3
   ZZ1 = Z1 - Z2
   ZZ2 = Z2 - Z3
   NORMAL(1) = YY1*ZZ2 - ZZ1*YY2
   NORMAL(2) = -(XX1 * ZZ2 - ZZ1 *XX2)
   NORMAL(3) = XX1 * YY2 - YY1 * XX2

   NOR = SQRT(NORMAL(1)**2 + NORMAL(2)**2 +NORMAL(3)**2)
   NORMAL(1) = NORMAL(1) / NOR
   NORMAL(2) = NORMAL(2) / NOR
   NORMAL(3) = NORMAL(3) / NOR
   x_centroid  = (X1 + X2 + X3) / 3.0_CGREAL
   y_centroid  = (Y1 + Y2 + Y3) / 3.0_CGREAL
   z_centroid  = (Z1 + Z2 + Z3) / 3.0_CGREAL
   XXC  = X4 - x_centroid
   YYC  = Y4 - y_centroid
   ZZC  = Z4 - z_centroid 

   S = NORMAL(1)* XXC + NORMAL(2) * YYC + NORMAL(3) * ZZC
   IF(S < 0 ) THEN
     DO I = 1, 3
       NORMAL(I) = - NORMAL(I)
     ENDDO
   ENDIF

ENDSUBROUTINE fea_compute_normal_3D 
!----------------------------------------------------------
SUBROUTINE fea_compute_normal_2D(X1,Y1, X2,Y2, X3,Y3, NORMAL)  

   USE global_parameters

   IMPLICIT NONE

   REAL(KIND=CGREAL) :: X1,Y1,X2,Y2,X3,Y3,NORMAL(3)
   REAL(KIND=CGREAL) :: XX1, YY1, S, NOR
   REAL(KIND=CGREAL) :: x_mid, y_mid
   REAL(KIND=CGREAL) :: XXC, YYC
   INTEGER           :: I

   XX1 = X1 - X2
   YY1 = Y1 - Y2
   NORMAL(1) = YY1
   NORMAL(2) = -XX1 
   NORMAL(3) = 0 

   NOR = SQRT(NORMAL(1)**2 + NORMAL(2)**2)
   NORMAL(1) = NORMAL(1) / NOR
   NORMAL(2) = NORMAL(2) / NOR
   x_mid  = (X1 + X2) / 2.0_CGREAL
   y_mid  = (Y1 + Y2) / 2.0_CGREAL
   XXC  = X3 - x_mid
   YYC  = Y3 - y_mid

   S = NORMAL(1)* XXC + NORMAL(2) * YYC 
   IF(S < 0 ) THEN
     DO I = 1, 2
       NORMAL(I) = - NORMAL(I)
     ENDDO
   ENDIF

ENDSUBROUTINE fea_compute_normal_2D
!----------------------------------------------------------

 SUBROUTINE fea_getpressure_SSM(xx, yy, zz, pma)                               !SER_TO_PAR. QX. CH29

   USE global_parameters
   USE grid_arrays
   USE boundary_arrays
   USE pressure_arrays
   USE flow_parameters

   IMPLICIT NONE

   REAL(KIND=CGREAL)  :: xx, yy, zz, pma
   INTEGER            :: i, ii, j, jj, k, kk, iflowc, fea_NPresStencil 
   INTEGER            :: icell(8), jcell(8), kcell(8)
   REAL(KIND=CGREAL)  :: totaldrv
   REAL(KIND=CGREAL)  :: pcell(8), drvcell(8)
   INTEGER            :: iil, jjl

   pcell   = 0.0_CGREAL
   drvcell = 0.0_CGREAL

   IF(xx <= x(1) .or. xx > x(nx_GLBL)) THEN
     ii = -1
   ELSE
     DO i  = 0, nx_GLBL-1
       IF(xx > xc(i) .and. xx <= xc(i+1)) then
         ii = i
         EXIT
       ENDIF
     ENDDO !i
   END IF

   IF(yy <= y(1) .or. yy > y(ny_GLBL)) THEN
     jj   = -1
   ELSE
     DO j  = 0, ny_GLBL-1
       IF( yy > yc(j) .and. yy <= yc(j+1)) then
         jj = j
         EXIT
       ENDIF
     ENDDO !j
   ENDIF

   IF (ndim == DIM_3D) THEN 

     IF(zz <= z(1) .or. zz > z(nz)) THEN
       kk   = -1
     ELSE
       DO k  = 0, nz-1
         IF( zz > zc(k) .and. zz <= zc(k+1)) THEN
           kk  = k
           EXIT
         ENDIF
       ENDDO !k
     ENDIF

   ELSE IF (ndim == DIM_2D) THEN
     kk = 1
   ENDIF

! Marker point is not inside the fluid region, pressure equal zero
   IF(ii < 0 .or. jj < 0 .or. kk < 0 ) THEN
     pma = 0.0_CGREAL
     RETURN
   ENDIF

   IF(x(myIs)> xx .or. x(myIe+1) <= xx .or. y(myJs) > yy .or. y(myJe+1) <= yy) THEN
    pma = 0.0_CGREAL
    RETURN
   ENDIF
   
   icell(1) = ii
   jcell(1) = jj
   kcell(1) = kk
   icell(2) = ii + 1
   jcell(2) = jj
   kcell(2) = kk
   icell(3) = ii + 1
   jcell(3) = jj + 1
   kcell(3) = kk
   icell(4) = ii
   jcell(4) = jj + 1
   kcell(4) = kk
   icell(5) = ii
   jcell(5) = jj
   kcell(5) = kk + 1
   icell(6) = ii + 1
   jcell(6) = jj
   kcell(6) = kk + 1
   icell(7) = ii + 1
   jcell(7) = jj + 1
   kcell(7) = kk + 1
   icell(8) = ii
   jcell(8) = jj + 1
   kcell(8) = kk + 1

   iflowc = 0
   fea_NPresStencil = 2**ndim
   
   DO i = 1, fea_NPresStencil
     iil =  G2LI(icell(i))
     jjl =  G2LJ(jcell(i)) 
     IF(iblank( iil, jjl, kcell(i)) == 0) THEN
       iflowc = iflowc + 1
       pcell(iflowc) = p(iil, jjl, kcell(i))
       drvcell(iflowc) = sqrt((xx- xc(icell(i)))**2 &
                       + (yy -yc(jcell(i)))**2 &
                       + (zz -zc(kcell(i)))**2)
       IF(drvcell(iflowc) == 0.0_CGREAL) THEN
         pma = pcell(iflowc)
         RETURN 
       ENDIF
       drvcell(iflowc) = 1.0_CGREAL / drvcell(iflowc)
     END IF
   END DO !i

   IF(iflowc == 0) THEN
     pma = 0.0_CGREAL
     RETURN
   END IF

   pma = 0.0_CGREAL
   totaldrv  =  0.0_CGREAL
   DO i = 1, iflowc
     totaldrv =  totaldrv + drvcell(i)
   END DO

   DO i = 1, iflowc
     pma =  pma + drvcell(i) * pcell(i) / totaldrv
   END DO

 END SUBROUTINE fea_getpressure_SSM
!----------------------------------------------------------

 SUBROUTINE fea_getpressure_GCM(xx, yy, zz, pma)                        !SER_TO_PAR. QX. CH30

   USE global_parameters
   USE grid_arrays
   USE boundary_arrays
   USE pressure_arrays
   USE flow_parameters
   USE flow_arrays

   IMPLICIT NONE

   REAL(KIND=CGREAL)  :: xx, yy, zz, pma
   INTEGER            :: i, ii, j, jj, k, kk, iflowc, fea_NPresStencil 
   INTEGER            :: icell(8), jcell(8), kcell(8)
   REAL(KIND=CGREAL)  :: totaldrv
   REAL(KIND=CGREAL)  :: pcell(8), drvcell(8)
   INTEGER            :: iil, jjl

   pcell   = 0.0_CGREAL
   drvcell = 0.0_CGREAL

   IF(xx <= x(1) .or. xx > x(nx_GLBL)) THEN
     ii = -1
   ELSE
     DO i  = 0, nx_GLBL-1
       IF(xx > xc(i) .and. xx <= xc(i+1)) then
         ii = i
          EXIT
       ENDIF
     ENDDO !i
   END IF

   IF(yy <= y(1) .or. yy > y(ny_GLBL)) THEN
     jj   = -1
   ELSE
     DO j  = 0, ny_GLBL-1
       IF( yy > yc(j) .and. yy <= yc(j+1)) then
         jj = j
         EXIT
       ENDIF
     ENDDO !j
   ENDIF
      
   IF (ndim == DIM_3D) THEN
     IF(zz <= z(1) .or. zz > z(nz)) THEN
       kk   = -1
     ELSE
       DO k  = 0, nz-1
         IF( zz > zc(k) .and. zz <= zc(k+1)) then
           kk  = k
           EXIT
         ENDIF
       ENDDO !k
     ENDIF
   ELSE IF (ndim == DIM_2D) THEN
     kk = 1
   ENDIF

! Marker point is not inside the fluid region, pressure equal zero
   IF(ii < 0 .or. jj < 0 .or. kk < 0 ) THEN
     pma = 0.0_CGREAL
     RETURN
   ENDIF

   IF(x(myIs)> xx .or. x(myIe+1) <= xx .or. y(myJs) > yy .or. y(myJe+1) <= yy) THEN
    pma = 0.0_CGREAL
    RETURN
   ENDIF

   icell(1) = ii
   jcell(1) = jj
   kcell(1) = kk
   icell(2) = ii + 1
   jcell(2) = jj
   kcell(2) = kk
   icell(3) = ii + 1
   jcell(3) = jj + 1
   kcell(3) = kk
   icell(4) = ii
   jcell(4) = jj + 1
   kcell(4) = kk
   icell(5) = ii
   jcell(5) = jj
   kcell(5) = kk + 1
   icell(6) = ii + 1
   jcell(6) = jj
   kcell(6) = kk + 1
   icell(7) = ii + 1
   jcell(7) = jj + 1
   kcell(7) = kk + 1
   icell(8) = ii
   jcell(8) = jj + 1
   kcell(8) = kk + 1

   iflowc = 0
   
   fea_NPresStencil = 2**ndim

   DO i = 1, fea_NPresStencil
     iil = G2LI(icell(i))
     jjl = G2LJ(jcell(i)) 
     IF(iblank( iil, jjl, kcell(i)) == 0 &
        .or. ghostCellMark(iil,jjl, kcell(i)) == 1) THEN
       iflowc = iflowc + 1
       pcell(iflowc) = p(iil, jjl, kcell(i))
       IF(ghostCellMark(iil,jjl, kcell(i)) ==1 )  pcell(iflowc) = pGhost(iil, jjl, kcell(i))
       drvcell(iflowc) = sqrt((xx- xc(icell(i)))**2 &
                       + (yy -yc(jcell(i)))**2 &
                       + (zz -zc(kcell(i)))**2 * REAL((ndim - DIM_2D), KIND=CGREAL))
       IF(drvcell(iflowc) == 0.0_CGREAL) THEN
         pma = pcell(iflowc)
         EXIT
       ENDIF
       drvcell(iflowc) = 1.0_CGREAL / drvcell(iflowc)
     END IF
   END DO !i

   IF(iflowc == 0) THEN
     pma = 0.0_CGREAL
     RETURN
   END IF

   pma = 0.0_CGREAL
   totaldrv  =  0.0_CGREAL
   DO i = 1, iflowc
     totaldrv =  totaldrv + drvcell(i)
   END DO
   DO i = 1, iflowc
     pma =  pma + drvcell(i) * pcell(i) / totaldrv
  END DO
  return

 END SUBROUTINE fea_getpressure_GCM
!---------------------------------------------------------

 SUBROUTINE fea_collisondetction(ibody)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE grid_arrays

   IMPLICIT NONE

   INTEGER           :: ibody
   REAL(KIND=CGREAL) :: pointlocation, contactlocation
   INTEGER           :: i, j

   fea_contactflag(ibody) = 0
   fea_ncontactpoint(ibody) = 0
   DO i = 1, fea_npt(ibody)
     SELECT CASE(fea_icontactdir(ibody))
       CASE(ICOORD)
         pointlocation = fea_cood(i, 1, ibody) &
                       + fea_d((i-1)*fea_nd + 1, ibody)
         contactlocation = x(fea_icontactplane(ibody))
         CASE(JCOORD)
         pointlocation = fea_cood(i, 2, ibody) &
                       + fea_d((i-1)*fea_nd + 2, ibody)
         contactlocation = y(fea_icontactplane(ibody))
       CASE(KCOORD)
         pointlocation = fea_cood(i, 3, ibody)&
                       + fea_d((i-1)*fea_nd + 3, ibody)
         contactlocation = z(fea_icontactplane(ibody))
     ENDSELECT


     IF(fea_icontactdirnorm(ibody) == 1) THEN
 !SURFACE NORMAL IS POSITIVE, THAT MEANS BODY PENETRATE SURFACE FROM UP
       IF(pointlocation < contactlocation) THEN
         fea_contactflag(ibody) = 1
         fea_ncontactpoint(ibody) = fea_ncontactpoint(ibody) + 1
         fea_icontactpoint(fea_ncontactpoint(ibody), ibody) = i
       ENDIF
     ELSE
       IF(pointlocation > contactlocation) THEN
         fea_contactflag(ibody) = 1
         fea_ncontactpoint(ibody) = fea_ncontactpoint(ibody) + 1
         fea_icontactpoint(fea_ncontactpoint(ibody), ibody) = i
       ENDIF
     ENDIF

   ENDDO  !i

 ENDSUBROUTINE fea_collisondetction
!---------------------------------

 SUBROUTINE fea_contact3D(ibody, rhsfea)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: i 
   REAL(KIND = CGREAL)      :: rhsfea(fea_npt(ibody)*fea_nd)

   CALL fea_contactsurf(ibody)

   DO i = 1, fea_ncontactsurf(ibody)
      CALL fea_contactforce3D(i, ibody)
   ENDDO !i

   DO i = 1, fea_npt(ibody)*fea_nd
      rhsfea(i) = fea_gp(i, ibody)
   END DO !i

 ENDSUBROUTINE
!-------------------------

 SUBROUTINE fea_contact2D(ibody, rhsfea)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
 
   IMPLICIT NONE

   INTEGER                  :: ibody
   INTEGER                  :: i
   REAL(KIND = CGREAL)      :: rhsfea(fea_npt(ibody)*fea_nd) 

   CALL fea_contactside(ibody)
   DO i = 1, fea_ncontactside(ibody)
     CALL fea_contactforce2D(i, ibody)
   ENDDO
  
   DO i = 1, fea_npt(ibody)*fea_nd
     rhsfea(i) = fea_gp(i, ibody)
   END DO !i

 ENDSUBROUTINE
!-------------------------------------------------------

 SUBROUTINE fea_contactsurf(ibody)

   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: i, j, k, ip, i1, i2, i3, i4

   fea_ncontactsurf(ibody) = 0

   DO j = 1, fea_nload(ibody) ! for every surface
     i1 = fea_iload(j, 1,ibody)
     i2 = fea_iload(j, 2,ibody)
     i3 = fea_iload(j, 3,ibody)
     DO i = 1, fea_ncontactpoint(ibody) ! for every contact point
       ip = fea_icontactpoint(i, ibody)
       IF(ip == i1 .OR. ip == i2 .OR. ip == i3)      THEN
         fea_ncontactsurf(ibody) = fea_ncontactsurf(ibody) + 1
         fea_icontactsurf(fea_ncontactsurf(ibody), ibody) = j
         EXIT
       ENDIF
     ENDDO !i
   ENDDO !j

 ENDSUBROUTINE fea_contactsurf
!----------------------------------------------------------

 SUBROUTINE fea_contactforce3D(isurf, ibody)
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE grid_arrays

   IMPLICIT NONE

   INTEGER           :: isurf, ibody

   INTEGER           :: n, i, j, k, idir
   INTEGER           :: ii, jj, kk, indx(4)
   INTEGER           :: j1, j2, j3
   REAL(KIND=CGREAL) :: xx(4), yy(4), zz(4)
   REAL(KIND=CGREAL) :: areaoftriangle, gs(3),  ps(3), peq(3), normal
   REAL(KIND=CGREAL) :: pointlocation, contactlocation
  
   n = fea_icontactsurf(isurf, ibody)
   idir = fea_icontactdir(ibody)

   SELECT CASE(fea_icontactdir(ibody))
     CASE(ICOORD)
       contactlocation = x(fea_icontactplane(ibody))
     CASE(JCOORD)
       contactlocation = y(fea_icontactplane(ibody))
     CASE(KCOORD)
       contactlocation = z(fea_icontactplane(ibody))
   ENDSELECT


   DO ii = 1, 4
     indx(ii) = fea_iload(n, ii, ibody)
     xx(ii) = fea_cood(indx(ii), 1, ibody) &
            + fea_d((indx(ii)-1)*fea_nd+1, ibody)
     yy(ii) = fea_cood(indx(ii), 2, ibody)&
            + fea_d((indx(ii)-1)*fea_nd+2, ibody)
     zz(ii) = fea_cood(indx(ii), 3, ibody) &
            + fea_d((indx(ii)-1)*fea_nd+3, ibody)
   ENDDO

   DO ii = 1, 3
     pointlocation = fea_cood(indx(ii), fea_icontactdir(ibody), ibody)&
                   + fea_d((indx(ii)-1)*fea_nd &
                   + fea_icontactdir(ibody), ibody)

     IF(fea_icontactdirnorm(ibody) == 1) THEN
       gs(ii) = contactlocation - pointlocation
     ELSE
       gs(ii) = pointlocation - contactlocation
     ENDIF

     IF(gs(ii) < 0) gs(ii) = 0
     ps(ii) = gs(ii) * fea_penaltycoeff(ibody)
   ENDDO

   CALL fea_compute_tri_area(xx(1), yy(1), zz(1), xx(2), yy(2),zz(2),&
                         xx(3), yy(3), zz(3), areaoftriangle)

   IF(fea_icontactdirnorm(ibody) == 1) THEN
    normal = 1.0_CGREAL
   ELSE
    normal = -1.0_CGREAL
   ENDIF

! Pei = A*(pi/6+pj/12+pk/12)
   DO j = 1, 3
     j1=  j
     j2 = mod(j+1, 3)
     IF (j2 == 0) j2=3
     j3 = mod(j+2, 3)
     IF (j3 == 0) j3=3
     peq(j) = areaoftriangle *(ps(j1) / 6.0_CGREAL + &
            ps(j2)/ 12.0_CGREAL + ps(j3)/ 12.0_CGREAL)
   ENDDO

   DO j = 1, 3
     fea_gp(fea_nd*(indx(j)-1)+idir, ibody) &
     = fea_gp(fea_nd*(indx(j)-1)+idir, ibody) +  peq(j) * normal
   ENDDO

 ENDSUBROUTINE fea_contactforce3D
!------------------------------------------------------------------------------------------------------

 SUBROUTINE fea_contactside(ibody)
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: i, j, ip, i1, i2

   fea_ncontactside(ibody) = 0

   DO j = 1, fea_nload(ibody) ! for every side 
     i1 = fea_iload(j, 1,ibody)
     i2 = fea_iload(j, 2,ibody)
     DO i = 1, fea_ncontactpoint(ibody) ! for every contact point
       ip = fea_icontactpoint(i, ibody)
       IF(ip == i1 .OR. ip == i2)      THEN
         fea_ncontactside(ibody) = fea_ncontactside(ibody) + 1
         fea_icontactside(fea_ncontactside(ibody), ibody) = j
         EXIT
       ENDIF
     ENDDO !i
   ENDDO !j

 ENDSUBROUTINE fea_contactside
!----------------------------------------------------------
 
 SUBROUTINE fea_contactforce2D(iside, ibody)
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE grid_arrays

   IMPLICIT NONE

   INTEGER           :: iside, ibody
   INTEGER           :: n, i, j, idir
   INTEGER           :: ii, jj, indx(3)
   REAL(KIND=CGREAL) :: xx(3), yy(3)
   REAL(KIND=CGREAL) :: lengthofside, gs(2),  ps(2), peq(2), normal
   REAL(KIND=CGREAL) :: pointlocation, contactlocation

   n = fea_icontactside(iside, ibody)
   idir = fea_icontactdir(ibody)

   SELECT CASE(fea_icontactdir(ibody))
     CASE(ICOORD)
       contactlocation = x(fea_icontactplane(ibody))
     CASE(JCOORD)
       contactlocation = y(fea_icontactplane(ibody))
   ENDSELECT


   DO ii = 1,2 
     indx(ii) = fea_iload(n, ii, ibody)
     xx(ii) = fea_cood(indx(ii), 1, ibody) &
            + fea_d((indx(ii)-1)*fea_nd+1, ibody)
     yy(ii) = fea_cood(indx(ii), 2, ibody)&
            + fea_d((indx(ii)-1)*fea_nd+2, ibody)
   ENDDO

   DO ii = 1,2 
     pointlocation = fea_cood(indx(ii), fea_icontactdir(ibody), ibody)&
                   + fea_d((indx(ii)-1)*fea_nd &
                   + fea_icontactdir(ibody), ibody)
     IF(fea_icontactdirnorm(ibody) == 1) THEN
       gs(ii) = contactlocation - pointlocation
     ELSE
       gs(ii) = pointlocation - contactlocation
     ENDIF
     IF(gs(ii) < 0) gs(ii) = 0
     ps(ii) = gs(ii) * fea_penaltycoeff(ibody)
   ENDDO
  
   lengthofside = sqrt((xx(1) - xx(2))**2.0_CGREAL &
                + (yy(1) - yy(2))**2.0_CGREAL)
   
   IF(fea_icontactdirnorm(ibody) == 1) THEN
     normal = 1.0_CGREAL
   ELSE  
     normal = -1.0_CGREAL
   ENDIF
         
! Pei = A*(pi/3+pj/6)
   peq(1) = lengthofside * (ps(1) / 3.0_CGREAL + ps(2) / 6.0_CGREAL)
   peq(2) = lengthofside * (ps(2) / 3.0_CGREAL + ps(1) / 6.0_CGREAL)
   
   DO j = 1, 2 
    fea_gp(fea_nd*(indx(j)-1)+idir, ibody) &
    = fea_gp(fea_nd*(indx(j)-1)+idir, ibody) +  peq(j) * normal
   ENDDO

 ENDSUBROUTINE fea_contactforce2D
!-----------------------------------------------------------------------------------------------

 SUBROUTINE fea_setnewmark(ibody)

   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays
  
   IMPLICIT NONE

   INTEGER           :: ibody
   
   INTEGER           :: i, j
   !REAL(KIND=CGREAL) :: keffecttempt(fea_npt(ibody)*fea_nd,fea_mband(ibody))

   fea_dt = dt / fea_dtratio
   keffecttempt = 0.0
   fea_nmc0 =  1.0_CGREAL / (fea_beta * fea_dt * fea_dt)
   fea_nmc1 = fea_gamma / (fea_beta * fea_dt)
   fea_nmc2 = 1.0_CGREAL / (fea_beta * fea_dt)
   fea_nmc3 = 0.5_CGREAL / fea_beta - 1.0_CGREAL
   fea_nmc4 = fea_gamma / fea_beta - 1.0_CGREAL
   fea_nmc5 = fea_dt / 2.0_CGREAL * (fea_gamma / fea_beta - 2.0_CGREAL)
   fea_nmc6 = fea_dt * (1.0_CGREAL - fea_gamma)
   fea_nmc7 = fea_gamma * fea_dt

 IF(fea_readkeffect == 0) THEN
! Form keffect for newmark scheme
   DO i = 1, fea_npt(ibody)*fea_nd 
   DO j = 1, fea_mband(ibody)
     keffecttempt(i, j) =  fea_gkm(i, j, ibody)&
                        + fea_nmc0 * fea_gmm(i, j, ibody)&
                        + fea_nmc1 * fea_gcm(i, j, ibody)
   END DO !j
   END DO !i 

! Triangle decompostion of keffect
   CALL fea_decomposition(fea_npt(ibody)*fea_nd,&
                          fea_mband(ibody),&
                          keffecttempt)
   
   DO i = 1, fea_npt(ibody)*fea_nd
   DO j = 1, fea_mband(ibody)
     fea_keffect(i, j, ibody) =  keffecttempt(i, j)
   END DO !j
   END DO !i    
  
  CALL fea_writekeffectfile(ibody)
ELSE
  CALL fea_readkeffectfile(ibody)
 ENDIF

 ENDSUBROUTINE fea_setnewmark
!----------------------------------------------------------

 SUBROUTINE fea_dynamicsolutionoutput(nbody, n)
   
   USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
  
   IMPLICIT NONE

   INTEGER          :: nbody, n
   INTEGER          :: i, j, ibody
   CHARACTER*19              :: fname1

   write(fname1,"('stress.',i7.7)") n

   OPEN(fea_stressfile, FILE=fname1,STATUS='UNKNOWN')

   IF (fea_nd == 3) THEN
     DO ibody = 1, nbody
      WRITE(fea_stressfile, *) 'VARIABLES = "X", "Y", "Z", "U", "V", "W"'
      WRITE(fea_stressfile, *) 'ZONE T = "DISP", F = FEPOINT, ET = TETRAHEDRON,N=', fea_npt(ibody), ',E=',  fea_nelement(ibody)
      DO i = 1, fea_npt(ibody)
        WRITE(fea_stressfile,*) (fea_cood(i, j, ibody) + fea_d(3*(i-1)+j, ibody), j = 1, fea_nd),(fea_V(3*(i-1)+j, ibody), j = 1, fea_nd)
      END DO
      DO i = 1, fea_nelement(ibody)
         WRITE(fea_stressfile,*) (fea_ielement(i, j, ibody), j = 6, 5+fea_nnode )
      END DO
     ENDDO
   ENDIF

   IF (fea_nd ==2) THEN
     DO ibody = 1, nbody
       WRITE (fea_stressfile, *) 'VARIABLES = "X", "Y", "U", "W"'
       WRITE (fea_stressfile, *)  'ZONE T = "DISP", F = FEPOINT, ET = TRIANGLE,N=', fea_npt(ibody), ',E=',  fea_nelement(ibody)
       DO i = 1, fea_npt(ibody)
         WRITE (fea_stressfile, *) (fea_cood(i, j, ibody) + fea_d(2*(i-1)+j, ibody), j = 1, fea_nd),(fea_V(2*(i-1)+j, ibody), j = 1, fea_nd)
       ENDDO
       DO i = 1, fea_nelement(ibody)
         WRITE(fea_stressfile,*) (fea_ielement(i, j, ibody), j = 5, 4+fea_nnode )
       ENDDO
     ENDDO
   ENDIF
 CLOSE(fea_stressfile)
 
  END SUBROUTINE fea_dynamicsolutionoutput
!----------------------------------------------------------

 SUBROUTINE fea_open_probe_files()

   USE fea_probe_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   CHARACTER*13          :: fea_probeFile
   CHARACTER*29         :: fea_inProbeFile

   INTEGER :: m, j

!  Open files for Probe

   DO m = 1, fea_nProbe
     fea_probeFile = TRIM("fea_probe_out")
     WRITE(fea_inProbeFile,101) fea_probeFile,m
     OPEN(UNIT=fea_ifuProbeOut+m-1,FILE=fea_inProbeFile,FORM='formatted',ACTION="WRITE")

!   Write Header
     WRITE(fea_ifuProbeOut+m-1,1000) &
           fea_iprobebody(m),fea_iprobenode(m),&
           (fea_cood(fea_iprobenode(m), j, fea_iprobebody(m)),j =1,fea_nd)
   END DO ! m

!   formats

101  FORMAT(a,'_',i3.3,'.dat')
1000 FORMAT('# probe data (time, u, v, w)',/, &
            '# ibody ',I5,', inode ',I5,/, &
            '# x=',E13.5,', y=',E13.5,', z=',E13.5)

 END SUBROUTINE fea_open_probe_files
!---------------------------------------------

 SUBROUTINE fea_read_probe_inputs()

   USE fea_probe_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER :: m

   OPEN(fea_probein,FILE='fea_probe_in.dat',STATUS='UNKNOWN')
   READ(fea_probein,*)fea_nProbe
   PRINT*,'   fea_nProbe = ',fea_nProbe

   ALLOCATE(fea_iprobebody(fea_nProbe))
   ALLOCATE(fea_iprobenode(fea_nProbe))

   PRINT*,'Reading fea_probe_in.dat File'
   DO m= 1, fea_nProbe
     READ(fea_probein,*)fea_iprobebody(m),fea_iprobenode(m)
   ENDDO ! m

 END SUBROUTINE fea_read_probe_inputs
!---------------------------------------------

 SUBROUTINE fea_write_probe_files()

    USE fea_probe_parameters
    USE flow_parameters
    USE finite_element_parameters
    USE finite_element_arrays

    IMPLICIT NONE

    INTEGER :: m, ibody, n
    REAL(KIND=CGREAL) :: uProbe,vProbe,wProbe, uuProbe, vvProbe, wwProbe

    DO m= 1, fea_nProbe
    !  write (*, *) 'fea_nprobe is ', fea_nprobe
      ibody = fea_iprobebody(m)
   !  write (*, *) 'ibody is ', ibody
      n     = fea_iprobenode(m)
   !  write (*, *) 'fea_iprobenode is', n

      IF (fea_nd == 3) THEN
        uProbe = fea_d(fea_nd * (n - 1) + 1, ibody)
        vProbe = fea_d(fea_nd * (n - 1) + 2, ibody)
        wProbe = fea_d(fea_nd * (n - 1) + 3, ibody)
        uuProbe = fea_v(fea_nd * (n - 1) + 1, ibody)
        vvProbe = fea_v(fea_nd * (n - 1) + 2, ibody)
        wwProbe = fea_v(fea_nd * (n - 1) + 3, ibody)
        WRITE(fea_ifuProbeOut+m-1,'(1PE14.7,6E18.7)') time,uProbe,vProbe,wProbe, uuProbe, vvProbe, wwProbe
      ELSE IF (fea_nd == 2) THEN
        uProbe = fea_d(fea_nd * (n - 1) + 1, ibody)
        vProbe = fea_d(fea_nd * (n - 1) + 2, ibody)
        uuProbe = fea_v(fea_nd * (n - 1) + 1, ibody)
        vvProbe = fea_v(fea_nd * (n - 1) + 2, ibody)
        WRITE(fea_ifuProbeOut+m-1,'(1PE14.7,4E18.7)') time,uProbe,vProbe, uuProbe, vvProbe
      ENDIF
    ENDDO ! m

 END SUBROUTINE fea_write_probe_files
!---------------------------------------------------------------------------------

 SUBROUTINE fea_readrestart
   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   INTEGER  :: ibody ,i ,j
   DO ibody = 1, fea_nbody
   DO i = 1, fea_npt(ibody)*fea_nd
     READ(fea_restartfilein) fea_d(i, ibody), fea_v(i, ibody),&
                             fea_a(i, ibody), fea_gp_old(i,ibody)
   END DO !i

 ! read the contact information
   IF(fea_contactprobflag==1 .and. fea_contact_model == GEO_CONTACT_MODEL) THEN
    READ(fea_restartfilein)fea_contactflag(ibody)
    IF(fea_contactflag(ibody) == 1) THEN
      READ(fea_restartfilein) fea_ncontactpoint(ibody)
      DO i = 1, fea_ncontactpoint(ibody)
       READ(fea_restartfilein) fea_icontactpoint(i, ibody)
      ENDDO
    ENDIF !(fea_contactflag(ibody) == 1
   ENDIF !fea_contactprobflag==1 .and. fea_contact_model == GEO_CONTACT_MODEL



   ENDDO !ibody

   CLOSE(fea_restartfilein)

 ENDSUBROUTINE fea_readrestart
!---------------------------------------------------------------------------------------

 SUBROUTINE fea_writerestart

   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   IMPLICIT NONE

   INTEGER  :: ibody ,i ,j
   SELECT CASE(idxRstrt)
     CASE(2)
       OPEN(fea_restartfileout, FILE='fea_restart_out1.dat',FORM='UNFORMATTED')
     CASE(1)
       OPEN(fea_restartfileout, FILE='fea_restart_out2.dat',FORM='UNFORMATTED')
   END SELECT ! indexRstrt

   DO ibody = 1, fea_nbody
   DO i = 1, fea_npt(ibody)*fea_nd
     WRITE(fea_restartfileout) fea_d(i, ibody), fea_v(i, ibody), &
                               fea_a(i, ibody), fea_gp(i,ibody)
   END DO !i
   
   ! record the contact information
   IF(fea_contactprobflag==1 .and. fea_contact_model == GEO_CONTACT_MODEL) THEN
    WRITE(fea_restartfileout)fea_contactflag(ibody)
    IF(fea_contactflag(ibody) == 1) THEN
      WRITE(fea_restartfileout) fea_ncontactpoint(ibody)
      DO i = 1, fea_ncontactpoint(ibody)
       WRITE(fea_restartfileout) fea_icontactpoint(i, ibody)
      ENDDO
    ENDIF !fea_contactflag(ibody) == 1
   ENDIF !fea_contactprobflag==1 .and. fea_contact_model == GEO_CONTACT_MODEL

   ENDDO !ibody

   CLOSE(fea_restartfileout)

 ENDSUBROUTINE fea_writerestart
!---------------------------------------------------------------------------------------
 SUBROUTINE fea_readMfile(ibody)
   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   
   IMPLICIT NONE

   INTEGER  :: ibody ,i ,j
   
    READ(fea_Mfile) 
    DO i = 1, fea_npt(ibody)*fea_nd
     DO j = 1, fea_mband(ibody)
      READ(fea_Mfile) fea_gmm(i,j,ibody),fea_gcm(i, j, ibody), fea_gkm(i, j, ibody) 
     END DO !j
    END DO !i
 

 ENDSUBROUTINE fea_readMfile
!---------------------------------------------------------------------------------------

 SUBROUTINE fea_writeMfile(ibody)
  USE global_parameters
  USE flow_parameters
  USE finite_element_parameters
  USE finite_element_arrays
  
  IMPLICIT NONE

  INTEGER  :: ibody ,i ,j
  IF(ImtheBOSS) THEN
    WRITE(fea_Mfile)ibody
    DO i = 1, fea_npt(ibody)*fea_nd
     DO j = 1, fea_mband(ibody)
      WRITE(fea_Mfile) fea_gmm(i,j,ibody),fea_gcm(i, j, ibody), fea_gkm(i, j, ibody)
     END DO !j
    END DO !i
  ENDIF 


 ENDSUBROUTINE fea_writeMfile
!---------------------------------------------------------------------------------------
 SUBROUTINE fea_readkeffectfile(ibody)
   USE global_parameters
   USE flow_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER  :: ibody ,i ,j

    READ(fea_keffectfile)
    DO i = 1, fea_npt(ibody)*fea_nd
     DO j = 1, fea_mband(ibody)
      READ(fea_keffectfile) fea_keffect(i,j,ibody)
     END DO !j
    END DO !i


 ENDSUBROUTINE fea_readkeffectfile
!---------------------------------------------------------------------------------------

 SUBROUTINE fea_writekeffectfile(ibody)
  USE global_parameters
  USE flow_parameters
  USE finite_element_parameters
  USE finite_element_arrays

  IMPLICIT NONE

  INTEGER  :: ibody ,i ,j
  IF(ImtheBOSS) THEN
    WRITE(fea_keffectfile)ibody
    DO i = 1, fea_npt(ibody)*fea_nd
     DO j = 1, fea_mband(ibody)
      WRITE(fea_keffectfile) fea_keffect(i,j,ibody)
     END DO !j
    END DO !i
  ENDIF


 ENDSUBROUTINE fea_writekeffectfile
!---------------------------------------------------------------------------------------

SUBROUTINE fea_body_marker_pressure
 USE global_parameters
 USE flow_parameters
 USE boundary_arrays 
 IMPLICIT NONE

 INTEGER :: ibody, ipt
 REAL(KIND=CGREAL) :: ptmp(nPtsMax)
 
 DO ibody = 1, nbody
  ptmp = 0.0_CGREAL
  IF(boundary_formulation == GCM_METHOD) THEN
   DO ipt = 1, nPtsBodyMarker(ibody)
     CALL fea_getpressure_GCM(xBodyMarker(ibody, ipt),yBodyMarker(ibody, ipt),zBodyMarker(ibody, ipt),ptmp(ipt))
   ENDDO
  ELSE
   DO ipt = 1, nPtsBodyMarker(ibody)
    CALL fea_getpressure_SSM(xBodyMarker(ibody, ipt),yBodyMarker(ibody, ipt),zBodyMarker(ibody, ipt),ptmp(ipt))
   ENDDO
  ENDIF
#   ifdef MPI
  CALL par_getSumRealArray(ptmp,  nPtsBodyMarker(ibody))
#   endif
  DO ipt = 1, nPtsBodyMarker(ibody) 
   pBodyMarker(ibody, ipt) = ptmp(ipt)
 ENDDO
ENDDO !ibody 

ENDSUBROUTINE fea_body_marker_pressure
!---------------------------------------------------------------------------------------

SUBROUTINE fea_compute_eqload(ibody)
  USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays
   USE boundary_arrays
   USE flow_parameters
   USE grid_arrays
   USE mpi

   IMPLICIT NONE

   INTEGER           :: ibody
   INTEGER           :: ie(3)
   INTEGER           :: i, j, k, j1, j2, j3, j4, ii, n, m
   REAL(KIND=CGREAL) :: xx(3), yy(3), zz(3), xx3, yy3,& 
                        xx4, yy4, zz4, A, ps(3), peq(3), normal(3), L

   DO i = 1, fea_nptmax*fea_nd
     fea_gp(i, ibody) = 0.0_CGREAL
   END DO !i

   DO i = 1, fea_nload(ibody)
     DO j = 1,fea_nd
       ie(j) = fea_iload(i, j, ibody)
    ENDDO

     DO j = 1, fea_nd
       n = featomarkerpointtable(ie(j), 1, ibody)
       m = featomarkerpointtable(ie(j), 2, ibody)
      IF(n <= 0 .or. m <= 0) THEN
        write(*,*) 'could not find the corresponding point in marker for:', ie(j), ibody
       stop
       ENDIF
       xx(j) = xBodyMarker(m, n)
       yy(j) = yBodyMarker(m, n)
       IF (fea_nd == DIM_3D) THEN
        zz(j) = zBodyMarker(m, n)
       ELSE
         zz(j) = zc(1)
       ENDIF
      ps(j) = pBodyMarker(m, n)
     ENDDO!j


     IF (fea_nd ==DIM_3D) THEN

       CALL fea_compute_tri_area(xx(1), yy(1), zz(1), &
                             xx(2), yy(2),zz(2), &
                             xx(3), yy(3), zz(3), A)
       j = fea_iload(i, 4, ibody)
       xx4 = fea_cood(j, 1, ibody) + fea_d(fea_nd * (j - 1) + 1, ibody)
       yy4 = fea_cood(j, 2, ibody) + fea_d(fea_nd * (j - 1) + 2, ibody)
       zz4 = fea_cood(j, 3, ibody) + fea_d(fea_nd * (j - 1) + 3, ibody)
! Pei = A(pi/6+pj/12+pk/12)
       DO j = 1, 3
         j1=  j
         j2 = mod(j+1, 3)
         IF (j2==0) j2=3
         j3 = mod(j+2, 3)
         IF (j3==0) j3=3
         peq(j) = A *(ps(j1) / 6.0_CGREAL + ps(j2)/ 12.0_CGREAL &
                + ps(j3)/ 12.0_CGREAL)
       ENDDO
       CALL fea_compute_normal_3D(xx(1), yy(1), zz(1), xx(2), yy(2),&
                            zz(2), xx(3), yy(3), zz(3), &
                            XX4,YY4,ZZ4,NORMAL)
       DO j = 1, 3
       DO k = 1, fea_nd
         fea_gp(fea_nd*(ie(j)-1)+k, ibody) &
          = fea_gp(fea_nd*(ie(j)-1)+k, ibody) +  peq(j) * normal(k)
       ENDDO
       ENDDO

     ELSE IF (fea_nd == DIM_2D) THEN

       L = sqrt((xx(1) - xx(2))**2.0_CGREAL + (yy(1) - yy(2))**2.0_CGREAL)
       
       j = fea_iload(i, 3, ibody)
       xx3 = fea_cood(j, 1, ibody) + fea_d(fea_nd * (j - 1) + 1, ibody)
       yy3 = fea_cood(j, 2, ibody) + fea_d(fea_nd * (j - 1) + 2, ibody)

       peq(1) = L * (ps(1)/3.0_CGREAL + ps(2)/6.0_CGREAL)
       peq(2) = L * (ps(2)/3.0_CGREAL + ps(1)/6.0_CGREAL)

       CALL fea_compute_normal_2D(xx(1), yy(1), xx(2), yy(2), xx3, yy3, NORMAL)
       DO j = 1,2
       DO k = 1, fea_nd
         fea_gp(fea_nd*(ie(j)-1)+k, ibody) &
         = fea_gp(fea_nd*(ie(j)-1)+k, ibody)&
         +  peq(j) * normal(k)
       ENDDO
       ENDDO
      

     ENDIF

     DO j = 1, fea_nd
     DO k = 1, fea_nfixpt(ibody)
       IF ( ie(j) == fea_ifixed(k,1,ibody) ) THEN
         DO ii = 1, fea_nd
           fea_gp(fea_nd *(ie(j)-1)+ii, ibody)=0.0_CGREAL
         ENDDO
         EXIT
       ENDIF
     ENDDO
     ENDDO

   END DO!i

   IF(ntime >= fea_nfilterstart .and. ntime <= fea_nfilterend) THEN
    IF(ntime == fea_nfilterstart)THEN
     DO i =  1, fea_npt(ibody)*fea_nd
      fea_gp_old(i, ibody) =  0.0_CGREAL !fea_gp(i, ibody)
     ENDDO
    ENDIF

    DO i =  1, fea_npt(ibody)*fea_nd
     fea_gp(i, ibody) = fea_gp(i, ibody) * fea_omegaload + (1.0_CGREAL - fea_omegaload) * fea_gp_old(i, ibody)
    ENDDO

    DO i =  1, fea_npt(ibody)*fea_nd
     fea_gp_old(i, ibody) = fea_gp(i, ibody)
    ENDDO
   ENDIF

ENDSUBROUTINE fea_compute_eqload
!------------------------------------------------------------------------------------------------------
SUBROUTINE fea_compute_initial_acceleration(ibody)
  USE global_parameters
   USE finite_element_parameters
   USE finite_element_arrays

   IMPLICIT NONE

   INTEGER           :: ibody

   INTEGER           :: i, j, ii, ik, n

   ALLOCATE(tempt(fea_npt(ibody)*fea_nd, fea_mband(ibody)),tempta(fea_npt(ibody)*fea_nd))

  CALL fea_compute_eqload(ibody)
! Compute initial acceleration
   DO i = 1, fea_npt(ibody)*fea_nd
     tempta(i) = fea_gp(i, ibody)
   END DO !i
   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i - fea_mband(ibody) + 1
     IF(ik < 1) ik = 1
     DO j = ik, i
       tempta(i) = tempta(i)&
                 - fea_gkm(j, i-j+1, ibody) * fea_d(j, ibody)&
                 - fea_gcm(j, i-j+1, ibody) * fea_v(j, ibody)
     END DO !j
   END DO !i

   DO i = 1, fea_npt(ibody)*fea_nd
     ik = i + fea_mband(ibody) - 1
     IF(ik > fea_npt(ibody)*fea_nd) ik = fea_npt(ibody) * fea_nd
     DO j = i+1, ik
       tempta(i) = tempta(i)&
                 - fea_gkm(i, j-i+1, ibody) * fea_d(j, ibody)&
                 - fea_gcm(i, j-i+1, ibody) * fea_v(j, ibody)
     END DO !j
   END DO !i

   DO i = 1, fea_npt(ibody)*fea_nd
   DO j = 1, fea_mband(ibody)
     tempt(i, j) = fea_gmm(i, j, ibody)
   END DO !j
   END DO !i

   CALL fea_decomposition(fea_npt(ibody)*fea_nd, fea_mband(ibody), tempt)
   CALL fea_backsubstitution(fea_npt(ibody)*fea_nd, &
                             fea_mband(ibody), tempt, tempta)

   DO i =1, fea_npt(ibody)*fea_nd
     fea_a(i, ibody) = tempta(i)
   END DO
  DEALLOCATE(tempt,tempta)

ENDSUBROUTINE fea_compute_initial_acceleration
!-------------------------------------------------------------------------------
SUBROUTINE fea_writeafile(ibody)
 USE global_parameters
  USE flow_parameters
  USE finite_element_parameters
  USE finite_element_arrays

  IMPLICIT NONE

  INTEGER  :: ibody ,i ,j
  IF(ImtheBOSS) THEN
    WRITE(fea_afile)ibody
    DO i = 1, fea_npt(ibody)*fea_nd
      WRITE(fea_afile) fea_a(i,ibody)
    END DO !i
  ENDIF

ENDSUBROUTINE fea_writeafile
!------------------------------------------------------------------------------ 
SUBROUTINE fea_readafile
  USE global_parameters
  USE flow_parameters 
  USE finite_element_parameters
  USE finite_element_arrays
   
  IMPLICIT NONE

  INTEGER  :: ibody ,i ,j
  IF(ImtheBOSS) THEN
    READ(fea_afile) 
    DO i = 1, fea_npt(ibody)*fea_nd
      READ(fea_afile) fea_a(i,ibody)
    END DO !i
  ENDIF 

ENDSUBROUTINE fea_readafile
!------------------------------------------------------------------------------- 
