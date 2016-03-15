!  Filename: CORE_CUTCELL3D.PAR.F90
!  Latest Modification: Dec, 22 2010 (ver. PAT 1.2.0)
!  Made by JHSeo

# define L2GI(i)      myIs+i-1
# define L2GJ(j)      myJs+j-1

# define G2LI(i)      i-(myIs-1)
# define G2LJ(j)      j-(myJs-1)


SUBROUTINE cutcell_rhs(nlu,nlv,nlw)
 USE global_parameters
 USE flow_parameters
 USE boundary_arrays
 USE grid_arrays
 USE unstructured_surface_arrays
 USE flow_arrays
 USE cutcell_arrays
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(nxc,nyc,nzc) :: nlu, nlv, nlw
 INTEGER i, j, k, iG, jG, kG
 REAL(KIND=CGREAL) fadvec(3), fdIFf(3), dt_adv, dt_diff
 
 
 DO i=1, nxc
 DO j=1, nyc
 DO k=1, nzc
 
  iG = L2GI(i)
  jG = L2GJ(j)
  kG = k
 
  CALL cutcell_flux_AD(i,j,k,fadvec,fdIFf)
  
  nlu(i,j,k) = nlu(i,j,k) + fadvec(1)*dxinv(iG)*dyinv(jG)*dzinv(k)
  nlv(i,j,k) = nlv(i,j,k) + fadvec(2)*dxinv(iG)*dyinv(jG)*dzinv(k)  
  nlw(i,j,k) = nlw(i,j,k) + fadvec(3)*dxinv(iG)*dyinv(jG)*dzinv(k)
 
 ENDDO
 ENDDO
 ENDDO
 
  
END SUBROUTINE cutcell_rhs


SUBROUTINE cutcell_flux_AD(i,j,k,fadvec,fdIFf)
 USE global_parameters
 USE flow_parameters
 USE boundary_arrays
 USE grid_arrays
 USE unstructured_surface_arrays
 USE cutcell_arrays
 USE flow_arrays
 USE GCM_arrays
 
 IMPLICIT NONE

 INTEGER i, j, k, cElem, iBody,iG,jG,kG
 REAL(KIND=CGREAL) xp, yp, zp, xBI, yBI, zBI, uBI, vBI, wBI, dist
 REAL(KIND=CGREAL) fadvec(3), fdiff(3)
 REAL(KIND=CGREAL) dudn, dvdn, dwdn, usigma, xBIN, yBIN, zBIN
 
 iG = L2GI(i)
 jG = L2GJ(j)
 kG = k
 
 IF(ghostcellsolid(i,j,k) == 1) THEN
  xp = xc(iG)
  yp = yc(jG)
  zp = zc(k)
 
  iBody = BodyNum(i,j,k)
  IF(iBody == 0) THEN
   iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
   IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
  ENDIF
 
  CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI, yBI, zBI, cElem ) 
  CALL GCM_calc_BIVelocity_Unstruc( i, j, k, xBI, yBI, zBI, cElem, uBI, vBI, wBI )
  
  xBIN = triElemNormx(iBody,cElem)
  yBIN = triElemNormy(iBody,cElem)
  zBIN = triElemNormz(iBody,cElem)  

!! dist = sqrt( (xp-xBI)**2 + (yp-yBI)**2 + (zp-zBI)**2 ) 

  usigma = xBIN*uBI + yBIN*vBI + zBIN*wBI
  
  fadvec(1) = uBI*usigma*cut_area(i,j,k) 
  fadvec(2) = vBI*usigma*cut_area(i,j,k)
  fadvec(3) = wBI*usigma*cut_area(i,j,k)
  
  dudn = 0.0!(uBI-u(i,j,k))/dist
  dvdn = 0.0!(vBI-v(i,j,k))/dist
  dwdn = 0.0!(wBI-w(i,j,k))/dist
  
  fdIFf(1) = viscTot(i,j,k)*dudn*cut_area(i,j,k)*dz(k)
  fdIFf(2) = viscTot(i,j,k)*dvdn*cut_area(i,j,k)*dz(k)  
  fdIFf(3) = viscTot(i,j,k)*dwdn*cut_area(i,j,k)*dz(k)
   
 ELSE

  fadvec = 0.0
  fdiff  = 0.0

 ENDIF
	
END SUBROUTINE cutcell_flux_AD


SUBROUTINE cutcell_dpdn(i,j,k,dpdn)
 USE global_parameters
 USE flow_parameters
 USE boundary_arrays
 USE grid_arrays
 USE unstructured_surface_arrays
 USE cutcell_arrays
 
 IMPLICIT NONE

 INTEGER i, j, k, cElem, iBody, iG, jG, kG
 REAL(KIND=CGREAL) xp, yp, zp, xBI, yBI, zBI, axBI, ayBI, azBI, xBIN, yBIN, zBIN, dpdn
 
 iG = L2GI(i)
 jG = L2GJ(j)
 kG = k
 
 IF(ghostcellsolid(i,j,k) == 1) THEN
  xp = xc(iG)
  yp = yc(jG)
  zp = zc(k)
 
  iBody = BodyNum(i,j,k)
  IF(iBody == 0) THEN
   iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
   IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
  ENDIF
 
  CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI, yBI, zBI, cElem ) 
!!  CALL GCM_calc_BIAccel_Unstruc( i, j, k, xBI, yBI, zBI, cElem, axBI, ayBI, azBI )
  
  xBIN = triElemNormx(iBody,cElem)
  yBIN = triElemNormy(iBody,cElem)
  zBIN = triElemNormz(iBody,cElem)  

  dpdn = 0.0 !!!  - (xBIN*axBI + yBIN*ayBI + zBIN*azBI)
 
 ELSE

  dpdn = 0.0

 ENDIF
	
END SUBROUTINE cutcell_dpdn


SUBROUTINE cutcell_norm_vel(i,j,k,usigma)
 USE global_parameters
 USE flow_parameters
 USE boundary_arrays
 USE grid_arrays
 USE unstructured_surface_arrays
 USE cutcell_arrays
 
 IMPLICIT NONE

 INTEGER i, j, k, cElem, iBody, iG, jG, kG
 REAL(KIND=CGREAL) xp, yp, zp, xBI, yBI, zBI, uBI, vBI, wBI, usigma, xBIN, yBIN, zBIN
 
 iG = L2GI(i)
 jG = L2GJ(j)
 kG = k
 
 IF(ghostcellsolid(i,j,k) == 1) THEN
  xp = xc(iG)
  yp = yc(jG)
  zp = zc(k)
 
  iBody = BodyNum(i,j,k)
  IF(iBody == 0) THEN
   iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
   IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
  ENDIF
 
  CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI, yBI, zBI, cElem ) 
  CALL GCM_calc_BIVelocity_Unstruc( i, j, k, xBI, yBI, zBI, cElem, uBI, vBI, wBI )

  xBIN = triElemNormx(iBody,cElem)
  yBIN = triElemNormy(iBody,cElem)
  zBIN = triElemNormz(iBody,cElem)  

  usigma = xBIN*uBI + yBIN*vBI + zBIN*wBI
 
 ELSE

  usigma = 0.0

 ENDIF
	
END SUBROUTINE cutcell_norm_vel


SUBROUTINE mix_conserve(q1)
 USE global_parameters
 USE flow_parameters
 USE flow_arrays
 USE CUTCELL_arrays
 USE GCM_arrays
 USE unstructured_surface_arrays
 USE boundary_arrays
 USE grid_arrays

 IMPLICIT NONE
 
 REAL(KIND=CGREAL) q1(nxc,nyc,nzc)
 REAL(KIND=CGREAL),DIMENSION(:,:,:),ALLOCATABLE :: q2, vol2, vola
 REAL(KIND=CGREAL),DIMENSION(:,:,:,:),ALLOCATABLE :: face_cor
 REAL(KIND=CGREAL) myTime, xBI, yBI, zBI, xBIN, yBIN, zBIN, xp, yp, zp, gamma
 INTEGER iBody, cElem
 INTEGER i_tgx,j_tgx,k_tgx,i_tgy,j_tgy,k_tgy,i_tgxy,j_tgxy,k_tgxy
 INTEGER i_tgz,j_tgz,k_tgz
 REAL(KIND=CGREAL) beta_x,beta_y,beta_xy,X_tgx,X_tgy,X_tgxy,beta_sum
 REAL(KIND=CGREAL) beta_z,X_tgz
 INTEGER i,j,k,ip,jp,kp
  
 ALLOCATE(q2(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(vola(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(vol2(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 
 q2   = 0.0_CGREAL
 vola = 0.0_CGREAL
 vol2 = 1.0_CGREAL
 
 gamma = 1.0
 
 
 DO i=1,nxc
 DO j=1,nyc
 DO k=1,nzc
 q2(i,j,k) = q1(i,j,k)*dx(L2GI(i))*dy(L2GJ(j))*dz(k)
 vola(i,j,k) = VOL_CELL(i,j,k)
 vol2(i,j,k) = VOL_CELL(i,j,k)*dx(L2GI(i))*dy(L2GJ(j))*dz(k)
 ENDDO
 ENDDO
 ENDDO
 
 CALL par_comm_var(q2,nxc,nyc,nzc,Ngl,myTime)
 CALL par_comm_var(vola,nxc,nyc,nzc,Ngl,myTime)
 CALL par_comm_var(vol2,nxc,nyc,nzc,Ngl,myTime)
 
 DO i=0,nxc+1
 DO j=0,nyc+1
 DO k=0,nzc+1
 
  IF(ghostcellmark(i,j,k).eq.1) THEN
  
   xp = xc(L2GI(i))
   yp = yc(L2GJ(j))
   zp = zc(k)
   
   iBody = BodyNum(i,j,k)
   IF(iBody == 0) THEN
    iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
    IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
   ENDIF
  
   CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI, yBI, zBI, cElem ) 

   xBIN = triElemNormx(iBody,cElem)
   yBIN = triElemNormy(iBody,cElem)
   zBIN = triElemNormz(iBody,cElem)  
   
   
   IF(xBIN .gt. 0.0) THEN
    ip = -1
   ELSE 
    ip = 1
   ENDIF
   
   IF(yBIN .gt. 0.0) THEN
    jp = -1
   ELSE
    jp = 1
   ENDIF
   
   IF(zBIN .gt. 0.0) THEN
    kp = -1
   ELSE
    kp = 1
   ENDIF
   
   i_tgx = i + ip
   j_tgx = j
   k_tgx = k
   
   i_tgy = i
   j_tgy = j + jp
   k_tgy = k
   
   i_tgz = i 
   j_tgz = j 
   k_tgz = k + kp
   
   beta_x =  xBIN**2*vola(i_tgx,j_tgx,k_tgx)**gamma * (1.0-ghostcellmark(i_tgx,j_tgx,k_tgx))
   beta_y =  yBIN**2*vola(i_tgy,j_tgy,k_tgy)**gamma * (1.0-ghostcellmark(i_tgy,j_tgy,k_tgy))
   beta_z =  zBIN**2*vola(i_tgz,j_tgz,k_tgz)**gamma * (1.0-ghostcellmark(i_tgz,j_tgz,k_tgz))
   
   IF(ndim .eq. 2) beta_z = 0.0_CGREAL
   
   beta_sum = beta_x + beta_y + beta_z
   
   beta_x = beta_x / beta_sum
   beta_y = beta_y / beta_sum
   beta_z = beta_z / beta_sum
   
   X_tgx =  beta_x *q2(i,j,k)   
   X_tgy =  beta_y *q2(i,j,k)  
   X_tgz =  beta_z *q2(i,j,k) 
   
   q2(i_tgx,j_tgx,k_tgx) = q2(i_tgx,j_tgx,k_tgx) + X_tgx
   q2(i_tgy,j_tgy,k_tgy) = q2(i_tgy,j_tgy,k_tgy) + X_tgy
   q2(i_tgz,j_tgz,k_tgz) = q2(i_tgz,j_tgz,k_tgz) + X_tgz
   
   vol2(i_tgx,j_tgx,k_tgx) = vol2(i_tgx,j_tgx,k_tgx) + beta_x*vola(i,j,k)*dx(L2GI(i))*dy(L2GJ(j))*dz(k)
   vol2(i_tgy,j_tgy,k_tgy) = vol2(i_tgy,j_tgy,k_tgy) + beta_y*vola(i,j,k)*dx(L2GI(i))*dy(L2GJ(j))*dz(k)
   vol2(i_tgz,j_tgz,k_tgz) = vol2(i_tgz,j_tgz,k_tgz) + beta_z*vola(i,j,k)*dx(L2GI(i))*dy(L2GJ(j))*dz(k)  
   
   q2(i,j,k) = 0.0
   
  ENDIF

 ENDDO
 ENDDO
 ENDDO
 
 DO i=1,nxc
 DO j=1,nyc
 DO k=1,nzc
 
  q1(i,j,k) = q2(i,j,k)*dxinv(L2GI(i))*dyinv(L2GJ(j))*dzinv(k)

!   IF(isolid(i,j,k).eq.0) THEN 
!    q1(i,j,k) = q2(i,j,k)/vol2(i,j,k)*vola(i,j,k)
!   ELSE
!    q1(i,j,k) = 0.0
!   ENDIF
  
 ENDDO
 ENDDO
 ENDDO
  
 DEALLOCATE(q2,vola,vol2) 
 
END SUBROUTINE mix_conserve

 
SUBROUTINE CutCell_Allocate_Memory
 USE global_parameters
 USE flow_parameters
 USE CUTCELL_arrays
 
 ALLOCATE(AREA_W(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(AREA_E(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(AREA_S(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(AREA_N(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(AREA_B(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(AREA_F(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(VOL_CELL(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(CUT_AREA(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(isolid(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(cent_x(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 ALLOCATE(cent_y(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1)) 
 ALLOCATE(cent_z(1-ngl:nxc+ngl,1-ngl:nyc+ngl,0:nzc+1))
 
END SUBROUTINE


SUBROUTINE init_CutCell
 USE global_parameters
 USE flow_parameters
 USE flow_arrays
 USE boundary_arrays
 USE grid_arrays
 USE gcm_arrays
 USE unstructured_surface_arrays
 USE CUTCELL_arrays

 IMPLICIT NONE
 
 REAL(KIND=CGREAL) xBI,yBI,zBI
 REAL(KIND=CGREAL) xBIN,yBIN,zBIN
 REAL(KIND=CGREAL) vertex(8,3)
 REAL(KIND=CGREAL) Area_frac(6)
 REAL(KIND=CGREAL) vol_frac
 REAL(KIND=CGREAL) Area_cut, cent(3)
 

 REAL(KIND=CGREAL) xBI1, yBI1, zBI1, xp, yp, zp
 REAL(KIND=CGREAL) xBIN1, yBIN1, zBIN1, myTime
 
 INTEGER cElem, iBody, iG, jG, kG, i, j, k
  
 Area_W = 1.0
 Area_E = 1.0
 Area_S = 1.0
 Area_N = 1.0
 Area_B = 1.0
 Area_F = 1.0
 VOL_CELL = 1.0
 CUT_AREA = 0.0
 isolid = 0
 
 IF(iCC .eq. 1) THEN
 
 DO k=1,nzc
 DO j=1,nyc
 DO i=1,nxc
 
  iG = L2GI(i)
  jG = L2GJ(j)
  kG = k
  
  cent_x(i,j,k) = xc(iG)
  cent_y(i,j,k) = yc(jG)
  cent_z(i,j,k) = zc(kG)
  
  IF(iblank_solid(i,j,k) .eq. 1) THEN
   
   Area_W(i,j,k) = 0.0
   Area_E(i,j,k) = 0.0
   Area_S(i,j,k) = 0.0
   Area_N(i,j,k) = 0.0
   Area_B(i,j,k) = 0.0
   Area_F(i,j,k) = 0.0
   
   VOL_CELL(i,j,k) = 0.0
   CUT_AREA(i,j,k) = 0.0
   isolid(i,j,k) = 1
   
  ENDIF
  
  IF(ghostcellsolid(i,j,k) .eq. 1) THEN
  
   xp = xc(iG)
   yp = yc(jG)
   zp = zc(kG)
 
   CALL GCM_Calc_BodyIntercept_Unstruc( i, j, k, xp, yp, zp, xBI1, yBI1, zBI1, cElem ) 
   
   xBI = xBI1
   yBI = yBI1
   zBI = zBI1
 
   iBody   = bodyNum(i,j,k)
	
   IF(iBody == 0) THEN
    iBody = max(bodyNum(i+1,j,k),bodyNum(i-1,j,k),bodyNum(i,j+1,k),bodyNum(i,j-1,k))
    IF( ndim == 3 ) iBody = max(iBody,bodyNum(i,j,k+1),bodyNum(i,j,k-1))
   ENDIF
   
   vertex(1,1) = x(iG)
   vertex(1,2) = y(jG+1)
   vertex(1,3) = z(kG)
   
   vertex(2,1) = x(iG)
   vertex(2,2) = y(jG)
   vertex(2,3) = z(kG)
   
   vertex(3,1) = x(iG+1)
   vertex(3,2) = y(jG)
   vertex(3,3) = z(kG)
   
   vertex(4,1) = x(iG+1)
   vertex(4,2) = y(jG+1)
   vertex(4,3) = z(kG)

   vertex(5,1) = x(iG)
   vertex(5,2) = y(jG+1)
   vertex(5,3) = z(kG+1)
   
   vertex(6,1) = x(iG)
   vertex(6,2) = y(jG)
   vertex(6,3) = z(kG+1)
   
   vertex(7,1) = x(iG+1)
   vertex(7,2) = y(jG)
   vertex(7,3) = z(kG+1)
   
   vertex(8,1) = x(iG+1)
   vertex(8,2) = y(jG+1)
   vertex(8,3) = z(kG+1)

   xBIN = triElemNormx(iBody,cElem)
   yBIN = triElemNormy(iBody,cElem)
   zBIN = triElemNormz(iBody,cElem)  
   
   CALL Cal_CutCell(vertex,xBI,yBI,zBI,xBIN,yBIN,zBIN,Area_frac,vol_frac,area_cut,cent)
   
   Area_W(i,j,k) = Area_frac(1)*dyinv(jG)*dzinv(k)
   Area_E(i,j,k) = Area_frac(2)*dyinv(jG)*dzinv(k)
   Area_S(i,j,k) = Area_frac(3)*dxinv(iG)*dzinv(k)
   Area_N(i,j,k) = Area_frac(4)*dxinv(iG)*dzinv(k)
   Area_B(i,j,k) = Area_frac(5)*dxinv(iG)*dyinv(jG)
   Area_F(i,j,k) = Area_frac(6)*dxinv(iG)*dyinv(jG)
      
   VOL_CELL(i,j,k) = vol_frac*dxinv(iG)*dyinv(jG)*dzinv(k)
   CUT_AREA(i,j,k) = area_cut
   
   IF(VOL_CELL(i,j,k) .gt. 0.0) THEN
     isolid(i,j,k) = 0
	 cent_x(i,j,k) = cent(1)
	 cent_y(i,j,k) = cent(2)
	 cent_z(i,j,k) = cent(3)
   ENDIF
   
  ENDIF
  
 ENDDO
 ENDDO
 ENDDO
 
! CALL par_comm_var(cent_x,nxc,nyc,nzc,Ngl,myTime)
! CALL par_comm_var(cent_y,nxc,nyc,nzc,Ngl,myTime)
! CALL par_comm_var(cent_z,nxc,nyc,nzc,Ngl,myTime)

 ENDIF !  iCC
 
END SUBROUTINE init_CutCell



SUBROUTINE Cal_CutCell(vertex_xyz,xBI,yBI,zBI,xBIN,yBIN,zBIN,Area_F,Vol_F,Cut_area,cent)
 USE global_parameters
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) vertex_xyz(8,3), vertex(4,3)
 REAL(KIND=CGREAL) xBI, yBI, zBI, xBIN, yBIN, zBIN
 REAL(KIND=CGREAL) Area_F(6)
 REAL(KIND=CGREAL) Vol_F, Cut_area, cal_cutArea2
 REAL(KIND=CGREAL) Area_Fraction, intercept_xyz(4,3)
 INTEGER vertex_to_face(6,4)
 INTEGER iPlane, nintercept, i, j, k
 REAL(KIND=CGREAL) edge_xyz1(6,3), edge_xyz2(6,3)
 INTEGER nedge, Sum_inout, F_inout, F_count
 REAL(KIND=CGREAL) cx, cy, cz
 REAL(KIND=CGREAL) cent(3), Face_cent(3)
 
 
 DATA vertex_to_face/6,7,3,4,1,5, &
                     2,3,2,1,2,6, &
					 1,4,6,5,3,7, &
					 5,8,7,8,4,8/
 
 
 nedge = 0
 SUM_inout = 0
 cent = 0.0
 F_count = 0
 
 DO i=1,6
 
 DO j=1,4
  vertex(j,1:3) = vertex_xyz(vertex_to_face(i,j),1:3)
 ENDDO
 
 iPlane = (i+1)/2
 CALL Cal_Face(vertex,xBI,yBI,zBI,xBIN,yBIN,zBIN,iPlane,Area_Fraction,intercept_xyz,nintercept,F_inout,Face_cent)
 SUM_inout = SUM_inout + F_inout
 
 IF(F_inout .lt. 4 ) THEN
  cent(1) = cent(1) + Face_cent(1)
  cent(2) = cent(2) + Face_cent(2)
  cent(3) = cent(3) + Face_cent(3)
  F_count = F_count + 1
 ENDIF
  
 
 IF(nintercept .gt. 0) THEN
 
 IF(nintercept .ne. 2) THEN
  WRITE(*,*) 'wrong intersections ', nintercept
  STOP
 ENDIF
 
 nedge = nedge + 1
 edge_xyz1(nedge,1:3) = intercept_xyz(1,1:3)
 edge_xyz2(nedge,1:3) = intercept_xyz(2,1:3)
 
 ENDIF
 
 Area_F(i) = Area_Fraction
  
 ENDDO
 
 IF(SUM_inout .eq. 0) THEN
  
  Cut_Area = 0.0
  Vol_F = abs(vertex_xyz(3,1)-vertex_xyz(2,1))*abs(vertex_xyz(4,2)-vertex_xyz(3,2))*abs(vertex_xyz(2,3)-vertex_xyz(6,3))
 
 ELSEIF(SUM_inout .eq. 24) THEN
 
  Cut_Area = 0.0
  Vol_F = 0.0
  
 ELSE
  
  Cut_Area = Cal_cutArea2(edge_xyz1,edge_xyz2,nedge,cx,cy,cz)
  
  Vol_F = Area_F(2)*vertex_xyz(3,1) - Area_F(1)*vertex_xyz(2,1) + &
          Area_F(4)*vertex_xyz(1,2) - Area_F(3)*vertex_xyz(2,2) + &
		  Area_F(6)*vertex_xyz(6,3) - Area_F(5)*vertex_xyz(2,3)
		  
  Vol_F = Vol_F + (cx*xBIN+cy*yBIN+cz*zBIN)*Cut_Area
  Vol_F = Vol_F/3.0
  
 ENDIF
 
 IF( Cut_Area .gt. 0.0 ) THEN
  F_count = F_count + 1
  cent(1) = cent(1) + cx
  cent(2) = cent(2) + cy
  cent(3) = cent(3) + cz
 ENDIF
 
 IF(F_count .gt. 0 ) THEN
  cent(1) = cent(1)/F_count
  cent(2) = cent(2)/F_count
  cent(3) = cent(3)/F_count
 ENDIF
 
END SUBROUTINE Cal_CutCell


SUBROUTINE Cal_Face(vertex_xyz,xBI,yBI,zBI,xBIN,yBIN,zBIN,iPlane,Area_Fraction,intercept_xyz,nintercept,SUM_inout,cent)
 USE global_parameters
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) vertex_xyz(4,3)
 INTEGER vertex_inout(4)
 INTEGER edge_int(4)
 REAL(KIND=CGREAL) Area_Fraction
 REAL(KIND=CGREAL) xBI,yBI,zBI
 REAL(KIND=CGREAL) xBIN,yBIN,zBIN
 INTEGER iPlane
 REAL(KIND=CGREAL) intercept_xyz(4,3)
 INTEGER vertex_to_edge(4,2) 
 INTEGER Check_inout
 REAL(KIND=CGREAL) Plane_intercept
 REAL(KIND=CGREAL) Cal_polyA
 INTEGER npoly
 INTEGER nintercept, SUM_inout, i, j, k, iCheck
 REAL(KIND=CGREAL) poly_x(8), poly_y(8), poly_z(8)
 REAL(KIND=CGREAL) x_int(4), y_int(4), z_int(4)
 REAL(KIND=CGREAL) cent(3), tol
 REAL(KIND=CGREAL) tt, xp, yp, zp, x1, y1, z1, x2, y2, z2, dx, dy
 
 tol = 1.e-10
 
 npoly = 0
 nintercept = 0
 SUM_inout = 0
 edge_int = 0
 cent = 0.0
 
 DATA vertex_to_edge/1,2,3,4,2,3,4,1/
 
 DO i=1,4
 vertex_inout(i) = Check_inout(xBI,yBI,zBI,xBIN,yBIN,zBIN,vertex_xyz(i,1),vertex_xyz(i,2),vertex_xyz(i,3))
 SUM_inout = SUM_inout + vertex_inout(i)
 ENDDO
 
 IF(SUM_inout .gt. 0 .and. SUM_inout .lt. 4) THEN  !! Cut-cell
 
 DO i=1,4
  iCheck = vertex_inout(vertex_to_edge(i,1)) + vertex_inout(vertex_to_edge(i,2))
  
  x1 = vertex_xyz(vertex_to_edge(i,1),1)
  y1 = vertex_xyz(vertex_to_edge(i,1),2)
  z1 = vertex_xyz(vertex_to_edge(i,1),3)
  
  x2 = vertex_xyz(vertex_to_edge(i,2),1)
  y2 = vertex_xyz(vertex_to_edge(i,2),2)
  z2 = vertex_xyz(vertex_to_edge(i,2),3)
  
  SELECT CASE(iCheck)
  
  CASE ( 0 )  ! entire edge is out of the body
  
  edge_int(i) = 0
    
  CASE ( 1 )  ! edge is cut by the body-surface
  
  tt = Plane_intercept(xBI,yBI,zBI,xBIN,yBIN,zBIN,x1,y1,z1,x2,y2,z2)
  
  IF( tt.ge.(0.0-tol) .and. tt.le.(1.0+tol) ) THEN
   
   xp = x1 + (x2-x1)*tt
   yp = y1 + (y2-y1)*tt
   zp = z1 + (z2-z1)*tt
   
   nintercept = nintercept + 1
   x_int(nintercept) = xp
   y_int(nintercept) = yp
   z_int(nintercept) = zp
   
   intercept_xyz(nintercept,1) = xp
   intercept_xyz(nintercept,2) = yp
   intercept_xyz(nintercept,3) = zp   
   
   edge_int(i) = nintercept
   
  ELSE
  
  WRITE(*,*) ' Cutcell : Geometric Error ! '
  WRITE(*,*) tt
  WRITE(*,*) x1, y1, z1
  WRITE(*,*) x2, y2, z2
  WRITE(*,*) xBI, yBI, zBI
  WRITE(*,*) xBIN, yBIN, zBIN
  STOP
  
  ENDIF
  
    
  CASE ( 2 )  ! entire edge is inside of the body
  
  edge_int(i) = 0
  
  END SELECT
  
 ENDDO
 
 !! Calc Surface Fraction
 !! Construct Polygon
 
 npoly = 0
 
 DO i=1,4
 
 IF(vertex_inout(i) .eq. 0) THEN
  npoly = npoly + 1
  poly_x(npoly) = vertex_xyz(i,1)
  poly_y(npoly) = vertex_xyz(i,2)
  poly_z(npoly) = vertex_xyz(i,3)
 ENDIF
 
 IF(edge_int(i) .ne. 0) THEN
  npoly = npoly + 1
  poly_x(npoly) = x_int(edge_int(i))
  poly_y(npoly) = y_int(edge_int(i))
  poly_z(npoly) = z_int(edge_int(i))
 ENDIF
 
 ENDDO
  
 ! Center of polygon
 
  DO i=1,npoly
  cent(1) = cent(1) + poly_x(i)/npoly
  cent(2) = cent(2) + poly_y(i)/npoly
  cent(3) = cent(3) + poly_z(i)/npoly
  ENDDO
 
 ! Calc area of polygon
 Area_Fraction = Cal_polyA(poly_x,poly_y,poly_z,npoly,iPlane)
  
 ELSEIF(sum_inout .eq. 0)  THEN !! Entire face is out of the body
 
  DO i=1,4
  cent(1) = cent(1) + 0.25*vertex_xyz(i,1)
  cent(2) = cent(2) + 0.25*vertex_xyz(i,2)
  cent(3) = cent(3) + 0.25*vertex_xyz(i,3)
  ENDDO
 
 SELECT CASE (iPlane)
 CASE ( 1 )
 dx = abs(vertex_xyz(3,2)-vertex_xyz(2,2))
 dy = abs(vertex_xyz(4,3)-vertex_xyz(3,3))
 CASE ( 2 )
 dx = abs(vertex_xyz(3,3)-vertex_xyz(2,3))
 dy = abs(vertex_xyz(4,1)-vertex_xyz(3,1)) 
 CASE ( 3 )
 dx = abs(vertex_xyz(3,1)-vertex_xyz(2,1))
 dy = abs(vertex_xyz(4,2)-vertex_xyz(3,2))
 END SELECT
 
 Area_Fraction = dx*dy
 
 ELSEIF(sum_inout .eq. 4) THEN  !! Entire face is inside of the body
  
 cent = 0.0
 Area_Fraction = 0.0
   
 ENDIF
  
 END SUBROUTINE Cal_Face
 

REAL(KIND=SELECTED_REAL_KIND(P=14,R=30)) FUNCTION Cal_cutArea2(xyz1,xyz2,nedge,cx,cy,cz)
  USE global_parameters

  IMPLICIT NONE
  
  REAL(KIND=CGREAL) xyz1(6,3), xyz2(6,3)
  INTEGER nedge, i
  REAL(KIND=CGREAL) Summa, tx1, ty1, tz1, tx2, ty2, tz2, nx, ny, nz, dS, cx, cy, cz, osign
  
  Summa = 0.0
  
  cx = 0.0 ; cy = 0.0 ; cz = 0.0
  
  DO i=1,nedge
  
  cx = cx + 0.5*(xyz2(i,1) + xyz1(i,1))/nedge
  cy = cy + 0.5*(xyz2(i,2) + xyz1(i,2))/nedge
  cz = cz + 0.5*(xyz2(i,3) + xyz1(i,3))/nedge
  
  ENDDO
    
  DO i=1,nedge
  
  tx1 = xyz1(i,1) - cx
  ty1 = xyz1(i,2) - cy
  tz1 = xyz1(i,3) - cz
  
  tx2 = xyz2(i,1) - cx
  ty2 = xyz2(i,2) - cy
  tz2 = xyz2(i,3) - cz
  
  nx = ty1*tz2 - tz1*ty2
  ny = tz1*tx2 - tx1*tz2
  nz = tx1*ty2 - ty1*tx2
  
  Summa = Summa + 0.5*sqrt(nx**2 + ny**2 + nz**2)
  
  ENDDO
  
  Cal_cutArea2 = Summa
  
END FUNCTION Cal_cutArea2

  
REAL(KIND=SELECTED_REAL_KIND(P=14,R=30)) FUNCTION Cal_polyA(xi,yi,zi,npoly,iPlane)
  USE global_parameters
  
  REAL(KIND=CGREAL) xi(8),x(0:8)
  REAL(KIND=CGREAL) yi(8),y(0:8)
  REAL(KIND=CGREAL) zi(8),z(0:8)
  INTEGER npoly, iPlane
  REAL(KIND=CGREAL) summa
  
  SELECT CASE ( iPlane )
  
  CASE ( 1 )  ! i-Plane (y-z)
  
  x(0) = yi(npoly)
  y(0) = zi(npoly)
  
  DO i=1,npoly
   x(i) = yi(i)
   y(i) = zi(i)
  ENDDO
  
  CASE ( 2 ) ! j-plane (z-x)
  
  x(0) = zi(npoly)
  y(0) = xi(npoly)
  
  DO i=1,npoly
   x(i) = zi(i)
   y(i) = xi(i)
  ENDDO
  
  CASE ( 3 ) ! k-plane (x-y)
  
  x(0) = xi(npoly)
  y(0) = yi(npoly)
   
  DO i=1,npoly
   x(i) = xi(i)
   y(i) = yi(i)
  ENDDO
  
  END SELECT
  
  summa = 0.0
  
  DO i=0,npoly-1
  summa = summa + ( x(i)*y(i+1) - x(i+1)*y(i) )
  ENDDO
  
  Cal_polyA = abs(0.5*summa)
    
 END FUNCTION Cal_polyA
 
 
!! FUNCTION to compute intersection btw the line and the plane
REAL(KIND=SELECTED_REAL_KIND(P=14,R=30)) FUNCTION Plane_Intercept(xBI,yBI,zBI,xBIN,yBIN,zBIN,x1,y1,z1,x2,y2,z2)
 USE global_parameters
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) xBI,yBI,zBI
 REAL(KIND=CGREAL) xBIN,yBIN,zBIN
 REAL(KIND=CGREAL) x1,y1,z1
 REAL(KIND=CGREAL) x2,y2,z2
 REAL(KIND=CGREAL) a,b,c,d,num,den,t
 
 d = -(xBI*XBIN + yBI*yBIN + zBI*zBIN)
 a = xBIN
 b = yBIN
 c = zBIN
 
 num = -d -a*x1 -b*y1 -c*z1
 den = a*(x2-x1) + b*(y2-y1) +c*(z2-z1)
 
 IF(abs(den) .lt. 1.e-14) THEN
  t = -255.0
 ELSE
  t = num/den
 ENDIF
 
 Plane_Intercept = t
 
END FUNCTION Plane_Intercept



INTEGER FUNCTION Check_inout(xBI,yBI,zBI,xBIN,yBIN,zBIN,xp,yp,zp)
 USE global_parameters
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) xBI,yBI,zBI
 REAL(KIND=CGREAL) xBIN,yBIN,zBIN
 REAL(KIND=CGREAL) xp,yp,zp
 REAL(KIND=CGREAL) VX,VY,VZ, cDot
 
 VX = xBI - xp
 VY = yBI - yp
 VZ = zBI - zp
 
 cDot = VX*xBIN + VY*yBIN + VZ*zBIN
 
 IF(cDot .gt. 0.0) THEN
  Check_inout = 0
 ELSE
  Check_inout = 1
 ENDIF
 
END FUNCTION Check_inout
  
  






