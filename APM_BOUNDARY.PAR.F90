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
!  Filename: APM_boundary.PAR.F90
!  Latest Modification: Sep, 11 2010 (ver. P2.0.0)
!  Made by Jung-Hee, Seo
! --------------------------------------------------------------------
! APM_boundary.F90
! Ver.  0.9.3.
! Date  08.17.2010
!
! Update History
! 2010.08.17: Revised normal vector computation
! 


MODULE APM_arrays
 INTEGER  Nmax, NCOL
 INTEGER, DIMENSION(:,:),   ALLOCATABLE :: i_data, j_data, k_data
 REAL(KIND= SELECTED_REAL_KIND(P=14,R=30)),    DIMENSION(:),     ALLOCATABLE :: Rmax, BIx, BIy, BIz, BInx, BIny, BInz
 REAL(KIND= SELECTED_REAL_KIND(P=14,R=30)),    DIMENSION(:,:,:), ALLOCATABLE :: APM_mat
END MODULE APM_arrays
!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE APM_allocate_memory(N_max,N_order)
 USE flow_parameters
 USE APM_arrays
 USE GCM_arrays
 
 IMPLICIT NONE
 
 INTEGER N_max, N_order, N_dim, i
 
 IF(ndim == 2) THEN
 SELECT CASE(N_order)
 CASE(1)
 NCOL = 4
 CASE(2)
 NCOL = 6
 CASE(3)
 NCOL = 10
 CASE(4)
 NCOL = 15
 CASE(5)
 NCOL = 21
 END SELECT
 ENDIF
 
 IF(ndim == 3) THEN
 SELECT CASE(N_order)
 CASE(1)
 NCOL = 8
 CASE(2)
 NCOL = 10
 CASE(3)
 NCOL = 20
 CASE(4)
 NCOL = 35
 END SELECT
 ENDIF
 
 Nmax = N_max
 
 ALLOCATE(i_data(nGhost,Nmax),j_data(nGhost,Nmax),k_data(nGhost,Nmax))
 ALLOCATE(Rmax(nGhost))
 ALLOCATE(BIx(nGhost),BIy(nGhost),BIz(nGhost))
 ALLOCATE(BInx(nGhost),BIny(nGhost),BInz(nGhost))
 ALLOCATE(APM_mat(nGhost,4,Nmax))

END SUBROUTINE APM_allocate_memory
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_deallocate
 USE APM_arrays
 
 IMPLICIT NONE
 
 DEALLOCATE(i_data,j_data,k_data)
 DEALLOCATE(Rmax, BIx, BIy, BIz)
 DEALLOCATE(BInx, BIny, BInz)
 DEALLOCATE(APM_mat)
 
END SUBROUTINE APM_deallocate
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_search_2D(gmap)
 USE flow_parameters
 USE APM_arrays
 USE GCM_arrays
 USE boundary_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 INTEGER gmap(mxa,mya,mza)
 INTEGER ISS, ISE, JSS, JSE, KSS, KSE
 INTEGER imin, imax, jmin, jmax, kmin, kmax
 INTEGER iis, iie, jjs, jje, kks, kke, iistep, jjstep, kkstep
 INTEGER iGC, jGC, kGC
 INTEGER Nnow
 INTEGER n, i, j, k
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: icheck
  
 REAL(KIND=CGREAL) xGC, yGC, zGC, dlength, Rnow, Rtemp, Rmin
 
 ALLOCATE(icheck(mxa,mya,mza))
 
 imin = ias
 jmin = jas
 kmin = kas
 
 imax = iae-1
 jmax = jae-1
 kmax = kae-1
  
 DO n=1, nGhost
 
  Nnow = 0
  Rnow = 0.0
  
  iGC = iGhost(n) + IFS
  jGC = jGhost(n) + JFS
  kGC = kGhost(n) + KFS
  
  xGC = xca(iGC)
  yGC = yca(jGC)
  zGC = zca(kGC)
   
  BIx(n) = xBodyIntercept(n)
  BIy(n) = yBodyIntercept(n)
  BIz(n) = zGC !zBodyIntercept(n)
  
  !! dlength = sqrt( (BIx(n) - xGC)**2 + (BIy(n) - yGC)**2 + (BIz(n) - zGC)**2 )
  
  BInx(n) =  xBodyInterceptNorm(n) !(BIx(n)-xGC)/dlength
  BIny(n) =  yBodyInterceptNorm(n) !(BIy(n)-yGC)/dlength
  BInz(n) =  zBodyInterceptNorm(n) !(BIz(n)-zGC)/dlength
  
  ISS = iGC -1 ; ISE = iGC +1
  JSS = jGC -1 ; JSE = jGC +1
  KSS = kGC -1 ; KSE = kGC +1
  
!  write(*,*) iGC, jGC, kGC
  Rmin = 0.0
  
  icheck = 0
  
  DO WHILE( Nnow .lt. Nmax-1 ) 
   
   ISS = max(ISS,imin)
   ISE = min(ISE,imax)
   JSS = max(JSS,jmin)
   JSE = min(JSE,jmax)
   KSS = max(KSS,kmin)
   KSE = min(KSE,kmax)
   
   ! for 2D
   
   k = 1+KFS
   
   IF(BIx(n) .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = 1
   ELSE
    iis = ISE ; iie = ISS ; iistep = -1
   ENDIF
   
   IF(BIy(n) .le. yGC) THEN
    jjs = JSS ; jje = JSE ; jjstep = JSE-JSS
   ELSE
    jjs = JSE ; jje = JSS ; jjstep = JSS-JSE
   ENDIF
   
   
   DO j = jjs, jje, jjstep
   DO i = iis, iie, iistep
    
	IF(icheck(i,j,k) .ne. 1) THEN
	
    Rtemp = sqrt( (xca(i)-BIx(n))**2 + (yca(j)-BIy(n))**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax-1) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
!     write(*,*) '*',i,j,k,Rtemp
     ENDIF
    ENDIF
    
	ENDIF
	
   ENDDO
   ENDDO
   
   IF(BIx(n) .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = ISE-ISS
   ELSE
    iis = ISE ; iie = ISS ; iistep = ISS-ISE
   ENDIF
   
   IF(BIy(n) .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF

   DO i = iis, iie, iistep  
   DO j = jjs, jje, jjstep
   
    IF( icheck(i,j,k) .ne. 1 ) THEN
   
    Rtemp = sqrt( (xca(i)-BIx(n))**2 + (yca(j)-BIy(n))**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax-1) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
!     write(*,*) '*',i,j,k,Rtemp
     ENDIF
    ENDIF
	
	ENDIF
   
   ENDDO
   ENDDO
   
   ISS = ISS -1 ; ISE = ISE +1
   JSS = JSS -1 ; JSE = JSE +1
   KSS = KSS -1 ; KSE = KSE +1
   
  ENDDO ! Nmax
  
  Rmax(n) = Rnow
  
!  write(*,*) ' Ghost cell # ',n, ' Rmax = ',Rnow
  
 ENDDO ! nGhost
 
 DEALLOCATE(icheck)

END SUBROUTINE APM_search_2D
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_search_3D(gmap)
 USE flow_parameters
 USE APM_arrays
 USE GCM_arrays
 USE boundary_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 INTEGER gmap(mxa,mya,mza)
 INTEGER ISS, ISE, JSS, JSE, KSS, KSE
 INTEGER imin, imax, jmin, jmax, kmin, kmax
 INTEGER iis, iie, jjs, jje, kks, kke, iistep, jjstep, kkstep
 INTEGER iGC, jGC, kGC
 INTEGER Nnow
 INTEGER n, i, j, k
 REAL(KIND=CGREAL) xGC, yGC, zGC, dlength, Rnow, Rtemp, Rmin
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: icheck
 
 ALLOCATE(icheck(mxa,mya,mza))
 
 imin = ias
 jmin = jas
 kmin = kas
 
 imax = iae-1
 jmax = jae-1
 kmax = kae-1
 
 DO n=1, nGhost
 
  Nnow = 0
  Rnow = 0.0
  
  iGC = iGhost(n) + IFS
  jGC = jGhost(n) + JFS
  kGC = kGhost(n) + KFS
  
  xGC = xca(iGC)
  yGC = yca(jGC)
  zGC = zca(kGC)
  
  BIx(n) = xBodyIntercept(n)
  BIy(n) = yBodyIntercept(n)
  BIz(n) = zBodyIntercept(n)
  
  !! dlength = sqrt( (BIx(n) - xGC)**2 + (BIy(n) - yGC)**2 + (BIz(n) - zGC)**2 )
  
  BInx(n) =  xBodyInterceptNorm(n) !(BIx(n)-xGC)/dlength
  BIny(n) =  yBodyInterceptNorm(n) !(BIy(n)-yGC)/dlength
  BInz(n) =  zBodyInterceptNorm(n) !(BIz(n)-zGC)/dlength
  
  ISS = iGC -1 ; ISE = iGC +1
  JSS = jGC -1 ; JSE = jGC +1
  KSS = kGC -1 ; KSE = kGC +1
  
  Rmin = 0.0
  
  icheck = 0
  
  DO WHILE( Nnow .lt. Nmax-1 ) 
   
   ISS = max(ISS,imin)
   ISE = min(ISE,imax)
   JSS = max(JSS,jmin)
   JSE = min(JSE,jmax)
   KSS = max(KSS,kmin)
   KSE = min(KSE,kmax)
   
   ! for 3D
   
  
   IF(BIx(n) .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = 1
   ELSE
    iis = ISE ; iie = ISS ; iistep = -1
   ENDIF
   
   IF(BIy(n) .le. yGC) THEN
    jjs = JSS ; jje = JSE ; jjstep = JSE-JSS
   ELSE
    jjs = JSE ; jje = JSS ; jjstep = JSS-JSE
   ENDIF
   
   IF(BIz(n) .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = 1
   ELSE
    kks = KSE ; kke = KSS ; kkstep = -1
   ENDIF
   
   
   DO j = jjs, jje, jjstep
   DO k = kks, kke, kkstep
   DO i = iis, iie, iistep
       
	IF( icheck(i,j,k) .ne. 1) THEN
    Rtemp = sqrt( (xca(i)-BIx(n))**2 + (yca(j)-BIy(n))**2 + (zca(k)-BIz(n))**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax-1) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
!     write(*,*) '*',i,j,k,Rtemp
     ENDIF
    ENDIF
	ENDIF
    
   ENDDO
   ENDDO
   ENDDO
   
   
   
   IF(BIx(n) .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = ISE-ISS
   ELSE
    iis = ISE ; iie = ISS ; iistep = ISS-ISE
   ENDIF
   
   IF(BIy(n) .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF
   
   IF(BIz(n) .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = 1
   ELSE
    kks = KSE ; kke = KSS ; kkstep = -1
   ENDIF

   DO i = iis, iie, iistep  
   DO k = kks, kke, kkstep
   DO j = jjs, jje, jjstep
   
    IF(icheck(i,j,k) .ne. 1) THEN
    Rtemp = sqrt( (xca(i)-BIx(n))**2 + (yca(j)-BIy(n))**2 + (zca(k)-BIz(n))**2 )     
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax-1) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
!     write(*,*) '*',i,j,k,Rtemp
     ENDIF
    ENDIF
	ENDIF
   
   ENDDO
   ENDDO
   ENDDO
   
   
   IF(BIx(n) .le. xGC) THEN
    iis = ISS+1 ; iie = ISE-1 ; iistep = 1
   ELSE
    iis = ISE-1 ; iie = ISS+1 ; iistep = -1
   ENDIF
   
   IF(BIy(n) .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF
   
   IF(BIz(n) .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = KSE-KSS
   ELSE
    kks = KSE ; kke = KSS ; kkstep = KSS-KSE
   ENDIF

   DO k = kks, kke, kkstep
   DO j = jjs, jje, jjstep
   DO i = iis, iie, iistep  
    
	IF(icheck(i,j,k) .ne. 1) THEN
    Rtemp = sqrt( (xca(i)-BIx(n))**2 + (yca(j)-BIy(n))**2 + (zca(k)-BIz(n))**2 )      
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax-1) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
!     write(*,*) '*',i,j,k,Rtemp
     ENDIF
    ENDIF
	ENDIF
   
   ENDDO
   ENDDO
   ENDDO
   
   ISS = ISS -1 ; ISE = ISE +1
   JSS = JSS -1 ; JSE = JSE +1
   KSS = KSS -1 ; KSE = KSE +1
   
  ENDDO ! Nmax
  
  Rmax(n) = Rnow
  
!  write(*,*) ' Ghost cell # ',n, ' Rmax = ',Rnow
  
 ENDDO ! nGhost
 
 DEALLOCATE(icheck)

END SUBROUTINE APM_search_3D
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_SolveMat(Norder,minRank,Condmax)
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE   :: xn, yn, zn, WORK1, sigma
 REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: VanMat, AMat, Umat, VTmat, Smat, Bmat
 REAL(KIND=CGREAL) xGC, yGC, zGC, wn, sum, SVDtol, S1, SRank, Condmax

 INTEGER n, l, m, j, INFO, Norder, Rank, minRank
 
 ALLOCATE(xn(Nmax), yn(Nmax), zn(Nmax))
 ALLOCATE(VanMat(Nmax,NCOL))
 ALLOCATE(Amat(Nmax,NCOL))
 ALLOCATE(Umat(Nmax,Nmax))
 ALLOCATE(VTmat(NCOL,NCOL))
 ALLOCATE(Smat(NCOL,Nmax))
 ALLOCATE(Bmat(NCOL,Nmax))
 ALLOCATE(sigma(NCOL))
 ALLOCATE(WORK1(5*Nmax))
 
 pi = acos(-1.)
 SVDtol = 1.e-10
 Condmax = 0.0
 minRank = NCOL
 
 DO n=1, nGhost
 
 xGC = xca(iGhost(n)+IFS)
 yGC = yca(jGhost(n)+JFS)
 zGC = zca(kGhost(n)+KFS)
 
 SELECT CASE (Ndim)
 
 CASE (2)
 
 xn(1) = (xGC - BIx(n))
 yn(1) = (yGC - BIy(n))
 zn(1) = 0.0
 
 DO l=2,Nmax
  xn(l) = (xca(i_data(n,l-1)) - BIx(n))
  yn(l) = (yca(j_data(n,l-1)) - BIy(n))
  zn(l) = 0.0
 ENDDO
 
 call VanMatrix2D(VanMat, xn, yn, Nmax, NCOL, Norder)
 
 CASE (3)
 
 xn(1) = (xGC - BIx(n))
 yn(1) = (yGC - BIy(n))
 zn(1) = (zGC - BIz(n))
 
 DO l=2,Nmax
  xn(l) = (xca(i_data(n,l-1)) - BIx(n))
  yn(l) = (yca(j_data(n,l-1)) - BIy(n))
  zn(l) = (zca(k_data(n,l-1)) - BIz(n))
 ENDDO
 
 call VanMatrix3D(VanMat, xn, yn, zn, Nmax, NCOL, Norder)
  
 END SELECT
 
 DO l=1,Nmax
 DO m=1,NCOL

 wn = 0.5*(1.+cos(pi*sqrt(xn(l)**2 + yn(l)**2 + zn(l)**2)/Rmax(n)))
 Amat(l,m) = wn*VanMat(l,m)
 
 ENDDO
 ENDDO
 
 call DGESVD('A','A',Nmax,NCOL,Amat,Nmax,sigma,Umat,Nmax,VTmat,NCOL,WORK1,5*Nmax,INFO)
 
 ! compute pesudo inverse
 
 ! sigma inverse
 
 Smat = 0.0
 S1 = sigma(1)
 
 DO l=1,NCOL
  IF(sigma(l).gt.SVDtol*S1) THEN
   Rank = l
   SRank = sigma(l)
   Smat(l,l) = 1./sigma(l)
  ENDIF
 ENDDO
 
 Condmax = MAX(Condmax,S1/SRank)
 minRank = MIN(minRank,Rank)
 
 ! pesudo inverse
 
 DO l=1,NCOL
 DO m=1,Nmax
  sum = 0.0
  DO j=1, NCOL
   sum = sum + VTmat(j,l)*Smat(j,m)
  ENDDO
  Bmat(l,m) = sum
 ENDDO
 ENDDO

 DO l=1,NCOL
 DO m=1,Nmax
  sum = 0.0
  DO j=1, Nmax
   sum = sum + Bmat(l,j)*Umat(m,j)
  ENDDO
  Smat(l,m) = sum
 ENDDO
 ENDDO
 
 
 DO l=1, 4
 DO m=1, Nmax
  wn = 0.5*(1.+cos(pi*sqrt(xn(m)**2 + yn(m)**2 + zn(m)**2)/Rmax(n)))
  APM_mat(n,l,m) = wn*Smat(l,m)
 ENDDO
 ENDDO

 ENDDO ! nGhost
 
 DEALLOCATE(xn,yn,zn,VanMat,AMat,Umat,VTmat,sigma,WORK1, Smat, Bmat)
 
END SUBROUTINE APM_SolveMat
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_setBC_DC(var,varBC)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) var(mxa,mya,mza)
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) varBC, sum
 
 
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS
  
  sum = 0.0
  
  DO j=2,Nmax
   sum = sum + APM_mat(n,1,j)*var(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
  ENDDO
  
  var(iGC,jGC,kGC) = (varBC - sum)/APM_mat(n,1,1)
 
 ENDDO

END SUBROUTINE APM_setBC_DC
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_setBC_NM(prs)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) prs(mxa,mya,mza)
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sum
  
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS
  
  sum = 0.0
  
  SELECT CASE(Ndim)
  
  CASE (2)
  
  ! c10
  DO j=2,Nmax
   sum = sum + APM_mat(n,2,j)*prs(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
  ENDDO
  
  ! c01
  DO j=2,Nmax
   sum = sum + APM_mat(n,3,j)*prs(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
  ENDDO
  
  prs(iGC,jGC,kGC) = (0.0 - sum)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n))
  
  CASE (3)
  
  ! c100
  DO j=2,Nmax
   sum = sum + APM_mat(n,2,j)*prs(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
  ENDDO
  
  ! c010
  DO j=2,Nmax
   sum = sum + APM_mat(n,3,j)*prs(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
  ENDDO
  
  ! c001
  DO j=2,Nmax
   sum = sum + APM_mat(n,4,j)*prs(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInz(n)
  ENDDO
  
  prs(iGC,jGC,kGC) = (0.0 - sum)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n)+APM_mat(n,4,1)*BInz(n))
  
  END SELECT

 ENDDO

END SUBROUTINE APM_setBC_NM
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_setBC_vel_SLIP(uvel,vvel,wvel)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(mxa,mya,mza) :: uvel, vvel, wvel
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sumU, sumV, sumW
 REAL(KIND=CGREAL) velnorm, utemp, vtemp, wtemp
 
  
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS
  
  ! Neumann velocity 
  
  sumU = 0.0
  sumV = 0.0
  sumW = 0.0
  
  SELECT CASE (Ndim)
  
  CASE (2)
  
  ! c10
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,2,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
   sumV = sumV + APM_mat(n,2,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
   sumW = sumW + APM_mat(n,2,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
  ENDDO
  
  ! c01
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,3,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
   sumV = sumV + APM_mat(n,3,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
   sumW = sumW + APM_mat(n,3,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
  ENDDO
  
  utemp = (0.0 - sumU)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n))
  vtemp = (0.0 - sumV)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n))
  wtemp = (0.0 - sumW)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n))
  
  CASE (3)
  
  ! c100
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,2,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
   sumV = sumV + APM_mat(n,2,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
   sumW = sumW + APM_mat(n,2,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInx(n)
  ENDDO
  
  ! c010
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,3,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
   sumV = sumV + APM_mat(n,3,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
   sumW = sumW + APM_mat(n,3,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BIny(n)
  ENDDO
  
  ! c001
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,4,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInz(n)
   sumV = sumV + APM_mat(n,4,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInz(n)
   sumW = sumW + APM_mat(n,4,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))*BInz(n)
  ENDDO
  
  utemp = (0.0 - sumU)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n)+APM_mat(n,4,1)*BInz(n))
  vtemp = (0.0 - sumV)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n)+APM_mat(n,4,1)*BInz(n))
  wtemp = (0.0 - sumW)/(APM_mat(n,2,1)*BInx(n)+APM_mat(n,3,1)*BIny(n)+APM_mat(n,4,1)*BInz(n))
  
  END SELECT
  
  velnorm = utemp*BInx(n) + vtemp*BIny(n) + wtemp*BInz(n)
  
  uvel(iGC,jGC,kGC) = utemp - velnorm*BInx(n)
  vvel(iGC,jGC,kGC) = vtemp - velnorm*BIny(n)
  wvel(iGC,jGC,kGC) = wtemp - velnorm*BInz(n)
  
  sumU = 0.0
  sumV = 0.0
  sumW = 0.0
  
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,1,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumV = sumV + APM_mat(n,1,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumW = sumW + APM_mat(n,1,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
  ENDDO
  
  utemp = (0.0 - sumU)/APM_mat(n,1,1)
  vtemp = (0.0 - sumV)/APM_mat(n,1,1)
  wtemp = (0.0 - sumW)/APM_mat(n,1,1)
  
  velnorm = utemp*BInx(n) + vtemp*BIny(n) + wtemp*BInz(n)
  
  uvel(iGC,jGC,kGC) = uvel(iGC,jGC,kGC) + velnorm*BInx(n)
  vvel(iGC,jGC,kGC) = vvel(iGC,jGC,kGC) + velnorm*BIny(n)
  wvel(iGC,jGC,kGC) = wvel(iGC,jGC,kGC) + velnorm*BInz(n)

 ENDDO !n

END SUBROUTINE APM_setBC_vel_SLIP
!-----------------------------------------------------------------------------------------------------------------------



SUBROUTINE APM_setBC_vel_NORM(uvel,vvel,wvel)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(mxa,mya,mza) :: uvel, vvel, wvel
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sumU, sumV, sumW
 REAL(KIND=CGREAL) velnorm, utemp, vtemp, wtemp
 
  
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS
  
  sumU = 0.0
  sumV = 0.0
  sumW = 0.0
  
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,1,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumV = sumV + APM_mat(n,1,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumW = sumW + APM_mat(n,1,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
  ENDDO
  
  utemp = (0.0 - sumU)/APM_mat(n,1,1)
  vtemp = (0.0 - sumV)/APM_mat(n,1,1)
  wtemp = (0.0 - sumW)/APM_mat(n,1,1)
  
  velnorm = utemp*BInx(n) + vtemp*BIny(n) + wtemp*BInz(n)
  
  uvel(iGC,jGC,kGC) =  velnorm*BInx(n)
  vvel(iGC,jGC,kGC) =  velnorm*BIny(n)
  wvel(iGC,jGC,kGC) =  velnorm*BInz(n)

 ENDDO !n

END SUBROUTINE APM_setBC_vel_NORM
!-----------------------------------------------------------------------------------------------------------------------



SUBROUTINE APM_setBC_vel_noslip(uvel,vvel,wvel)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(mxa,mya,mza) :: uvel, vvel, wvel
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sumU, sumV, sumW
 REAL(KIND=CGREAL) velnorm, utemp, vtemp, wtemp
 
! uGC = 0.0
! vGC = 0.0
! wGC = 0.0
  
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS
  
  sumU = 0.0
  sumV = 0.0
  sumW = 0.0
  
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,1,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumV = sumV + APM_mat(n,1,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumW = sumW + APM_mat(n,1,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
  ENDDO
  
  uvel(iGC,jGC,kGC) = (0.0 - sumU)/APM_mat(n,1,1)
  vvel(iGC,jGC,kGC) = (0.0 - sumV)/APM_mat(n,1,1)
  wvel(iGC,jGC,kGC) = (0.0 - sumW)/APM_mat(n,1,1)

 ENDDO !n

END SUBROUTINE APM_setBC_vel_noslip
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_setBC_vel0(uvel,vvel,wvel,fsmach)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL), DIMENSION(mxa,mya,mza) :: uvel, vvel, wvel
 REAL(KIND=CGREAL) fsmach
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sumU, sumV, sumW
 REAL(KIND=CGREAL) velnorm, utemp, vtemp, wtemp
 REAL(KIND=CGREAL) uBI, vBI, wBI
 
! uGC = 0.0
! vGC = 0.0
! wGC = 0.0
  
 DO n=1,nGhost
 
  iGC = iGhost(n)+IFS
  jGC = jGhost(n)+JFS
  kGC = kGhost(n)+KFS

  uBI = uBodyIntercept(n)*fsmach
  vBI = vBodyIntercept(n)*fsmach
  wBI = wBodyIntercept(n)*fsmach
  
  sumU = 0.0
  sumV = 0.0
  sumW = 0.0
  
  DO j=2,Nmax
   sumU = sumU + APM_mat(n,1,j)*uvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumV = sumV + APM_mat(n,1,j)*vvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
   sumW = sumW + APM_mat(n,1,j)*wvel(i_data(n,j-1),j_data(n,j-1),k_data(n,j-1))
  ENDDO
  
  uvel(iGC,jGC,kGC) = (uBI - sumU)/APM_mat(n,1,1)
  vvel(iGC,jGC,kGC) = (vBI - sumV)/APM_mat(n,1,1)
  wvel(iGC,jGC,kGC) = (wBI - sumW)/APM_mat(n,1,1)

 ENDDO !n

END SUBROUTINE APM_setBC_vel0
!-----------------------------------------------------------------------------------------------------------------------



SUBROUTINE VanMatrix2D(VanMat, xn, yn, Nmax, NCOL, Norder)

 IMPLICIT NONE

 INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=14,R=30)  
 INTEGER Nmax, NCOL, Norder
 REAL(KIND=CGREAL), DIMENSION(Nmax, NCOL) ::  VanMat
 REAL(KIND=CGREAL), DIMENSION(Nmax) :: xn, yn
 INTEGER m
 
 SELECT CASE (Norder)
 
 CASE (1)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = xn(m)*yn(m)
 ENDDO
 
 
 CASE (2)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = xn(m)**2
  VanMat(m,5) = yn(m)**2
  VanMat(m,6) = xn(m)*yn(m)
 ENDDO
 
 CASE (3)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = xn(m)**2
  VanMat(m,5) = yn(m)**2
  VanMat(m,6) = xn(m)**3
  VanMat(m,7) = yn(m)**3
  VanMat(m,8) = xn(m)*yn(m)
  VanMat(m,9) = xn(m)**2*yn(m)
  VanMat(m,10) = xn(m)*yn(m)**2
 ENDDO
 
 CASE (4)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = xn(m)**2
  VanMat(m,5) = yn(m)**2
  VanMat(m,6) = xn(m)**3
  VanMat(m,7) = yn(m)**3
  VanMat(m,8) = xn(m)**4
  VanMat(m,9) = yn(m)**4
  VanMat(m,10) = xn(m)*yn(m)
  VanMat(m,11) = xn(m)*yn(m)**2
  VanMat(m,12) = xn(m)**2*yn(m)
  VanMat(m,13) = xn(m)**2*yn(m)**2
  VanMat(m,14) = xn(m)*yn(m)**3
  VanMat(m,15) = xn(m)**3*yn(m)
 ENDDO
 
 CASE (5)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = xn(m)**2
  VanMat(m,5) = yn(m)**2
  VanMat(m,6) = xn(m)**3
  VanMat(m,7) = yn(m)**3
  VanMat(m,8) = xn(m)**4
  VanMat(m,9) = yn(m)**4
  VanMat(m,10) = xn(m)**5
  VanMat(m,11) = yn(m)**5
  VanMat(m,12) = xn(m)*yn(m)
  VanMat(m,13) = xn(m)*yn(m)**2
  VanMat(m,14) = xn(m)**2*yn(m)
  VanMat(m,15) = xn(m)**2*yn(m)**2
  VanMat(m,16) = xn(m)*yn(m)**3
  VanMat(m,17) = xn(m)**3*yn(m)
  VanMat(m,18) = xn(m)*yn(m)**4
  VanMat(m,19) = xn(m)**2*yn(m)**3
  VanMat(m,20) = xn(m)**3*yn(m)**2
  VanMat(m,21) = xn(m)**4*yn(m)
 ENDDO
 
 END SELECT
 
END SUBROUTINE VanMatrix2D 
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE VanMatrix3D(VanMat, xn, yn, zn, Nmax, NCOL, Norder)
 IMPLICIT NONE
 
 INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=14,R=30) 
 INTEGER Nmax, NCOL, Norder
 REAL(KIND=CGREAL), DIMENSION(Nmax, NCOL) ::  VanMat
 REAL(KIND=CGREAL), DIMENSION(Nmax) :: xn, yn, zn
 INTEGER m
 
 SELECT CASE (Norder)
 
 CASE (1)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = zn(m)
  VanMat(m,5) = xn(m)*yn(m)
  VanMat(m,6) = xn(m)*zn(m)
  VanMat(m,7) = yn(m)*zn(m)
  VanMat(m,8) = xn(m)*yn(m)*zn(m)
 ENDDO
 
 CASE (2)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = zn(m)
  VanMat(m,5) = xn(m)**2
  VanMat(m,6) = yn(m)**2
  VanMat(m,7) = zn(m)**2
  VanMat(m,8) = xn(m)*yn(m)
  VanMat(m,9) = xn(m)*zn(m)
  VanMat(m,10) = yn(m)*zn(m)
 ENDDO
 
 
 CASE (3)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = zn(m)
  VanMat(m,5) = xn(m)**2
  VanMat(m,6) = yn(m)**2
  VanMat(m,7) = zn(m)**2
  VanMat(m,8) = xn(m)**3
  VanMat(m,9) = yn(m)**3
  VanMat(m,10) = zn(m)**3
  VanMat(m,11) = xn(m)*yn(m)
  VanMat(m,12) = xn(m)*zn(m)
  VanMat(m,13) = yn(m)*zn(m)
  VanMat(m,14) = xn(m)*yn(m)*zn(m)
  VanMat(m,15) = xn(m)**2*yn(m)
  VanMat(m,16) = xn(m)**2*zn(m)
  VanMat(m,17) = xn(m)*yn(m)**2
  VanMat(m,18) = xn(m)*zn(m)**2
  VanMat(m,19) = yn(m)*zn(m)**2
  VanMat(m,20) = yn(m)**2*zn(m)
 ENDDO
 
 CASE (4)
 
 DO m = 1, Nmax
  VanMat(m,1) = 1
  VanMat(m,2) = xn(m)
  VanMat(m,3) = yn(m)
  VanMat(m,4) = zn(m)
  VanMat(m,5) = xn(m)**2
  VanMat(m,6) = yn(m)**2
  VanMat(m,7) = zn(m)**2
  VanMat(m,8) = xn(m)**3
  VanMat(m,9) = yn(m)**3
  VanMat(m,10) = zn(m)**3
  VanMat(m,11) = xn(m)**4
  VanMat(m,12) = yn(m)**4
  VanMat(m,13) = zn(m)**4
  VanMat(m,14) = xn(m)*yn(m)
  VanMat(m,15) = xn(m)*zn(m)
  VanMat(m,16) = yn(m)*zn(m)
  VanMat(m,17) = xn(m)*yn(m)*zn(m)
  VanMat(m,18) = xn(m)**2*yn(m)
  VanMat(m,19) = xn(m)**2*zn(m)
  VanMat(m,20) = xn(m)**2*yn(m)*zn(m)
  VanMat(m,21) = xn(m)*yn(m)**2
  VanMat(m,22) = yn(m)**2*zn(m)
  VanMat(m,23) = xn(m)*yn(m)**2*zn(m)
  VanMat(m,24) = xn(m)*zn(m)**2
  VanMat(m,25) = yn(m)*zn(m)**2
  VanMat(m,26) = xn(m)*yn(m)*zn(m)**2
  VanMat(m,27) = xn(m)**2*yn(m)**2
  VanMat(m,28) = xn(m)**2*zn(m)**2
  VanMat(m,29) = yn(m)**2*zn(m)**2
  VanMat(m,30) = xn(m)**3*yn(m)
  VanMat(m,31) = xn(m)**3*zn(m)
  VanMat(m,32) = xn(m)*yn(m)**3
  VanMat(m,33) = yn(m)**3*zn(m)
  VanMat(m,34) = xn(m)*zn(m)**3
  VanMat(m,35) = yn(m)*zn(m)**3
 ENDDO
 
 END SELECT
 
END SUBROUTINE VanMatrix3D 
!-----------------------------------------------------------------------------------------------------------------------

