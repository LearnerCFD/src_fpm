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
!  Filename: APM_fresh_arrays.PAR.F90
!  Latest Modification: Sep, 11 2010 (ver. P2.0.0)
!  Made by Jung-Hee, Seo
! --------------------------------------------------------------------
MODULE APM_fresh_arrays
 INTEGER  Nmax, NCOL
 INTEGER, DIMENSION(:,:),   ALLOCATABLE :: i_data, j_data, k_data
 REAL(KIND= SELECTED_REAL_KIND(P=14,R=30)),    DIMENSION(:),     ALLOCATABLE :: Rmax
 REAL(KIND= SELECTED_REAL_KIND(P=14,R=30)),    DIMENSION(:,:,:), ALLOCATABLE :: APM_mat
END MODULE APM_fresh_arrays
!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE APM_fresh_allocate_memory(N_max,N_order)
 USE flow_parameters
 USE APM_fresh_arrays
 USE GCM_arrays
 USE boundary_arrays
 
 IMPLICIT NONE
 
 INTEGER N_max, N_order, N_dim

 
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
 
 ALLOCATE(i_data(nFresh,Nmax),j_data(nFresh,Nmax),k_data(nFresh,Nmax))
 ALLOCATE(Rmax(nFresh))
 ALLOCATE(APM_mat(nFresh,1,Nmax))
 
END SUBROUTINE APM_fresh_allocate_memory
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_fresh_deallocate
 USE APM_fresh_arrays
 
 IMPLICIT NONE
 
 DEALLOCATE(i_data,j_data,k_data)
 DEALLOCATE(Rmax)
 DEALLOCATE(APM_mat)
 
END SUBROUTINE APM_fresh_deallocate
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_fresh_search_2D(gmap, fmap)
 USE flow_parameters
 USE APM_fresh_arrays
 USE GCM_arrays
 USE boundary_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 INTEGER gmap(mxa,mya,mza), fmap(mxa,mya,mza)
 INTEGER ISS, ISE, JSS, JSE, KSS, KSE
 INTEGER imin, imax, jmin, jmax, kmin, kmax
 INTEGER iis, iie, jjs, jje, kks, kke, iistep, jjstep, kkstep
 INTEGER iGC, jGC, kGC
 INTEGER Nnow
 INTEGER n, i, j, k
 REAL(KIND=CGREAL) xGC, yGC, zGC, dlength, Rnow, Rtemp, Rmin, xBI, yBI, zBI
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: icheck
 
 ALLOCATE(icheck(mxa,mya,mza))
 
 imin = ias
 jmin = jas
 kmin = kas
 
 imax = iae-1
 jmax = jae-1
 kmax = kae-1
 
  
 DO n=1, nFresh
 
  Nnow = 0
  Rnow = 0.0
  
  iGC = iFresh(n)+IFS
  jGC = jFresh(n)+JFS
  kGC = kFresh(n)+KFS
  
  xGC = xca(iGC)
  yGC = yca(jGC)
  zGC = zca(kGC)
  
  xBI = xBodyInterceptFresh(n)
  yBI = yBodyInterceptFresh(n)
  zBI = zBodyInterceptFresh(n)
  
  
  ISS = iGC -1 ; ISE = iGC +1
  JSS = jGC -1 ; JSE = jGC +1
  KSS = kGC -1 ; KSE = kGC +1
  
!  WRITE(*,*) iGC, jGC, kGC
  Rmin = 0.0
  
  icheck = 0
  
  DO WHILE( Nnow .lt. Nmax ) 
   
   ISS = max(ISS,imin)
   ISE = min(ISE,imax)
   JSS = max(JSS,jmin)
   JSE = min(JSE,jmax)
   KSS = max(KSS,kmin)
   KSE = min(KSE,kmax)
   
   ! for 2D
   
   k = 1+KFS
   
   IF(xBI .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = 1
   ELSE
    iis = ISE ; iie = ISS ; iistep = -1
   ENDIF
   
   IF(yBI .le. yGC) THEN
    jjs = JSS ; jje = JSE ; jjstep = JSE-JSS
   ELSE
    jjs = JSE ; jje = JSS ; jjstep = JSS-JSE
   ENDIF
   
   
   DO j = jjs, jje, jjstep
   DO i = iis, iie, iistep
    
	IF( icheck(i,j,k) .ne. 1 ) THEN
    Rtemp = sqrt( (xca(i)-xGC)**2 + (yca(j)-yGC)**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) .and. (fmap(i,j,k) .ne. 1) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
     ENDIF
    ENDIF
	ENDIF
    
   ENDDO
   ENDDO
   
   IF(xBI .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = ISE-ISS
   ELSE
    iis = ISE ; iie = ISS ; iistep = ISS-ISE
   ENDIF
   
   IF(yBI .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF

   DO j = jjs, jje, jjstep
   DO i = iis, iie, iistep
    
	IF( icheck(i,j,k) .ne. 1 ) THEN
    Rtemp = sqrt( (xca(i)-xGC)**2 + (yca(j)-yGC)**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) .and. (fmap(i,j,k) .ne. 1) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
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
  
 ENDDO ! nGhost
 
 DEALLOCATE(icheck)

END SUBROUTINE APM_fresh_search_2D
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_fresh_search_3D(gmap, fmap)
 USE flow_parameters
 USE APM_fresh_arrays
 USE GCM_arrays
 USE boundary_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 INTEGER gmap(mxa,mya,mza), fmap(mxa,mya,mza)
 INTEGER ISS, ISE, JSS, JSE, KSS, KSE
 INTEGER imin, imax, jmin, jmax, kmin, kmax
 INTEGER iis, iie, jjs, jje, kks, kke, iistep, jjstep, kkstep
 INTEGER iGC, jGC, kGC
 INTEGER Nnow
 INTEGER n, i, j, k
 REAL(KIND=CGREAL) xGC, yGC, zGC, dlength, Rnow, Rtemp, Rmin, xBI, yBI, zBI
 INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: icheck
 
 ALLOCATE(icheck(mxa,mya,mza))
 
 imin = ias
 jmin = jas
 kmin = kas
 
 imax = iae-1
 jmax = jae-1
 kmax = kae-1
 
  
 DO n=1, nFresh
 
  Nnow = 0
  Rnow = 0.0
  
  iGC = iFresh(n)+IFS
  jGC = jFresh(n)+JFS
  kGC = kFresh(n)+KFS
  
  xGC = xca(iGC)
  yGC = yca(jGC)
  zGC = zca(kGC)
  
  xBI = xBodyInterceptFresh(n)
  yBI = yBodyInterceptFresh(n)
  zBI = zBodyInterceptFresh(n)
  
  ISS = iGC -1 ; ISE = iGC +1
  JSS = jGC -1 ; JSE = jGC +1
  KSS = kGC -1 ; KSE = kGC +1
  
  Rmin = 0.0
  icheck = 0
  
  DO WHILE( Nnow .lt. Nmax ) 
   
   ISS = max(ISS,imin)
   ISE = min(ISE,imax)
   JSS = max(JSS,jmin)
   JSE = min(JSE,jmax)
   KSS = max(KSS,kmin)
   KSE = min(KSE,kmax)
   
   ! for 3D
   
  
   IF(xBI .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = 1
   ELSE
    iis = ISE ; iie = ISS ; iistep = -1
   ENDIF
   
   IF(yBI .le. yGC) THEN
    jjs = JSS ; jje = JSE ; jjstep = JSE-JSS
   ELSE
    jjs = JSE ; jje = JSS ; jjstep = JSS-JSE
   ENDIF
   
   IF(zBI .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = 1
   ELSE
    kks = KSE ; kke = KSS ; kkstep = -1
   ENDIF
   
   
   DO j = jjs, jje, jjstep
   DO k = kks, kke, kkstep
   DO i = iis, iie, iistep
    
	IF( icheck(i,j,k) .ne. 1 ) THEN
    Rtemp = sqrt( (xca(i)-xGC)**2 + (yca(j)-yGC)**2 + (zca(k)-zGC)**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) .and. (fmap(i,j,k) .ne. 1) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
     ENDIF
    ENDIF
	ENDIF
    
   ENDDO
   ENDDO
   ENDDO
   
   
   
   IF(xBI .le. xGC) THEN
    iis = ISS ; iie = ISE ; iistep = ISE-ISS
   ELSE
    iis = ISE ; iie = ISS ; iistep = ISS-ISE
   ENDIF
   
   IF(yBI .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF
   
   IF(zBI .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = 1
   ELSE
    kks = KSE ; kke = KSS ; kkstep = -1
   ENDIF

   DO j = jjs, jje, jjstep
   DO k = kks, kke, kkstep
   DO i = iis, iie, iistep
    
	IF( icheck(i,j,k) .ne. 1 ) THEN
    Rtemp = sqrt( (xca(i)-xGC)**2 + (yca(j)-yGC)**2 + (zca(k)-zGC)**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) .and. (fmap(i,j,k) .ne. 1) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
     ENDIF
    ENDIF
	ENDIF
    
   ENDDO
   ENDDO
   ENDDO
   
   
   IF(xBI .le. xGC) THEN
    iis = ISS+1 ; iie = ISE-1 ; iistep = 1
   ELSE
    iis = ISE-1 ; iie = ISS+1 ; iistep = -1
   ENDIF
   
   IF(yBI .le. yGC) THEN
    jjs = JSS+1 ; jje = JSE-1 ; jjstep = 1
   ELSE
    jjs = JSE-1 ; jje = JSS+1 ; jjstep = -1
   ENDIF
   
   IF(zBI .le. zGC) THEN
    kks = KSS ; kke = KSE ; kkstep = KSE-KSS
   ELSE
    kks = KSE ; kke = KSS ; kkstep = KSS-KSE
   ENDIF

   DO j = jjs, jje, jjstep
   DO k = kks, kke, kkstep
   DO i = iis, iie, iistep
    
	IF( icheck(i,j,k) .ne. 1 ) THEN
    Rtemp = sqrt( (xca(i)-xGC)**2 + (yca(j)-yGC)**2 + (zca(k)-zGC)**2 )    
    IF((gmap(i,j,k) .ne. 1) .and. (Rtemp .ge. Rmin) .and. (fmap(i,j,k) .ne. 1) ) THEN
     Nnow = Nnow + 1
     IF(Nnow .le. Nmax) THEN
     i_data(n,Nnow) = i ; j_data(n,Nnow) = j ; k_data(n,Nnow) = k
     Rnow = max(Rnow, Rtemp)
	 icheck(i,j,k) = 1
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
  
 ENDDO ! nGhost
 
 DEALLOCATE(icheck)

END SUBROUTINE APM_fresh_search_3D
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_fresh_SolveMat(Norder,minRank,Condmax)
 USE flow_parameters
 USE GCM_arrays
 USE APM_fresh_arrays
 USE boundary_arrays
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
 
 DO n=1, nFresh
 
 xGC = xca(iFresh(n)+IFS)
 yGC = yca(jFresh(n)+JFS)
 zGC = zca(kFresh(n)+KFS)
 
 SELECT CASE (Ndim)
 
 CASE (2)
 
 DO l=1,Nmax
  xn(l) = (xca(i_data(n,l)) - xGC)
  yn(l) = (yca(j_data(n,l)) - yGC)
  zn(l) = 0.0
 ENDDO
 
 call VanMatrix2D(VanMat, xn, yn, Nmax, NCOL, Norder)
 
 CASE (3)
 
 DO l=1,Nmax
  xn(l) = (xca(i_data(n,l)) - xGC)
  yn(l) = (yca(j_data(n,l)) - yGC)
  zn(l) = (zca(k_data(n,l)) - zGC)
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
 
 
 DO l=1, 1
 DO m=1, Nmax
  wn = 0.5*(1.+cos(pi*sqrt(xn(m)**2 + yn(m)**2 + zn(m)**2)/Rmax(n)))
  APM_mat(n,l,m) = wn*Smat(l,m)
 ENDDO
 ENDDO

 ENDDO ! nGhost
 
 DEALLOCATE(xn,yn,zn,VanMat,AMat,Umat,VTmat,sigma,WORK1, Smat, Bmat)
 
END SUBROUTINE APM_fresh_SolveMat
!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE APM_set_Fresh(var)
 USE global_parameters
 USE flow_parameters
 USE GCM_arrays
 USE APM_fresh_arrays
 USE boundary_arrays
 USE Agrid
 
 IMPLICIT NONE
 
 REAL(KIND=CGREAL) var(mxa,mya,mza)
 INTEGER n, j, iGC, jGC, kGC
 REAL(KIND=CGREAL) sum
 
 
 DO n=1,nFresh
 
  iGC = iFresh(n)+IFS
  jGC = jFresh(n)+JFS
  kGC = kFresh(n)+KFS
  
  sum = 0.0
  
  DO j=1,Nmax
   sum = sum + APM_mat(n,1,j)*var(i_data(n,j),j_data(n,j),k_data(n,j))
!   WRITE(*,*) i_data(n,j), j_data(n,j), k_data(n,j), var(i_data(n,j),j_data(n,j),k_data(n,j))
  ENDDO
  
  var(iGC,jGC,kGC) = sum
!  WRITE(*,*) iGC,jGC,kGC,' Fresh cell value = ', var(iGC,jGC,kGC)

  
 ENDDO

END SUBROUTINE
!-----------------------------------------------------------------------------------------------------------------------

