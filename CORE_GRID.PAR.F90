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
!  Filename: CORE_GRID.PAR.F90
!  Latest Modification: December 10, 2008 (ver. 2.0.1)
!  Made by S. A. Mohsen Karimian
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     make_grid()
!     metrics()
! --------------------------------------------------------------------

!---------------------------------------------------------------------
! grid convention
!
!      ny+1---------------------------------------
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!       ny +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!     ^    +--++--+--+--+--+---+--+--+--+--+--++--+
!   dy|    |  ||  |  |  |  | * |  |  |  |  |  ||  |
!     -  j +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        2 +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        1 +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        0 ---------------------------------------
!          0  1   2        i                  nx  nx+1
!                          <--->
!                           dx
!
!---------------------------------------------------------------------



SUBROUTINE make_grid()

! ----------------------------------------------------------
!  This subroutine generates and/or reads a Cartesian grid.
! ----------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL):: delta
    INTEGER :: i,junk
 
 
 
    IF (xgrid_unif==UNIFORM_GRID) then
      delta  = xout/REAL(nxc_GLBL,KIND=CGREAL)
       DO i=1,nx_GLBL
        x(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=10,FILE='xgrid.dat')

      DO i=1,nx_GLBL
        read(10,*) junk,x(i)
      ENDDO

      CLOSE(UNIT=10)
    ENDIF
    x(0)         = -x(2)
    x(nx_GLBL+1) =  x(nx_GLBL) + x(nx_GLBL) - x(nx_GLBL-1)
    
    IF (ygrid_unif==UNIFORM_GRID) then
      delta  = yout/REAL(nyc_GLBL,KIND=CGREAL)
      DO i=1,ny_GLBL
        y(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=11,FILE='ygrid.dat')

      DO i=1,ny_GLBL
        read(11,*) junk,y(i)
      ENDDO

      close(UNIT=11)
    ENDIF
    y(0)         = -y(2)
    y(ny_GLBL+1) =  y(ny_GLBL) + y(ny_GLBL) - y(ny_GLBL-1)

    IF (zgrid_unif==UNIFORM_GRID) then
      delta  = zout/REAL(nzc,KIND=CGREAL)
      DO i=1,nz
        z(i) = REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=12,FILE='zgrid.dat')

      DO i=1,nz
        read(12,*) junk,z(i)
      ENDDO

      close(UNIT=12)
    ENDIF
    z(0)    = -z(2)
    z(nz+1) =  z(nz) + z(nz) - z(nz-1)



    xc     = 0.0_CGREAL
    yc     = 0.0_CGREAL
    zc     = 0.0_CGREAL

    dx     = 0.0_CGREAL
    dy     = 0.0_CGREAL
    dz     = 0.0_CGREAL

    dxc    = 0.0_CGREAL
    dyc    = 0.0_CGREAL
    dzc    = 0.0_CGREAL
 
    dxinv  = 0.0_CGREAL
    dyinv  = 0.0_CGREAL
    dzinv  = 0.0_CGREAL

    dxcinv = 0.0_CGREAL
    dycinv = 0.0_CGREAL
    dzcinv = 0.0_CGREAL
 
!   Loop on cells in X direction
!   ----------------------------
    DO i=0,nxc_GLBL+1
      xc(i)    = 0.5_CGREAL*(x(i+1)+x(i)) 
      dx(i)    =             x(i+1)-x(i)  
      dxinv(i) = 1.0_CGREAL/dx(i)
    ENDDO
   
!   Loop on cells in Y direction
!   ----------------------------
    DO i=0,nyc_GLBL+1
      yc(i)    = 0.5_CGREAL*(y(i+1)+y(i)) 
      dy(i)    =             y(i+1)-y(i)  
      dyinv(i) = 1.0_CGREAL/dy(i)
    ENDDO
   
!   Loop on cells in Z direction
!   ----------------------------
    DO i=0,nzc+1
      zc(i)    = 0.5_CGREAL*(z(i+1)+z(i)) 
      dz(i)    =             z(i+1)-z(i)  
      dzinv(i) = 1.0_CGREAL/dz(i)
    ENDDO

!   Loop on cells in X direction
!   ----------------------------
    DO i=1,nxc_GLBL+1
      dxc(i)    = xc(i)-xc(i-1)  
      dxcinv(i) = 1.0_CGREAL/dxc(i)
    ENDDO
   
!   Loop on cells in Y direction
!   ----------------------------
    DO i=1,nyc_GLBL+1
      dyc(i)    = yc(i)-yc(i-1)  
      dycinv(i) = 1.0_CGREAL/dyc(i)
    ENDDO
   
!   Loop on cells in Z direction
!   ----------------------------
    DO i=1,nzc+1
      dzc(i)    = zc(i)-zc(i-1)  
      dzcinv(i) = 1.0_CGREAL/dzc(i)
    ENDDO

END SUBROUTINE make_grid  
!---------------------------------------------------------------------



SUBROUTINE metrics()

! ----------------------------------------------------------
!  This subroutine set the interpolation metrics as follow:
!
!            |--*--|--*--|--*--|-----|
!              i-1 w  i  e i+1 
!
!        g(w) = fx(i)g(i) + (1-fx(i))*g(i-1)
! ----------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k
    
    REAL(KIND=CGREAL) :: Xmid, Dext, rxt



    DO i=1,nxc_GLBL+1
      fx(i) = ( x(i) - xc(i-1) )/( xc(i) - xc(i-1) )
    ENDDO

    DO j=1,nyc_GLBL+1
      fy(j) = ( y(j) - yc(j-1) )/( yc(j) - yc(j-1) )
    ENDDO

    DO k=1,nzc+1
      fz(k) = ( z(k) - zc(k-1) )/( zc(k) - zc(k-1) )
    ENDDO
    
!   --------------------------------------
!    SAMK: Adding an optional buffer zone
!   --------------------------------------
    Xmid=5.d-1*(Xext+Xout)
    Dext=Xout-Xext
    
    IF (extended_outflow==1) THEN
      DO i=0,nxc_GLBL+1
        IF (xc(i)<=Xext) THEN
          damper(i)=1.d0
        ELSE
!!          damper(i)=Xout-Dext*sin(5.d-1*(xc(i)    -Xext)*pi/Dext+5.d-1*pi)
!!          damper(i)=          sin(5.d-1*(damper(i)-Xext)*pi/Dext+5.d-1*pi)

          damper(i)=          sin(5.d-1*(xc(i)-Xext)*pi/Dext+5.d-1*pi)
          damper(i)=1.d0+(DampFact-1.d0)*(1.d0-damper(i))
        END IF
      END DO
    ELSE
      damper(:)=1.d0
    END IF

END SUBROUTINE metrics
!------------------------------------------------------------------------------
