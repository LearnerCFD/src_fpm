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
!  Filename: CORE_POISSON.PAR.F90
!  Latest Modification: October 29, 2008 (ver. P2.0.0)
!  Made by S. A. Mohsen Karimian

!  Latest Modification: Oct, 20 2010 (ver. PAT 1.2.0)
!  by JHSeo
! --------------------------------------------------------------------

! --------------------------------------------------------------------
!  This file contains the following subroutines:
!     rhs_poisson(sum)
!     solve_poisson()
!     update_pressure_freshcells(var)
!     calc_residual(var,r,resm)
!     itsolv(var,r)
!     subtract_mean_pressure()
! --------------------------------------------------------------------



! Compile-time function definitions
! ---------------------------------
# define L2GI(i)       myIs+i-1
# define L2GJ(j)       myJs+j-1


SUBROUTINE rhs_poisson(sum) 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
	USE multiuse_arrays
	USE unstructured_surface_arrays
	USE cutcell_arrays
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)
	USE scalar	
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END)
    IMPLICIT NONE
    
    REAL(KIND=CGREAL), INTENT(OUT)    :: sum
    
    INTEGER :: i, j, k, iG, jG, iKill
	REAL(KIND=CGREAL) usigma


!    rhs = [ d(U)/dx + d(V)/dy + d(W)/dz ] / dt
!   --------------------------------------------
    sum = 0.0_CGREAL

    DO k = 1,nzc
    DO j = 1,nyc
      jG=L2GJ(j)

      DO i = 1,nxc
        iG=L2GI(i)
		
		IF(iCC) call Cutcell_norm_vel(i,j,k,usigma)
		
		iKill = isolid(i,j,k)*iCC + iblank(i,j,k)*(1-iCC)
		
        nlu(i,j,k)=((face_ue(i,j,k)*AREA_E(i,j,k)-face_uw(i,j,k)*AREA_W(i,j,k))*dxinv(iG) &
                   +(face_vn(i,j,k)*AREA_N(i,j,k)-face_vs(i,j,k)*AREA_S(i,j,k))*dyinv(jG) &
                   +(face_wf(i,j,k)*AREA_F(i,j,k)-face_wb(i,j,k)*AREA_B(i,j,k))*dzinv(k ) &
                   + CUT_AREA(i,j,k)*usigma*dxinv(iG)*dyinv(jG)*dzinv(k) )*REAL(1-iKill,KIND=CGREAL)                  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(START)                    
        zero(i,j,k) = 0.0_CGREAL				  
!     ----------------------------------------------NEWLY ADDED BY ZHANG CHAO(END) 

        sum = sum + nlu(i,j,k)*dx(iG)*dy(jG)*dz(k)
      ENDDO
    ENDDO
    ENDDO

#ifdef MPI
    CALL par_GetSumReal(sum)
#endif

 IF(iCC) THEN

 IF(monitorON) write(*,*) 'Sum before : ', sum
 
 call mix_Conserve(nlu)
 
 ENDIF ! iCC
 
 sum = 0.0_CGREAL
 
 do i=1,nxc
 do j=1,nyc
 do k=1,nzc
  iG=L2GI(i)
  jG=L2GJ(j)
  sum = sum + nlu(i,j,k)*dx(iG)*dy(jG)*dz(k)*(1.0-iblank_solid(i,j,k))
  nlu(i,j,k) = nlu(i,j,k)*dtinv  
 enddo
 enddo
 enddo
 
#ifdef MPI
    CALL par_GetSumReal(sum)
#endif

 IF(monitorON) write(*,*) 'Sum RHS : ', sum
 

END SUBROUTINE rhs_poisson
!-------------------------------------------------------------------------------




SUBROUTINE solve_poisson(totComm)

    USE global_parameters
    USE flow_parameters
    USE pressure_arrays
    USE grid_arrays
    USE multiuse_arrays


#   ifdef MPI
      use par_parameters
#   endif

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT(OUT) :: totComm



!   Local variables
!   ---------------
    INTEGER           :: iter,loc(3)
    REAL(KIND=CGREAL) :: maxres
    real(kind=CGREAL) :: compTime,commTime,compPerc,commPerc

    integer:: i,j,k
    INTEGER :: ierr


    IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime > 1) &
      CALL update_pressure_freshcells(pPrime)!...................COMPLETE(SAMK)

    IF (pbcx1 == PBC_DIRICHLET .OR. &   !
        pbcy1 == PBC_DIRICHLET .OR. &   ! Added by H. Luo
        pbcz1 == PBC_DIRICHLET .OR. &   !
        pbcx2 == PBC_DIRICHLET .OR. &   ! Modified by SAMK
        pbcy2 == PBC_DIRICHLET .OR. &   !
        pbcz2 == PBC_DIRICHLET) THEN    !

      CALL set_pressure_dirichlet_bc()!..........................COMPLETE(SAMK)
      CALL update_pressure_dirichlet_bc()
    END IF

    iter   = 0
    maxres = 1.0E10_CGREAL 
    compTime=0.
    commTime=0.
    compPerc=0.
    commPerc=0.
    totComm =0.

    CALL calc_residual(pPrime,nlu,maxres,loc)!...................COMPLETE(SAMK)


    IF (monitorON) THEN
      write(STDOUT,'(5X,A)') '------------------------------------------------------------------------'
      write(STDOUT,'(5X,A)') 'iter.     Residual        Position     Comp. Time(sec.) Comm. Time(sec.)'
      write(STDOUT,'(5X,A)') '----- --------------- ---------------- ---------------- ----------------'
      WRITE(STDOUT,'(5X,I5.5,X,E15.8,X,A,I4.4,A,I4.4,A,I4.4,A)') &
      iter,maxres,'(',loc(1),',',loc(2),',',loc(3),')'
    END IF ! ntime 

    SELECT CASE(it_solver_type)

    CASE (IT_SOLVER_TYPE_LSOR)

      DO WHILE ((iter .LT. iterMaxPoisson) .AND. (maxres .GT. restolPoisson))

        CALL itsolv(pPrime,nlu,compTime)

#       ifdef MPI
          call par_comm_POISSON(commTime)
          totComm=totComm+commTime
#       endif
!!        CALL update_pressure_dirichlet_bc()
!!        IF (boundary_formulation == GCM_METHOD) THEN
!!          CALL set_outer_ghost_pres(pPrime)
!!          CALL GCM_p_set_bc_internal(pPrime)
!!          CALL GCM_enforce_p_compatibility(pPrime)
!!        ENDIF

!!        IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
!!            bcy1 == BC_TYPE_PERIODIC .OR. &
!!            bcz1 == BC_TYPE_PERIODIC) THEN 
!!          CALL enforce_p_periodic(pPrime)
!!        END IF

        IF (bcx1 == BC_TYPE_PERIODIC .OR. &         !Ehsan added for Periodic BC
            bcy1 == BC_TYPE_PERIODIC .OR. &         
            bcz1 == BC_TYPE_PERIODIC) THEN        
          do i=0,nxc+1                               
            do j=0,nyc+1                              
              pPrime(i,j,0) = pPrime(i,j,nzc)       
              pPrime(i,j,nzc+1) = pPrime(i,j,1)          
            enddo                                   
          enddo                                

        END IF                             

        CALL calc_residual(pPrime,nlu,maxres,loc)

        iter = iter + 1
        
        compPerc=100.*compTime/(compTime+commTime)
        commPerc=100.*commTime/(compTime+commTime)
        IF (monitorON) &
          WRITE(STDOUT,'(5X,I5.5,X,E15.8,X,A,I4.4,A,I4.4,A,I4.4,A,2(X,1F7.3,A,F6.1,A))') &
            iter,maxres,'(',loc(1),',',loc(2),',',loc(3),')',compTime,'(',compPerc,'%)',commTime,'(',commPerc,'%)'

      ENDDO

    CASE (IT_SOLVER_TYPE_MG) 
            
      DO WHILE ((iter .LT. iterMaxPoisson) .AND. (maxres .GT. restolPoisson))


        CALL mg_solver(pPrime,nlu,compTime,commTime)!.................ONLY_PER(SAMK)
       

        CALL update_pressure_dirichlet_bc() 
        IF (boundary_formulation == GCM_METHOD) THEN
          CALL set_outer_ghost_pres(pPrime)!..........................ONLY_PER(SAMK)
          CALL GCM_p_set_bc_internal(pPrime)!....................COMPLETE(SAMK)
          CALL GCM_enforce_p_compatibility(pPrime)!..............COMPLETE(SAMK)
        ENDIF ! boundary_formulation
  
 
!!        IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
!!            bcy1 == BC_TYPE_PERIODIC .OR. &
!!            bcz1 == BC_TYPE_PERIODIC) THEN 
!!          CALL enforce_p_periodic(pPrime)
!!        END IF

        IF (bcx1 == BC_TYPE_PERIODIC .OR. &                                  !Ehsan added for Periodic BC
            bcy1 == BC_TYPE_PERIODIC .OR. &                              
            bcz1 == BC_TYPE_PERIODIC) THEN                              
          do i=0,nxc+1                             
            do j=0,nyc+1                             
              pPrime(i,j,0) = pPrime(i,j,nzc)            
              pPrime(i,j,nzc+1) = pPrime(i,j,1)        
            enddo                                  
          enddo                                  
        END IF                               


        CALL calc_residual(pPrime,nlu,maxres,loc)!...............COMPLETE(SAMK)

        iter = iter + 1
        
        compPerc=100.*compTime/(compTime+commTime)
        commPerc=100.*commTime/(compTime+commTime)
        IF (monitorON) &
          WRITE(STDOUT,'(5X,I5.5,X,E15.8,X,A,I4.4,A,I4.4,A,I4.4,A,2(X,1F7.3,A,F6.1,A))') &
            iter,maxres,'(',loc(1),',',loc(2),',',loc(3),')',compTime,'(',compPerc,'%)',commTime,'(',commPerc,'%)'

!DEBUG
!        IF (monitorON) write(299,*) iter,maxres
!DEBUG
      ENDDO

    END SELECT

!   Add the ups and ums back for periodic condition in certain direction
!   --------------------------------------------------------------------
!!    IF (bcx1 == BC_TYPE_PERIODIC .OR. & 
!!        bcy1 == BC_TYPE_PERIODIC .OR. &
!!        bcz1 == BC_TYPE_PERIODIC) THEN 
!!      CALL add_up_um 
!!    END IF

!   Add the ups and ums back for Dirichlet BC
!   -----------------------------------------
    IF (pbcx1 == PBC_DIRICHLET .OR. &   !
        pbcy1 == PBC_DIRICHLET .OR. &   ! Added by H. Luo
        pbcz1 == PBC_DIRICHLET .OR. &   !
        pbcx2 == PBC_DIRICHLET .OR. &   ! Modified by SAMK
        pbcy2 == PBC_DIRICHLET .OR. &   !
        pbcz2 == PBC_DIRICHLET) THEN    !
 
      CALL unset_pressure_dirichlet_bc()!........................COMPLETE(SAMK)
    END IF

    IF (monitorON) THEN
      WRITE(STDOUT,'(5X,A)') '------------------------------------------------------------------------'

      IF (maxres > restolPoisson ) THEN
        WRITE(STDOUT,'(5X,A,I5,A)') 'WARNING: Pressure did not converge in ',itermaxPoisson,' iterations!'
        WRITE(STDOUT,*)
      END IF
    ENDIF




END SUBROUTINE solve_poisson
!---------------------------------------------------------------------



SUBROUTINE update_pressure_freshcells(var)

! -------------------------------------------------------------------------------------------------
!  A fresh-cell was a ghost-cell in the previous time step and its value is stored in pGhost array
!  At current time-step we need to transfer this value to pPrime so it can be correctly used
!  as initial guess for current time-step. POC RM
! -------------------------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE
    
    REAL(KIND=CGREAL), INTENT (INOUT) ::var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)

    REAL(KIND=CGREAL) :: denom, myTime

    INTEGER :: i,j,k,ii,jj,kk,n,iBody
    INTEGER :: iG,jG,kG,ib,ibG,ibF
    INTEGER :: iBeg,iEnd,jBeg,jEnd,kBeg,kEnd



    IF (boundary_formulation == GCM_METHOD) THEN

!     estimate for fresh cell pressure value
!     --------------------------------------
      DO k=1,nzc
      DO j=myJLL,myJUL
      DO i=myILL,myIUL
        IF (fresh_cell(i,j,k) == 1) THEN
           var(i,j,k)    = pGhost(i,j,k) 
           pGhost(i,j,k) = 0.0_CGREAL
        ENDIF
      ENDDO
      ENDDO
      ENDDO

!     Now obtain estimate of pressure for newly created ghost cell.
!     -------------------------------------------------------------
      DO n=1,nGhost
        iG = iGhost(n)
        jG = jGhost(n)
        kG = kGhost(n)
        
        IF (iG<1 .OR. iG>nxc .OR. &
            jG<1 .OR. jG>nyc .OR. &
            kG<1 .OR. kG>nzc) THEN
          IF (ImTheBOSS) WRITE(*,*) 'Error in ghost cell number ',n,'. Location out of domain. Terminating...'
          CALL flow_stop
          STOP
        END IF

        iBody = bodyNum(iG,jG,kG)

        IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
            unstruc_surface_type(iBody) /= MEMBRANE) THEN
          IF (ImTheBOSS) WRITE(*,*) 'A ghost cell is assigned to a non-solid, non-membrane body! Terminating ...'
          CALL flow_stop
          stop
        END IF

!       this identifies newly created ghost-cell:
!       -----------------------------------------
        IF ( ABS(pGhost(iG,jG,kG)) <= 1.0E-10 ) THEN   
          
          denom = 0.0_CGREAL
          
          iBeg = MAX(myILL,iG-1)
          iEnd = MIN(myIUL,iG+1)
          jBeg = MAX(myJLL,jG-1)
          jEnd = MIN(myJUL,jG+1)
          kBeg = MAX(1    ,kG-1)
          kEnd = MIN(nzc  ,kG+1)

          DO kk=kBeg,kEnd
          DO jj=jBeg,jEnd
          DO ii=iBeg,iEnd
            IF ( unstruc_surface_type(iBody) == SOLID_BODY) THEN
              ibG = iblank_solid(iG,jG,kG) 
              ib  = iblank_solid(ii,jj,kk) 
            ELSE
              IF (ALLOCATED(iblank_memb)) THEN
                ibG = iblank_memb(iG,jG,kG,iBody-nBody_solid)
                ib  = iblank_memb(ii,jj,kk,iBody-nBody_solid)
              ENDIF
            ENDIF
            IF (ibG /= ib) THEN  

!             pGhost is obtained as mean of adjacent fluid nodal values
!             ---------------------------------------------------------
              pGhost(iG,jG,kG) =  pGhost(iG,jG,kG) + var(ii,jj,kk) 
              denom            = denom + 1.0_CGREAL
            ENDIF
          ENDDO
          ENDDO
          ENDDO

          pGhost(iG,jG,kG) =  pGhost(iG,jG,kG)/denom
   
        ENDIF ! pGhost

      ENDDO ! nGhost

    ELSE  ! SSM_METHOD

!     estimate for fresh cell pressure value
!     --------------------------------------
      DO k=1,nzc
      DO j=myJmin,myJmax
      DO i=myImin,myImax

        IF (fresh_cell(i,j,k) == 1) THEN
          iBody = bodyNum(i,j,k)
          IF (unstruc_surface_type(iBody) /= SOLID_BODY .AND. &
              unstruc_surface_type(iBody) /= MEMBRANE) THEN
            IF (ImTheBOSS) WRITE(*,*) 'A ghost cell is assigned to a non-solid, non-membrane body at', i, j, k,' Terminating ...'
            IF (ImTheBOSS) WRITE(*,*) 'fresh_cell(i, j, k) =', fresh_cell(i, j, k)
            CALL write_dump()
            CALL flow_stop
            stop
          END IF

          denom      = 0.0_CGREAL
          var(i,j,k) = 0.0_CGREAL

          iBeg = MAX(myILL,i-1)
          iEnd = MIN(myIUL,i+1)
          jBeg = MAX(myJLL,j-1)
          jEnd = MIN(myJUL,j+1)
          kBeg = MAX(1    ,k-1)
          kEnd = MIN(nzc  ,k+1)

          DO kk=kBeg,kEnd
          DO jj=jBeg,jEnd
          DO ii=iBeg,iEnd
            IF (unstruc_surface_type(iBody) == SOLID_BODY) THEN
              ibF = iblank_solid(i,j,k)
              ib  = iblank_solid(ii,jj,kk)
            ELSE
              IF ( ALLOCATED(iblank_memb) .EQV. .TRUE. ) THEN
                ibF = iblank_memb(i,j,k,iBody-nBody_solid)
                ib  = iblank_memb(ii,jj,kk,iBody-nBody_solid)
              ENDIF
            ENDIF
            IF (ibF == ib .AND. fresh_cell(ii,jj,kk) /= 1) THEN  
!             fresh cell value is mean of adjacent fluid nodal values.
!             --------------------------------------------------------
              var(i,j,k) =  var(i,j,k) + var(ii,jj,kk)
              denom      = denom + 1.0_CGREAL
            ENDIF
          ENDDO
          ENDDO
          ENDDO

          var(i,j,k) =  var(i,j,k)/denom

        ENDIF ! fresh_cell

      ENDDO
      ENDDO
      ENDDO

#     ifdef MPI
!       Communicating the outmost layer of the subdomain.
!       -------------------------------------------------
        CALL par_comm_outermost_int (fresh_cell,myTime)
        CALL par_comm_outermost_real(var       ,myTime)
#     endif

    ENDIF ! boundary_formulation

END SUBROUTINE update_pressure_freshcells
!---------------------------------------------------------------------



SUBROUTINE calc_residual(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
	USE Cutcell_arrays

#ifdef MPI
    use mpi
    use par_parameters
#endif

    IMPLICIT NONE

!   Parameters
!   ----------
    REAL(KIND=CGREAL), INTENT (IN)  ::var(1-Ngl:nxc+Ngl,1-Ngl:nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL), INTENT (IN)  ::r(nxc,nyc,nzc)
    REAL(KIND=CGREAL), INTENT (OUT) ::resm

    INTEGER,           INTENT (OUT) ::loc(3)

!   Local Variables
!   ---------------
    INTEGER           :: i, j, k, iG, jG, kG

    REAL(KIND=CGREAL) :: res, dpdn

    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL) :: bmy,bpy,bcy
    REAL(KIND=CGREAL) :: bmz,bpz,bcz

    REAL(KIND=CGREAL) :: kill4GCMxm(nxc_GLBL), kill4GCMxp(nxc_GLBL)
    REAL(KIND=CGREAL) :: kill4GCMym(nyc_GLBL), kill4GCMyp(nyc_GLBL)
    REAL(KIND=CGREAL) :: kill4GCMzm(nzc), kill4GCMzp(nzc)



#ifdef MPI
    integer           :: imax, ierr

    DOUBLE PRECISION :: res_LOCL(4)
    DOUBLE PRECISION :: res_GLBL(4,0:maxRank-1)
#endif



!   Initializing Neumann/Dirichlet BC flags (SAMK/RM)
!   -------------------------------------------------
    Kill4GCMxm(:       ) = 1.0_CGREAL - gcmFlag
    Kill4GCMxp(:       ) = 1.0_CGREAL - gcmFlag
    Kill4GCMxm(1       ) = 1.0_CGREAL
    Kill4GCMxp(nxc_GLBL) = 1.0_CGREAL

    Kill4GCMym(:       ) = 1.0_CGREAL - gcmFlag
    Kill4GCMyp(:       ) = 1.0_CGREAL - gcmFlag
    Kill4GCMym(1       ) = 1.0_CGREAL
    Kill4GCMyp(nyc_GLBL) = 1.0_CGREAL

    Kill4GCMzm(:  ) = 1.0_CGREAL - gcmFlag
    Kill4GCMzp(:  ) = 1.0_CGREAL - gcmFlag
    Kill4GCMzm(1  ) = 1.0_CGREAL
    Kill4GCMzp(nzc) = 1.0_CGREAL

    resm = 0.0_CGREAL

    DO k=1,nzc
      kG=k

      DO j=1,nyc
        jG=L2GJ(j)   

        DO i=1,nxc
          iG=L2GI(i)

          bmx =   dxcinv(iG)  *dxinv(iG)*(1.0_CGREAL - ium(i,j,k)*Kill4GCMxm(iG))*AREA_W(i,j,k)*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)
          bpx =   dxcinv(iG+1)*dxinv(iG)*(1.0_CGREAL - iup(i,j,k)*Kill4GCMxp(iG))*AREA_E(i,j,k)*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)
          bcx = - ( bmx + bpx )
  
          bmy =   dycinv(jG)  *dyinv(jG)*(1.0_CGREAL - jum(i,j,k)*Kill4GCMym(jG))*AREA_S(i,j,k)*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)
          bpy =   dycinv(jG+1)*dyinv(jG)*(1.0_CGREAL - jup(i,j,k)*Kill4GCMyp(jG))*AREA_N(i,j,k)*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)
          bcy = - ( bmy + bpy )

          bmz =   dzcinv(kG)  *dzinv(kG)*(1.0_CGREAL - kum(i,j,k)*Kill4GCMzm(kG))*KillFor2D*AREA_B(i,j,k)*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)
          bpz =   dzcinv(kG+1)*dzinv(kG)*(1.0_CGREAL - kup(i,j,k)*Kill4GCMzp(kG))*KillFor2D*AREA_F(i,j,k)*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC)
          bcz = - ( bmz + bpz )

          bmx = bmx*REAL(1-iblank(i,j,k),KIND=CGREAL)
          bpx = bpx*REAL(1-iblank(i,j,k),KIND=CGREAL)

          bmy = bmy*REAL(1-iblank(i,j,k),KIND=CGREAL)
          bpy = bpy*REAL(1-iblank(i,j,k),KIND=CGREAL)

          bmz = bmz*REAL(1-iblank(i,j,k),KIND=CGREAL)
          bpz = bpz*REAL(1-iblank(i,j,k),KIND=CGREAL)

          bc =   (bcx+bcy+bcz)*REAL(1-iblank(i,j,k),KIND=CGREAL)  &
                              +REAL(  iblank(i,j,k),KIND=CGREAL)      

		!!  CALL cutcell_dpdn(i,j,k,dpdn)
		!!  dpdn = dpdn*cut_area(i,j,k)*dxinv(iG)*dyinv(jG)*dzinv(kG)

          res     = r(i,j,k) - var(i,j,k)*bc                                    &
                             - var(i-1,j,k)*bmx*(1.0_CGREAL - ium(i,j,k) )      &
                             - var(i+1,j,k)*bpx*(1.0_CGREAL - iup(i,j,k) )      &
                             - var(i,j-1,k)*bmy*(1.0_CGREAL - jum(i,j,k) )      &
                             - var(i,j+1,k)*bpy*(1.0_CGREAL - jup(i,j,k) )      &
                             - var(i,j,k-1)*bmz*(1.0_CGREAL - kum(i,j,k) )      &
                             - var(i,j,k+1)*bpz*(1.0_CGREAL - kup(i,j,k) )      &
                             - pGhost(i-1,j,k)*bmx*ium(i,j,k)*gcmFlag           &
                             - pGhost(i+1,j,k)*bpx*iup(i,j,k)*gcmFlag           &
                             - pGhost(i,j-1,k)*bmy*jum(i,j,k)*gcmFlag           &
                             - pGhost(i,j+1,k)*bpy*jup(i,j,k)*gcmFlag           &
                             - pGhost(i,j,k-1)*bmz*kum(i,j,k)*gcmFlag           &
                             - pGhost(i,j,k+1)*bpz*kup(i,j,k)*gcmFlag           & 
                             + boundPresSource(i-1,j,k)                         &
                             + boundPresSource(i+1,j,k)                         &
                             + boundPresSource(i,j-1,k)                         &
                             + boundPresSource(i,j+1,k)                         &
                             + boundPresSource(i,j,k-1)                         &
                             + boundPresSource(i,j,k+1)                          

          res = res*REAL(1-iblank(i,j,k),KIND=CGREAL)

          IF (ABS(res) > resm ) THEN
            resm = ABS(res)
            loc(1) = iG
            loc(2) = jG
            loc(3) = kG
          ENDIF

        ENDDO
      ENDDO
    ENDDO

#ifdef MPI
    res_LOCL(1)=resm
    res_LOCL(2)=real(loc(1))+.01
    res_LOCL(3)=real(loc(2))+.01
    res_LOCL(4)=real(loc(3))+.01
    
    res_GLBL=0.


    call mpi_Allgather(res_LOCL,4,MPI_DOUBLE_PRECISION, &
                       res_GLBL,4,MPI_DOUBLE_PRECISION, &
                       mpi_comm_world, ierr)




    res=0.0_CGREAL
    imax = 0
    do i=0,maxRank-1
      if (res_GLBL(1,i)>res) then
        imax=i
        res=REAL(res_GLBL(1,i), KIND=CGREAL)
      end if
    end do
   
    resm  =    REAL(res_GLBL(1,imax), KIND=CGREAL)
    loc(1)=int(res_GLBL(2,imax))
    loc(2)=int(res_GLBL(3,imax))
    loc(3)=int(res_GLBL(4,imax))
#endif

END SUBROUTINE calc_residual
!-------------------------------------------------------------------------------  



SUBROUTINE itsolv(var,r,totalTime)

! ----------------------------------------
!  Line SOR with Gauss Siedel as smoother
! ----------------------------------------

    USE global_parameters
    USE flow_parameters
!!    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
!!    USE multiuse_arrays
    USE solver_arrays
!!    USE GCM_arrays
	USE cutcell_arrays
#ifdef MPI
    use mpi
#endif

    IMPLICIT NONE

!   parameters
!   ----------
    REAL(KIND=CGREAL), INTENT (INOUT) :: var(-(Ngl-1):nxc+Ngl,-(Ngl-1):nyc+Ngl,0:nzc+1)
    REAL(KIND=CGREAL), INTENT (IN)    :: r(nxc,nyc,nzc)
    real(kind=CGREAL), intent (out)   :: totalTime  

!   Local variables
!   ---------------
    INTEGER :: i , j , k
    INTEGER :: iG, jG, kG
!!    INTEGER :: iBody, iRow, iG,jG, iNode, jNode, n
!!    INTEGER :: loc(3)
!!    REAL(KIND=CGREAL) :: maxres 

    INTEGER :: clock1, clock2, clock_rate

    REAL(KIND=CGREAL) :: startTime, endTime
    
!!    REAL(KIND=CGREAL) :: VolCell,killForGCM
    REAL(KIND=CGREAL) :: killForGCM
    real(kind=CGREAL) :: KillItIm, KillItIp
    real(kind=CGREAL) :: KillItJm, KillItJp
    real(kind=CGREAL) :: KillItKm, KillItKp



#ifdef MPI
    startTime = MPI_WTIME()
#else
    call system_clock(clock1)
#endif

!   Note: homogeneous Neumann BC and Dirichlet BC have been hardcoded here.
!   -----------------------------------------------------------------------
    KillForGCM = 1.0d0-gcmFlag 
  
!!    CALL enforce_p_periodic(var)
    IF ( ndim == DIM_3D           .AND. &
         bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN
         do i=0,nxc+1                                    !Ehsan added for periodic BC
           do j=0,nyc+1                
             var(i,j,0) = var(i,j,nzc)     
             var(i,j,nzc+1) = var(i,j,1)       
           enddo       
         enddo           
    END IF ! bcz1

!   Line solve in the x-direction
!   -----------------------------
    DO k=1,nzc
      kG=k

      KillItKm=KillForGCM
      KillItKp=KillForGCM
      if (kG==1) then
        KillItKm=1.0d0
      elseif (kG==nzc) then
        KillItKp=1.0d0
      endif
            
      DO j=1,nyc
        jG=L2GJ(j)

        KillItJm=KillForGCM
        KillItJp=KillForGCM
        if (jG==1) then
          KillItJm=1.0d0
        elseif (jG==nyc_GLBL) then
          KillItJp=1.0d0
        endif

        DO i=1,nxc
          iG=L2GI(i)

          KillItIm=KillForGCM
          KillItIp=KillForGCM
          if (iG==1) then
            KillItIm=1.0d0
          elseif (iG==nxc_GLBL) then
            KillItIp=1.0d0
          endif

          amx(i) =   dxcinv(iG)  *dxinv(iG)*(1.0_CGREAL - ium(i,j,k)*KillItIm)*AREA_W(i,j,k)*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)
          apx(i) =   dxcinv(iG+1)*dxinv(iG)*(1.0_CGREAL - iup(i,j,k)*KillItIp)*AREA_E(i,j,k)*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)
          acx(i) = - ( amx(i) + apx(i) )
          
          amy(j) =   dycinv(jG)  *dyinv(jG)*(1.0_CGREAL - jum(i,j,k)*KillItJm)*AREA_S(i,j,k)*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)
          apy(j) =   dycinv(jG+1)*dyinv(jG)*(1.0_CGREAL - jup(i,j,k)*KillItJp)*AREA_N(i,j,k)*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)
          acy(j) = - ( amy(j) + apy(j) )

          amz(k) =   dzcinv(kG)  *dzinv(kG)*(1.0_CGREAL - kum(i,j,k)*KillItKm)*KillFor2D*AREA_B(i,j,k)*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)
          apz(k) =   dzcinv(kG+1)*dzinv(kG)*(1.0_CGREAL - kup(i,j,k)*KillItKp)*KillFor2D*AREA_F(i,j,k)*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC)
          acz(k) = - ( amz(k) + apz(k) )

          rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum(i,j,k))&
                            - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup(i,j,k))&
                            - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum(i,j,k))&
                            - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup(i,j,k))!!&
!!                          - pGhost(i-1,j,k)*amx(i)*ium(i,j,k)*gcmFlag  &
!!                          - pGhost(i+1,j,k)*apx(i)*iup(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j-1,k)*amy(j)*jum(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j+1,k)*apy(j)*jup(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j,k-1)*amz(k)*kum(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j,k+1)*apz(k)*kup(i,j,k)*gcmFlag  &
!!                          + boundPresSource(i-1,j,k) &
!!                          + boundPresSource(i+1,j,k) &
!!                          + boundPresSource(i,j-1,k) &
!!                          + boundPresSource(i,j+1,k) &
!!                          + boundPresSource(i,j,k-1) &
!!                          + boundPresSource(i,j,k+1)  

          amx(i) = amx(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL-ium(i,j,k))
          apx(i) = apx(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL-iup(i,j,k))

          acx(i) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
                                         +REAL(  iblank(i,j,k),KIND=CGREAL)                     

          rhs(i) = rhs(i)*REAL(1-iblank(i,j,k),KIND=CGREAL)

      ENDDO ! i 

!!      IF (bcx1 == BC_TYPE_PERIODIC .AND. bcx2 == BC_TYPE_PERIODIC) THEN  !!!!!!!!!probably broken
!!        rhs(1)    = rhs(1)    - var(nx-1,j,k)*amx(1)
!!        rhs(nx-1) = rhs(nx-1) - var(1,j,k)*apx(nx-1)
!!      END IF

      rhs(1)   = rhs(1)  -var(0,j,k)    *amx(1)  *(1.0_CGREAL - ium(1  ,j,k))*REAL(1-iblank(1  ,j,k),KIND=CGREAL)   ! Added by SAMK
      rhs(nxc) = rhs(nxc)-var(nxc+1,j,k)*apx(nxc)*(1.0_CGREAL - iup(nxc,j,k))*REAL(1-iblank(nxc,j,k),KIND=CGREAL)   ! Added by SAMK

      CALL tdma(amx,acx,apx,rhs(1:nxc),dummy(1:nxc),1,nxc)
      DO i=1,nxc
        var(i,j,k) = var(i,j,k) + omega*(dummy(i)-var(i,j,k))
      ENDDO

    ENDDO ! j 
    ENDDO ! k

!!    CALL enforce_p_periodic(var)

    IF ( ndim == DIM_3D           .AND. &
         bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN
         do i=0,nxc+1                                    !Ehsan added for periodic BC
           do j=0,nyc+1
             var(i,j,0) = var(i,j,nzc)
             var(i,j,nzc+1) = var(i,j,1)
           enddo
         enddo
    END IF ! bcz1

!   Line solve in the y-direction
!   -----------------------------
    DO k=1,nzc
      kG=k

      KillItKm=KillForGCM
      KillItKp=KillForGCM
      if (kG==1) then
        KillItKm=1.0d0
      elseif (kG==nzc) then
        KillItKp=1.0d0
      endif
            
      DO i=1,nxc
        iG=L2GI(i)

        KillItIm=KillForGCM
        KillItIp=KillForGCM
        if (iG==1) then
          KillItIm=1.0d0
        elseif (iG==nxc_GLBL) then
          KillItIp=1.0d0
        endif

        DO j=1,nyc
          jG=L2GJ(j)

          KillItJm=KillForGCM
          KillItJp=KillForGCM
          if (jG==1) then
            KillItJm=1.0d0
          elseif (jG==nyc_GLBL) then
            KillItJp=1.0d0
          endif

          amx(i) =   dxcinv(iG)  *dxinv(iG)*(1.0_CGREAL - ium(i,j,k)*KillItIm)*AREA_W(i,j,k)*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)
          apx(i) =   dxcinv(iG+1)*dxinv(iG)*(1.0_CGREAL - iup(i,j,k)*KillItIp)*AREA_E(i,j,k)*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)
          acx(i) = - ( amx(i) + apx(i) )
          
          amy(j) =   dycinv(jG)  *dyinv(jG)*(1.0_CGREAL - jum(i,j,k)*KillItJm)*AREA_S(i,j,k)*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)
          apy(j) =   dycinv(jG+1)*dyinv(jG)*(1.0_CGREAL - jup(i,j,k)*KillItJp)*AREA_N(i,j,k)*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)
          acy(j) = - ( amy(j) + apy(j) )

          amz(k) =   dzcinv(kG)  *dzinv(kG)*(1.0_CGREAL - kum(i,j,k)*KillItKm)*KillFor2D*AREA_B(i,j,k)*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)
          apz(k) =   dzcinv(kG+1)*dzinv(kG)*(1.0_CGREAL - kup(i,j,k)*KillItKp)*KillFor2D*AREA_F(i,j,k)*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC)
          acz(k) = - ( amz(k) + apz(k) )


          rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium(i,j,k) )  &  
                            - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup(i,j,k) )  & 
                            - var(i,j,k-1)*amz(k)*(1.0_CGREAL - kum(i,j,k) )  &
                            - var(i,j,k+1)*apz(k)*(1.0_CGREAL - kup(i,j,k) )  !!&
!!                          - pGhost(i-1,j,k)*amx(i)*ium(i,j,k)*gcmFlag  &
!!                          - pGhost(i+1,j,k)*apx(i)*iup(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j-1,k)*amy(j)*jum(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j+1,k)*apy(j)*jup(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j,k-1)*amz(k)*kum(i,j,k)*gcmFlag  &
!!                          - pGhost(i,j,k+1)*apz(k)*kup(i,j,k)*gcmFlag  &
!!                          + boundPresSource(i-1,j,k) &
!!                          + boundPresSource(i+1,j,k) &
!!                          + boundPresSource(i,j-1,k) &
!!                          + boundPresSource(i,j+1,k) &
!!                          + boundPresSource(i,j,k-1) &
!!                          + boundPresSource(i,j,k+1) 

!         Modify rhs and coefficients for Ghost cells in GCM
!         --------------------------------------------------
          amy(j) = amy(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL - jum(i,j,k))
          apy(j) = apy(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL - jup(i,j,k))
          acy(j) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
                                         +REAL(  iblank(i,j,k),KIND=CGREAL)       

          rhs(j) = rhs(j)*REAL(1-iblank(i,j,k),KIND=CGREAL)
        ENDDO

!!      IF (bcy1 == BC_TYPE_PERIODIC.AND.bcy2 == BC_TYPE_PERIODIC) THEN
!!        rhs(1)    = rhs(1) - var(i,ny-1,k)*amy(1)
!!        rhs(ny-1) = rhs(ny-1) - var(i,1,k)*apy(ny-1)
!!      END IF

        rhs(1)  =rhs(1)  -var(i,0,k)    *amy(1)  *(1.0_CGREAL - jum(i,1  ,k))*REAL(1-iblank(i,1  ,k),KIND=CGREAL)   ! Added by SAMK
        rhs(nyc)=rhs(nyc)-var(i,nyc+1,k)*apy(nyc)*(1.0_CGREAL - jup(i,nyc,k))*REAL(1-iblank(i,nyc,k),KIND=CGREAL)   ! Added by SAMK

        CALL tdma(amy,acy,apy,rhs(1:nyc),dummy(1:nyc),1,nyc)

        DO j=1,nyc
          var(i,j,k) = var(i,j,k) + omega*(dummy(j)-var(i,j,k))
        ENDDO

      ENDDO
    ENDDO

!   Line solver in the z-direction
!   ------------------------------
    IF (ndim == DIM_3D) THEN
!!      CALL enforce_p_periodic(var)

    IF (    bcz1 == BC_TYPE_PERIODIC .AND. &
         bcz2 == BC_TYPE_PERIODIC       ) THEN
         do i=0,nxc+1                                    !Ehsan added for periodic BC
           do j=0,nyc+1
             var(i,j,0) = var(i,j,nzc)
             var(i,j,nzc+1) = var(i,j,1)
           enddo
         enddo
    END IF ! bcz1

      DO j=1,nyc
        jG=L2GJ(j)

        KillItJm=KillForGCM
        KillItJp=KillForGCM
        if (jG==1) then
          KillItJm=1.0d0
        elseif (jG==nyc_GLBL) then
          KillItJp=1.0d0
        endif

        DO i=1,nxc
          iG=L2GI(i)

          KillItIm=KillForGCM
          KillItIp=KillForGCM
          if (iG==1) then
            KillItIm=1.0d0
          elseif (iG==nxc_GLBL) then
            KillItIp=1.0d0
          endif

          DO k=1,nzc
            kG=k

            KillItKm=KillForGCM
            KillItKp=KillForGCM
            if (kG==1) then
              KillItKm=1.0d0
            elseif (kG==nzc) then
              KillItKp=1.0d0
            endif
            
			amx(i) =   dxcinv(iG)  *dxinv(iG)*(1.0_CGREAL - ium(i,j,k)*KillItIm)*AREA_W(i,j,k)*(1.-ghostcellmark(i-1,j,k)*iSSMP*iCC)
			apx(i) =   dxcinv(iG+1)*dxinv(iG)*(1.0_CGREAL - iup(i,j,k)*KillItIp)*AREA_E(i,j,k)*(1.-ghostcellmark(i+1,j,k)*iSSMP*iCC)
			acx(i) = - ( amx(i) + apx(i) )
          
			amy(j) =   dycinv(jG)  *dyinv(jG)*(1.0_CGREAL - jum(i,j,k)*KillItJm)*AREA_S(i,j,k)*(1.-ghostcellmark(i,j-1,k)*iSSMP*iCC)
			apy(j) =   dycinv(jG+1)*dyinv(jG)*(1.0_CGREAL - jup(i,j,k)*KillItJp)*AREA_N(i,j,k)*(1.-ghostcellmark(i,j+1,k)*iSSMP*iCC)
			acy(j) = - ( amy(j) + apy(j) )

			amz(k) =   dzcinv(kG)  *dzinv(kG)*(1.0_CGREAL - kum(i,j,k)*KillItKm)*KillFor2D*AREA_B(i,j,k)*(1.-ghostcellmark(i,j,k-1)*iSSMP*iCC)
			apz(k) =   dzcinv(kG+1)*dzinv(kG)*(1.0_CGREAL - kup(i,j,k)*KillItKp)*KillFor2D*AREA_F(i,j,k)*(1.-ghostcellmark(i,j,k+1)*iSSMP*iCC)
			acz(k) = - ( amz(k) + apz(k) )

            rhs(k) = r(i,j,k) - var(i,j-1,k)*amy(j)*(1.0_CGREAL - jum(i,j,k) ) &
                              - var(i,j+1,k)*apy(j)*(1.0_CGREAL - jup(i,j,k) ) &
                              - var(i-1,j,k)*amx(i)*(1.0_CGREAL - ium(i,j,k) ) &
                              - var(i+1,j,k)*apx(i)*(1.0_CGREAL - iup(i,j,k) ) !!&
!!                            - pGhost(i-1,j,k)*amx(i)*ium(i,j,k)*gcmFlag  &
!!                            - pGhost(i+1,j,k)*apx(i)*iup(i,j,k)*gcmFlag  &
!!                            - pGhost(i,j-1,k)*amy(j)*jum(i,j,k)*gcmFlag  &
!!                            - pGhost(i,j+1,k)*apy(j)*jup(i,j,k)*gcmFlag  &
!!                            - pGhost(i,j,k-1)*amz(k)*kum(i,j,k)*gcmFlag  &
!!                            - pGhost(i,j,k+1)*apz(k)*kup(i,j,k)*gcmFlag  &
!!                            + boundPresSource(i-1,j,k) &
!!                            + boundPresSource(i+1,j,k) &
!!                            + boundPresSource(i,j-1,k) &
!!                            + boundPresSource(i,j+1,k) &
!!                            + boundPresSource(i,j,k-1) &
!!                            + boundPresSource(i,j,k+1) 

            amz(k) = amz(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL - kum(i,j,k))
            apz(k) = apz(k)*REAL(1-iblank(i,j,k),KIND=CGREAL)*(1.0_CGREAL - kup(i,j,k))
            acz(k) = (acx(i)+acy(j)+acz(k))*REAL(1-iblank(i,j,k),KIND=CGREAL) &
                                           +REAL(  iblank(i,j,k),KIND=CGREAL)

            rhs(k) = rhs(k)*(1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL) ) 
          ENDDO ! k
      
!!        IF (bcz1 == BC_TYPE_PERIODIC.AND.bcz2 == BC_TYPE_PERIODIC) THEN
!!          rhs(1)    = rhs(1)    - var(i,j,nz-1)*amz(1)
!!          rhs(nz-1) = rhs(nz-1) - var(i,j,1)*apz(nz-1)
!!        END IF
  
          rhs(1)  =rhs(1)  -var(i,j,0)    *amz(1)  *(1.0_CGREAL - kum(i,j,1)  )*REAL(1-iblank(i,j,1  ),KIND=CGREAL)   ! Added by SAMK
          rhs(nzc)=rhs(nzc)-var(i,j,nzc+1)*apz(nzc)*(1.0_CGREAL - kup(i,j,nzc))*REAL(1-iblank(i,j,nzc),KIND=CGREAL)   ! Added by SAMK

          IF (bcz1 == BC_TYPE_PERIODIC.AND.bcz2 == BC_TYPE_PERIODIC) THEN               !Ehsan added for Periodic BC
            rhs(1)    = rhs(1)    - var(i,j,nzc)*amz(1)                       
            rhs(nzc) = rhs(nzc) - var(i,j,1)*apz(nzc)                     
          END IF                                                   

          CALL tdma(amz,acz,apz,rhs(1:nzc),dummy(1:nzc),1,nzc)

          DO k=1,nzc
            var(i,j,k) = var(i,j,k) + omega*(dummy(k)-var(i,j,k))
          ENDDO ! k

        ENDDO ! i 
      ENDDO ! j
    ENDIF ! ndim

#ifdef MPI
    endTime = MPI_WTIME()
    totalTime = endTime-startTime
#else
    call system_clock(clock2, clock_rate)
    totalTime= REAL(clock2-clock1)/REAL(clock_rate)
#endif

END SUBROUTINE itsolv
!---------------------------------------------------------------------



SUBROUTINE subtract_mean_pressure()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE pressure_arrays
    USE grid_arrays
    USE boundary_arrays
	USE CUTCELL_arrays

    IMPLICIT NONE

!   Local variables
!   ---------------
    INTEGER :: i, j, k, iG, jG

    REAL(KIND=CGREAL) :: pTotal, volTotal, pAverage



    IF (pbcx1 /= PBC_DIRICHLET .AND. &
        pbcy1 /= PBC_DIRICHLET .AND. &
        pbcz1 /= PBC_DIRICHLET .AND. &
        pbcx2 /= PBC_DIRICHLET .AND. &
        pbcy2 /= PBC_DIRICHLET .AND. &
        pbcz2 /= PBC_DIRICHLET) THEN  

!     subtracting average pressure
!     ----------------------------
      pTotal   = 0.0_CGREAL
      volTotal = 0.0_CGREAL

      DO k = 1, nzc
      DO j = 1, nyc
        jG=L2GJ(j)

        DO i = 1, nxc
          iG=L2GI(i)

          pTotal   = pTotal &
                   + pPrime(i,j,k)*dx(iG)*dy(jG)*dz(k)*REAL(1-iblank(i,j,k)     ,KIND=CGREAL)*VOL_CELL(i,j,k)  &
                   + pGhost(i,j,k)*dx(iG)*dy(jG)*dz(k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)*VOL_CELL(i,j,k)

          volTotal = volTotal &
                   + dx(iG)*dy(jG)*dz(k)*REAL(1-iblank(i,j,k)     ,KIND=CGREAL)*VOL_CELL(i,j,k) &
                   + dx(iG)*dy(jG)*dz(k)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)*VOL_CELL(i,j,k)
        END DO
      END DO
      END DO
    
	  
#ifdef MPI
      CALL par_getSumReal(pTotal)
      CALL par_getSumReal(volTotal)
#endif
 
      pAverage = pTotal/volTotal
	  
	  IF(monitorON) THEN
	  write(*,*) 'VolTotal = ', volTotal
	  write(*,*) 'pAverage = ', pAverage
	  ENDIF
	  
      DO k = 1, nzc
      DO j = myJLL, myJUL
      DO i = myILL, myIUL
        pPrime(i,j,k) = (pPrime(i,j,k) - pAverage)*REAL(1-iblank(i,j,k)     ,KIND=CGREAL)
        pGhost(i,j,k) = (pGhost(i,j,k) - pAverage)*REAL(ghostCellMark(i,j,k),KIND=CGREAL)
      END DO
      END DO
      END DO

    END IF

END SUBROUTINE subtract_mean_pressure
