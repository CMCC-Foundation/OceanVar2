!======================================================================
!
! This file is part of Oceanvar.
!
!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version (GPL-3.0-or-later).
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!======================================================================
!-----------------------------------------------------------------------
!                                                                      !
!> Initialize the diffusion filter                                       
!!
!! It computes the diffusion coeffients besed on correlation radii.
!! It computes the main diagonal, upper diagonal and lower diagonal 
!! of the matrix to solve the diffusion equation using an implicit
!! scheme.
!!
!                                                                      !
! Version 1: Mario Adani and Francesco Carere 2023                     !
!-----------------------------------------------------------------------
SUBROUTINE ini_dflt

   USE grd_str
   USE dfl_str
   USE mpi_str
   USE drv_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

! Constants
   REAL(r8)             :: r_earth
! Indices
   INTEGER(i4)          :: i,j,k
! Local variables

! 0D
   REAL(r8)             :: mu, eps
   REAL(r8)             :: r_earth2
   REAL(r8)             :: g
   INTEGER(i4)          :: nx,ny
! 1D
   REAL(r8),POINTER     :: dk(:),b(:),c(:)
   REAL(r8),POINTER     :: a(:),e(:),f(:)
   REAL(r8),POINTER     :: la(:),lo(:)
! 2D
   REAL(r8),POINTER     :: gdx(:,:),gdy(:,:)
! 3D
   REAL(r8),POINTER     :: msk_eps(:,:,:)
   INTEGER(i4),POINTER  :: msk(:,:,:)

!---
!- check stability
   IF ( dfl%nt .LE. 2 ) THEN
      WRITE (drv%dia,*),'NT can must be greater THEN 2'
      STOP
   ENDIF
!---
!- square earth radius
   r_earth  = phy%re*1000._r8
   r_earth2 = r_earth**2
!---
!- degree to radians 
   ALLOCATE ( grd%lat_rad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae), &
              grd%lon_rad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   grd%lon_rad = grd%lon
   grd%lat_rad = grd%lat
   WHERE ( grd%lon_rad  .GT. 180 ) grd%lon_rad = grd%lon_rad - 360.0_r8
   grd%lon_rad=grd%lon_rad*phy%pi/180._r8
   grd%lat_rad=grd%lat_rad*phy%pi/180._r8
!---
!- compute lon,lat difference
   ALLOCATE ( grd%dlon(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae), &
              grd%dlat(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

   DO j = 1-grd%jas,grd%jm+grd%jae
      DO i = 1-grd%ias,grd%im+grd%iae-1
         grd%dlon(i,j) = ABS(grd%lon_rad(i,j)-grd%lon_rad(i+1,j  ))
      ENDDO
   ENDDO
   grd%dlon(grd%im+grd%iae,:) = grd%dlon(grd%im+grd%iae-1,:)

   DO j = 1-grd%jas,grd%jm+grd%jae-1
      DO i = 1-grd%ias,grd%im+grd%iae
         grd%dlat(i,j) = ABS(grd%lat_rad(i,j)-grd%lat_rad(i  ,j+1))
      ENDDO
   ENDDO
   grd%dlat(:,grd%jm+grd%jae) = grd%dlat(:,grd%jm+grd%jae-1)
   WHERE ( grd%dlat == 0 ) grd%dlat=MINVAL(grd%dlat,mask=(grd%dlat .NE. 0))

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%dlat)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%dlon)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%lat_rad)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%lon_rad)
   ENDIF

!---
!-Correlation radius 
   ALLOCATE (dfl%rx3d (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
             dfl%ry3d (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   IF ( dfl%rd_corr ) THEN
      CALL rdcrl
   ELSE
      dfl%rx3d(:,:,:) = dfl%rx ! meters
      dfl%ry3d(:,:,:) = dfl%ry ! meters
   ENDIF

!---
!- decrease correlation radius close to the coast
   IF ( dfl%USE_cst ) THEN
      DO k = 1, grd%km
         WHERE (grd%distc3d(:,:,k) .LT. dfl%cst_dst )
            dfl%rx3d(:,:,k)  = MAX( (dfl%rx3d(:,:,k)*grd%distc3d(:,:,k))/dfl%cst_dst,&
                                     r_earth*grd%dlon )
            dfl%ry3d(:,:,k)  = MAX( (dfl%ry3d(:,:,k)*grd%distc3d(:,:,k))/dfl%cst_dst,  &
                                     r_earth*grd%dlat*COS(grd%lat_rad) )
         END WHERE
      ENDDO
   ENDIF

!---
!- from correlation radis to diffusion coefficient 
   ALLOCATE ( dfl%kx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%ky(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   dfl%kx(:,:,:) = dfl%rx3d(:,:,:)**2/DBLE(2*(dfl%nt)) ! meter**2/s
   dfl%ky(:,:,:) = dfl%ry3d(:,:,:)**2/DBLE(2*(dfl%nt)) ! meter**2/s

!---
!- define boundary condition
   IF ( dfl%USE_bc ) THEN
      ALLOCATE (msk_eps(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
      IF (dfl%bc_type == 'DIRICHLET') mu = 1
      IF (dfl%bc_type == 'NEUMANN')   mu = 0
      eps = (2*mu)/(mu+2*(1-mu))
      msk_eps(:,:,:) = 1._r8
      WHERE ( grd%msk(:,:,:) .EQ. 0 ) msk_eps(:,:,:) = eps
      dfl%kx(:,:,:) = dfl%kx(:,:,:) * msk_eps(:,:,:)
      dfl%ky(:,:,:) = dfl%ky(:,:,:) * msk_eps(:,:,:)
      DEALLOCATE ( msk_eps )
   ENDIF

!---
!- read/compute  weights
   IF ( dfl%rd_wgh ) THEN
      CALL rdwgh
   ELSE
! analitical weights
      ALLOCATE ( gdx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae), &
                 gdy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae), &
                 msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
                 dfl%wgh(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km)  )
      gdx(:,:) = r_earth*COS(grd%lat_rad)*grd%dlat
      gdy(:,:) = r_earth*grd%dlon
      DO k = 1, grd%km
         dfl%wgh(:,:,k) = ( 2._r8 * DBLE(dfl%nt-1) * DSQRT(dfl%kx(:,:,k) * dfl%ky(:,:,k)) ) / &
            (gdx(:,:)*gdy(:,:) ) * grd%msk(:,:,k)
      ENDDO
      msk = 0
      WHERE ( grd%msk(:,:,:)     .EQ.  1._r8                .AND. &
             (grd%distc3d(:,:,:) .LT. 2._r8*dfl%rx3d(:,:,:) .OR.  &
              grd%distc3d(:,:,:) .LT. 2._r8*dfl%ry3d(:,:,:)) ) msk(:,:,:) = 1
      DO k = 1,grd%km
         DO ny = 1,grd%jm
            DO nx = 1,grd%im
               IF ( msk(nx,ny,k) .EQ. 1 ) THEN
                  g = EXP(-(2._r8*grd%distc3d(nx,ny,k))**2/(2*dfl%rx3d(nx,ny,k)*dfl%ry3d(nx,ny,k)))/(phy%pi*2._r8)
                  IF ( dfl%bc_type == 'NEUMANN' )  dfl%wgh(nx,ny,k) = dfl%wgh(nx,ny,k) / (1._r8+g)
                  IF ( dfl%bc_type == 'DIRICHLET' ) dfl%wgh(nx,ny,k) = dfl%wgh(nx,ny,k) /(1._r8-g)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE ( gdx,gdy,msk )
   ENDIF

!---
!- compute coefficients
   ALLOCATE ( dfl%Ax(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%Ay(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%Lx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%Ly(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%Ux(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km), &
              dfl%Uy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
!longitude
   ALLOCATE ( dk(1-grd%ias:grd%im+grd%iae),&
              a (1-grd%ias:grd%im+grd%iae),&
              b (1-grd%ias:grd%im+grd%iae),&
              c (1-grd%ias:grd%im+grd%iae),&
              e (1-grd%ias:grd%im+grd%iae),&
              f (1-grd%ias:grd%im+grd%iae),&
              la(1-grd%ias:grd%im+grd%iae),&
              lo(1-grd%ias:grd%im+grd%iae))

   nx = SIZE(dk)
   DO k = 1,grd%km
      DO j = 1-grd%jas,grd%jm+grd%jae
         la(:) = grd%lat_rad(:,j)
         lo(:) = grd%dlon(:,j)
         dk(:) = dfl%kx(:,j,k)/(r_earth2*COS(la)**2)

         CALL cmp_coefficients(nx,dk(:),lo(:),        &  !input
                               a(:),b(:),c(:))           !output

         CALL tridiagLU(1,nx,a(:),b(:),c(:),          &  !input
                          e(:),f(:))                     !output
         dfl%Ax(:,j,k) = a
         dfl%Ux(:,j,k) = e
         dfl%Lx(:,j,k) = f
      ENDDO
   ENDDO
   DEALLOCATE( dk, a, b, c,  e, f, lo, la )


!latitude
   ALLOCATE ( dk(1-grd%jas:grd%jm+grd%jae),&
              a (1-grd%jas:grd%jm+grd%jae),&
              b (1-grd%jas:grd%jm+grd%jae),&
              c (1-grd%jas:grd%jm+grd%jae),&
              e (1-grd%jas:grd%jm+grd%jae),&
              f (1-grd%jas:grd%jm+grd%jae),&
              la(1-grd%jas:grd%jm+grd%jae))
   ny = size(dk)
   DO k = 1,grd%km
      DO i = 1-grd%ias,grd%im+grd%iae
         dk(:) = dfl%ky(i,:,k)/r_earth2
         la(:) = grd%dlat(i,:)
         CALL cmp_coefficients(ny,dk,la(:),      &  !input
                               a(:),b(:),c(:))      !output

         CALL tridiagLU(0,ny,a(:),b(:),c(:),     &  !input
                        e(:),f(:))                  !output

         dfl%Ay(i,:,k) = a
         dfl%Uy(i,:,k) = e
         dfl%Ly(i,:,k) = f
      ENDDO
   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km , dfl%Ay)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km,  dfl%Uy)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km,  dfl%Ly)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km , dfl%Ax)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km,  dfl%Ux)
      CALL exo_mpi( 1_i4, 1_i4,                                         &
                   -1_i4, 1_i4, grd%im+1,                               &
                   -1_i4, 1_i4, grd%jm+1,                               &
                   1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,  grd%km,  dfl%Lx)
   ENDIF

   DEALLOCATE ( dk, a, b, c,  e, f, la )

END SUBROUTINE ini_dflt
!================================================================
SUBROUTINE cmp_coefficients(n,alpha,dx,a,b,c)

!-----------------------------------------------------------------------
!                                                                      !
!>  Algorithm to compute the implicit coefficient                                      
!!
!!
!!
!
! Version 1: Mario Adani  2023                                         !
!-----------------------------------------------------------------------


   USE set_knd

   IMPLICIT NONE

   INTEGER(i4),INTENT(IN) :: n
   REAL(r8),   INTENT(IN) :: alpha(n)
   REAL(r8),   INTENT(IN) :: dx(n)
   REAL(r8),   INTENT(OUT):: a(n),b(n),c(n)

   INTEGER(i4) :: i
   REAL(r8)    :: coefficient(n)

   coefficient(:) = alpha(:)/dx(:)**2
   DO i=1,n-1
      a(i) = -coefficient(i)
      c(i) = -coefficient(i+1)
      b(i) = 1+coefficient(i+1)+coefficient(i)
   ENDDO
   a(n) = a(n-1)
   b(n) = b(n-1)
   c(n) = c(n-1)

END SUBROUTINE cmp_coefficients
!------------------------------------------
SUBROUTINE tridiagLU(md,n,a,b,c,e,f)

!----------------------------------------------------------------------------------------
!> tridiagLU Obtain the LU factorization of a tridiagonal matrix
!!
!! Synopsis: [e,f] = tridiag(a,b,c)
!!
!! Input: a,b,c = vectors defining the tridiagonal matrix. a is the
!! subdiagonal, b is the main diagonal, and c is the superdiagonal
!!
!! Output: e,f = vectors defining the L and U factors of the tridiagonal matrix
!!
!
! Francesco Carere 2024:  Much more efficient parallelisations 
!                         of tridiagonal LU decomposition exist
!----------------------------------------------------------------------------------------

   USE set_knd
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4),INTENT(IN) :: md, n
   REAL(r8),INTENT(IN)    :: a(n),b(n),c(n)
   REAL(r8),INTENT(OUT)   :: e(n),f(n)
   INTEGER(i4)            :: i
   INTEGER                :: mpi_ijr, mpi_rcv, mpi_snd, mpi_ijrm, ierr, mpistat(mpi_status_size)

   IF ( md .EQ. 1 ) THEN !modus "i" loop
      mpi_ijr  = mpi%ir
      mpi_rcv  = mpi%lft
      mpi_snd  = mpi%rgh
      mpi_ijrm = mpi%irm
   ELSE
      mpi_ijr  = mpi%jr
      mpi_rcv  = mpi%bot
      mpi_snd  = mpi%top
      mpi_ijrm = mpi%jrm
   ENDIF

   e(:) = 0._r8
   f(:) = 0._r8
   IF ( mpi_ijr .EQ. 1 ) THEN
      e(1) = b(1)
      f(1) = c(1)/b(1)
   ELSE
      CALL MPI_RECV(f(1),1,mpi%r8,mpi_rcv,1,mpi%comm,mpistat,ierr)
   ENDIF

   DO i = 2,n
      e(i) = b(i) - a(i)*f(i-1)
      f(i) = c(i) / e(i)
   ENDDO
   IF ( mpi%nproc .GT. 1 ) THEN
      IF ( mpi_ijr .NE. mpi_ijrm ) CALL MPI_SEND(f(n-1),1,mpi%r8,mpi_snd,1,mpi%comm,ierr)
   ENDIF


END SUBROUTINE tridiagLU
!================================================================
SUBROUTINE tridiagLUSolve(md,n,d,a,e,f,v)

!----------------------------------------------------------------------------------------
!> tridiagLUsolve Solve (LU)*v = d where L and U are LU factors of a tridiagonal matric
!!
!! Synopsis: v = tridiagLUsolve(d,e,f)
!! v = tridiagLUsolve(d,e,f,v)
!!
!! Input: d = right hand side vector of the system of equatoins
!! e,f = vectors defining the L and U factors of the tridiagonal matrix.
!! e and f are obtained with the tridiagLU FUNCTION
!! v = solution vector. If v is supplied, the elements of v are over-
!! written (thereby saving the memory allocation step). If v is not
!! supplied, it is created. v is USEd as a scratch vector in the
!! forward solve.
!!
!! Output: v = solution vector
!
! Francesco Carere 2024:  Much more efficient parallelisations 
!                         of tridiagonal LU decomposition exist
!----------------------------------------------------------------------------------------

   USE set_knd
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4),INTENT(IN)    :: md, n
   REAL(r8),   INTENT(IN)    :: d(n),e(n),f(n),a(n)
   REAL(r8),   INTENT(INOUT) :: v(n)
   INTEGER(i4)               :: i
   INTEGER                   :: mpi_ijr, mpi_bck, mpi_fwd, mpi_ijrm, ierr, mpistat(mpi_status_size)

   IF (md.EQ.1) THEN !modus "i" loop
      mpi_ijr  = mpi%ir
      mpi_bck  = mpi%lft
      mpi_fwd  = mpi%rgh
      mpi_ijrm = mpi%irm
   ELSE
      mpi_ijr  = mpi%jr
      mpi_bck  = mpi%bot
      mpi_fwd  = mpi%top
      mpi_ijrm = mpi%jrm
   ENDIF

! --- Forward substitution to solve L*w = d
   IF ( mpi_ijr .EQ. 1 ) THEN
      v(1) = d(1) / e(1)
   ELSE
      CALL MPI_RECV(v(1),1,mpi%r8,mpi_bck,1,mpi%comm,mpistat,ierr)
   ENDIF
   DO i = 2,n
      v(i) = ( d(i) - a(i)*v(i-1) ) / e(i)
   ENDDO
   IF ( mpi_ijr .NE. mpi_ijrm ) THEN
      CALL MPI_SEND(v(n-1),1,mpi%r8,mpi_fwd,1,mpi%comm,ierr)
      CALL MPI_RECV(v(n)  ,1,mpi%r8,mpi_fwd,1,mpi%comm,mpistat,ierr )
   ENDIF
! --- Backward substitution to solve U*v = w
   DO i = n-1,1,-1
      v(i) = v(i) - f(i) * v(i+1)
   ENDDO
   IF ( mpi_ijr .NE. 1 ) CALL MPI_SEND(v(2),1,mpi%r8,mpi_bck,1,mpi%comm,ierr)
END SUBROUTINE tridiagLUsolve
!================================================================
SUBROUTINE tridiagLUsolve_adj(md,n, db, a, e, f, vb)

!----------------------------------------------------------------------------------------
!> tridiagLUsolve Solve adjoint
!!
!!
!!
! Francesco Carere 2024:  Much more efficient parallelisations 
!                         of tridiagonal LU decomposition exist
!----------------------------------------------------------------------------------------

   USE set_knd
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4), INTENT(IN) :: md,n
   REAL(r8), INTENT(IN)    :: e(n), f(n), a(n)
   REAL(r8)                :: db(n)
   REAL(r8), INTENT(INOUT) :: vb(n)
   INTEGER(i4)             :: i
   REAL(r8)                :: tempb
   INTEGER                 :: mpi_ijr, mpi_bck, mpi_fwd, mpi_ijrm, ierr, mpistat(mpi_status_size)

   IF ( md .EQ. 1 ) THEN !modus "i" loop
      mpi_ijr  = mpi%ir
      mpi_bck  =mpi%lft
      mpi_fwd  = mpi%rgh
      mpi_ijrm = mpi%irm
   ELSE
      mpi_ijr  = mpi%jr
      mpi_bck  = mpi%bot
      mpi_fwd  = mpi%top
      mpi_ijrm = mpi%jrm
   ENDIF

   IF ( mpi_ijr .NE. 1 ) CALL MPI_RECV(vb(1),1,mpi%r8,mpi_bck,1,mpi%comm,mpistat,ierr)
   DO i = 1,n-1,1
      vb(i+1) = vb(i+1) - f(i)*vb(i)
   ENDDO
   IF ( mpi_ijr .NE. mpi_ijrm ) THEN
      CALL MPI_SEND(vb(n-1),1,mpi%r8,mpi_fwd,1,mpi%comm,ierr)
      CALL MPI_RECV(vb(n-1)  ,1,mpi%r8,mpi_fwd,1,mpi%comm,mpistat,ierr )
   ELSE
      tempb   = vb(n)/e(n)
      vb(n)   = 0.0_r8
      db(n)   = db(n) + tempb
      vb(n-1) = vb(n-1) - a(n)*tempb
   ENDIF
   DO i  =  n-1,2,-1
      tempb   = vb(i)/e(i)
      vb(i)   = 0.0_r8
      db(i)   = db(i) + tempb
      vb(i-1) = vb(i-1) - a(i)*tempb
   ENDDO
   IF (  mpi_ijr  .EQ.  1  ) THEN
      db(1) = db(1) + vb(1)/e(1)
      vb(1) = 0.0_r8
   ELSE
      CALL MPI_SEND(vb(1),1,mpi%r8,mpi_bck,1,mpi%comm,ierr)
   ENDIF

END SUBROUTINE tridiagLUsolve_adj
