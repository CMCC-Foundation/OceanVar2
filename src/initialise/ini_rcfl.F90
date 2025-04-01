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
! Define filter constants                                              !
!                                                                      !
! Version 1: S.Dobricic 2006
! Version 2: S.Dobricic and R.Farina 2013                              !
!-----------------------------------------------------------------------
SUBROUTINE ini_rcfl

   USE set_knd
   USE drv_str
   USE grd_str
   USE eof_str
   USE cns_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   REAL(r8)                    :: mpimx, mpimn
   INTEGER(i4)                 :: k, nspl, i, j, kk, ik, jk, iter,my_i
   INTEGER(i4)                 :: iw, kst, ken, kln, kln0, ierr
   REAL(r8)                    :: E, dst
   REAL(r8)    , ALLOCATABLE   :: sfct(:), al(:), bt(:), sc(:)
   REAL(r8)    , ALLOCATABLE   :: sc2(:,:),sc2a(:,:)
   INTEGER(i4) , ALLOCATABLE   :: jnxx(:)
   REAL(r8)                    :: rone, L(0:2)
   REAL(r8)                    :: sigma,q,b0,b1,b2,b3,b
   REAL(r8)                    :: mat_bc_vect(9)
   INTEGER                     :: irecvi, isENDi
   INTEGER,        ALLOCATABLE :: istatus(:)

   ALLOCATE ( istatus(mpi_status_size) )

   rone = 1.0

   IF ( rcf%loc .GT. rcf%L ) THEN
      rcf%L = rcf%L * rcf%loc / sqrt( rcf%loc**2 - rcf%L**2 )
   ENDIF

   L(0) = rcf%L
   L(1) = rcf%L*0.9
   L(2) = rcf%L*0.1

   mat_bc_vect(:) = 0.0_r8

! ---
! Recursive filter constants
!
! 1) Create table

   DO iter=1,2

      rcf%L = L(iter)

      nspl = MAX(grd%jmg,grd%img)
      ALLOCATE ( sfct(nspl), jnxx(nspl), al(nspl), bt(nspl) )

      sfct(:) = 0.0_r8
      jnxx(:) = 0.0_r8
      al(:)   = 0.0_r8
      bt(:)   = 0.0_r8

      rcf%ntb = MIN(20,MIN(grd%jmg,grd%img))

      ALLOCATE ( rcf%al(rcf%ntb) )
      ALLOCATE ( rcf%sc(rcf%ntb) )

      rcf%dsmn =  1.e20_r8
      rcf%dsmx = -1.e20_r8
      DO j = 1,grd%jm
         DO i = 1,grd%im
            rcf%dsmn = MIN(rcf%dsmn,MIN(grd%dx(i,j),grd%dy(i,j)))
            rcf%dsmx = MAX(rcf%dsmx,MAX(grd%dx(i,j),grd%dy(i,j)))
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL MPI_REDUCE( rcf%dsmx, mpimx, 1, mpi%r8, MPI_MAX, 0, mpi%comm, ierr)
         IF ( mpi%myrank==0 ) THEN
            rcf%dsmx = mpimx
         ENDIF
         CALL MPI_BCAST( rcf%dsmx, 1, mpi%r8, 0, mpi%comm, ierr)
         CALL MPI_REDUCE( rcf%dsmn, mpimn, 1, mpi%r8, mpi_MIN, 0, mpi%comm, ierr)
         IF ( mpi%myrank==0 ) THEN
            rcf%dsmn = mpimn
         ENDIF
         CALL MPI_BCAST( rcf%dsmn, 1, mpi%r8, 0, mpi%comm, ierr)
      ENDIF

      rcf%dsmx = rcf%dsmx + MAX(rone,(rcf%dsmx-rcf%dsmn)/(rcf%ntb-2.))

      rcf%dsl = (rcf%dsmx-rcf%dsmn) / (rcf%ntb-1.)

      iw   = rcf%ntb/mpi%nproc + 1
      kst  = MIN( 1, rcf%ntb + 1)
      ken  = MIN( kst + iw - 1, rcf%ntb)
      kln0 = (ken-kst+1)
      kst  = MIN( mpi%myrank * iw + 1, rcf%ntb + 1)
      ken  = MIN( kst + iw - 1, rcf%ntb)
      kln  = ken - kst + 1
      ALLOCATE ( sc(kln))
      ALLOCATE ( sc2(kln0,mpi%nproc))
      ALLOCATE ( sc2a(kln0,mpi%nproc))
      kk = 0

      DO k = kst,ken
         dst = rcf%dsmn + (k-1.) * rcf%dsl
         sfct(:) = 0.
         sfct(nspl/2+1) = 1.
         sigma=rcf%L/dst
         CALL def_coef(sigma,b1,b2,b3,b)
         CALL rcfl_2(nspl,sfct,b1,b2,b3,b)
         CALL rcfl_2_ad(nspl,sfct,b1,b2,b3,b)
         kk = kk + 1
         sc(kk) = sfct(nspl/2+1)
      ENDDO

      DO k=1,mpi%nproc
         sc2(1:kln,k) = sc(1:kln)
      ENDDO
      CALL MPI_ALLTOALL( sc2, kln0, mpi%r8, sc2a, kln0, mpi%r8, mpi%comm, ierr)

      DO k = 1,mpi%nproc
         kst = MIN( (k-1) * iw + 1, rcf%ntb + 1)
         ken = MIN( kst + iw - 1, rcf%ntb)
         kln = ken - kst + 1
         rcf%sc(kst:ken) = sc2a(1:kln,k)
      ENDDO

      DEALLOCATE ( sc, sc2, sc2a )
      DEALLOCATE ( sfct, jnxx, al, bt )

      DO j = 1,grd%jm
         DO i = 1,grd%im
            IF ( grd%msk(i,j,1) .EQ. 1 ) THEN
               dst = ( grd%dx(i,j) - rcf%dsmn )/rcf%dsl
               k   = int(dst) + 1
               dst = dst - REAL(k-1)
               grd%scx(i,j,iter) = sqrt( 1./ (rcf%sc(k)*(1.-dst) + rcf%sc(k+1)*dst) )
               dst = ( grd%dy(i,j) - rcf%dsmn )/rcf%dsl
               k   = int(dst) + 1
               dst = dst - REAL(k-1)
               grd%scy(i,j,iter) = sqrt( 1./ (rcf%sc(k)*(1.-dst) + rcf%sc(k+1)*dst) )
            ELSE
               grd%scx(i,j,iter) = 0._r8
               grd%scy(i,j,iter) = 0._r8

            ENDIF
         ENDDO
      ENDDO

      DO j=1,grd%jm
         DO i=1,grd%im
            IF ( grd%msk(i,j,1) .EQ. 1 ) THEN
               dst =grd%dx(i,j)
               sigma= rcf%L/dst
               CALL  def_coef(sigma,grd%alx(i,j,iter),grd%btx(i,j,iter),grd%gmx(i,j,iter),grd%dlx(i,j,iter))
               CALL  def_coef_bc(mat_bc_vect,grd%alx(i,j,iter), grd%btx(i,j,iter), grd%gmx(i,j,iter))
               grd%mat_bc_x(:,i,j,iter)=mat_bc_vect(:)
            ELSE
               grd%mat_bc_x(:,i,j,iter)=0._r8
            ENDIF
         END DO
      END DO




      DO j=1,grd%im
         DO i=1,grd%jm
            IF ( grd%msk(j,i,1) .EQ. 1 ) THEN
               dst =grd%dy(j,i)
               sigma= rcf%L/dst
               CALL  def_coef(sigma,grd%aly(i,j,iter),grd%bty(i,j,iter),grd%gmy(i,j,iter),grd%dly(i,j,iter))
               CALL  def_coef_bc(mat_bc_vect,grd%aly(i,j,iter), grd%bty(i,j,iter), grd%gmy(i,j,iter))
               grd%mat_bc_y(:,i,j,iter)=mat_bc_vect(:)
            ELSE
               grd%mat_bc_y(:,i,j,iter)=0._r8
            ENDIF
         END DO
      END DO


   ENDDO

! 2) Mask Selection in x direction
   grd%tmx(:) = 0
   DO k = 1,grd%km
      DO j = 1,grd%jm
         DO i = 1-3,grd%im+3
            grd%tmx((k-1)*grd%jm+j) = grd%tmx((k-1)*grd%jm+j) + int(grd%msr(i,j,k))
         END DO
         IF ( grd%tmx((k-1)*grd%jm+j) .EQ. (grd%im+6) ) THEN
            grd%tmx((k-1)*grd%jm+j)=1   !!!! all Ocean
         ELSEIF ( grd%tmx((k-1)*grd%jm+j) .GT. 0 .AND. grd%tmx((k-1)*grd%jm+j).LT.(grd%im+6) ) THEN
            grd%tmx((k-1)*grd%jm+j)=2   !!!Ocean and Costs
         ENDIF
      END DO
   END DO

! 3) Mask Selection in y direction
   grd%tmy(:) = 0
   DO k = 1,grd%km
      DO j = 1,grd%im
         DO i = 1-3,grd%jm+3
            grd%tmy((k-1)*grd%im+j) = grd%tmy((k-1)*grd%im+j) + int(grd%msr(j,i,k))
         END DO
         IF ( grd%tmy((k-1)*grd%im+j) .EQ. (grd%jm+6) ) THEN
            grd%tmy((k-1)*grd%im+j)=1   !!!! all Ocean
         ELSEIF ( grd%tmy((k-1)*grd%im+j) .GT. 0 .AND. grd%tmy((k-1)*grd%im+j) .LT. (grd%jm+6) ) THEN
            grd%tmy((k-1)*grd%im+j)=2   !!!Ocean and Costs
         ENDIF
      END DO
   END DO

! 4) Apply scaling
   grd%scx(:,:,1) = grd%scx(:,:,1)*0.9
   grd%scx(:,:,2) = 1.0
   grd%scy(:,:,1) = grd%scy(:,:,1)*0.9
   grd%scy(:,:,2) = 1.0

   rcf%L = L(0)

   DEALLOCATE ( istatus )

END SUBROUTINE ini_rcfl
