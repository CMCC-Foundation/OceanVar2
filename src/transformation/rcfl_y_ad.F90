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
!> Recursive filter in y direction - adjoint                            
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE rcfl_y_ad( im, jm, km, fld, alp, bta, gam, del, mat_bc_y)

   USE set_knd
   USE cns_str
   USE mpi_str
   USE grd_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)               :: im, jm, km, jas, jae
   REAL(r8)                  :: fld(im,jm,km)
   REAL(r8)                  :: alp(jm,im), bta(jm,im)
   REAL(r8)                  :: gam(jm,im), del(jm,im)
   REAL(r8)                  :: mat_bc_y(9,jm,im)
   REAL(r8)                  :: coef
   INTEGER(i4)               :: i,j,k, ktr, llr, lr, js, je, is, ik,ii
   INTEGER                   :: irecvi, isENDi, ierr, npnt
   INTEGER(i4),  ALLOCATABLE :: jp(:)
   REAL(r8),     ALLOCATABLE :: a(:,:), b(:,:), c(:,:),c_r(:,:), msk(:,:)
   REAL(r8),     ALLOCATABLE :: bffr(:), bffs(:)
   INTEGER,      ALLOCATABLE :: istatus(:)
   REAL(r8),     ALLOCATABLE :: v(:),u(:),y(:)

   npnt=3
   jas=3
   jae=3
   is=0

   ALLOCATE ( a(jm,im*km), b(-2:jm,im*km), c(1:jm+jae,im*km), c_r(3,im*km), msk(1-jas:jm+jae,im*km))
   ALLOCATE ( bffr(im*km*npnt), bffs(im*km*npnt) )
   ALLOCATE ( istatus(mpi_status_size) )
   ALLOCATE ( jp(im*km))
   ALLOCATE ( v(3), u(3), y(3) )

   ik = MIN(2_i4,km)

   a(:,:)   = 0.0_r8
   b(:,:)   = 0.0_r8
   c(:,:)   = 0.0_r8
   c_r(:,:) = 0.0_r8
   msk(:,:) = 0.0_r8
   v(:)     = 0.0_r8
   u(:)     = 0.0_r8
   is       = 0

   DO k = 1,km
      DO j = 1,im
         jp((k-1)*im+j) = j
      ENDDO
   ENDDO

   DO k = 1,km
      DO j = 1,im
         IF ( grd%tmy((k-1)*grd%im+j) .EQ. 1 .OR. grd%tmy((k-1)*grd%im+j) .EQ. 2 ) THEN
            DO i = 1,jm
               c(i,(k-1)*im+j) = fld(j,i,k)
            ENDDO
         ENDIF
         IF (grd%tmy((k-1)*grd%im+j) .EQ. 2) THEN
            DO i = 1-jas,jm+jae
               msk(i,(k-1)*im+j) = grd%msr(j,i,k)
            ENDDO
         ENDIF
      ENDDO
   ENDDO

   ! positive direction
   DO lr = 1,mpi%thi(ik)

      js = grd%irs(lr,ik)
      je = grd%ire(lr,ik)
! mpi receive
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .GT. 1 ) THEN
         CALL MPI_IRECV( bffr,npnt*grd%imr(lr,ik), mpi%r8, mpi%bot, 1, mpi%comm, irecvi, ierr)
         CALL MPI_WAIT( irecvi, istatus, ierr)
         DO j = js,je
            c_r(1,j) =  bffr(j-js+1)
            c_r(2,j) =  bffr(j-js+1+grd%imr(lr,ik))
            c_r(3,j) =  bffr(j-js+1+2*grd%imr(lr,ik))
         ENDDO
      ENDIF

      DO j = js,je
         IF ( grd%tmy(j).EQ.1) THEN

            c(1:3,j) = c(1:3,j) + c_r(1:3,j)
            DO i=1,jm
               c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
               c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
               c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
               b(i, j)  = b(i  ,j) + del(i,jp(j))*c(i,j)
            END DO

         ELSEIF ( grd%tmy(j).EQ.2) THEN

            v(:) = 0.0_r8
            IF ( msk(1,j) .EQ. 1 .AND. msk(2,j) .EQ. 1.0 .AND. msk(3,j) .EQ. 0.0_r8 ) THEN
               c(1,j) = c(1,j) + c_r(1,j)
               c(2,j) = c(2,j) + c_r(2,j)
               v(2)   = c_r(3,j)
            ELSEIF ( msk(1,j) .EQ. 1.0 .AND. msk(2,j) .EQ. 0.0_r8 ) THEN
               c(1,j) = c(1,j) + c_r(1,j)
               v(2)   = c_r(2,j)
               v(3)   = c_r(3,j)
            ELSE
               c(1:3,j) = c(1:3,j) + msk(1:3,j)*c_r(1:3,j)
            ENDIF

            DO i=1,jm
               IF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 0.0_r8 ) THEN
                  v(1)       = v(1) + c(i,j)
                  u(:)       = 0.0_r8
                  u(1)       =  mat_bc_y(1,i,jp(j))*v(1) + mat_bc_y(4,i,jp(j))*v(2) + mat_bc_y(7,i,jp(j))*v(3)
                  u(2)       =  mat_bc_y(2,i,jp(j))*v(1) + mat_bc_y(5,i,jp(j))*v(2) + mat_bc_y(8,i,jp(j))*v(3)
                  u(3)       =  mat_bc_y(3,i,jp(j))*v(1) + mat_bc_y(6,i,jp(j))*v(2) + mat_bc_y(9,i,jp(j))*v(3)
                  b(i,    j) = b(i,    j) + del(i,jp(j))*u(1)
                  b(i-1,  j) = b(i-1,  j) + del(i,jp(j))*u(2)
                  b(i-2,  j) = b(i-2,  j) + del(i,jp(j))*u(3)
                  v(:)=0.0_r8
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 0.0_r8 ) THEN
                  c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
                  v(2)     = v(2)     + bta(i,jp(j))*c(i,j)
                  v(3)     = v(3)     + gam(i,jp(j))*c(i,j)
                  b(i, j)  = b(i  ,j) + del(i,jp(j))*c(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 1.0 .AND. msk(i+3,j) .EQ. 0.0_r8 ) THEN
                  c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
                  c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
                  v(2)     = v(2)     + gam(i,jp(j))*c(i,j)
                  b(i, j)  = b(i,j)   + del(i,jp(j))*c(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 1.0 .AND. msk(i+3,j) .EQ. 1.0 ) THEN
                  c(i+1,j) = c(i+1,j) + alp(i,jp(j))*c(i,j)
                  c(i+2,j) = c(i+2,j) + bta(i,jp(j))*c(i,j)
                  c(i+3,j) = c(i+3,j) + gam(i,jp(j))*c(i,j)
                  b(i, j)  =  b(i ,j) + del(i,jp(j))*c(i,j)
               ENDIF
            END DO

            IF ( msk(jm+1,j) .EQ. 1.0 .AND. msk(jm+2,j) .EQ. 1.0 .AND. msk(jm+3,j) .EQ. 0.0_r8 ) THEN
               c(jm+3,j) = v(2)
            ELSEIF ( msk(jm+1,j) .EQ. 1.0 .AND. msk(jm+2,j) .EQ. 0.0_r8 ) THEN
               c(jm+2,j) = v(2)
               c(jm+3,j) = v(3)
            ENDIF

         ENDIF
      ENDDO

! mpi send
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .LT. mpi%jrm ) THEN
         DO j = js,je
            bffs(j-js+1)                   =  c(jm+1,j)
            bffs(j-js+1+grd%imr(lr,ik))    =  c(jm+2,j)
            bffs(j-js+1+2* grd%imr(lr,ik)) =  c(jm+3,j)
         END DO
         CALL mpi_ISEND( bffs, npnt*grd%imr(lr,ik), mpi%r8, mpi%top, 1, mpi%comm, isENDi, ierr)
         CALL MPI_WAIT( isENDi, istatus, ierr)
      ENDIF

   ENDDO

   ! negative direction
   DO lr = 1,mpi%thi(ik)

      js = grd%irs(lr,ik)
      je = grd%ire(lr,ik)
! mpi receive
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .LT. mpi%jrm ) THEN
         CALL MPI_IRECV( bffr,npnt*grd%imr(lr,ik), mpi%r8, mpi%top, 1, mpi%comm, irecvi, ierr)
         CALL MPI_WAIT( irecvi, istatus, ierr)
         DO j = js,je
            b(jm,j)   =  b(jm,j)   + bffr(j-js+1)
            b(jm-1,j) =  b(jm-1,j) + bffr(j-js+1+grd%imr(lr,ik))
            b(jm-2,j) =  b(jm-2,j) + bffr(j-js+1+2*grd%imr(lr,ik))
         ENDDO
      ENDIF

      DO j = js,je
         IF ( grd%tmy(j) .EQ. 1 ) THEN
            DO i=jm,1+is,-1
               b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
               b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
               b(i-3,j) = b(i-3,j) + gam(i,jp(j))*b(i,j)
               a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
            ENDDO
         ELSEIF ( grd%tmy(j) .EQ. 2 ) THEN
            DO i = jm,1+is,-1
               IF ( msk(i,j) .EQ. 1.0 .AND. msk(i-1,j) .EQ. 0.0_r8 ) THEN
                  a(i,j)  = a(i,j) + del(i,jp(j))*b(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i-1,j) .EQ. 1.0 .AND. msk(i-2,j) .EQ. 0) THEN
                  b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
                  a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i-1,j) .EQ. 1.0 .AND. msk(i-2,j) .EQ. 1.0 .AND. msk(i-3,j) .EQ. 0 ) THEN
                  b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
                  b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
                  a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i-1,j) .EQ. 1.0 .AND. msk(i-2,j) .EQ. 1.0 .AND. msk(i-3,j) .EQ. 1 ) THEN
                  b(i-1,j) = b(i-1,j) + alp(i,jp(j))*b(i,j)
                  b(i-2,j) = b(i-2,j) + bta(i,jp(j))*b(i,j)
                  b(i-3,j) = b(i-3,j) + gam(i,jp(j))*b(i,j)
                  a(i  ,j) = a(i  ,j) + del(i,jp(j))*b(i,j)
               ENDIF
            ENDDO
         ENDIF
      END DO

! mpi send
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .GT. 1 ) THEN
         DO j = js,je
            bffs(j-js+1)                  = b(0,j)
            bffs(j-js+1+grd%imr(lr,ik))   = b(-1,j)
            bffs(j-js+1+2*grd%imr(lr,ik)) = b(-2,j)
         ENDDO
         CALL mpi_ISEND( bffs, npnt*grd%imr(lr,ik), mpi%r8, mpi%bot, 1, mpi%comm, isENDi, ierr)
         CALL MPI_WAIT( isENDi, istatus, ierr)
      ENDIF

   ENDDO

   DO k = 1,km
      DO j = 1,im
         IF ( grd%tmy((k-1)*grd%im+j) .EQ. 1 .OR. grd%tmy((k-1)*grd%im+j) .EQ. 2 ) THEN
            DO i = 1,jm
               fld(j,i,k) = a(i,(k-1)*im+j)
            ENDDO
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE ( jp )
   DEALLOCATE ( istatus )
   DEALLOCATE ( bffr, bffs )
   DEALLOCATE ( v, u, y )
   DEALLOCATE ( a, b, c, c_r, msk )

END SUBROUTINE rcfl_y_ad
