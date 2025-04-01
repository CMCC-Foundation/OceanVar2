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
!> Recursive filter in x direction                                      !
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE rcfl_x( im, jm, km, fld, alp, bta, gam, del, mat_bc_x)

   USE set_knd
   USE cns_str
   USE mpi_str
   USE grd_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)               :: im, jm, km, ias, iae, ias_new, iae_new
   REAL(r8)                  :: fld(im,jm,km)
   REAL(r8)                  :: alp(im,jm), bta(im,jm)
   REAL(r8)                  :: gam(im,jm), del(im,jm)
   REAL(r8)                  :: mat_bc_x(9,im,jm)
   REAL(r8)                  :: coef
   INTEGER(i4)               :: i,j,k, ktr, llr, lr, js, je, is, ik,ii
   INTEGER                   :: irecvi, isENDi, ierr, npnt
   REAL(r8),     ALLOCATABLE :: v(:), u(:), y(:)
   INTEGER(i4),  ALLOCATABLE :: jp(:)
   REAL(r8),     ALLOCATABLE :: a(:,:), b(:,:), c(:,:), c_s(:,:) , msk(:,:)
   REAL(r8),     ALLOCATABLE :: bffr(:), bffs(:)
   INTEGER,      ALLOCATABLE :: istatus(:)

   npnt = 3
   ias  = 3
   iae  = 3

   ALLOCATE ( a(im,jm*km), b(1-npnt:im,jm*km), c(im+npnt,jm*km), c_s(3,jm*km), msk(1-ias:im+iae,jm*km) )
   ALLOCATE ( v(3), u(3), y(3) )
   aLLOCATE ( bffr(jm*km*npnt), bffs(jm*km*npnt) )
   ALLOCATE ( istatus(mpi_status_size) )
   ALLOCATE ( jp(jm*km))

   ik = MIN(2_i4,km)

   a(:,:)   = 0.0_r8
   b(:,:)   = 0.0_r8
   c(:,:)   = 0.0_r8
   msk(:,:) = 0.0_r8
   v(:)     = 0.0_r8
   u(:)     = 0.0_r8
   c_s(:,:) = 0.0_r8

   DO k = 1,km
      DO j = 1,jm
         jp((k-1)*jm+j) = j
      ENDDO
   ENDDO

   DO k = 1,km
      DO j = 1,jm
         IF ( grd%tmx((k-1)*grd%jm+j) .EQ. 1 .OR. grd%tmx((k-1)*grd%jm+j) .EQ. 2 ) THEN
            DO i = 1,im
               a(i,(k-1)*jm+j) = fld(i,j,k)
            ENDDO
         ENDIF
         IF ( grd%tmx((k-1)*grd%jm+j) .EQ. 2 ) THEN
            DO i = 1-ias,im+iae
               msk(i,(k-1)*jm+j) =grd%msr(i,j,k)
            ENDDO
         ELSEIF ( grd%tmx((k-1)*grd%jm+j) .EQ. 1 ) THEN
            DO i=1-ias,im+iae
               msk(i,(k-1)*jm+j) = 1.0 
            ENDDO
         ENDIF
      ENDDO
   ENDDO

! positive direction
   DO lr = 1,mpi%thj(ik)

      js = grd%jrs(lr,ik)
      je = grd%jre(lr,ik)
      ! MPI receive
      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .GT. 1) THEN
         CALL MPI_IRECV( bffr, npnt*grd%jmr(lr,ik), mpi%r8, mpi%lft, 1, mpi%comm, irecvi, ierr)
         CALL MPI_WAIT( irecvi, istatus, ierr)
         DO j = js,je
            b(-2,j) = bffr(j-js+1                    )
            b(-1,j) = bffr(j-js+1 + grd%jmr(lr,ik)   )
            b( 0,j) = bffr(j-js+1 + grd%jmr(lr,ik)*2 )
         ENDDO
      ENDIF

      IF ( mpi%ir .EQ. 1 ) THEN
         is = 1
         DO j=js,je
            b( 1,j) = del(1,jp(j))*a(1,j)*msk(1,j)
         END DO
      ELSE
         is = 0
      ENDIF

      DO j = js,je
         IF ( grd%tmx(j) .EQ. 1 ) THEN
            DO i = 1+is,im
               b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ gam(i,jp(j))*b(i-3,j)+del(i,jp(j))*a(i,j)
            ENDDO
         ELSEIF ( grd%tmx(j) .EQ. 2 ) THEN
            DO i = 1+is,im
               IF ( msk(i,j).EQ.1.0 .AND. msk(i-1,j).EQ.0.0_r8 ) THEN
                  b(i,j) = del(i,jp(j))*a(i,j)
               ELSEIF ( msk(i,j).EQ.1.0 .AND. msk(i-1,j) .EQ. 1.0 .AND. msk(i-2,j).EQ.0) THEN
                  b(i,j) = alp(i,jp(j))*b(i-1,j)+ del(i,jp(j))*a(i,j)
               ELSEIF ( msk(i,j).EQ.1.0 .AND. msk(i-1,j).EQ.1.0 .AND. msk(i-2,j).EQ.1.0 .AND. msk(i-3,j).EQ.0 ) THEN
                  b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ del(i,jp(j))*a(i,j)
               ELSEIF ( msk(i,j).EQ.1.0 .AND. msk(i-1,j).EQ.1.0 .AND. msk(i-2,j).EQ. 1.0 .AND. msk(i-3,j).EQ.1 ) THEN
                  b(i,j) = alp(i,jp(j))*b(i-1,j) + bta(i,jp(j))*b(i-2,j)+ gam(i,jp(j))*b(i-3,j)+del(i,jp(j))*a(i,j)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

! mpi send
      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .LT. mpi%irm ) THEN
         DO j = js,je
            bffs(j-js+1                    ) = b(im-2,j)
            bffs(j-js+1 + grd%jmr(lr,ik)   ) = b(im-1,j)
            bffs(j-js+1 + grd%jmr(lr,ik)*2 ) = b(im  ,j)
         ENDDO
         CALL mpi_ISEND( bffs, npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, isENDi, ierr)
         CALL MPI_WAIT( isENDi, istatus, ierr)
      ENDIF

   ENDDO

   !negative direction
   DO lr = 1,mpi%thj(ik)

      js = grd%jrs(lr,ik)
      je = grd%jre(lr,ik)
      ! MPI receive
      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .LT. mpi%irm ) THEN
         CALL MPI_IRECV( bffr, npnt*grd%jmr(lr,ik), mpi%r8, mpi%rgh, 1, mpi%comm, irecvi, ierr)
         CALL MPI_WAIT( irecvi, istatus, ierr)
         DO j = js,je
            c(im+1,j) = bffr(j-js+1                    )
            c(im+2,j) = bffr(j-js+1 + grd%jmr(lr,ik)   )
            c(im+3,j) = bffr(j-js+1 + grd%jmr(lr,ik)*2 )
         ENDDO
      ENDIF

      is=0

      DO j = js,je
         IF ( grd%tmx(j) .EQ. 1 ) THEN

            DO i=im-is,1,-1
               c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*c(i+3,j)+del(i,jp(j))*b(i,j)
            ENDDO
            c_s(1:3,j) = c(1:3,j)

         ELSEIF ( grd%tmx(j) .EQ. 2 ) THEN

            IF ( msk(im+1,j).EQ.1.0 .AND. msk(im+2,j).EQ.1.0 .AND. msk(im+3,j).EQ.0.0_r8 ) THEN
               v(2) = c(im+3,j)
            ELSEIF ( msk(im+1,j) .EQ. 1.0 .AND. msk(im+2,j) .EQ. 0.0_r8 ) THEN
               v(2) = c(im+2,j)
               v(3) = c(im+3,j)
            ENDIF
            DO i = im-is,1,-1
               IF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 0.0_r8 ) THEN
                  coef = del(i,jp(j))
                  DO ii=1,3
                     u(ii) = b(i+1-ii,  j)*coef
                     v(ii) = 0.0_r8
                  END DO
                  v(1)   = v(1)+mat_bc_x(1,i,jp(j))*u(1)+mat_bc_x(2,i,jp(j))*u(2)+mat_bc_x(3,i,jp(j))*u(3)
                  v(2)   = v(2)+mat_bc_x(4,i,jp(j))*u(1)+mat_bc_x(5,i,jp(j))*u(2)+mat_bc_x(6,i,jp(j))*u(3)
                  v(3)   = v(3)+mat_bc_x(7,i,jp(j))*u(1)+mat_bc_x(8,i,jp(j))*u(2)+mat_bc_x(9,i,jp(j))*u(3)
                  c(i,j) = v(1)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 0.0_r8 ) THEN
                  c(i,j) =  alp(i,jp(j))*C(i+1,j) + bta(i,jp(j))*v(2)+ gam(i,jp(j))*v(3)+del(i,jp(j))*b(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 1.0 .AND. msk(i+3,j) .EQ. 0.0_r8 ) THEN
                  c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*v(2)+del(i,jp(j))*b(i,j)
               ELSEIF ( msk(i,j) .EQ. 1.0 .AND. msk(i+1,j) .EQ. 1.0 .AND. msk(i+2,j) .EQ. 1.0 .AND. msk(i+3,j) .EQ. 1.0 ) THEN
                  c(i,j) = alp(i,jp(j))*c(i+1,j) + bta(i,jp(j))*c(i+2,j)+ gam(i,jp(j))*c(i+3,j)+del(i,jp(j))*b(i,j)
               ENDIF
            ENDDO
            c_s(1:3,j) = c(1:3,j)
            IF ( msk(1,j) .EQ. 1.0 .AND. msk(2,j) .EQ. 1.0 .AND. msk(3,j) .EQ. 0.0_r8 ) THEN
               c_s(3,j) = v(2)
            ELSEIF ( msk(1,j) .EQ. 1.0 .AND. msk(2,j) .EQ. 0.0_r8 ) THEN
               c_s(3,j) = v(3)
               c_s(2,j) = v(2)
            ENDIF

         ENDIF
      ENDDO


! MPI send
      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .GT. 1 ) THEN
         DO j = js,je
            bffs(j-js+1                   ) = c_s(1,j)
            bffs(j-js+1 + grd%jmr(lr,ik)  ) = c_s(2,j)
            bffs(j-js+1 + grd%jmr(lr,ik)*2) = c_s(3,j)
         ENDDO
         CALL mpi_ISEND( bffs, npnt*grd%jmr(lr,ik), mpi%r8, mpi%lft, 1, mpi%comm, isENDi, ierr)
         CALL MPI_WAIT( isENDi, istatus, ierr)
      ENDIF

   ENDDO

   DO k = 1,km
      DO j = 1,jm
         IF ( grd%tmx((k-1)*grd%jm+j) .EQ. 1 .OR. grd%tmx((k-1)*grd%jm+j) .EQ. 2 ) THEN
            DO i = 1,im
               fld(i,j,k) = c(i,(k-1)*jm+j) 
            ENDDO
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE ( jp )
   DEALLOCATE ( istatus )
   DEALLOCATE ( bffr, bffs )
   DEALLOCATE ( v, u, y )
   DEALLOCATE ( a, b, c, c_s, msk )

END SUBROUTINE rcfl_x
