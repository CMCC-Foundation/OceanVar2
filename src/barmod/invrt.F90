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
!> Implicit solver - overrelaxation                                     
!!
!! It solves:
!!  dn(t+1)-dn(t-1) = (du(t)/dx+ dv(t)/dy)
!!
! Version 1  : Srdjan Dobricic 2007                                    !
! Version 1.1: Paolo  Oddo     2014                                    !
!-----------------------------------------------------------------------
SUBROUTINE invrt( kstp )

   USE set_knd
   USE bmd_str
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   REAL(r8), ALLOCATABLE  :: res(:,:)
   REAL(r8)       :: reser, resem, reserp
   INTEGER(i4)    :: i, j, icnt, kstp
   INTEGER        :: ierr

   ALLOCATE ( res(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

! MPI exchange
   IF ( mpi%nproc .GT. 1 )                                                                  THEN
      CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
   ENDIF

   reser = 1.e20_r8

   bmd%itr(kstp) = 0_i4

   DO icnt=1,bmd%ncnt      ! DO CONVERGENCE LOOP

      IF ( reser .GT. bmd%resem ) THEN  ! IF CHECK RESIDUAL

         bmd%itr(kstp) = bmd%itr(kstp) + 1
         ! ------------------------------------------------------------------------------
         DO j = grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
            DO i = grd%is2,grd%im-1+grd%iae,2       ! 2,jm-1,2
               res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                          bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                          bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
               bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
            ENDDO
         ENDDO

         DO j = grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
            DO i = grd%is3,grd%im-1+grd%iae,2       ! 3,im-1,2
               res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                          bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                          bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
               bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
            ENDDO
         ENDDO
         ! ------------------------------------------------------------------------------

! MPI exchange
         IF ( mpi%nproc .GT. 1 ) THEN
            CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,   &
                          1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
         ENDIF

         ! ------------------------------------------------------------------------------
         DO j = grd%js3,grd%jm-1+grd%jae,2          ! 3,jm-1,2
            DO i = grd%is2,grd%im-1+grd%iae,2       ! 2,im-1,2
               res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                          bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                          bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
               bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
            ENDDO
         ENDDO

         DO j = grd%js2,grd%jm-1+grd%jae,2          ! 2,jm-1,2
            DO i = grd%is3,grd%im-1+grd%iae,2       ! 3,im-1,2
               res(i,j) = bmd%a1(i,j)*bmd%eta(i+1,j  )+bmd%a2(i,j)*bmd%eta(i-1,j  )+        &
                          bmd%a3(i,j)*bmd%eta(i  ,j+1)+bmd%a4(i,j)*bmd%eta(i  ,j-1)-        &
                          bmd%a0(i,j)*bmd%eta(i  ,j  ) - bmd%rgh(i,j)
               bmd%eta(i,j) = bmd%eta(i,j)+bmd%ovr*res(i,j)/bmd%a0(i,j)*bmd%mst(i,j)
            ENDDO
         ENDDO
         ! ------------------------------------------------------------------------------

! MPI exchange
         IF ( mpi%nproc .GT. 1 ) THEN
            CALL exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                          1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
         ENDIF

         reser = 0.0_r8
         DO j = 2-grd%jas,grd%jm-1+grd%jae                  ! 2,jm-1
            DO i = 2-grd%ias,grd%im-1+grd%iae               ! 2,im-1
               reser = reser + ABS(res(i,j))
            ENDDO
         ENDDO
         IF ( bmd%bnm .GT. 0.0_r8 ) reser = reser/bmd%bnm

         IF ( mpi%nproc .GT. 1 ) THEN
            reserp = reser
            CALL MPI_ALLREDUCE( reserp, reser, 1, mpi%r8, MPI_SUM, mpi%comm, ierr)
         ENDIF

      ENDIF                  ! IF CHECK RESIDUAL

   ENDDO                   ! DO CONVERGENCE LOOP

   DEALLOCATE ( res )

END SUBROUTINE invrt
