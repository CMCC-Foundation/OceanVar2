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
!> Divergence damping adjoint                                           
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!-----------------------------------------------------------------------
SUBROUTINE div_dmp_ad

   USE set_knd
   USE drv_str
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)             :: k, kdiv, i, j
   REAL(r8), ALLOCATABLE   :: div(:,:,:)

   ALLOCATE ( div(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )

   DO kdiv = 1,100
      div(:,:,:) = 0.0_r8

      DO k = grd%km,1,-1
         DO j = 2-grd%jas,grd%jm-1+grd%jae
            DO i = 1,grd%im-1+grd%iae
               div(i,j  ,k) = div(i,j  ,k) + grd%vvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
               div(i,j-1,k) = div(i,j-1,k) + grd%vvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
            ENDDO
         ENDDO
         DO j = 1,grd%jm-1+grd%jae
            DO i = 2-grd%ias,grd%im-1+grd%iae
               div(i  ,j,k) = div(i  ,j,k) + grd%uvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
               div(i-1,j,k) = div(i-1,j,k) + grd%uvl(i,j,k)*0.2 * grd%adxdy**2 / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            ENDDO
         ENDDO
      ENDDO

      IF ( mpi%nproc.GT.1 ) THEN
         CALL exo_mpi( 0_i4, 0_i4,                                       &
                      -1_i4, 0_i4, grd%im,                               &
                      -1_i4, 0_i4, grd%jm,                               &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, div)
      ENDIF

      DO k = grd%km,1,-1
         IF ( mpi%nproc .GT. 1 .AND. mpi%ir .LT. mpi%irm ) grd%uvl(grd%im+1,:,:) = 0.0_r8
         IF ( mpi%nproc .GT. 1 .AND. mpi%jr .LT. mpi%jrm ) grd%vvl(:,grd%jm+1,:) = 0.0_r8
         DO j = 1,grd%jm-1+grd%jae
            DO i = 1,grd%im-1+grd%iae
               grd%uvl(i  ,j  ,k) = grd%uvl(i  ,j  ,k) - div(i,j,k) * grd%dy(i  ,j  ) / grd%dx(i,j) / grd%dy(i,j)
               grd%uvl(i+1,j  ,k) = grd%uvl(i+1,j  ,k) + div(i,j,k) * grd%dy(i+1,j  ) / grd%dx(i,j) / grd%dy(i,j)
               grd%vvl(i  ,j  ,k) = grd%vvl(i  ,j  ,k) - div(i,j,k) * grd%dx(i  ,j  ) / grd%dy(i,j) / grd%dx(i,j)
               grd%vvl(i  ,j+1,k) = grd%vvl(i  ,j+1,k) + div(i,j,k) * grd%dx(i  ,j+1) / grd%dy(i,j) / grd%dx(i,j)
            ENDDO
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL exo_mpi( 0_i4, 1_i4,                                       &
                       1_i4, grd%im+1, 1_i4,                             &
                       0_i4, 1_i4, grd%jm+1,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%uvl)
         CALL exo_mpi( 0_i4, 1_i4,                                       &
                       0_i4, 1_i4, grd%im+1,                             &
                       1_i4, grd%jm+1, 1_i4,                             &
                       1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, grd%vvl)
      ENDIF

   ENDDO

   DEALLOCATE ( div )

END SUBROUTINE div_dmp_ad
