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
!> Calculate horizontal velocity from geostrophic formula (adjoint)     
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Bug correction 21.04.2009  thanks to Andrea Storto                   !
!            Mario Adani 2023 (Limit geostrophy (+-2degrees))          !
!-----------------------------------------------------------------------
SUBROUTINE get_vel_ad

   USE set_knd
   USE grd_str
   USE bmd_str
   USE mpi_str
   USE drv_str
   USE bal_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   REAL(r8), ALLOCATABLE, DIMENSION (:,:,:)  :: ud, vd
   INTEGER(i4)                               :: k, i, j, lvl

   ALLOCATE ( ud(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( vd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )

   IF ( drv%bal(drv%ktr) .EQ. 1 ) THEN
      lvl = bal%nlevs
   ELSE
      lvl = grd%km
   ENDIF

   ud(:,:,:) = 0.0_r8
   vd(:,:,:) = 0.0_r8

   DO k = lvl,1,-1

      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .LT. mpi%irm) vd(grd%im+1,:,k) = 0.0_r8
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .GT. 1      ) vd(:,0,k)        = 0.0_r8
      DO j = 2-grd%jas,grd%jm      ! 2:grd%jm
         DO i = 1,grd%im-1+grd%iae   ! 1:grd%im-1
            vd(i  ,j  ,k) = vd(i  ,j  ,k) + grd%vvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
            vd(i  ,j-1,k) = vd(i  ,j-1,k) + grd%vvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
            vd(i+1,j  ,k) = vd(i+1,j  ,k) + grd%vvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
            vd(i+1,j-1,k) = vd(i+1,j-1,k) + grd%vvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i,j-1,k)
            grd%vvl_ad(i,j,k) = 0.0_r8
         ENDDO
      ENDDO

      IF ( mpi%nproc .GT. 1 .AND. mpi%ir .GT. 1       ) ud(0,:       ,k) = 0.0_r8
      IF ( mpi%nproc .GT. 1 .AND. mpi%jr .LT. mpi%jrm ) ud(:,grd%jm+1,k) = 0.0_r8
      DO j = 1,grd%jm-1+grd%jae
         DO i = 2-grd%ias,grd%im
            ud(i  ,j  ,k) = ud(i  ,j  ,k) + grd%uvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            ud(i-1,j  ,k) = ud(i-1,j  ,k) + grd%uvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            ud(i  ,j+1,k) = ud(i  ,j+1,k) + grd%uvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            ud(i-1,j+1,k) = ud(i-1,j+1,k) + grd%uvl_ad(i,j,k)*0.25_r8 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            grd%uvl_ad(i,j,k) = 0.0_r8
         ENDDO
      ENDDO

   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                   -1_i4, 0_i4,   grd%im,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, ud )
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                   -1_i4, 0_i4,   grd%jm,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, grd%km, vd )
   ENDIF

   IF ( mpi%nproc .GT. 1 .AND. mpi%jr .GT. 1 ) grd%eta_ad(:,0) = 0.0_r8
   IF ( mpi%nproc .GT. 1 .AND. mpi%ir .GT. 1 ) grd%eta_ad(0,:) = 0.0_r8
   DO k = 1,lvl
      DO j = 2-grd%jas,grd%jm      ! 2:grd%jm
       DO i = 1,grd%im               ! 1:grd%im
            IF ( ABS(grd%lat(i,j)) .GT. 2._r8 )  THEN
               grd%b_y(i,j,k)    = grd%b_y(i,j,k)    - ud(i,j,k)         / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
               grd%eta_ad(i,j  ) = grd%eta_ad(i,j  ) - ud(i,j,k)*phy%g / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
               grd%eta_ad(i,j-1) = grd%eta_ad(i,j-1) + ud(i,j,k)*phy%g / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
            ENDIF
         ENDDO
      ENDDO
      DO j = 1,grd%jm
       DO i = 2-grd%ias,grd%im
            IF ( ABS(grd%lat(i,j)) .GT. 2._r8 )  THEN
               grd%b_x(i,j,k)    = grd%b_x(i,j,k)    + vd(i,j,k)         / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
               grd%eta_ad(i-1,j) = grd%eta_ad(i-1,j) - vd(i,j,k)*phy%g / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
               grd%eta_ad(i  ,j) = grd%eta_ad(i  ,j) + vd(i,j,k)*phy%g / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k) / grd%f(i,j)
            ENDIF
         ENDDO
      ENDDO
      ud(:,:,k) = 0.0_r8
      vd(:,:,k) = 0.0_r8
   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 0_i4, 0_i4,                                         &
                   -1_i4, 0_i4,   grd%im,                               &
                   -1_i4, 0_i4,   grd%jm,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%eta_ad )
   ENDIF

   grd%uvl_ad(:,:,:) = 0.0_r8
   grd%vvl_ad(:,:,:) = 0.0_r8

   grd%b_y(:,2-grd%jas:grd%jm,1) = grd%b_y(:,2-grd%jas:grd%jm,1) * 0.5_r8
   grd%b_x(2-grd%ias:grd%im,:,1) = grd%b_x(2-grd%ias:grd%im,:,1) * 0.5_r8

   DO k = 2,grd%km
      grd%b_y(:,2-grd%jas:grd%jm,k-1) = grd%b_y(:,2-grd%jas:grd%jm,k-1) + grd%b_y(:,2-grd%jas:grd%jm,k) * 0.5_r8
      grd%b_y(:,2-grd%jas:grd%jm,k)   = grd%b_y(:,2-grd%jas:grd%jm,k) * 0.5
      grd%b_x(2-grd%ias:grd%im,:,k-1) = grd%b_x(2-grd%ias:grd%im,:,k-1) + grd%b_x(2-grd%ias:grd%im,:,k) * 0.5_r8
      grd%b_x(2-grd%ias:grd%im,:,k)   = grd%b_x(2-grd%ias:grd%im,:,k) * 0.5
   ENDDO

   DEALLOCATE ( ud, vd )

END SUBROUTINE get_vel_ad
