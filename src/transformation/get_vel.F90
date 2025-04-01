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
!> Calculate horizontal velocity from geostrophic formula               
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
! Bug correction 21.04.2009  thanks to Andrea Storto                   !
!            Mario Adani 2023 (Limit geostrophy (+-2degrees))          !
!-----------------------------------------------------------------------
SUBROUTINE get_vel

   USE set_knd
   USE drv_str
   USE grd_str
   USE mpi_str
   USE bal_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   REAL(r8), ALLOCATABLE, DIMENSION (:,:,:)  :: ud, vd

   INTEGER(i4)                               :: k, i, j, lr, lvl


   ALLOCATE ( ud(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( vd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )

   IF ( drv%bal(drv%ktr) .EQ. 1 ) THEN
      lvl = bal%nlevs
   ELSE
      lvl = grd%km
   ENDIF

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 0_i4, 1_i4,                                       &
                    1_i4, grd%im, 0_i4,                               &
                    1_i4, grd%jm, 0_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, grd%eta)
   ENDIF

   DO k = grd%km,2,-1
      grd%b_x(:,:,k) = ( grd%b_x(:,:,k) + grd%b_x(:,:,k-1) ) * 0.5
      grd%b_y(:,:,k) = ( grd%b_y(:,:,k) + grd%b_y(:,:,k-1) ) * 0.5
   ENDDO

   grd%b_x(:,:,1) = grd%b_x(:,:,1) * 0.5
   grd%b_y(:,:,1) = grd%b_y(:,:,1) * 0.5

   grd%uvl(:,:,:) = 0.0_r8
   grd%vvl(:,:,:) = 0.0_r8

   ud(:,:,:) = 0.0_r8
   vd(:,:,:) = 0.0_r8

   DO k=1,lvl 

      DO j = 1,grd%jm
       DO i = 2-grd%ias,grd%im
            IF ( ABS(grd%lat(i,j)) .GT. 2._r8 )  THEN
               vd(i,j,k) =   ( (grd%eta(i,j)-grd%eta(i-1,j))*phy%g + grd%b_x(i,j,k) )               &
                  / grd%dx(i,j) * grd%msk(i,j,k)*grd%msk(i-1,j,k)  / grd%f(i,j)
            ENDIF
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm
       DO i = 1,grd%im
            IF ( ABS(grd%lat(i,j)) .GT. 2._r8 ) THEN
               ud(i,j,k) = - ( (grd%eta(i,j)-grd%eta(i,j-1))*phy%g + grd%b_y(i,j,k) )               &
                  / grd%dy(i,j) * grd%msk(i,j,k)*grd%msk(i,j-1,k) / grd%f(i,j)
            ENDIF
         ENDDO
      ENDDO

   ENDDO

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 1_i4,                                       &
                    1_i4, grd%im, 0_i4,                               &
                   -1_i4, 1_i4, grd%jm+1,                             &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, ud)
      CALL exo_mpi( 1_i4, 1_i4,                                       &
                   -1_i4, 1_i4, grd%im+1,                             &
                    1_i4, grd%jm, 0_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, vd)
   ENDIF

   DO k = 1,lvl
      DO j = 1,grd%jm-1+grd%jae
         DO i = 2-grd%ias,grd%im
            grd%uvl(i,j,k) = ( ud(i,j,k)+ud(i-1,j,k)+ud(i,j+1,k)+ud(i-1,j+1,k) )*0.25 * grd%msk(i,j,k)*grd%msk(i-1,j,k)
         ENDDO
      ENDDO
      DO j = 2-grd%jas,grd%jm      ! 2:grd%jm
         DO i = 1,grd%im-1+grd%iae   ! 1:grd%im-1
            grd%vvl(i,j,k) = ( vd(i,j,k)+vd(i,j-1,k)+vd(i+1,j,k)+vd(i+1,j-1,k) )*0.25  * grd%msk(i,j,k)*grd%msk(i,j-1,k)
         ENDDO
      ENDDO
   ENDDO

   DEALLOCATE ( ud, vd )

END SUBROUTINE get_vel
