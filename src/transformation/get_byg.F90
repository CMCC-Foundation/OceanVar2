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
!> Calculate vertical integral of bouyancy gradient                     
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
! Version 1.1: P.Oddo 2014                                             !
! Version 2:   Mario Adani 2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE get_byg

   USE set_knd
   USE drv_str
   USE grd_str
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INCLUDE 'mpif.h'

   INTEGER(i4)           :: i, j, k
   REAL   (r8)           :: rho_unescotl
   REAL   (r8)           :: grho0

   grho0 = phy%g / phy%rho0

! ---
! Assume that density is a simple linear function of temperature and salinity

   IF ( drv%nneos(drv%ktr) .EQ. 1 ) THEN

      grd%dns(:,:,:) = (-grd%alpha*grd%tem(:,:,:) +          &
                         grd%beta*grd%sal(:,:,:))            &
                         * grd%msk(:,:,:) * grho0

   ELSEIF ( drv%nneos(drv%ktr) .EQ. 2 ) THEN

      grd%dns(:,:,:) = (-grd%alpha3d(:,:,:)*grd%tem(:,:,:) +  &
                         grd%beta3d (:,:,:)*grd%sal(:,:,:))   &
                         * grd%msk(:,:,:) * grho0

   ELSEIF ( drv%nneos(drv%ktr) .EQ. 3 ) THEN

      DO k = 1,grd%km
         DO j = 1-grd%jas,grd%jm
            DO i = 1-grd%ias,grd%im
               grd%dns(i,j,k) = rho_unescotl(grd%salb(i,j,k),grd%temb(i,j,k),&
                                             grd%sal(i,j,k),grd%tem(i,j,k))  &
                                 * grd%msk(i,j,k) * grho0
            ENDDO
         ENDDO
      ENDDO

   ELSE

      WRITE (drv%dia,*)' -----------------------------------------------'
      WRITE (drv%dia,*)' get_byg: Unsupported  equation of state option '
      WRITE (drv%dia,*)' Please choose drv%nneos from 1 to 3.           '
      WRITE (drv%dia,*)' -----------------------------------------------'
      CALL FLUSH(drv%dia)
      CALL abort

   ENDIF

   IF ( mpi%nproc .GT. 1 ) THEN

      CALL exo_mpi( 0_i4, 1_i4,                                       &
                    1_i4, grd%im, 0_i4,                               &
                    1_i4, grd%jm, 0_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%dns)

   ENDIF

! ---
! Bouyancy force

   grd%b_x(:,:,:) = 0.0_r8
   grd%b_y(:,:,:) = 0.0_r8

   DO j = 2-grd%jas,grd%jm
      DO i = 2-grd%ias,grd%im
         grd%b_x(i,j,1) = (grd%dns(i,j,1)-grd%dns(i-1,j,1))*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i-1,j,1)
         grd%b_y(i,j,1) = (grd%dns(i,j,1)-grd%dns(i,j-1,1))*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i,j-1,1)
      ENDDO
   ENDDO

   DO k = 2,grd%km
      DO j = 2-grd%jas,grd%jm
         DO i = 2-grd%ias,grd%im
            grd%b_x(i,j,k) = grd%b_x(i,j,k-1) + (grd%dns(i,j,k)-grd%dns(i-1,j,k))*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i-1,j,k)
            grd%b_y(i,j,k) = grd%b_y(i,j,k-1) + (grd%dns(i,j,k)-grd%dns(i,j-1,k))*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i,j-1,k)
         ENDDO
      ENDDO
   ENDDO

! ---
! Verical integral of bouyancy force
   grd%bx(:,:) = 0.0_r8
   grd%by(:,:) = 0.0_r8
   DO k = 1,grd%km
      DO j = 2-grd%jas,grd%jm
         DO i = 2-grd%ias,grd%im
            grd%bx(i,j) = grd%bx(i,j) + grd%b_x(i,j,k)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            grd%by(i,j) = grd%by(i,j) + grd%b_y(i,j,k)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
         ENDDO
      ENDDO
   ENDDO

   DO j = 2-grd%jas,grd%jm
      DO i = 2-grd%ias,grd%im
         grd%bx(i,j) = grd%bx(i,j)/grd%dx(i,j)
         grd%by(i,j) = grd%by(i,j)/grd%dy(i,j)
      ENDDO
   ENDDO

END SUBROUTINE get_byg
