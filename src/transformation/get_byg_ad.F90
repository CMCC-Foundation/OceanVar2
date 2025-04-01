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
!> Calculate vertical integral of bouyancy gradient (adjoint)           
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
! Version 2: Mario Adani and Francesco Carere 2024                     !
!-----------------------------------------------------------------------
SUBROUTINE get_byg_ad

   USE set_knd
   USE drv_str
   USE grd_str
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)    :: i, j, k
   REAL(r8)       :: tad,sad
   REAL(r8)       :: grho0

   grho0 = phy%g / phy%rho0

   grd%dns(:,:,:) = 0.0_r8

   DO j = 2-grd%jas,grd%jm
      DO i = 2-grd%ias,grd%im
         grd%bx(i,j) = grd%bx(i,j)/grd%dx(i,j)
         grd%by(i,j) = grd%by(i,j)/grd%dy(i,j)
      ENDDO
   ENDDO

   DO k = 1,grd%km
      DO j = 2-grd%jas,grd%jm
         DO i = 2-grd%ias,grd%im
            grd%b_x(i,j,k) = grd%b_x(i,j,k) + grd%bx(i,j)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i-1,j,k)
            grd%b_y(i,j,k) = grd%b_y(i,j,k) + grd%by(i,j)*grd%dz(k) * grd%msk(i,j,k)*grd%msk(i,j-1,k)
         ENDDO
      ENDDO
   ENDDO

! j --------------------
   DO k = grd%km,2,-1
      DO j = 2-grd%jas,grd%jm
         DO i = 2-grd%ias,grd%im
            grd%dns(i,j  ,k) = grd%dns(i,j  ,k) +               &
                               grd%b_y(i,j,k)*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i,j-1,k)
            grd%dns(i,j-1,k) = grd%dns(i,j-1,k) -               &
                               grd%b_y(i,j,k)*grd%dz(k)*grd%msk(i,j,k)*grd%msk(i,j-1,k)
            grd%b_y(i,j,k-1) = grd%b_y(i,j,k-1) + grd%b_y(i,j,k)
            grd%b_y(i,j,k)   = 0.0_r8
         ENDDO
      ENDDO
   ENDDO
   DO j = 2-grd%jas,grd%jm
      DO i = 2-grd%ias,grd%im
         grd%dns(i,j  ,1) = grd%dns(i,j  ,1) +               &
                            grd%b_y(i,j,1)*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i,j-1,1)
         grd%dns(i,j-1,1) = grd%dns(i,j-1,1) -               &
                            grd%b_y(i,j,1)*grd%dz(1)*grd%msk(i,j,1)*grd%msk(i,j-1,1)
      ENDDO
   ENDDO
   IF (mpi%nproc.GT.1) THEN

      CALL exo_mpi( 0_i4, 0_i4,                                        &
                    0_i4, 0_i4, grd%im,                               &
                   -1_i4, 0_i4, grd%jm,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%dns)

   ENDIF
! i --------------------

   DO k = grd%km,2,-1
      DO j = 2-grd%jas,grd%jm
         DO i = 2-grd%ias,grd%im
            grd%dns(i  ,j,k) = grd%dns(i  ,j,k) +               &
                               grd%b_x(i  ,j,k)*grd%dz(k)*grd%msk(i  ,j,k)*grd%msk(i-1,j,k)
            grd%dns(i-1,j,k) = grd%dns(i-1,j,k) -               &
                                grd%b_x(i  ,j,k)*grd%dz(k)*grd%msk(i  ,j,k)*grd%msk(i-1,j,k)
            grd%b_x(i  ,j,k-1) = grd%b_x(i  ,j,k-1) + grd%b_x(i  ,j,k)
            grd%b_x(i  ,j,k)   = 0.0_r8
         ENDDO
      ENDDO
   ENDDO
   DO j = 2-grd%jas,grd%jm
      DO i = 2-grd%ias,grd%im
         grd%dns(i  ,j,1) = grd%dns(i  ,j,1) +               &
                            grd%b_x(i  ,j,1)*grd%dz(k)*grd%msk(i  ,j,1)*grd%msk(i-1,j,1)
         grd%dns(i-1,j,1) = grd%dns(i-1,j,1) -               &
                            grd%b_x(i  ,j,1)*grd%dz(k)*grd%msk(i  ,j,1)*grd%msk(i-1,j,1)
      ENDDO
   ENDDO
!--------------------

   IF ( mpi%nproc .GT. 1 ) THEN

      CALL exo_mpi( 0_i4, 0_i4,                                        &
                   -1_i4, 0_i4, grd%im,                               &
                    0_i4, 0_i4, grd%jm,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%dns)

   ENDIF

   IF ( drv%nneos(drv%ktr) .EQ. 1 ) THEN

      grd%tem_ad(:,:,:) = grd%tem_ad(:,:,:) - grd%alpha*grd%dns(:,:,:)  * grd%msk(:,:,:) * grho0
      grd%sal_ad(:,:,:) = grd%sal_ad(:,:,:) + grd%beta *grd%dns(:,:,:)  * grd%msk(:,:,:) * grho0

   ELSEIF ( drv%nneos(drv%ktr) .EQ. 2 ) THEN

      grd%tem_ad(:,:,:) = grd%tem_ad(:,:,:) - grd%alpha3d(:,:,:)*grd%dns(:,:,:) * grd%msk(:,:,:) * grho0
      grd%sal_ad(:,:,:) = grd%sal_ad(:,:,:) + grd%beta3d (:,:,:)*grd%dns(:,:,:) * grd%msk(:,:,:) * grho0

   ELSEIF ( drv%nneos(drv%ktr) .EQ. 3 ) THEN

      DO k = 1,grd%km
         DO j = 1-grd%jas,grd%jm
            DO i = 1-grd%ias,grd%im
               CALL rho_unescoad(grd%dns(i,j,k),grd%salb(i,j,k),grd%temb(i,j,k),sad, tad)
               grd%tem_ad(i,j,k) = grd%tem_ad(i,j,k) + tad*grd%msk(i,j,k) * grho0
               grd%sal_ad(i,j,k) = grd%sal_ad(i,j,k) + sad*grd%msk(i,j,k) * grho0
            ENDDO
         ENDDO
      ENDDO
   ELSE

      WRITE (drv%dia,*)' --------------------------------------------------'
      WRITE (drv%dia,*)' get_byg_ad: Unsupported  equation of state option '
      WRITE (drv%dia,*)' Please choose drv%nneos from 1 to 3.              '
      WRITE (drv%dia,*)' --------------------------------------------------'
      CALL FLUSH(drv%dia)
      CALL abort

   ENDIF

END SUBROUTINE get_byg_ad
