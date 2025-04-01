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
!>  Adjoint of balance operator based on dynamic height
!!
!! It computes density anomalies from sea surface elevation
!! and temperture and salinity anomalies from density using adjoint
!! version of equation of state. 
!!
!                                                                      !
! Version 1: Andrea Storto 2021                                        !
!            Mario  Adani  2024                                        !
!-----------------------------------------------------------------------
SUBROUTINE bal_op_ad

   USE set_knd
   USE grd_str
   USE drv_str
   USE bal_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)    :: k,j,i
   REAL(r8)       :: tad,sad

!-----------------------------
!Compute dynamic height adjoint
   grd%dns = 0._r8
   DO k = bal%nlevs,1,-1
      grd%dns(:,:,k) = - bal%dhdz(k)*grd%eta_ad(:,:) / phy%rho0
   ENDDO

!-----------------------------
!Compute dynamic height adjoint
   IF ( drv%nneos(drv%ktr) .EQ. 1 ) THEN
      DO k = 1,grd%km
         grd%tem_ad(:,:,k) = grd%tem_ad(:,:,k) - grd%alpha*grd%dns(:,:,k) * grd%msk(:,:,k)
         grd%sal_ad(:,:,k) = grd%sal_ad(:,:,k) + grd%beta *grd%dns(:,:,k) * grd%msk(:,:,k)
      ENDDO
   ELSEIF ( drv%nneos(drv%ktr) .EQ. 2 ) THEN
      DO k = 1,grd%km
         grd%tem_ad(:,:,k) = grd%tem_ad(:,:,k) - grd%alpha3d(:,:,k)*grd%dns(:,:,k) * grd%msk(:,:,k)
         grd%sal_ad(:,:,k) = grd%sal_ad(:,:,k) + grd%beta3d (:,:,k)*grd%dns(:,:,k) * grd%msk(:,:,k)
      ENDDO
   ELSEIF ( drv%nneos(drv%ktr) .EQ. 3 ) THEN
      DO k = 1,grd%km
         DO j = 1-grd%jas,grd%jm+grd%jae
            DO i = 1-grd%ias,grd%im+grd%iae
               CALL rho_unescoad(grd%dns(i,j,k),grd%salb(i,j,k),grd%temb(i,j,k),sad, tad)
               grd%tem_ad(i,j,k) = grd%tem_ad(i,j,k) + tad*grd%msk(i,j,k)
               grd%sal_ad(i,j,k) = grd%sal_ad(i,j,k) + sad*grd%msk(i,j,k)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      write(drv%dia,*)' -------------------------------------------------'
      write(drv%dia,*)' bal_op_ad: Unsupported  equation of state option.'
      write(drv%dia,*)' Please choose drv%nneos from 1 to 3.             '
      write(drv%dia,*)' -------------------------------------------------'
      CALL abort
   ENDIF

END SUBROUTINE bal_op_ad
