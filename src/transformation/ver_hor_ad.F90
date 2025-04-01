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
!> Apply transformation from physical to control space                        !
!
!! - divergence dumping (adjoint)
!! - velocity  (adjoint)
!! - barotropic model/dynamic height (adjoint)
!! - buoyancy (adjoint)
!! - recursive/diffusive filter (adjoint)
!! - EOFs (adjoint)
!                                                                      !
! Version 1: Srdjan Dobricic                     2006                  !
! Version 2: Srdjan Dobricic                     2007                  !
! Version 3: Srdjan Dobricic and R. Farina       2013                  !
!     Symmetric calculation in presence of coastal boundaries          !
!     eta_ad, tem_ad, and sal_ad are here temporary arrays             !
! Version 4: Mario Adani and  Francesco Carere   2023                  !
!     Diffusion Filter and Balance model based on Dynamic Height       !
!-----------------------------------------------------------------------
SUBROUTINE ver_hor_ad

   USE set_knd
   USE grd_str
   USE eof_str
   USE cns_str
   USE drv_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)    :: k, ione, ierr, i, j, iter,init, niter
   LOGICAL        :: eta_condition

   ione = 1

   eta_condition = ( (drv%bmd(drv%ktr)+drv%bal(drv%ktr) .EQ. 0)                      .OR.  &
                     ((drv%bal(drv%ktr) .EQ. 1) .AND. (drv%ssh_unbalanced(drv%ktr))) )

! ---
! Divergence damping
   IF ( drv%dda(drv%ktr) .EQ. 1 ) THEN
      CALL div_dmp_ad
   ENDIF

! ---
! Velocity
   CALL get_vel_ad

! ---
! Simplified balance operator
   IF ( drv%bal(drv%ktr) .EQ. 1 ) THEN
      CALL bal_op_ad
   ENDIF

! ---
! Barotropic model
   IF (drv%bmd(drv%ktr) .EQ. 1) THEN
      CALL bar_mod_ad
   ENDIF

! ---
! Bouyancy force
      CALL get_byg_ad

! ---
! Horizontal localization
   IF (eta_condition) THEN
      grd%eta_ad(1:grd%im,1:grd%jm) = grd%eta_ad(1:grd%im,1:grd%jm)   * grd%loc(1:grd%im,1:grd%jm)
   ENDIF
   DO k=1,grd%km
      grd%tem_ad(1:grd%im,1:grd%jm,k) = grd%tem_ad(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
      grd%sal_ad(1:grd%im,1:grd%jm,k) = grd%sal_ad(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
   ENDDO

! ---
! Horizontal filter
   IF (drv%filter .EQ. 1 ) THEN
! Recursive Filter Adj
      CALL recursive_filter_ad

   ELSEIF ( drv%filter .EQ. 2 ) THEN
! Diffusion Filter Adj
      CALL diffusive_filter_ad
   ENDIF

! ---
! Vertical EOFs
   CALL veof_ad

END SUBROUTINE ver_hor_ad
