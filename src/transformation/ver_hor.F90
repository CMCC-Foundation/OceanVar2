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
!> Apply trasformations from control to physical space
!!
!! - EOFs
!! - recursive/diffusive filter
!! - buoyancy
!! - barotropic model/dynamic height
!! - velocity 
!! - divergence dumping
!                                                                      !
! Version 1: Srdjan Dobricic                     2006                  !
! Version 2: Srdjan Dobricic                     2007                  !
! Version 3: Srdjan Dobricic and R. Farina       2013                  !
!     Symmetric calculation in presence of coastal boundaries          !
!     eta_ad, tem_ad, and sal_ad are here temporary arrays             !
! Version 4: Mario Adani and  Francesco Carere   2023                  !
!     Diffusion Filter and Balance model based on Dynamic Height       !
!-----------------------------------------------------------------------
SUBROUTINE ver_hor

   USE set_knd
   USE grd_str
   USE eof_str
   USE drv_str
   USE mpi_str
   USE adjck_str

   IMPLICIT NONE

   INTEGER(i4)    :: k
   LOGICAL        :: eta_condition

   eta_condition = ( (drv%bmd(drv%ktr)+drv%bal(drv%ktr) .EQ. 0)                      .OR.  &
                     ((drv%bal(drv%ktr) .EQ. 1) .AND. (drv%ssh_unbalanced(drv%ktr))) )

! ---
! Vertical EOFs
   CALL veof

!---
! Apply Recursive filter
   IF ( drv%filter .EQ. 1 ) THEN
      CALL recursive_filter

!---
! Apply Diffusive filter
   ELSEIF ( drv%filter .EQ. 2 ) THEN
      IF ( (adjck%dfl) .AND. (mpi%irm .EQ. 1) .AND. (mpi%jrm .EQ. 1) ) THEN
         CALL adjck_dfl
      ENDIF

      CALL diffusive_filter
   ENDIF

! ---
! Horizontal localization
   IF (eta_condition) THEN
      grd%eta(1:grd%im,1:grd%jm) = grd%eta(1:grd%im,1:grd%jm)   * grd%loc(1:grd%im,1:grd%jm)
   ENDIF
   DO k=1,grd%km
      grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
      grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
   ENDDO

! ---
! Get Bouyancy force
   CALL get_byg

! ---
! Simplified balance operator
   IF (drv%bal(drv%ktr) .EQ. 1) THEN
      IF ( (adjck%bal) .AND. (mpi%irm .EQ. 1) .AND. (mpi%jrm .EQ. 1) ) THEN
         CALL adjck_balop ! Only serial for now
      ENDIF
      CALL bal_op
   ENDIF

! ---
! Barotropic model
   IF ( drv%bmd(drv%ktr) .EQ. 1 ) THEN
      IF ( (adjck%byg) .AND. (mpi%irm .EQ. 1) .AND. (mpi%jrm .EQ. 1) ) THEN
         CALL adjck_byg ! Only serial for now
      ENDIF
      IF ( (adjck%bmd) .AND. (mpi%irm .EQ. 1) .AND. (mpi%jrm .EQ. 1) ) THEN
         CALL adjck_bmd ! Only serial for now
      ENDIF
      CALL bar_mod
   ENDIF

! ---
! Geostrophic velocity
   CALL get_vel

! ---
! Divergence damping
   IF ( drv%dda(drv%ktr) .EQ. 1 .OR. (drv%ddi(drv%ktr) .EQ. 1 .AND. drv%ktr .EQ. drv%ntr) ) THEN
      CALL div_dmp
   ENDIF

END SUBROUTINE ver_hor
