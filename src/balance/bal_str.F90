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
!----------------------------------------------------------------------------
!                                                                           !
!>  Structure for simplified balance operator based on dynamic height
!                                                                           !
! Version 1: Andrea Storto 2021                                             !
!            Mario  Adani  2024                                             !
!----------------------------------------------------------------------------
MODULE bal_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Variable related to dynamic height operator
   TYPE bal_t

      REAL(r8), ALLOCATABLE :: dhdz(:)         !< Width of the layer
      INTEGER(i4)           :: nlevs           !< Model level of no motion
      REAL(r8)              :: lnm             !< Level of no motion depth

   END TYPE bal_t

   TYPE (bal_t)             :: bal             !< initialize derived typ

END MODULE bal_str
