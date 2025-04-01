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
!> Data structure for debugging
!!
!! Cost function, control vector and optimisation arrays                
!! Some variables for debuging                                       
!!
!                                                                      !
! Version 1: Francesco Carere 2024                                     !
!-----------------------------------------------------------------------
MODULE rpr_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

!---
!> Variable for reporducibility
   TYPE rpr_t

      INTEGER(i8)               ::  obs_nog          !< Number of global observations
      REAL(r8), ALLOCATABLE     ::  obs_amo(:)       !< Number of threads in i direction (recursive filter)
      INTEGER(i4), POINTER      ::  idxmap(:)        !< Map of the grid index in serial execution

   END TYPE rpr_t

   TYPE (rpr_t)                  :: rpr              !< initialize derived type

END MODULE rpr_str
