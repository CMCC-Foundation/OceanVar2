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
!>  Structure of adjoint checks
!                                                                      
! Version 1: Mario Adani    2023                                       !
!-----------------------------------------------------------------------
MODULE adjck_str

   IMPLICIT NONE

!>  Logical variable to activate the adjoint checks
   TYPE adjck_t

      LOGICAL    ::      bal  !< adjoint check for dynamic height
      LOGICAL    ::      bmd  !< adjoint check for barotropic model
      LOGICAL    ::      byg  !< adjoint check for buoyancy force
      LOGICAL    ::      dfl  !< adjoint check for diffusive filter

   END TYPE adjck_t

   TYPE(adjck_t) :: adjck !< initialize derived type

END MODULE adjck_str
