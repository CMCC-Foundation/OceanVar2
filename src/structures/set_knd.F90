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
!> The precision of real and integer 
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
MODULE set_knd

   IMPLICIT NONE

   PUBLIC

   INTEGER, PARAMETER ::                &
      r4 = SELECTED_REAL_KIND( 6, 37),  &  !< REAL*4
      r8 = SELECTED_REAL_KIND(12,307)      !< REAL*8

   INTEGER, PARAMETER ::                &
      i4 = SELECTED_INT_KIND(9) ,       &  !< INTEGER*4
      i8 = SELECTED_INT_KIND(14)           !< INTEGER*8

END MODULE set_knd
