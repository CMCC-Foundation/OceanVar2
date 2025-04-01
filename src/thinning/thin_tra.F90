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
!> Apply  thinning to ARGO  trajectory                                  
!!
!! Not implemented
!!
!                                                                      !
! Version 1: Name Surname Year                                         !
!-----------------------------------------------------------------------
SUBROUTINE thin_tra

   USE drv_str
   USE obs_str, ONLY : thin, tra

   IMPLICIT NONE

   WRITE (drv%dia,*) 'WARNING: Thinning of ARGO  trajectory not implemented'

END SUBROUTINE thin_tra

