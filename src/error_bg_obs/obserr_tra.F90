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
!> Argo trajectories observational error                                !
!!
!! Not implemented! It uses the values read in the file tra_mis.dat
!!
!                                                                      !
! Version 1: Name Surname Year                                         !
!-----------------------------------------------------------------------
SUBROUTINE obserr_tra

   USE set_knd
   USE obs_str, ONLY : tra
   USE drv_str

   IMPLICIT NONE

   WRITE (drv%dia,*) ' --- Argo trajectories error not implemented yet.    '
   WRITE (drv%dia,*) ' ...applying the error read from observational file. '

END SUBROUTINE obserr_tra

