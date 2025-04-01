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
!> Convert from control to v                                          
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE cnv_ctv

   USE set_knd
   USE grd_str
   USE ctl_str
   USE eof_str

   IMPLICIT NONE

   INTEGER(i4)   :: i,j,k, kk

   kk = 0
   DO k = 1,ros%neof
      DO j = 1,grd%jm
         DO i = 1,grd%im
            kk = kk+1
            grd%ro(i,j,k) = ctl%x_c(kk)
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE cnv_ctv
