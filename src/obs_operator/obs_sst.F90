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
!> Apply observational operator for SST                                
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Jenny Pistoia   2013                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_sst

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, kk, k

   DO kk = 1,sst%no
      IF ( sst%flc(kk) .EQ. 1 ) THEN
         i = sst%ib(kk)
         j = sst%jb(kk)
         k = sst%kb(kk)
         sst%inc(kk) = sst%pq1(kk) * grd%tem(i  ,j ,k ) +       &
                       sst%pq2(kk) * grd%tem(i+1,j ,k ) +       &
                       sst%pq3(kk) * grd%tem(i  ,j+1,k) +       &
                       sst%pq4(kk) * grd%tem(i+1,j+1,k)
      ENDIF
   ENDDO

END SUBROUTINE obs_sst
