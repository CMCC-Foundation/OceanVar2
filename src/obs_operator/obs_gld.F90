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
!> Apply observational operator for gliders                            
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_gld

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,gld%no
      IF ( (gld%flc(kk) .EQ. 1 .OR. gld%fls(kk) .EQ. 1) .AND. gld%par(kk) .EQ. 1 ) THEN
         i = gld%ib(kk)
         j = gld%jb(kk)
         k = gld%kb(kk)
         gld%inc(kk) = gld%pq1(kk) * grd%tem(i  ,j  ,k  ) +       &
                       gld%pq2(kk) * grd%tem(i+1,j  ,k  ) +       &
                       gld%pq3(kk) * grd%tem(i  ,j+1,k  ) +       &
                       gld%pq4(kk) * grd%tem(i+1,j+1,k  ) +       &
                       gld%pq5(kk) * grd%tem(i  ,j  ,k+1) +       &
                       gld%pq6(kk) * grd%tem(i+1,j  ,k+1) +       &
                       gld%pq7(kk) * grd%tem(i  ,j+1,k+1) +       &
                       gld%pq8(kk) * grd%tem(i+1,j+1,k+1)
      ELSEIF ( (gld%flc(kk) .EQ. 1 .OR. gld%fls(kk) .EQ. 1) .AND. gld%par(kk) .EQ. 2 ) THEN
         i = gld%ib(kk)
         j = gld%jb(kk)
         k = gld%kb(kk)
         gld%inc(kk) = gld%pq1(kk) * grd%sal(i  ,j  ,k  ) +       &
                       gld%pq2(kk) * grd%sal(i+1,j  ,k  ) +       &
                       gld%pq3(kk) * grd%sal(i  ,j+1,k  ) +       &
                       gld%pq4(kk) * grd%sal(i+1,j+1,k  ) +       &
                       gld%pq5(kk) * grd%sal(i  ,j  ,k+1) +       &
                       gld%pq6(kk) * grd%sal(i+1,j  ,k+1) +       &
                       gld%pq7(kk) * grd%sal(i  ,j+1,k+1) +       &
                       gld%pq8(kk) * grd%sal(i+1,j+1,k+1)
      ENDIF
   ENDDO

END SUBROUTINE obs_gld
