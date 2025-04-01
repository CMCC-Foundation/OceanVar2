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
!> Apply observational operator for velocities from drifters           
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_vdr

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,vdr%no
      IF ( (vdr%flc(kk) .EQ. 1 .OR. vdr%fls(kk) .EQ. 1) .AND. vdr%par(kk) .EQ. 1 ) THEN
         i = vdr%ib(kk)
         j = vdr%jb(kk)
         k = vdr%kb(kk)
         vdr%inc(kk) = vdr%pq1(kk) * grd%uvl(i  ,j  ,k  ) +       &
                       vdr%pq2(kk) * grd%uvl(i+1,j  ,k  ) +       &
                       vdr%pq3(kk) * grd%uvl(i  ,j+1,k  ) +       &
                       vdr%pq4(kk) * grd%uvl(i+1,j+1,k  ) +       &
                       vdr%pq5(kk) * grd%uvl(i  ,j  ,k+1) +       &
                       vdr%pq6(kk) * grd%uvl(i+1,j  ,k+1) +       &
                       vdr%pq7(kk) * grd%uvl(i  ,j+1,k+1) +       &
                       vdr%pq8(kk) * grd%uvl(i+1,j+1,k+1)
      ELSEIF ( (vdr%flc(kk) .EQ. 1 .OR. vdr%fls(kk) .EQ. 1) .AND. vdr%par(kk) .EQ. 2 ) THEN
         i = vdr%ib(kk)
         j = vdr%jb(kk)
         k = vdr%kb(kk)
         vdr%inc(kk) = vdr%pq1(kk) * grd%vvl(i  ,j  ,k  ) +       &
                       vdr%pq2(kk) * grd%vvl(i+1,j  ,k  ) +       &
                       vdr%pq3(kk) * grd%vvl(i  ,j+1,k  ) +       &
                       vdr%pq4(kk) * grd%vvl(i+1,j+1,k  ) +       &
                       vdr%pq5(kk) * grd%vvl(i  ,j  ,k+1) +       &
                       vdr%pq6(kk) * grd%vvl(i+1,j  ,k+1) +       &
                       vdr%pq7(kk) * grd%vvl(i  ,j+1,k+1) +       &
                       vdr%pq8(kk) * grd%vvl(i+1,j+1,k+1)
      ENDIF
   ENDDO

END SUBROUTINE obs_vdr
