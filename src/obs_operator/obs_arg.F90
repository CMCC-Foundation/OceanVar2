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
!> Apply observational operator for ARGO floats                         
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_arg

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,arg%no
      IF ( ( arg%flc(kk) .EQ. 1 .OR. arg%fls(kk) .EQ. 1 ) .AND. arg%par(kk) .EQ. 1 ) THEN
         i = arg%ib(kk)
         j = arg%jb(kk)
         k = arg%kb(kk)
         arg%inc(kk) = arg%pq1(kk) * grd%tem(i  ,j  ,k  ) +       &
                       arg%pq2(kk) * grd%tem(i+1,j  ,k  ) +       &
                       arg%pq3(kk) * grd%tem(i  ,j+1,k  ) +       &
                       arg%pq4(kk) * grd%tem(i+1,j+1,k  ) +       &
                       arg%pq5(kk) * grd%tem(i  ,j  ,k+1) +       &
                       arg%pq6(kk) * grd%tem(i+1,j  ,k+1) +       &
                       arg%pq7(kk) * grd%tem(i  ,j+1,k+1) +       &
                       arg%pq8(kk) * grd%tem(i+1,j+1,k+1)
      ELSEIF ( (arg%flc(kk) .EQ. 1 .OR. arg%fls(kk) .EQ. 1) .AND. arg%par(kk) .EQ. 2 ) THEN
         i = arg%ib(kk)
         j = arg%jb(kk)
         k = arg%kb(kk)
         arg%inc(kk) = arg%pq1(kk) * grd%sal(i  ,j  ,k  ) +       &
                       arg%pq2(kk) * grd%sal(i+1,j  ,k  ) +       &
                       arg%pq3(kk) * grd%sal(i  ,j+1,k  ) +       &
                       arg%pq4(kk) * grd%sal(i+1,j+1,k  ) +       &
                       arg%pq5(kk) * grd%sal(i  ,j  ,k+1) +       &
                       arg%pq6(kk) * grd%sal(i+1,j  ,k+1) +       &
                       arg%pq7(kk) * grd%sal(i  ,j+1,k+1) +       &
                       arg%pq8(kk) * grd%sal(i+1,j+1,k+1)
      ENDIF
   ENDDO

END SUBROUTINE obs_arg
