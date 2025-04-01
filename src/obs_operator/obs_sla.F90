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
!> Apply observational operator for SLA                               
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_sla

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k

   DO k=1,sla%no
      IF ( sla%flc(k) .EQ. 1 .OR. sla%fls(k) .EQ. 1 ) THEN
         i = sla%ib(k)
         j = sla%jb(k)
         sla%inc(k) = sla%pq1(k) * grd%eta(i  ,j  ) +       &
                      sla%pq2(k) * grd%eta(i+1,j  ) +       &
                      sla%pq3(k) * grd%eta(i  ,j+1) +       &
                      sla%pq4(k) * grd%eta(i+1,j+1)
      ENDIF
   ENDDO

END SUBROUTINE obs_sla
