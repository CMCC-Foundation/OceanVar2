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
!> Apply observational operator for XBT floats (adjoint)                
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!            Mario Adani     2024  (reproducibility)                   !
!-----------------------------------------------------------------------
SUBROUTINE obs_xbt_ad

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,xbt%no
      IF ( xbt%flc(kk) .EQ. 1 .AND. xbt%par(kk) .EQ. 1 ) THEN
         obs%k = obs%k + 1
         i = xbt%ib(kk)
         j = xbt%jb(kk)
         k = xbt%kb(kk)
#ifdef REPRO
         grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + DBLE(SNGL(xbt%pq1(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + DBLE(SNGL(xbt%pq2(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + DBLE(SNGL(xbt%pq3(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + DBLE(SNGL(xbt%pq4(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + DBLE(SNGL(xbt%pq5(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + DBLE(SNGL(xbt%pq6(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + DBLE(SNGL(xbt%pq7(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + DBLE(SNGL(xbt%pq8(kk) * obs%gra(obs%k)))
#else
         grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + xbt%pq1(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + xbt%pq2(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + xbt%pq3(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + xbt%pq4(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + xbt%pq5(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + xbt%pq6(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + xbt%pq7(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + xbt%pq8(kk) * obs%gra(obs%k)
#endif
      ENDIF
   ENDDO

END SUBROUTINE obs_xbt_ad
