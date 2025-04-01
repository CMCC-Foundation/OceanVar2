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
!> Apply observational operator for SST (adjoint)                       
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!            Mario Adani     2024  (reproducibility)                   !
!-----------------------------------------------------------------------
SUBROUTINE obs_sst_ad

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, kk, k

   DO kk = 1,sst%no
      IF ( sst%flc(kk) .EQ. 1 ) THEN
         obs%k = obs%k + 1
         i = sst%ib(kk)
         j = sst%jb(kk)
         k = sst%kb(kk)
#ifdef REPRO
         grd%tem_ad(i  ,j ,k ) = grd%tem_ad(i  ,j ,k ) + DBLE(SNGL(sst%pq1(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j ,k ) = grd%tem_ad(i+1,j ,k ) + DBLE(SNGL(sst%pq2(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j+1,k) = grd%tem_ad(i  ,j+1,k) + DBLE(SNGL(sst%pq3(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j+1,k) = grd%tem_ad(i+1,j+1,k) + DBLE(SNGL(sst%pq4(kk) * obs%gra(obs%k)))
#else
         grd%tem_ad(i  ,j ,k ) = grd%tem_ad(i  ,j ,k ) + sst%pq1(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j ,k ) = grd%tem_ad(i+1,j ,k ) + sst%pq2(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j+1,k) = grd%tem_ad(i  ,j+1,k) + sst%pq3(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j+1,k) = grd%tem_ad(i+1,j+1,k) + sst%pq4(kk) * obs%gra(obs%k)
#endif
      ENDIF
   ENDDO

END SUBROUTINE obs_sst_ad
