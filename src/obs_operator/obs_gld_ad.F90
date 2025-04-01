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
!> Apply observational operator for gliders (adjoint)                  
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!            Mario Adani     2024  (reproducibility)                   !
!-----------------------------------------------------------------------
SUBROUTINE obs_gld_ad

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,gld%no
      IF ( gld%flc(kk) .EQ. 1 .AND. gld%par(kk) .EQ. 1 ) THEN
         obs%k = obs%k + 1
         i = gld%ib(kk)
         j = gld%jb(kk)
         k = gld%kb(kk)
#ifdef REPRO
         grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + DBLE(SNGL(gld%pq1(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + DBLE(SNGL(gld%pq2(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + DBLE(SNGL(gld%pq3(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + DBLE(SNGL(gld%pq4(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + DBLE(SNGL(gld%pq5(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + DBLE(SNGL(gld%pq6(kk) * obs%gra(obs%k)))
         grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + DBLE(SNGL(gld%pq7(kk) * obs%gra(obs%k)))
         grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + DBLE(SNGL(gld%pq8(kk) * obs%gra(obs%k)))
#else
         grd%tem_ad(i  ,j  ,k  ) = grd%tem_ad(i  ,j  ,k  ) + gld%pq1(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j  ,k  ) = grd%tem_ad(i+1,j  ,k  ) + gld%pq2(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j+1,k  ) = grd%tem_ad(i  ,j+1,k  ) + gld%pq3(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j+1,k  ) = grd%tem_ad(i+1,j+1,k  ) + gld%pq4(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j  ,k+1) = grd%tem_ad(i  ,j  ,k+1) + gld%pq5(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j  ,k+1) = grd%tem_ad(i+1,j  ,k+1) + gld%pq6(kk) * obs%gra(obs%k)
         grd%tem_ad(i  ,j+1,k+1) = grd%tem_ad(i  ,j+1,k+1) + gld%pq7(kk) * obs%gra(obs%k)
         grd%tem_ad(i+1,j+1,k+1) = grd%tem_ad(i+1,j+1,k+1) + gld%pq8(kk) * obs%gra(obs%k)
#endif
      ELSEIF ( gld%flc(kk) .EQ. 1 .AND. gld%par(kk) .EQ. 2 ) THEN
         obs%k = obs%k + 1
         i = gld%ib(kk)
         j = gld%jb(kk)
         k = gld%kb(kk)
#ifdef REPRO
         grd%sal_ad(i  ,j  ,k  ) = grd%sal_ad(i  ,j  ,k  ) + DBLE(SNGL(gld%pq1(kk) * obs%gra(obs%k)))
         grd%sal_ad(i+1,j  ,k  ) = grd%sal_ad(i+1,j  ,k  ) + DBLE(SNGL(gld%pq2(kk) * obs%gra(obs%k)))
         grd%sal_ad(i  ,j+1,k  ) = grd%sal_ad(i  ,j+1,k  ) + DBLE(SNGL(gld%pq3(kk) * obs%gra(obs%k)))
         grd%sal_ad(i+1,j+1,k  ) = grd%sal_ad(i+1,j+1,k  ) + DBLE(SNGL(gld%pq4(kk) * obs%gra(obs%k)))
         grd%sal_ad(i  ,j  ,k+1) = grd%sal_ad(i  ,j  ,k+1) + DBLE(SNGL(gld%pq5(kk) * obs%gra(obs%k)))
         grd%sal_ad(i+1,j  ,k+1) = grd%sal_ad(i+1,j  ,k+1) + DBLE(SNGL(gld%pq6(kk) * obs%gra(obs%k)))
         grd%sal_ad(i  ,j+1,k+1) = grd%sal_ad(i  ,j+1,k+1) + DBLE(SNGL(gld%pq7(kk) * obs%gra(obs%k)))
         grd%sal_ad(i+1,j+1,k+1) = grd%sal_ad(i+1,j+1,k+1) + DBLE(SNGL(gld%pq8(kk) * obs%gra(obs%k)))
#else
         grd%sal_ad(i  ,j  ,k  ) = grd%sal_ad(i  ,j  ,k  ) + gld%pq1(kk) * obs%gra(obs%k)
         grd%sal_ad(i+1,j  ,k  ) = grd%sal_ad(i+1,j  ,k  ) + gld%pq2(kk) * obs%gra(obs%k)
         grd%sal_ad(i  ,j+1,k  ) = grd%sal_ad(i  ,j+1,k  ) + gld%pq3(kk) * obs%gra(obs%k)
         grd%sal_ad(i+1,j+1,k  ) = grd%sal_ad(i+1,j+1,k  ) + gld%pq4(kk) * obs%gra(obs%k)
         grd%sal_ad(i  ,j  ,k+1) = grd%sal_ad(i  ,j  ,k+1) + gld%pq5(kk) * obs%gra(obs%k)
         grd%sal_ad(i+1,j  ,k+1) = grd%sal_ad(i+1,j  ,k+1) + gld%pq6(kk) * obs%gra(obs%k)
         grd%sal_ad(i  ,j+1,k+1) = grd%sal_ad(i  ,j+1,k+1) + gld%pq7(kk) * obs%gra(obs%k)
         grd%sal_ad(i+1,j+1,k+1) = grd%sal_ad(i+1,j+1,k+1) + gld%pq8(kk) * obs%gra(obs%k)
#endif
      ENDIF
   ENDDO

END SUBROUTINE obs_gld_ad
