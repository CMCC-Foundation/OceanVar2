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
! Apply observational operator for velocities from drifters (adjoint)  !
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!            Mario Adani     2024  (reproducibility)                   !
!-----------------------------------------------------------------------
SUBROUTINE obs_vdr_ad

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,vdr%no
      IF ( vdr%flc(kk) .EQ. 1 .AND. vdr%par(kk) .EQ. 1 ) THEN
         obs%k = obs%k + 1
         i = vdr%ib(kk)
         j = vdr%jb(kk)
         k = vdr%kb(kk)
#ifdef REPRO
         grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + DBLE(SNGL(vdr%pq1(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + DBLE(SNGL(vdr%pq2(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + DBLE(SNGL(vdr%pq3(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + DBLE(SNGL(vdr%pq4(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i  ,j  ,k+1) = grd%uvl_ad(i  ,j  ,k+1) + DBLE(SNGL(vdr%pq5(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i+1,j  ,k+1) = grd%uvl_ad(i+1,j  ,k+1) + DBLE(SNGL(vdr%pq6(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i  ,j+1,k+1) = grd%uvl_ad(i  ,j+1,k+1) + DBLE(SNGL(vdr%pq7(kk) * obs%gra(obs%k)))
         grd%uvl_ad(i+1,j+1,k+1) = grd%uvl_ad(i+1,j+1,k+1) + DBLE(SNGL(vdr%pq8(kk) * obs%gra(obs%k)))
#else
         grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + vdr%pq1(kk) * obs%gra(obs%k)
         grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + vdr%pq2(kk) * obs%gra(obs%k)
         grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + vdr%pq3(kk) * obs%gra(obs%k)
         grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + vdr%pq4(kk) * obs%gra(obs%k)
         grd%uvl_ad(i  ,j  ,k+1) = grd%uvl_ad(i  ,j  ,k+1) + vdr%pq5(kk) * obs%gra(obs%k)
         grd%uvl_ad(i+1,j  ,k+1) = grd%uvl_ad(i+1,j  ,k+1) + vdr%pq6(kk) * obs%gra(obs%k)
         grd%uvl_ad(i  ,j+1,k+1) = grd%uvl_ad(i  ,j+1,k+1) + vdr%pq7(kk) * obs%gra(obs%k)
         grd%uvl_ad(i+1,j+1,k+1) = grd%uvl_ad(i+1,j+1,k+1) + vdr%pq8(kk) * obs%gra(obs%k)
#endif
      ELSEIF ( vdr%flc(kk) .EQ. 1 .AND. vdr%par(kk) .EQ. 2 ) THEN
         obs%k = obs%k + 1
         i = vdr%ib(kk)
         j = vdr%jb(kk)
         k = vdr%kb(kk)
#ifdef REPRO
         grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + DBLE(SNGL(vdr%pq1(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + DBLE(SNGL(vdr%pq2(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + DBLE(SNGL(vdr%pq3(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + DBLE(SNGL(vdr%pq4(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i  ,j  ,k+1) = grd%vvl_ad(i  ,j  ,k+1) + DBLE(SNGL(vdr%pq5(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i+1,j  ,k+1) = grd%vvl_ad(i+1,j  ,k+1) + DBLE(SNGL(vdr%pq6(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i  ,j+1,k+1) = grd%vvl_ad(i  ,j+1,k+1) + DBLE(SNGL(vdr%pq7(kk) * obs%gra(obs%k)))
         grd%vvl_ad(i+1,j+1,k+1) = grd%vvl_ad(i+1,j+1,k+1) + DBLE(SNGL(vdr%pq8(kk) * obs%gra(obs%k)))
#else
         grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + vdr%pq1(kk) * obs%gra(obs%k)
         grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + vdr%pq2(kk) * obs%gra(obs%k)
         grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + vdr%pq3(kk) * obs%gra(obs%k)
         grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + vdr%pq4(kk) * obs%gra(obs%k)
         grd%vvl_ad(i  ,j  ,k+1) = grd%vvl_ad(i  ,j  ,k+1) + vdr%pq5(kk) * obs%gra(obs%k)
         grd%vvl_ad(i+1,j  ,k+1) = grd%vvl_ad(i+1,j  ,k+1) + vdr%pq6(kk) * obs%gra(obs%k)
         grd%vvl_ad(i  ,j+1,k+1) = grd%vvl_ad(i  ,j+1,k+1) + vdr%pq7(kk) * obs%gra(obs%k)
         grd%vvl_ad(i+1,j+1,k+1) = grd%vvl_ad(i+1,j+1,k+1) + vdr%pq8(kk) * obs%gra(obs%k)
#endif
      ENDIF
   ENDDO

END SUBROUTINE obs_vdr_ad
