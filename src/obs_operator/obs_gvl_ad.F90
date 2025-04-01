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
!> Apply observational operator for velocities from gliders  (adjoint)  
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!            Mario Adani     2024  (reproducibility)                   !
!-----------------------------------------------------------------------
SUBROUTINE obs_gvl_ad

   USE set_knd
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k, kk

   DO kk = 1,gvl%no
      IF ( gvl%flc(kk) .EQ. 1 .AND. gvl%par(kk) .EQ. 1 ) THEN
         obs%k = obs%k + 1
         i = gvl%ib(kk)
         j = gvl%jb(kk)
         DO k = 1,gvl%kb(kk)+1
#ifdef REPRO
            grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + DBLE(SNGL(gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + DBLE(SNGL(gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + DBLE(SNGL(gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + DBLE(SNGL(gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
#else
            grd%uvl_ad(i  ,j  ,k  ) = grd%uvl_ad(i  ,j  ,k  ) + gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%uvl_ad(i+1,j  ,k  ) = grd%uvl_ad(i+1,j  ,k  ) + gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%uvl_ad(i  ,j+1,k  ) = grd%uvl_ad(i  ,j+1,k  ) + gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%uvl_ad(i+1,j+1,k  ) = grd%uvl_ad(i+1,j+1,k  ) + gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
#endif
         ENDDO
      ELSEIF ( gvl%flc(kk) .EQ. 1 .AND. gvl%par(kk) .EQ. 2 ) THEN
         obs%k = obs%k + 1
         i = gvl%ib(kk)
         j = gvl%jb(kk)
         DO k = 1,gvl%kb(kk)+1
#ifdef REPRO
            grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + DBLE(SNGL(gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + DBLE(SNGL(gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + DBLE(SNGL(gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
            grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + DBLE(SNGL(gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)))
#else
            grd%vvl_ad(i  ,j  ,k  ) = grd%vvl_ad(i  ,j  ,k  ) + gvl%pq1(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%vvl_ad(i+1,j  ,k  ) = grd%vvl_ad(i+1,j  ,k  ) + gvl%pq2(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%vvl_ad(i  ,j+1,k  ) = grd%vvl_ad(i  ,j+1,k  ) + gvl%pq3(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
            grd%vvl_ad(i+1,j+1,k  ) = grd%vvl_ad(i+1,j+1,k  ) + gvl%pq4(kk) * gvl%dzr(k,kk) * obs%gra(obs%k)
#endif
         ENDDO
      ENDIF
   ENDDO

END SUBROUTINE obs_gvl_ad
