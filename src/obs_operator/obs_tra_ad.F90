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
SUBROUTINE obs_tra_ad

!-----------------------------------------------------------------------
!                                                                      !
!> Call Argo trajectory model (adjoint)                                
!!
!!
!!
!                                                                      !
! Version 1: Vincent Taillandier, Srdjan Dobricic 2007                 !
!                                                                      !
!-----------------------------------------------------------------------
! ADANI: Not sure if it is reproducible in mpi. Still not tested !!!
!        At the moment we'll leave it as it is
!-----------------------------------------------------------------------

   USE set_knd
   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, img, jmg, k, km

   IF ( tra%ncc .GT. 0 ) THEN

      IF (mpi%myrank.EQ.0) THEN

         DO k=1,tra%no
            IF ( tra%flc(k) .EQ. 1 ) THEN
               obs%k = obs%k + 1
               tra%xtl_ad(k) = obs%gra(obs%k)
            ENDIF
         ENDDO

         DO k = 1,tra%no
            IF ( tra%flc(k) .EQ. 1 ) THEN
               obs%k = obs%k + 1
               tra%ytl_ad(k) = obs%gra(obs%k)
            ENDIF
         ENDDO

         tra%uvl_ad(:,:) = 0.0_r8
         tra%vvl_ad(:,:) = 0.0_r8

         CALL mod_trj_ad( tra%im,tra%jm,tra%umn,tra%vmn,tra%dx,tra%dy,tra%flc, tra%fls, &
                          tra%nt,tra%no,tra%xmn,tra%ymn,tra%dtm,                        &
                          tra%uvl_ad,tra%vvl_ad,tra%xtl_ad,tra%ytl_ad )

      ENDIF

      IF ( mpi%myrank .EQ. 0 ) THEN
         img = grd%img
         jmg = grd%jmg
      ELSE
         img = 1
         jmg = 1
      ENDIF

      k   = tra%lev
      km  = grd%km

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL gta_mpi( 0_i4, img, jmg, k, km, grd%uvl_ad, tra%uvl_ad)
         CALL gta_mpi( 0_i4, img, jmg, k, km, grd%vvl_ad, tra%vvl_ad)
      ELSE
         grd%uvl_ad(:,:,k) = grd%uvl_ad(:,:,k) + tra%uvl_ad(:,:)
         grd%vvl_ad(:,:,k) = grd%vvl_ad(:,:,k) + tra%vvl_ad(:,:)
      ENDIF

   ENDIF

END SUBROUTINE obs_tra_ad
