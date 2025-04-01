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
!> Initialize for adjoint calculations                                 
!!
!! It computes the gradient of the observational component
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Andrea Storto   2021                                      !
!            Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE res_inc

   USE set_knd
   USE grd_str
   USE obs_str
   USE bmd_str
   USE ctl_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4) :: jstart, jend

   grd%eta_ad(:,:  ) = 0.0_r8
   grd%tem_ad(:,:,:) = 0.0_r8
   grd%sal_ad(:,:,:) = 0.0_r8
   grd%uvl_ad(:,:,:) = 0.0_r8
   grd%vvl_ad(:,:,:) = 0.0_r8

   grd%b_x(:,:,:) = 0.0_r8
   grd%b_y(:,:,:) = 0.0_r8
   grd%bx(:,:) = 0.0_r8
   grd%by(:,:) = 0.0_r8

! ---
!First compute without considering huber norm
   obs%gra(:) = obs%amo(:) / obs%err(:)

! ---
! If necessary overwrite obs%gra considering huber norm

! SLA observations
   jstart = 1
   jend   = sla%nc
   IF ( huberqc%sla .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! ARGO observations
   jstart = sla%nc + 1
   jend   = sla%nc + arg%nc
   IF ( huberqc%arg .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! XBT observations
   jstart = sla%nc + arg%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc
   IF ( huberqc%xbt ) CALL huber_resinc(jstart,jend)
! Glider observations
   jstart = sla%nc + arg%nc + xbt%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc
   IF ( huberqc%gld .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! Argo trajectory observations
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc
   IF ( huberqc%tra .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! Trajectory observations of surface drIFters
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc
   IF ( huberqc%trd .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! Observations of drIFter velocities
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc
   IF ( huberqc%vdr .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! Observations of glider velocities
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc
   IF ( huberqc%gvl .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

! SST observations
   jstart = sla%nc + arg%nc + xbt%nc + gld%nc +   &
            2 * tra%nc + 2 * trd%nc + vdr%nc + gvl%nc + 1
   jend   = sla%nc + arg%nc + xbt%nc + gld%nc +     &
            2 * tra%nc + 2 * trd%nc  + vdr%nc + gvl%nc + sst%nc
   IF ( huberqc%sst .AND. (ctl%isave(30) .GT. huberqc%iter) ) CALL huber_resinc(jstart,jend)

END SUBROUTINE res_inc
!===========================================
!======================================================================
SUBROUTINE huber_resinc(jstart,jend)

   USE set_knd
   USE obs_str, ONLY : obs

   IMPLICIT NONE

   INTEGER(i4) :: jo, jstart, jend
   REAL(r8)    :: zj

   DO jo  =  jstart,jend

      IF ( obs%amo(jo)  .GT.  ABS(obs%ahub(jo,1)) .AND. obs%amo(jo)  .LE.  ABS(obs%ahub2(jo,1)) ) THEN

         obs%gra(jo) = obs%ahub(jo,1) / obs%err(jo)

      ELSEIF ( obs%amo(jo)  .LT.  -ABS(obs%ahub(jo,2)) .AND. obs%amo(jo)  .GE.  -ABS(obs%ahub2(jo,2)) ) THEN

         obs%gra(jo) = -obs%ahub(jo,2) / obs%err(jo)

      ELSEIF ( obs%amo(jo)  .GT.  ABS(obs%ahub2(jo,1)) ) THEN

         zj          = 2._r8*obs%ahub(jo,1)*DSQRT(obs%ahub2(jo,1))
         obs%gra(jo) = 0.5_r8*(zj/DSQRT(obs%amo(jo))) / obs%err(jo)

      ELSEIF ( obs%amo(jo)  .LT.  -ABS(obs%ahub2(jo,2)) ) THEN

         zj          = -2._r8*obs%ahub(jo,2)*DSQRT(obs%ahub2(jo,2))
         obs%gra(jo) = 0.5_r8*(zj/DSQRT(ABS(obs%amo(jo)))) / obs%err(jo)

      ELSE

         obs%gra(jo) = obs%amo(jo) / obs%err(jo)

      ENDIF

   ENDDO

END SUBROUTINE

