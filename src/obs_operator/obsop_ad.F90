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
!> Apply observational operators - adjoint                              
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE obsop_ad

   USE set_knd
   USE obs_str
   USE grd_str
   USE mpi_str

   IMPLICIT NONE

! ---
   obs%k = 0

! ---
! Satellite observations of SLA
   CALL obs_sla_ad

! ---
! ARGO observations
   CALL obs_arg_ad

! ---
! XBT observations
   CALL obs_xbt_ad

! ---
! Glider observations
   CALL obs_gld_ad

! ---
! Observations of Argo trajectories
   IF (tra%no.GT.0) CALL obs_tra_ad

! ---
! Observations of trajectories of surface drIFters
   IF (trd%no.GT.0) CALL obs_trd_ad

! ---
! Observations of velocity from drIFters
   CALL obs_vdr_ad

! ---
! Observations of velocity from gliders
   CALL obs_gvl_ad

! ---
! Observations of SST
   CALL obs_sst_ad

! ---
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae ,   1_i4, grd%eta_ad)
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%tem_ad)
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%sal_ad)
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%uvl_ad)
      CALL exo_mpi( 1_i4, 0_i4,                                         &
                    1_i4, grd%im+1, 1_i4,                               &
                    1_i4, grd%jm+1, 1_i4,                               &
                    1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , grd%km, grd%vvl_ad)
   ENDIF

END SUBROUTINE obsop_ad
