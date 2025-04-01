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
!> Barotropic model operator data structure
!
! Version 1: Srdjan Dobricic 2007                                      !
!-----------------------------------------------------------------------
MODULE bmd_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Variable related to barotropic model
   TYPE bmd_t
      INTEGER(i4)              ::  ncnt         !< Maximum number of iterations in the implicit solver
      REAL(r8)                 ::  ovr          !< Over-relaxation factor
      REAL(r8)                 ::  resem        !< Stopping criteria
      REAL(r8)                 ::  bnm          !< Number of sea points

      REAL(r8)                 ::  dt           !< Time step
      INTEGER(i4)              ::  nstp         !< Number of time steps per day
      REAL(r8)                 ::  ndy          !< Number of simulation days
      REAL(r8)                 ::  ady          !< Number of averaging days
      INTEGER(i4)              ::  nstps        !< Number of time steps of the main loop
      INTEGER(i4)              ::  nstpa        !< Number of time steps for averaging
      REAL(r8)                 ::  alp1         !< Weighting factor in the trapezoidal scheme
      REAL(r8)                 ::  alp2         !< Weighting factor in the trapezoidal scheme
      REAL(r8)                 ::  fc1          !< Friction intensity (1/dt) [1/s] applyed to t-1
      REAL(r8)                 ::  fc2          !< Friction intensity (1/dt) [1/s] applyed to t+1
      REAL(r8)                 ::  df1          !< Friction intensity ( fc1 * (dx*dy)**2 ) [m**2/s]
      REAL(r8)                 ::  df2          !< Friction intensity ( fc2 * (dx*dy)**2 ) [m**2/s]

      INTEGER(i4), POINTER     ::  itr(:)       !< Number of iterations in the solver
      REAL(r8),    POINTER     ::  mst(:,:)     !< Sea-land mask on t points
      REAL(r8),    POINTER     ::  msu(:,:)     !< Sea-land mask on u points
      REAL(r8),    POINTER     ::  msv(:,:)     !< Sea-land mask on v points
      REAL(r8),    POINTER     ::  hgt(:,:)     !< Depth on t points
      REAL(r8),    POINTER     ::  hgu(:,:)     !< Depth on u points
      REAL(r8),    POINTER     ::  hgv(:,:)     !< Depth on v points
      REAL(r8),    POINTER     ::  dxu(:,:)     !< DX on u points
      REAL(r8),    POINTER     ::  dyu(:,:)     !< DY on u points
      REAL(r8),    POINTER     ::  dxv(:,:)     !< DX on v points
      REAL(r8),    POINTER     ::  dyv(:,:)     !< DY on v points
      REAL(r8),    POINTER     ::   a1(:,:)     !< Constant for over-relaxation
      REAL(r8),    POINTER     ::   a2(:,:)     !< Constant for over-relaxation
      REAL(r8),    POINTER     ::   a3(:,:)     !< Constant for over-relaxation
      REAL(r8),    POINTER     ::   a4(:,:)     !< Constant for over-relaxation
      REAL(r8),    POINTER     ::   a0(:,:)     !< Constant for over-relaxation
      REAL(r8),    POINTER     ::  bxby(:,:)    !< Laplacian of buoyancy force
      REAL(r8),    POINTER     ::   rgh(:,:)    !< Right hand side of equation for eta
      REAL(r8),    POINTER     ::   etb(:,:)    !< Eta at t-1
      REAL(r8),    POINTER     ::    ub(:,:)    !< U at t-1
      REAL(r8),    POINTER     ::    vb(:,:)    !< V at t-1
      REAL(r8),    POINTER     ::    un(:,:)    !< U at t
      REAL(r8),    POINTER     ::    vn(:,:)    !< V at t
      REAL(r8),    POINTER     ::   eta(:,:)    !< Eta at t+1
      REAL(r8),    POINTER     ::    ua(:,:)    !< U at t+1
      REAL(r8),    POINTER     ::    va(:,:)    !< V at t+1
      REAL(r8),    POINTER     ::   etm(:,:)    !< Averaged eta
      REAL(r8),    POINTER     ::    um(:,:)    !< Averaged u
      REAL(r8),    POINTER     ::    vm(:,:)    !< Averaged v
      REAL(r8),    POINTER     ::   div(:,:)    !< Divergence at t-1
      REAL(r8),    POINTER     ::    cu(:,:)    !< Coriolis term on u points
      REAL(r8),    POINTER     ::    cv(:,:)    !< Coriolis term on v points
      REAL(r8),    POINTER     ::   dux(:,:)    !< Friction on U in x
      REAL(r8),    POINTER     ::   duy(:,:)    !< Friction on U in y
      REAL(r8),    POINTER     ::   dvx(:,:)    !< Friction on V in x
      REAL(r8),    POINTER     ::   dvy(:,:)    !< Friction on V in y
      REAL(r8),    POINTER     ::   etx(:,:)    !< Free surface gradient at t-1
      REAL(r8),    POINTER     ::   ety(:,:)    !< Free surface gradient at t-1

   END TYPE bmd_t

   TYPE (bmd_t)                 :: bmd    !< initialize derived type

END MODULE bmd_str
