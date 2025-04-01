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
!> Structure of cost function, control vector and optimisation arrays
!                                                                      !
! Version 1: Srdjan Dobricic   2006                                    !
!            Francesco Carere  2024                                    !
!-----------------------------------------------------------------------
MODULE ctl_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Structure for lbfgs
   TYPE lbfgs_t

      INTEGER(i4)                     ::  n               !< size of the optimisation vector
      INTEGER(i4)                     ::  m               !< number of copies to be saved
      CHARACTER*60                    ::  task, csave     !< working string
      LOGICAL, DIMENSION(4)           ::  lsave           !< working array
      INTEGER(i4), DIMENSION(44)      ::  isave           !< working array
      INTEGER(i4), POINTER            ::  nbd(:), iwa(:)  !< working array
      INTEGER(i4)                     ::  iprint          !< working array
      REAL(r8)                        ::  f_b             !< The background cost function
      REAL(r8)                        ::  f_o             !< The observational cost function
      DOUBLE PRECISION                ::  f_c, factr      !< The cost function, accuracy
      DOUBLE PRECISION                ::  pgtol, pgper    !< Stopping criteria, percentage of initial gradient
      DOUBLE PRECISION, DIMENSION(29) :: dsave            !< working array
      DOUBLE PRECISION, POINTER       ::  x_c(:)          !< The control vector (background - analyses)
      DOUBLE PRECISION, POINTER       ::  g_c(:)          !< The gradient of f_c
      DOUBLE PRECISION, POINTER       ::  l_c(:), u_c(:)  !< lower and upper bound of x_c
      DOUBLE PRECISION, POINTER       ::  wa(:), sg(:), sgo(:), yg(:), ygo(:),      &
                                          ws(:,:), wy(:,:), sy(:,:), ss(:,:),       &
                                          yy(:,:), wt(:,:), wn(:,:), snd(:,:),      &
                                          z_c(:), r_c(:), d_c(:), t_c(:)        !< working arrays

   END TYPE lbfgs_t

   TYPE (lbfgs_t)                     :: ctl         !< initialize derived type

! ---
!> Structure for lbfgs global version 
   TYPE lbfgs_tg 

      INTEGER(i4)                     ::  n                      !< size of the optimisation vector
      INTEGER(i4), ALLOCATABLE        ::  n_g(:), n_cum(:)       !
      DOUBLE PRECISION,  ALLOCATABLE  ::  x_c(:), g_c(:)         ! The control vector (background - analyses),The gradient of f_c
      DOUBLE PRECISION,  POINTER      ::  l_c(:),  u_c(:)        !< lower and upper bound of x_c
      INTEGER,           POINTER      ::  nbd(:), iwa(:)         !< working array
      DOUBLE PRECISION,  POINTER      ::  ws(:,:), wy(:,:),  &
                                          z_c(:), r_c(:), d_c(:), t_c(:)        !< working arrays
   END TYPE lbfgs_tg

   TYPE(lbfgs_tg)                 :: ctl_glob       !<initialize derived type

END MODULE ctl_str
