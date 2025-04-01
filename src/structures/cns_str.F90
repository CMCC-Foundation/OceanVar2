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
!> Structure of constants 
!
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Mario Adani 2024                                          !
!-----------------------------------------------------------------------
MODULE cns_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

!---
!> recursive filter
   TYPE rcf_t

      INTEGER(i4)          ::  ntr     !< No. of iterations (half of)
      REAL(r8)             ::  dx      !< Grid resolution (m)
      REAL(r8)             ::  L       !< Correlation radius
      REAL(r8)             ::  loc     !< Localization (=1 localize)
      REAL(r8)             ::  E       !< Norm
      REAL(r8)             ::  alp     !< Filter weight
      INTEGER(i4)          ::  ntb     !< Number of points in the table
      REAL(r8)             ::  dsmn    !< Minimum distance
      REAL(r8)             ::  dsmx    !< Maximum distance
      REAL(r8)             ::  dsl     !< Table increment
      REAL(r8), POINTER    ::  al(:)   !< Filter weights in the table
      REAL(r8), POINTER    ::  sc(:)   !< Filter scaling factors in the table
      REAL(r8)             ::  scl     !< Scaling factor
      INTEGER(i4)          ::  kstp    !< Step for extended points

   END TYPE rcf_t

   TYPE (rcf_t)              :: rcf    !< initialize derived type

!---
!> physical constant
   TYPE phy_t

      REAL(r8)             :: g       = 9.80665_r8                    !< Graviational acceleration
      REAL(r8)             :: rho0    = 1025._r8                      !< Reference density
      REAL(r8)             :: re      = 6371.229_r8                   !< Earth Radius 
      REAL(r8)             :: pi      = 3.141592653589793_r8          !< Pi
      REAL(r8)             :: d2r     = 3.141592653589793_r8/180._r8  !< Degree to Radians
      REAL(r8)             :: rsmall  = EPSILON( 1.0_r8 )             !< smallest real computer value
      REAL(r8)             :: omega2  = 0.00014584_r8                 !< 2*rotation rate of the Earth 

   END TYPE phy_t

   TYPE (phy_t)              :: phy !<initialize derived type

END MODULE cns_str
