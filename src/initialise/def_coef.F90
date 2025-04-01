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
!> Compute Recursive Gaussian Filter Coefficient 
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE def_coef(sigma,b1,b2,b3,b)

   IMPLICIT NONE

   REAL*8  sigma,b0,b1,b2,b3,b,q

   IF (sigma .GE.  2.5) THEN
      q = 0.98711 * sigma - 0.96330
   ELSE
      q = 3.97156 - 4.14554 * sqrt(1 - 0.26891 * sigma);
   ENDIF

   b0 = 1.57825 + (2.44413 * q) + (1.4281 * q**2) + (0.422205 * q**3)
   b1 =( (2.44413 * q)  + (2.85619 * q**2) + (1.26661 * q**3))/b0
   b2 =( -(1.4281 * q**2) - (1.26661 * q**3))/b0
   b3 =( 0.422205 * q**3)/b0
   b  = 1 -  ( b1 + b2 + b3 );

END SUBROUTINE def_coef




