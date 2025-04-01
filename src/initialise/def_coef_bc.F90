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
SUBROUTINE def_coef_bc(vect,a1,a2,a3)

   IMPLICIT NONE

   REAL*8   ::   vect(9)
   REAL*8   ::   a1,a2,a3,m_c

!---
! Triggs - Sdika Backward conditions
   m_c=1./(1.+a1-a2+a3)
   m_c=m_c/(1.-a1-a2-a3)
   m_c=m_c/(1.+a2+(a1-a3)*a3)

   vect(1) =m_c*(-a3*a1+1.-a3**2-a2)
   vect(2) =m_c*(a3+a1)*(a2+a3*a1)
   vect(3) =m_c*a3*(a1+a3*a2)
   vect(4) =m_c*(a1+a3*a2)
   vect(5) =m_c*(-(a2-1.)*(a2+a3*a1))
   vect(6) =m_c*(-a3*(a3*a1+a3**2+a2-1))
   vect(7) =m_c*(a3*a1+a2+a1**2-a2**2)
   vect(8) =m_c*((a1*a2)+(a3*a2**2)-(a1*a3**2)-(a3**3)-(a3*a2)+a3)
   vect(9) =m_c*(a3*(a1+a3*a2))

END SUBROUTINE def_coef_bc
