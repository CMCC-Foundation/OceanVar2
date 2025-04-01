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
!>  Gaussian Filter                                                     
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE rcfl_2(jpj,G,b1,b2,b3,b)

  INTEGER  c,jpj
  REAL*8   W(jpj)
  REAL*8   G(jpj)
  REAL*8   b1,b2,b3,b
 
  ! ---
  ! Forward differences
  W(1) = b * G(1)
  W(2) = b * G(2) + b1 * W(1)
  W(3) = b * G(3) + b1 * W(2) + b2 * W(1)
  DO  c = 4,jpj
     W(c) = b * G(c) +b1 * W(c-1) + b2 * W(c-2) + b3 * W(c-3)
  ENDDO

  ! ---
  ! Backward  differences
  G(jpj)   = b * W(jpj)
  G(jpj-1) = b * W(jpj-1) + b1 * G(jpj)
  G(jpj-2) = b * W(jpj-2) + b1 * G(jpj-1) + b2 * G(jpj)
  DO  c = jpj-3,1,-1
     G(c) = b * W(c) +  b1 * G(c+1) + b2 * G(c+2) + b3 * G(c+3)
  ENDDO

END SUBROUTINE rcfl_2
