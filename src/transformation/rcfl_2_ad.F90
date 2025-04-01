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
!>  Gaussian Filter adjoint                                             
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE rcfl_2_ad(jpj,G,b1,b2,b3,b)

   INTEGER       c,jpj
   REAL*8        W(jpj)
   REAL*8        G(jpj)
   REAL*8        b1,b2,b3,b
   REAL*8        A(jpj)

   W(:)=0
   A(:)=0
   !---
   ! Forward  Differences Adjoint 
   DO  c =1,jpj-3
      G(c+1)=G(c+1)+b1*G(c)
      G(c+2)=G(c+2)+b2*G(c)
      G(c+3)=G(c+3)+b3*G(c)
      W(c  )=W(c  )+b *G(c)
   ENDDO

   W(jpj-2) = W(jpj-2) + b*G(jpj-2)
   G(jpj-1) = G(jpj-1) + b1*G(jpj-2)
   G(jpj)   = G(jpj)   + b2*G(jpj-2)

   W(jpj-1) = W(jpj-1) + b*G(jpj-1)
   G(jpj  ) = G(jpj  ) + b1*G(jpj-1)

   W(jpj)   = W(jpj)   + b*G(jpj)

   !---
   ! Backward Differences in Adjoint
   G(:) = 0.0d0
   DO c =jpj,4,-1
      G(c  )=G(c  )+b *W(c)
      W(c-1)=W(c-1)+b1*W(c)
      W(c-2)=W(c-2)+b2*W(c)
      W(c-3)=W(c-3)+b3*W(c)
   END DO

   G(3)=G(3)+b *W(3)
   W(2)=W(2)+b1*W(3)
   W(1)=W(1)+b2*W(3)

   G(2)=G(2)+ b*W(2)
   W(1)=W(1)+b1*W(2)

   G(1) =G(1)+b * W(1)

END SUBROUTINE rcfl_2_ad
