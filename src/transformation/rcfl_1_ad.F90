subroutine rcfl_1_ad( jm, alp, bta, fld)

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!--------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!                                                                      !
! Recursive filter in y direction - adjoint
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use cns_str

 implicit none

 INTEGER(i4)    :: jm

 REAL(r8)       :: fld(jm)
 REAL(r8)       :: a(jm), b(jm), c(jm)
 REAL(r8)       :: alp(jm), bta(jm)

 INTEGER(i4)    :: j, ktr


        a(:) = 0.0
        b(:) = 0.0
        c(:) = fld(:)

       do ktr = 1,rcf%ntr

! negative direction 
        b(:) = 0.0

         do j=1,jm-1
          c(j+1) = c(j+1) + bta(j)*c(j) 
          b(j)   = (1.-bta(j))*c(j)
         enddo


         if( ktr.eq.1 )then
           b(jm) = b(jm) + c(jm) / (1.+bta(jm))
         else
           b(jm  ) = b(jm  ) + (1.-bta(jm)) * c(jm) / (1.-bta(jm)**2)**2
           b(jm-1) = b(jm-1) - (1.-bta(jm)) * bta(jm)**3 * c(jm) / (1.-bta(jm)**2)**2
         endif

! positive direction 
        a(:) = 0.0

         do j=jm,2,-1
          b(j-1) = b(j-1) + alp(j)*b(j)
          a(j) = a(j) + (1.-alp(j))*b(j)
         enddo


         if( ktr.eq.1 )then
           a(1) = a(1) + (1.-alp(1)) * b(1)
         elseif( ktr.eq.2 )then
           a(1) = a(1) + b(1) / (1.+alp(1))
         else
           a(1) = a(1) + (1.-alp(1)) * b(1) / (1.-alp(1)**2)**2
           a(2) = a(2) - (1.-alp(1)) * alp(1)**3 * b(1) / (1.-alp(1)**2)**2
         endif


         c(:) = a(:)

       enddo

          fld(:) = c(:) 

end subroutine rcfl_1_ad
