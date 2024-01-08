subroutine rcfl_1( jm, alp, bta, fld)

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
! Recursive filter in y direction
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


          a(:) = fld(:)
          b(:) = 0.0
          c(:) = 0.0

       do ktr = 1,rcf%ntr

! positive direction
         if( ktr.eq.1 )then
           b(1) = (1.-alp(1)) * a(1)
         elseif( ktr.eq.2 )then
            b(1) = a(1) / (1.+alp(1))
         else
            b(1) = (1.-alp(1)) * (a(1)-alp(1)**3 * a(2)) / (1.-alp(1)**2)**2
         endif

        do j=2,jm
             b(j) = alp(j)*b(j-1) + (1.-alp(j))*a(j)
        enddo

! negative direction
         if( ktr.eq.1 )then
           c(jm) = b(jm) / (1.+bta(jm))
         else
           c(jm) = (1.-bta(jm)) * (b(jm)-bta(jm)**3 * b(jm-1)) / (1.-bta(jm)**2)**2
         endif

         do j=jm-1,1,-1
          c(j) = bta(j)*c(j+1) + (1.-bta(j))*b(j)
         enddo

         a(:) = c(:)

       enddo

         fld(:) = a(:) 

end subroutine rcfl_1
