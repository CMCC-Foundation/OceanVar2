subroutine obs_sst_ad

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
! Apply observational operator for SLA (adjoint)                       !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str

 implicit none

 INTEGER(i4)   ::  i, j, kk, k

 do kk=1,sst%no
  if(sst%flc(kk).eq.1)then

    obs%k = obs%k + 1
    i=sst%ib(kk)
    j=sst%jb(kk)
    k=sst%kb(kk)
    grd%tem_ad(i  ,j ,k ) = grd%tem_ad(i  ,j ,k ) + sst%pq1(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j ,k ) = grd%tem_ad(i+1,j ,k ) + sst%pq2(kk) * obs%gra(obs%k)
    grd%tem_ad(i  ,j+1,k) = grd%tem_ad(i  ,j+1,k) + sst%pq3(kk) * obs%gra(obs%k)
    grd%tem_ad(i+1,j+1,k) = grd%tem_ad(i+1,j+1,k) + sst%pq4(kk) * obs%gra(obs%k)
  endif

 enddo


end subroutine obs_sst_ad
