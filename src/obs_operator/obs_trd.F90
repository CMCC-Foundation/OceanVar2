subroutine obs_trd

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
!                                                                          !
!    This file is part of OceanVar.                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Call drifter trajectory model                                        !
!                                                                      !
! Version 1: V. Taillandier, S. Dobricic 2007                          !
!                                                                      !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use obs_str
 use mpi_str

 implicit none

 INTEGER(i4)   ::  i, j, img, jmg, k, km

 if(trd%ncc.gt.0 .or. trd%ncs.gt.0)then

   if(mpi%myrank.eq.0) then
    img = grd%img
    jmg = grd%jmg
   else
    img = 1
    jmg = 1
   endif

    k   = trd%lev
    km  = grd%km

   if(mpi%nproc.gt.1)then
      call gth_mpi( img, jmg, k, km, grd%uvl, trd%uvl)
      call gth_mpi( img, jmg, k, km, grd%vvl, trd%vvl)
   else
      trd%uvl(:,:) = grd%uvl(:,:,k) 
      trd%vvl(:,:) = grd%vvl(:,:,k)
   endif


 if(mpi%myrank.eq.0) then

  call mod_trj_tl( trd%im,trd%jm,trd%umn,trd%vmn,trd%dx,trd%dy,trd%flc, trd%fls, &
                   trd%nt,trd%no,trd%xmn,trd%ymn,trd%dtm,                        &
                   trd%uvl,trd%vvl,trd%xtl,trd%ytl )

   do k=1,trd%no

    if(trd%flc(k).eq.1 .or. trd%fls(k).eq.1)then

      print*,' === >> ',trd%flc(k),trd%fls(k)

      trd%inx(k) = trd%xtl(k)
      trd%iny(k) = trd%ytl(k)

      i=int(trd%xmn(trd%nt+1,k)+trd%xtl(k))
      j=int(trd%ymn(trd%nt+1,k)+trd%ytl(k))
!      trd%loa(k) = trd%lon(i,j) +     &
!                   (trd%xmn(trd%nt+1,k)+trd%xtl(k)-i)*(trd%lon(i+1,j)-trd%lon(i,j))
!      trd%laa(k) = trd%lat(i,j) +     &
!                   (trd%ymn(trd%nt+1,k)+trd%ytl(k)-j)*(trd%lat(i,j+1)-trd%lat(i,j))

   endif

  enddo


  endif


 endif

end subroutine obs_trd
