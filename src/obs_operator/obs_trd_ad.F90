subroutine obs_trd_ad

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna, and                    !
!                   Vincent Taillandier, Locean, Paris                     !
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
! Call drifter trajectory model (adjoint)                              !
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

 if(trd%ncc.gt.0)then

  if(mpi%myrank.eq.0) then

   do k=1,trd%no

    if(trd%flc(k).eq.1 )then

       obs%k = obs%k + 1

       trd%xtl_ad(k) = obs%gra(obs%k)

    endif
  
   enddo

   do k=1,trd%no

    if(trd%flc(k).eq.1 )then

       obs%k = obs%k + 1

       trd%ytl_ad(k) = obs%gra(obs%k)

    endif

   enddo

  call mod_trj_ad( trd%im,trd%jm,trd%umn,trd%vmn,trd%dx,trd%dy,trd%flc, trd%fls, &
                   trd%nt,trd%no,trd%xmn,trd%ymn,trd%dtm,                        &
                   trd%uvl_ad,trd%vvl_ad,trd%xtl_ad,trd%ytl_ad )


  endif

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
      call gta_mpi( 0_i4, img, jmg, k, km, grd%uvl_ad, trd%uvl_ad)
      call gta_mpi( 0_i4, img, jmg, k, km, grd%vvl_ad, trd%vvl_ad)
   else
      grd%uvl_ad(:,:,k) = grd%uvl_ad(:,:,k) + trd%uvl_ad(:,:)
      grd%vvl_ad(:,:,k) = grd%vvl_ad(:,:,k) + trd%vvl_ad(:,:)
   endif

 endif


end subroutine obs_trd_ad
