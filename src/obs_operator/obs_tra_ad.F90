subroutine obs_tra_ad

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
! Call Argo trajectory model (adjoint)                                 !
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

 if(tra%ncc.gt.0)then

  if(mpi%myrank.eq.0) then

   do k=1,tra%no

    if(tra%flc(k).eq.1 )then

       obs%k = obs%k + 1

       tra%xtl_ad(k) = obs%gra(obs%k)

    endif

   enddo

   do k=1,tra%no

    if(tra%flc(k).eq.1 )then

       obs%k = obs%k + 1

       tra%ytl_ad(k) = obs%gra(obs%k)

    endif

   enddo

    tra%uvl_ad(:,:) = 0.0
    tra%vvl_ad(:,:) = 0.0

    call mod_trj_ad( tra%im,tra%jm,tra%umn,tra%vmn,tra%dx,tra%dy,tra%flc, tra%fls, &
                     tra%nt,tra%no,tra%xmn,tra%ymn,tra%dtm,                        &
                     tra%uvl_ad,tra%vvl_ad,tra%xtl_ad,tra%ytl_ad )


  endif


   if(mpi%myrank.eq.0) then
    img = grd%img
    jmg = grd%jmg
   else
    img = 1
    jmg = 1
   endif

    k   = tra%lev
    km  = grd%km

   if(mpi%nproc.gt.1)then
      call gta_mpi( 0_i4, img, jmg, k, km, grd%uvl_ad, tra%uvl_ad)
      call gta_mpi( 0_i4, img, jmg, k, km, grd%vvl_ad, tra%vvl_ad)
   else
      grd%uvl_ad(:,:,k) = grd%uvl_ad(:,:,k) + tra%uvl_ad(:,:)
      grd%vvl_ad(:,:,k) = grd%vvl_ad(:,:,k) + tra%vvl_ad(:,:)
   endif

 endif


end subroutine obs_tra_ad
