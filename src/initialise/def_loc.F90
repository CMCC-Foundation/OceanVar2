subroutine def_loc

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
! READ parameters of the MFS_16_72 grid                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use set_knd
  use cns_str
  use grd_str
  use obs_str
  use mpi_str
  use drv_str
  use netcdf

  implicit none

  integer               :: ierr, i, j, k, iter, ii, jj, iext
  integer(i4)           :: img, jmg
  real(r4), allocatable :: loct(:,:)
  real(r4), allocatable :: locs(:,:)
  real(r4), allocatable ::  lon(:,:)
  real(r4), allocatable ::  lat(:,:)
  real(r4), allocatable ::  msk(:,:)
  real(r8)              :: dxx, dyy, dst


 if(rcf%loc.gt.rcf%L)then

   if(mpi%myrank.eq.0) then
    img = grd%img
    jmg = grd%jmg
   else
    img = 1
    jmg = 1
   endif


  allocate ( loct(img,jmg) )
  allocate ( locs(img,jmg) )
  allocate (  lon(img,jmg) )
  allocate (  lat(img,jmg) )
  allocate (  msk(img,jmg) )


! --------
! Ajoint of bservational operators
! ---

     grd%eta_ad(:,:  ) = 0.0
     grd%tem_ad(:,:,:) = 0.0
     grd%sal_ad(:,:,:) = 0.0
     grd%uvl_ad(:,:,:) = 0.0
     grd%vvl_ad(:,:,:) = 0.0
     obs%gra(:) = obs%res(:)

     call obsop_ad

!--

     grd%loc(:,:) = 0.0

     do k=1,grd%km
      do j=1,grd%jm
       do i=1,grd%im
        if(grd%tem_ad(i,j,k).ne.0.0) grd%loc(i,j) = 1.0
        if(grd%sal_ad(i,j,k).ne.0.0) grd%loc(i,j) = 1.0
        if(grd%uvl_ad(i,j,k).ne.0.0) grd%loc(i,j) = 1.0
        if(grd%vvl_ad(i,j,k).ne.0.0) grd%loc(i,j) = 1.0
       enddo
      enddo
     enddo
      do j=1,grd%jm
       do i=1,grd%im
        if(grd%eta_ad(i,j).ne.0.0) grd%loc(i,j) = 1.0
       enddo
      enddo

   if(mpi%nproc.gt.1)then
      call gth_mpi( img, jmg, 1_i4, 1_i4, grd%loc, loct)
      call gth_mpi( img, jmg, 1_i4, 1_i4, grd%lon,  lon)
      call gth_mpi( img, jmg, 1_i4, 1_i4, grd%lat,  lat)
      call gth_mpi( img, jmg, 1_i4, 1_i4,                &
           grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,1),  msk)
   else
      loct(:,:) = grd%loc(:,:)
       lon(:,:) = grd%lon(:,:)
       lat(:,:) = grd%lat(:,:)
       msk(:,:) = grd%msk(:,:,1)
   endif


  if(mpi%myrank.eq.0) then

write(drv%dia,*) 'heavy computation'

     locs(:,:) = loct(:,:)
    do j=1,jmg-1
     do i=1,img-1
      if(locs(i,j).eq.0.0 .and. msk(i,j).eq.1.0) then
          dxx = 6371.*3.14/180. * (lon(i,j)-lon(i+1,j  )) * cos(lat(i,j)*3.14/180.)
          dyy = 6371.*3.14/180. * (lat(i,j)-lat(i+1,j  ))
          dst = sqrt(dxx**2+dyy**2)
          dxx = 6371.*3.14/180. * (lon(i,j)-lon(i  ,j+1)) * cos(lat(i,j)*3.14/180.)
          dyy = 6371.*3.14/180. * (lat(i,j)-lat(i  ,j+1))
          dst = min(dst,sqrt(dxx**2+dyy**2))
          iext = int(4.*rcf%loc*0.001) / dst
         do jj=max(j-iext,1),min(j+iext,jmg)
          do ii=max(i-iext,1),min(i+iext,img)
           if(loct(ii,jj).eq.1.0)then
              dxx = 6371.*3.14/180. * (lon(i,j)-lon(ii,jj)) * cos(lat(i,j)*3.14/180.)
              dyy = 6371.*3.14/180. * (lat(i,j)-lat(ii,jj))
              dst = sqrt((dxx**2+dyy**2)/(2.*(rcf%loc*0.001)**2))
              locs(i,j) = max(locs(i,j),real(exp(-dst)))
           endif
          enddo
         enddo
      endif
     enddo
    enddo
     loct(:,:) = locs(:,:)

  open(121,form='unformatted',file='locs.dat')
  write(121) loct
  close(121)

   endif

   if(mpi%nproc.gt.1)then
      call gta_mpi( 1_i4, img, jmg, 1_4, 1_i4, grd%loc, loct)
   else
      grd%loc(:,:) = loct(:,:)
   endif

    deallocate ( loct )
    deallocate ( locs )
    deallocate ( lon  )
    deallocate ( lat  )
    deallocate ( msk  )

    do j=1-grd%jas,grd%jm+grd%jae
     do i=1-grd%ias,grd%im+grd%iae
      grd%loc(i,j) = min(grd%loc(i,j),1.0)
     enddo
    enddo

 else

      grd%loc(:,:) = 1.0

 endif

end subroutine def_loc
