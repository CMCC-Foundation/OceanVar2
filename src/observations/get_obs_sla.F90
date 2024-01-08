subroutine get_obs_sla

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
! Load SLA observations                                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! Modify the track recognition P. Oddo
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, i, iter
  REAL(r8)      ::  sumt, sumi, timp, dxx, dyy, dsm, satino

  sla%no = 0
  sla%nc = 0

   
  open(511,file='sla_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations 

   read(511) sla%no
   write(drv%dia,*) ' --- No of SLA obs: ',  sla%no, obs%sla


   if(sla%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( sla%ino(sla%no), sla%flg(sla%no), sla%flc(sla%no))
   ALLOCATE ( sla%lon(sla%no), sla%lat(sla%no), sla%tim(sla%no))
   ALLOCATE ( sla%val(sla%no), sla%bac(sla%no), sla%inc(sla%no))
   ALLOCATE ( sla%bia(sla%no), sla%err(sla%no))
   ALLOCATE ( sla%res(sla%no), sla%b_a(sla%no))
   ALLOCATE ( sla%ib(sla%no), sla%jb(sla%no))
   ALLOCATE ( sla%pb(sla%no), sla%qb(sla%no))
   ALLOCATE ( sla%pq1(sla%no), sla%pq2(sla%no), sla%pq3(sla%no), sla%pq4(sla%no))
   ALLOCATE ( sla%dpt(sla%no))
   ALLOCATE ( sla%dtm(sla%no))
   ALLOCATE ( sla%rss(sla%no), sla%ins(sla%no))
   ALLOCATE ( sla%fls(sla%no))

   sla%rss(:) = 0.0
   sla%ins(:) = 0.0
   sla%fls(:) = 0.0

! ---
! Initialise quality flag
   sla%flc(:) = 1

! ---
! Level corresponding to the minimum depth
   sla%kdp=grd%km
  do k=grd%km, 1, -1
   if(grd%dep(k).ge.sla%dep) sla%kdp = k
  enddo

       read (511)                                              &
        sla%ino(1:sla%no), sla%flg(1:sla%no)                   &
       ,sla%lon(1:sla%no), sla%lat(1:sla%no), sla%tim(1:sla%no)&
       ,sla%val(1:sla%no), sla%bac(1:sla%no)                   &
       ,sla%err(1:sla%no), sla%res(1:sla%no)                   &
       ,sla%ib(1:sla%no), sla%jb(1:sla%no)                     &
       ,sla%pb(1:sla%no), sla%qb(1:sla%no), sla%dtm(1:sla%no)
    close(511)
!
!  do k=1,sla%no
!       write(*,*)                                              &
!        sla%ino(k), sla%flg(k)                   &
!       ,sla%lon(k), sla%lat(k), sla%tim(k)&
!       ,sla%val(k), sla%bac(k)                   &
!       ,sla%err(k), sla%res(k)                   &
!       ,sla%ib(k), sla%jb(k)                     &
!       ,sla%pb(k), sla%qb(k), sla%dtm(k)
!
!  enddo


   if(obs%sla.eq.0) sla%flg(:) = -1


   ! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flg(k).eq.1)then
    sla%nc = sla%nc + 1
   endif
  enddo
   write(drv%dia,*) 'Number of good SLA observations(before removing observations, sla%flg): ',  sla%nc
  !do k=1,sla%no
  ! if(sla%flg(k).eq.1)then
  ! write(*,*) sla%res(k), 'Paolo1',k
  ! endif
  ! enddo
! ---
! Remove bias along each track and obseravtaions with large residuals

 do iter = 1,3

!bias
   sla%bia(:) = 0.0
   timp = sla%tim(1)
   satino=sla%ino(1)
   dsm = 100.
   i1 = 1
 do k=2,sla%no

      dxx = 6371.*3.14/180. * (sla%lon(k)-sla%lon(k-1)) * cos(sla%lat(k)*3.14/180.)
      dyy = 6371.*3.14/180. * (sla%lat(k)-sla%lat(k-1))

  !if((sla%tim(k).ne.timp .or. sqrt(dxx**2+dyy**2).gt.dsm) .and. k.gt.i1)then
  if((sla%ino(k).ne.satino .or. sqrt(dxx**2+dyy**2).gt.dsm) .and. k.gt.i1)then
     sumt = 0.0
     sumi = 0.0

    do i=i1,k-1
     if(sla%flg(i).eq.1)then
      sumt = sumt + sla%res(i)
      sumi = sumi + 1.0
     endif
    enddo

     !if(sumi.gt.0.) write(*,*)'dentro primo if con',sla%ino(k),satino,sqrt(dxx**2+dyy**2),k,i1,iter
     !if(sumi.gt.0.) write(*,*)'e quindi ottengo sumt=',sumt,k,iter

     if(sumi.gt.0.) sumt = sumt/sumi
     !write(*,*)'che diviso sumt=',sumt,k,iter

    do i=i1,k-1
     sla%res(i) = sla%res(i) - sumt
     sla%bia(i) = sumt
    enddo

     timp = sla%tim(k)
     satino=sla%ino(k)
     i1 = k
  else if(k.eq.sla%no .and. k.ge.i1)then           !-!

     !write(*,*)'Sono in else if con k',k,sla%no,i1
     sumt = 0.0                                      !
     sumi = 0.0                                      !
    do i=i1,k                            !-!         !
     if(sla%flg(i).eq.1)then              !!         !
      sumt = sumt + sla%res(i)             !         !
      sumi = sumi + 1.0                    !         !
     endif                                !!         !
    enddo                                !-!         !
     if(sumi.gt.0.) sumt = sumt/sumi                 !
    do i=i1,k                            !-!         !
     sla%res(i) = sla%res(i) - sumt        !         !
     sla%bia(i) = sumt                     !         !
    enddo                                !-!         !
  endif                                            !-!
 enddo

  !do k=1,sla%no
  ! if(sla%flg(k).eq.1)then
  ! write(*,*) sla%res(k), 'Residuals at iter=',iter,k
  ! endif
  ! enddo

 enddo ! iter

! residual check
  do k=1,sla%no
   if(abs(sla%res(k)).gt.0.3) sla%flg(k) = 0
  enddo


! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flc(k).eq.1) then
    sla%nc = sla%nc + 1
   endif
  enddo

   write(drv%dia,*) 'Number of good SLA observations(after removing observations in sla%flc): ',  sla%nc 

! ---
! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flg(k).eq.1)then
    sla%nc = sla%nc + 1
   else
    sla%flc(k) = 0
!    sla%bia(k) = 0.
!    sla%res(k) = 0.
    sla%inc(k) = 0.
    sla%b_a(k) = 0.
    sla%pq1(k) = 0.
    sla%pq2(k) = 0.
    sla%pq3(k) = 0.
    sla%pq4(k) = 0.
   endif
  enddo

  sla%flc(:) = sla%flg(:)

  write(drv%dia,*) 'Number of good SLA observations(after removing observations in sla%flg): ',  sla%nc
 ! do k=1,sla%no
 !  if(sla%flg(k).eq.1)then
 !  write(*,*) sla%res(k), 'Paolo'
 !  endif
 !  enddo
1111 continue


end subroutine get_obs_sla

subroutine int_par_sla

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str

 implicit none

  integer(i4)   ::  k, ierr
  integer(i4)   ::  i1, kk, i, j1, j
  integer(i8)   ::  klev
  real(r8)      ::  p1, q1
  real(r8)      ::  msk4, div_x, div_y, rmn, dst, dstm
  real(r8)      ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
  real(r8)      ::  lat_lt_rot, lon_rb_rot

  rmn = 1.e-6

 if(sla%no.gt.0) then

       sla%flc(:) = sla%flg(:)

       sla%dpt(:) = 0.0

! ---
! Adjust longitudes
      if(grd%bwst.gt.180.)then
        do k=1,sla%no
         if( sla%lon(k).lt.0.0)then
            sla%lon(k) = sla%lon(k) + 360.
         endif
        enddo
      endif

! ---
! Interpolation parameters

  call int_obs_hor ( sla%no, sla%lat, sla%lon, sla%flc, sla%ib, sla%jb, sla%pb, sla%qb)


    do kk = 1,sla%no
     if(sla%flc(kk).eq.1)then
       i1 = sla%ib(kk)
       j1 = sla%jb(kk)
       sla%dpt(kk) = max(grd%hgt(i1,j1),max(grd%hgt(i1+1,j1),max(grd%hgt(i1,j1+1),grd%hgt(i1+1,j1+1))))
       msk4 = grd%msk(i1,j1,sla%kdp) + grd%msk(i1+1,j1,sla%kdp) + grd%msk(i1,j1+1,sla%kdp) + grd%msk(i1+1,j1+1,sla%kdp)
      if(msk4.lt.4.0) sla%flc(kk) = 0
     endif
    enddo



! ---
! Horizontal interpolation parameters for each masked grid
          klev = 1
       do k = 1,sla%no
        if(sla%flc(k) .eq. 1) then
          call int_obs_pq( sla%ib(k), sla%jb(k), klev, sla%pb(k), sla%qb(k),  &
                           sla%pq1(k), sla%pq2(k), sla%pq3(k), sla%pq4(k))
        endif
       enddo

! ---
! Count good observations
    sla%nc = 0
  do k=1,sla%no
   if(sla%flc(k).eq.1)then
    sla%nc = sla%nc + 1
   endif
  enddo

   write(drv%dia,*) 'Number of good SLA observations: ',  sla%nc

  sla%inc(:) = 0.0

 endif



end subroutine int_par_sla
