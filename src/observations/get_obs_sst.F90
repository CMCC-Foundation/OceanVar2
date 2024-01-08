subroutine get_obs_sst

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
! Load SST observations                                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           ! 
! Version 2: J.Pistoia 2013                                            !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str
 use mpi_str

 implicit none

  INTEGER       ::  ierr
  INTEGER(i4)   ::  k
  INTEGER(i4)   ::  i1, i, iter
  REAL(r8)      ::  sumt, sumi, timp, dxx, dyy, dsm

  sst%no = 0
  sst%nc = 0

   
  open(511,file='sst_mis.dat',form='unformatted',status='old',err=1111)

 read(511) sst%no
   write(drv%dia,*) 'Number of SST observations: ',  sst%no, obs%sst

   if(sst%no.eq.0)then
      close(511)
      return
   endif

! ---
! Allocate memory for observations 
   allocate ( sst%ino(sst%no),sst%par(sst%no), sst%flg(sst%no),sst%flc(sst%no))
   allocate ( sst%lon(sst%no), sst%lat(sst%no), sst%dpt(sst%no), sst%tim(sst%no))
   allocate ( sst%val(sst%no), sst%bac(sst%no), sst%inc(sst%no))
   allocate ( sst%bia(sst%no), sst%err(sst%no), sst%res(sst%no),sst%b_a(sst%no))
   allocate ( sst%ib(sst%no),sst%pb(sst%no), sst%jb(sst%no), sst%qb(sst%no))
   allocate ( sst%kb(sst%no), sst%rb(sst%no))
   allocate ( sst%pq1(sst%no), sst%pq2(sst%no), sst%pq3(sst%no), sst%pq4(sst%no)) 
   allocate ( sst%fls(sst%no))


! ---
! Initialise quality flag
   sst%flc(:) = 1

     read (511)                                              &
     sst%ino(1:sst%no), sst%flg(1:sst%no), sst%par(1:sst%no) &
    ,sst%lon(1:sst%no), sst%lat(1:sst%no)                    &
    ,sst%dpt(1:sst%no), sst%tim(1:sst%no)                    &
    ,sst%val(1:sst%no), sst%bac(1:sst%no)                    &
    ,sst%err(1:sst%no), sst%res(1:sst%no)                    &
    ,sst%ib(1:sst%no), sst%jb(1:sst%no), sst%kb(1:sst%no)    &
    ,sst%pb(1:sst%no), sst%qb(1:sst%no), sst%rb(1:sst%no)

     close(511)


    sst%kb(:) = 1

!print*, sst%flg(1:sst%no)

   if(obs%sst.eq.0) sst%flg(:) = -1

! ---
! Remove bias along each track and obseravtaions with large residuals

! residual check
  do k=1,sst%no
   if(abs(sst%res(k)).gt.0.4) sst%flg(k) = 0
  enddo


! ---
! Count good observations
    sst%nc = 0
  do k=1,sst%no
   if(sst%flg(k).eq.1)then
    sst%nc = sst%nc + 1
   else
    sst%inc(k) = 0.
    sst%b_a(k) = 0.
    sst%pq1(k) = 0.
    sst%pq2(k) = 0.
    sst%pq3(k) = 0.
    sst%pq4(k) = 0.
   endif
  enddo

  sst%flc(:) = sst%flg(:)


1111 continue


end subroutine get_obs_sst

subroutine int_par_sst

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
  integer(i4)   ::  i1, kk, i, j1, j, idep
  integer(i8)   ::  klev
  real(r8)      ::  p1, q1
  real(r8)      ::  msk4, div_x, div_y, rmn, dst, dstm
  real(r8)      ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
  real(r8)      ::  lat_lt_rot, lon_rb_rot

  rmn = 1.e-6

 if(sst%no.gt.0) then

       sst%flc(:) = sst%flg(:)

       sst%dpt(:) = 0.0

! ---
! Adjust longitudes
      if(grd%bwst.gt.180.)then
        do k=1,sst%no
         if( sst%lon(k).lt.0.0)then
            sst%lon(k) = sst%lon(k) + 360.
         endif
        enddo
      endif

! ---
! Interpolation parameters

  call int_obs_hor ( sst%no, sst%lat, sst%lon, sst%flc, sst%ib, sst%jb, sst%pb, sst%qb)


    do kk = 1,sst%no
     if(sst%flc(kk).eq.1)then
       i1 = sst%ib(kk)
       j1 = sst%jb(kk)
       idep = sst%kb(kk)
       msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
      if(msk4.lt.1.0) sst%flc(kk) = 0
     endif
    enddo



! ---
! Horizontal interpolation parameters for each masked grid
          klev = 1
       do k = 1,sst%no
        if(sst%flc(k) .eq. 1) then
          call int_obs_pq( sst%ib(k), sst%jb(k), klev, sst%pb(k), sst%qb(k),  &
                           sst%pq1(k), sst%pq2(k), sst%pq3(k), sst%pq4(k))
        endif
       enddo

! ---
! Count good observations
    sst%nc = 0
  do k=1,sst%no
   if(sst%flc(k).eq.1)then
    sst%nc = sst%nc + 1
   endif
  enddo

   write(drv%dia,*) 'Number of good SST observations: ',  sst%nc


 endif




end subroutine int_par_sst
