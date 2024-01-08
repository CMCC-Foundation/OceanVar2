subroutine get_obs_arg

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
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
! Load ARGO observations                                               !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str
 use mpi_str

 implicit none

  INTEGER(i4)   ::  k, ierr
  INTEGER(i4)   ::  i1, kk, i
  REAL          ::  cnti


    arg%no = 0
    arg%nc = 0

   
  open(511,file='arg_mis.dat',form='unformatted',status='old',err=1111)

! ---
! Allocate memory for observations 

   read(511) arg%no

   write(drv%dia,*)' --- No of ARGO obs: ',arg%no

   if(arg%no.eq.0)then
      close(511)
      return
   endif

   ALLOCATE ( arg%ino(arg%no), arg%flg(arg%no), arg%flc(arg%no), arg%par(arg%no))
   ALLOCATE ( arg%lon(arg%no), arg%lat(arg%no), arg%dpt(arg%no), arg%tim(arg%no))
   ALLOCATE ( arg%val(arg%no), arg%bac(arg%no), arg%inc(arg%no))
   ALLOCATE ( arg%bia(arg%no), arg%err(arg%no))
   ALLOCATE ( arg%res(arg%no), arg%b_a(arg%no))
   ALLOCATE ( arg%ib(arg%no), arg%jb(arg%no), arg%kb(arg%no))
   ALLOCATE ( arg%pb(arg%no), arg%qb(arg%no), arg%rb(arg%no))
   ALLOCATE ( arg%pq1(arg%no), arg%pq2(arg%no), arg%pq3(arg%no), arg%pq4(arg%no))
   ALLOCATE ( arg%pq5(arg%no), arg%pq6(arg%no), arg%pq7(arg%no), arg%pq8(arg%no))
   ALLOCATE ( arg%rss(arg%no), arg%ins(arg%no))
   ALLOCATE ( arg%fls(arg%no))



   arg%rss(:) = 0.0
   arg%ins(:) = 0.0
   arg%fls(:) = 0.0
   arg%bia(:) = 0.0

       read (511)                                               &
        arg%ino(1:arg%no), arg%flg(1:arg%no), arg%par(1:arg%no) &
       ,arg%lon(1:arg%no), arg%lat(1:arg%no)                    &
       ,arg%dpt(1:arg%no), arg%tim(1:arg%no)                    &
       ,arg%val(1:arg%no), arg%bac(1:arg%no)                    &
       ,arg%err(1:arg%no), arg%res(1:arg%no)                    &
       ,arg%ib(1:arg%no), arg%jb(1:arg%no), arg%kb(1:arg%no)    &
       ,arg%pb(1:arg%no), arg%qb(1:arg%no), arg%rb(1:arg%no)
    close(511)

    arg%nc = 0
  do k=1,arg%no
   if(arg%flg(k).eq.1 )arg%nc = arg%nc+1
  enddo
!  do k=1,arg%no
!      write(*,*)                                               &
!        arg%ino(k), arg%flg(k), arg%par(k) &
!       ,arg%lon(k), arg%lat(k)                    &
!       ,arg%dpt(k), arg%tim(k)                    &
!       ,arg%val(k), arg%bac(k)                    &
!       ,arg%err(k), arg%res(k)                    &
!       ,arg%ib(k), arg%jb(k), arg%kb(k)    &
!       ,arg%pb(k), arg%qb(k), arg%rb(k)
!  enddo

! ---
! Initialise quality flag
   if(obs%arg.eq.0) arg%flg(:) = -1
   if(obs%arg.eq.0) write(drv%dia,*)'Bad quality flag ',obs%arg

! ---
! Vertical interpolation parameters
    do k = 1,arg%no
     if(arg%flg(k).eq.1)then
       arg%kb(k) = grd%km-1
     do kk = 1,grd%km-1
      if( arg%dpt(k).ge.grd%dep(kk) .and. arg%dpt(k).lt.grd%dep(kk+1) ) then
       arg%kb(k) = kk
       arg%rb(k) = (arg%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
      endif
     enddo
     endif
    enddo


! residual check
  do k=1,arg%no
   if(arg%par(k).eq.1 .and. abs(arg%res(k)).gt.5.0) arg%flg(k) = 0
   if(arg%par(k).eq.2 .and. abs(arg%res(k)).gt.2.0) arg%flg(k) = 0
!   if(arg%par(k).eq.1 .and. abs(arg%res(k)).gt.5.0) write(drv%dia,*)'ARGOn ',k, arg%res(k),'Temp'
!   if(arg%par(k).eq.2 .and. abs(arg%res(k)).gt.2.0) write(drv%dia,*)'ARGOn ',k, arg%res(k),'Salt'
!   if(arg%par(k).eq.1) arg%err(k) = 0.1
!   if(arg%par(k).eq.2) arg%err(k) = 0.1
!   if(arg%tim(k).gt.0.4) arg%flg(k) = 0
!!   if(arg%lat(k).gt.44) arg%flg(k) = 0
  enddo

    go to 1000
! ---
! Thin observations
    do k = 1,arg%no-1
     if(arg%flg(k).eq.1)then
         kk = k + 1
         cnti = 1.
       do while(kk.le.arg%no .and. arg%kb(k).eq.arg%kb(min(kk,arg%no)) .and.  arg%flg(min(kk,arg%no)).eq.1)
         arg%val(k) = arg%val(k) + arg%val(kk)
         arg%bac(k) = arg%bac(k) + arg%bac(kk)
         arg%res(k) = arg%res(k) + arg%res(kk)
         arg%dpt(k) = arg%dpt(k) + arg%dpt(kk)
         arg%flg(kk) = 0
         cnti = cnti + 1.
         kk = kk + 1
       enddo
         arg%val(k) = arg%val(k)/cnti
         arg%bac(k) = arg%bac(k)/cnti
         arg%res(k) = arg%res(k)/cnti
         arg%dpt(k) = arg%dpt(k)/cnti
         arg%rb(k) = (arg%dpt(k) - grd%dep(arg%kb(k))) / (grd%dep(arg%kb(k)+1) - grd%dep(arg%kb(k)))
     endif
    enddo

1000 continue

! ---
! Count good observations
    arg%nc = 0
  do k=1,arg%no
   if(arg%flg(k).eq.1)then
    arg%nc = arg%nc + 1
   else
    arg%bia(k) = 0.
    arg%res(k) = 0.
    arg%inc(k) = 0.
    arg%b_a(k) = 0.
    arg%pq1(k) = 0.
    arg%pq2(k) = 0.
    arg%pq3(k) = 0.
    arg%pq4(k) = 0.
    arg%pq5(k) = 0.
    arg%pq6(k) = 0.
    arg%pq7(k) = 0.
    arg%pq8(k) = 0.
   endif
  enddo

  arg%flc(:) = arg%flg(:)


1111 continue

end subroutine get_obs_arg



subroutine int_par_arg

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
 use mpi_str

 implicit none

  INTEGER(i4)   ::  i, k,kk, ierr
  INTEGER(i4)   ::  i1, j1, k1, idep
  INTEGER(i8)   ::  klev
  REAL(r8)      ::  p1, q1, r1
  REAL(r8)      ::  msk4, div_x, div_y, rmn


  rmn = 1.e-6

 if(arg%no.gt.0) then

   arg%flc(:) = arg%flg(:)

! ---
! Adjust longitudes
      if(grd%bwst.gt.180.)then
        do k=1,arg%no
         if( arg%lon(k).lt.0.0)then
            arg%lon(k) = arg%lon(k) + 360.
         endif
        enddo
      endif

! ---
! Horizontal interpolation parameters

  call int_obs_hor ( arg%no, arg%lat, arg%lon, arg%flc, arg%ib, arg%jb, arg%pb, arg%qb)

! ---
! Undefine masked for multigrid
    do k = 1,arg%no
     if(arg%flc(k).eq.1)then
      i1 = arg%ib(k)
      j1 = arg%jb(k)
      idep = arg%kb(k)+1
      msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
      if(msk4.lt.1.) arg%flc(k) = 0
     endif
    enddo

! ---
! Horizontal interpolation parameters for each masked grid
       do k = 1,arg%no
        if(arg%flc(k) .eq. 1) then

          klev = arg%kb(k)
          call int_obs_pq( arg%ib(k), arg%jb(k), klev, arg%pb(k), arg%qb(k),  &
                           arg%pq1(k), arg%pq2(k), arg%pq3(k), arg%pq4(k))
          klev = arg%kb(k) + 1
          call int_obs_pq( arg%ib(k), arg%jb(k), klev, arg%pb(k), arg%qb(k),  &
                           arg%pq5(k), arg%pq6(k), arg%pq7(k), arg%pq8(k))

         r1=arg%rb(k)
          arg%pq1(k) = (1.-r1) * arg%pq1(k)
          arg%pq2(k) = (1.-r1) * arg%pq2(k)
          arg%pq3(k) = (1.-r1) * arg%pq3(k)
          arg%pq4(k) = (1.-r1) * arg%pq4(k)
          arg%pq5(k) =     r1  * arg%pq5(k)
          arg%pq6(k) =     r1  * arg%pq6(k)
          arg%pq7(k) =     r1  * arg%pq7(k)
          arg%pq8(k) =     r1  * arg%pq8(k)

        endif
       enddo

! ---
! Count good observations
    arg%nc = 0
  do k=1,arg%no
   if(arg%flc(k).eq.1)then
    arg%nc = arg%nc + 1
   endif
  enddo

  arg%inc(:) = 0.0

 endif

end subroutine int_par_arg
