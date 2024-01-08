subroutine wrt_dia

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
! Write outputs and diagnostics                                        !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! Version 1.1: P. Oddo  2014                                           !
!-----------------------------------------------------------------------


 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str
 use bmd_str
 use mpi_str
 use netcdf

 implicit none

 include "mpif.h"
  INTEGER                    :: status
  INTEGER                    :: ierr
  INTEGER                    :: fh
  INTEGER                    :: filetype
  INTEGER                    :: one
  INTEGER                    :: bffsz
  INTEGER(mpi_offset_kind)   :: zero
  INTEGER, allocatable       ::  map(:)

 INTEGER(i4)                 :: i, j, k, kk
 CHARACTER                   :: fgrd

 INTEGER(i8)               :: maxno
 INTEGER(i8), allocatable  ::  allflg(:)
 REAL(r8), allocatable     ::  allinc(:), allb_a(:)

 INTEGER   :: nobs
 INTEGER   :: stat, ncid
 INTEGER   :: idlon, idlat, iddep
 INTEGER   :: ideta, idtem, idsal, iduvl, idvvl, idmsk

   write(drv%dia,*) ' --------------------------------------------------'
   write(drv%dia,*) ' writes corrections and statistics '
   write(drv%dia,*) ' --------------------------------------------------'

   write(fgrd,'(i1)')drv%ktr

! ---
! Innovations

   call write_inn(   1_i4, grd%eta, 'eta' )
   call write_inn( grd%km, grd%tem, 'tem' )
   call write_inn( grd%km, grd%sal, 'sal' )
   call write_inn( grd%km, grd%uvl, 'uvl' )
   call write_inn( grd%km, grd%vvl, 'vvl' )

! ---
! Eigenvalues

!   open(101,file=drv%flag//drv%date//'eiv_'//fgrd//'.dat',form='unformatted')
!   open(101,file='eiv_'//fgrd//'.dat',form='unformatted')
!    write(101) real(grd%msk(:,:,1))
!   do k=1,ros%neof
!    write(101) grd%ro(:,:,k)
!   enddo
!   close(101)

! ---
! Observations

!  open(215,file=drv%flag//drv%date//'obs_'//fgrd//'.dat',form='unformatted')
!  open(215,file='obs_'//fgrd//'.dat',form='unformatted')

   if(mpi%myrank.eq.0) open(215,file='obs.dat',form='unformatted')

! Oddo diagnostic
   if(sla%no.gt.0)then
   if(mpi%myrank.eq.0) open(216,file='sla_stat.txt',form='formatted')
   endif
   if(arg%no.gt.0)then
   if(mpi%myrank.eq.0) open(217,file='arg_stat.txt',form='formatted')
   endif
   if(xbt%no.gt.0)then
   if(mpi%myrank.eq.0) open(218,file='xbt_stat.txt',form='formatted')
   endif
   if(gld%no.gt.0)then
   if(mpi%myrank.eq.0) open(219,file='gld_stat.txt',form='formatted')
   endif
! end Oddo statistic

   maxno = max(sst%no,max(sla%no,max(arg%no,max(xbt%no,gld%no))))
   maxno = max(gvl%no,max(vdr%no,maxno))

   if(maxno.gt.0) ALLOCATE ( allflg(maxno), allinc(maxno), allb_a(maxno) )

 if(sla%no.gt.0)then
    allflg(:) = 0.0
    allinc(:) = 0.0
    allb_a(:) = 0.0
    nobs = sla%no
   call mpi_reduce(  sla%flg, allflg, nobs, mpi%i8  ,   &
                     mpi_max, 0, mpi%comm, ierr)
   call mpi_reduce(  sla%inc, allinc, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
   call mpi_reduce(  sla%b_a, allb_a, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
   write(215) sla%no
   if(sla%no.gt.0) then
    sla%flg(1:sla%no) = allflg(1:sla%no)
    sla%inc(1:sla%no) = allinc(1:sla%no)
    sla%b_a(1:sla%no) = allb_a(1:sla%no)
         write(215)                                 &
        sla%ino(1:sla%no), sla%flg(1:sla%no),sla%lon(1:sla%no) &
       ,sla%lat(1:sla%no), sla%tim(1:sla%no),sla%val(1:sla%no) &
       ,sla%bac(1:sla%no), sla%err(1:sla%no),sla%res(1:sla%no) &
       ,sla%bia(1:sla%no), sla%inc(1:sla%no),sla%b_a(1:sla%no) &
       ,sla%dpt(1:sla%no), sla%dtm(1:sla%no)
       write(216,"(i8)") sla%no
       do kk=1,sla%no
       write (216,"(2I8,12F10.5)")                 &
        sla%ino(kk),sla%flg(kk),sla%lon(kk),sla%lat(kk),sla%tim(kk) &
       ,sla%val(kk),sla%bac(kk),sla%err(kk),sla%res(kk),sla%bia(kk) &
       ,sla%inc(kk),sla%b_a(kk),sla%dpt(kk),sla%dtm(kk)
       enddo
   endif
 endif

 if(arg%no.gt.0)then
    allflg(:) = 0.0
    allinc(:) = 0.0
    allb_a(:) = 0.0
    nobs = arg%no
   call mpi_reduce(  arg%flg, allflg, nobs, mpi%i8  ,   &
                     mpi_max, 0, mpi%comm, ierr)
   call mpi_reduce(  arg%inc, allinc, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
   call mpi_reduce(  arg%b_a, allb_a, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
   write(215) arg%no
   if(arg%no.ne.0)  then
    arg%flg(1:arg%no) = allflg(1:arg%no)
    arg%inc(1:arg%no) = allinc(1:arg%no)
    arg%b_a(1:arg%no) = allb_a(1:arg%no)
       write (215)                                  &
        arg%ino(1:arg%no),arg%flg(1:arg%no),arg%par(1:arg%no) &
       ,arg%lon(1:arg%no),arg%lat(1:arg%no),arg%dpt(1:arg%no) &
       ,arg%tim(1:arg%no),arg%val(1:arg%no),arg%bac(1:arg%no) &
       ,arg%err(1:arg%no),arg%res(1:arg%no),arg%bia(1:arg%no) &
       ,arg%inc(1:arg%no),arg%b_a(1:arg%no)
       write(217,"(i8)") arg%no
       do kk=1,arg%no
       write (217,"(3I8,10F10.5)")                  &
        arg%ino(kk),arg%flg(kk),arg%par(kk)         &
       ,arg%lon(kk),arg%lat(kk),arg%dpt(kk)         &
       ,arg%tim(kk),arg%val(kk),arg%bac(kk)         &
       ,arg%err(kk),arg%res(kk),arg%bia(kk)         &
       ,arg%inc(kk)
       enddo
   endif
 endif

 if(xbt%no.gt.0)then
    allflg(:) = 0.0
    allinc(:) = 0.0
    allb_a(:) = 0.0
    nobs = xbt%no
   call mpi_reduce(  xbt%flg, allflg, nobs, mpi%i8  ,   &
                     mpi_max, 0, mpi%comm, ierr)
   call mpi_reduce(  xbt%inc, allinc, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
   call mpi_reduce(  xbt%b_a, allb_a, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
   write(215) xbt%no
   if(xbt%no.ne.0) then
    xbt%flg(1:xbt%no) = allflg(1:xbt%no)
    xbt%inc(1:xbt%no) = allinc(1:xbt%no)
    xbt%b_a(1:xbt%no) = allb_a(1:xbt%no)
         write (215)                                  &
        xbt%ino(1:xbt%no),xbt%flg(1:xbt%no),xbt%par(1:xbt%no) &
       ,xbt%lon(1:xbt%no),xbt%lat(1:xbt%no),xbt%dpt(1:xbt%no) &
       ,xbt%tim(1:xbt%no),xbt%val(1:xbt%no),xbt%bac(1:xbt%no) &
       ,xbt%err(1:xbt%no),xbt%res(1:xbt%no),xbt%bia(1:xbt%no) &
       ,xbt%inc(1:xbt%no),xbt%b_a(1:xbt%no) 
       write(218,"(i8)") xbt%no
       do kk=1,xbt%no
       write (218,"(3I8,11F10.5)")                  &
        xbt%ino(kk),xbt%flg(kk),xbt%par(kk) &
       ,xbt%lon(kk),xbt%lat(kk),xbt%dpt(kk) &
       ,xbt%tim(kk),xbt%val(kk),xbt%bac(kk) &
       ,xbt%err(kk),xbt%res(kk),xbt%bia(kk) &
       ,xbt%inc(kk),xbt%b_a(kk)   
       enddo
   endif
 endif

 if(gld%no.gt.0)then
    allflg(:) = 0.0
    allinc(:) = 0.0
    allb_a(:) = 0.0
    nobs = gld%no
   call mpi_reduce(  gld%flg, allflg, nobs, mpi%i8  ,   &
                     mpi_max, 0, mpi%comm, ierr)
   call mpi_reduce(  gld%inc, allinc, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
   call mpi_reduce(  gld%b_a, allb_a, nobs, mpi%r8  ,   &
                     mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
   write(215) gld%no
   if(gld%no.ne.0) then
    gld%flg(1:gld%no) = allflg(1:gld%no)
    gld%inc(1:gld%no) = allinc(1:gld%no)
    gld%b_a(1:gld%no) = allb_a(1:gld%no)
         write (215)                                  &
        gld%ino(1:gld%no),gld%flg(1:gld%no),gld%par(1:gld%no) &
       ,gld%lon(1:gld%no),gld%lat(1:gld%no),gld%dpt(1:gld%no) &
       ,gld%tim(1:gld%no),gld%val(1:gld%no),gld%bac(1:gld%no) &
       ,gld%err(1:gld%no),gld%res(1:gld%no),gld%err(1:gld%no) &
       ,gld%inc(1:gld%no),gld%b_a(1:gld%no)
       write(219,"(i8)") gld%no
       do kk=1,gld%no
       write (219,"(3I8,11F10.5)")                  &
        gld%ino(kk),gld%flg(kk),gld%par(kk) &
       ,gld%lon(kk),gld%lat(kk),gld%dpt(kk) &
       ,gld%tim(kk),gld%val(kk),gld%bac(kk) &
       ,gld%err(kk),gld%res(kk),gld%err(kk) &
       ,gld%inc(kk),gld%b_a(kk)   
       enddo
   endif
 endif

 if(sst%no.gt.0)then
   allflg(1:sst%no) = 0
   allinc(1:sst%no) = 0.0
   allb_a(1:sst%no) = 0.0
    nobs = sst%no
    call mpi_reduce(  sst%flg, allflg, nobs, mpi%i8  ,   &
                      mpi_max, 0, mpi%comm, ierr)
    call mpi_reduce(  sst%inc, allinc, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
    call mpi_reduce(  sst%b_a, allb_a, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
   write(215) sst%no
   if(sst%no.ne.0) then
    sst%flg(1:sst%no) = allflg(1:sst%no)
    sst%inc(1:sst%no) = allinc(1:sst%no)
    sst%b_a(1:sst%no) = allb_a(1:sst%no)
         write (215)                                  &
      sst%ino(1:sst%no), sst%flg(1:sst%no), sst%par(1:sst%no) &
     ,sst%lon(1:sst%no), sst%lat(1:sst%no)                    &
     ,sst%dpt(1:sst%no), sst%tim(1:sst%no)                    &
     ,sst%val(1:sst%no), sst%bac(1:sst%no)                    &
     ,sst%err(1:sst%no), sst%res(1:sst%no)                    &
     ,sst%ib(1:sst%no), sst%jb(1:sst%no), sst%kb(1:sst%no)    &
     ,sst%pb(1:sst%no), sst%qb(1:sst%no), sst%rb(1:sst%no)
   endif
 endif

 if(vdr%no.ne.0) then
   allflg(1:vdr%no) = 0
   allinc(1:vdr%no) = 0.0
   allb_a(1:vdr%no) = 0.0
    nobs = vdr%no
    call mpi_reduce(  vdr%flg, allflg, nobs, mpi%i8  ,   &
                      mpi_max, 0, mpi%comm, ierr)
    call mpi_reduce(  vdr%inc, allinc, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
    call mpi_reduce(  vdr%b_a, allb_a, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
    write(215) vdr%no
   if(vdr%no.ne.0) then
    vdr%flg(1:vdr%no) = allflg(1:vdr%no)
    vdr%inc(1:vdr%no) = allinc(1:vdr%no)
    vdr%b_a(1:vdr%no) = allb_a(1:vdr%no)
         write (215)                                  &
      vdr%ino(1:vdr%no), vdr%flg(1:vdr%no), vdr%par(1:vdr%no) &
     ,vdr%lon(1:vdr%no), vdr%lat(1:vdr%no)                    &
     ,vdr%dpt(1:vdr%no), vdr%tim(1:vdr%no)                    &
     ,vdr%tms(1:vdr%no), vdr%tme(1:vdr%no)                    &
     ,vdr%val(1:vdr%no), vdr%bac(1:vdr%no)                    &
     ,vdr%err(1:vdr%no), vdr%res(1:vdr%no)                    &
     ,vdr%bia(1:vdr%no), vdr%inc(1:vdr%no)                    &
     ,vdr%b_a(1:vdr%no)
   endif
 endif

 if(gvl%no.ne.0) then
   allflg(1:gvl%no) = 0
   allinc(1:gvl%no) = 0.0
   allb_a(1:gvl%no) = 0.0
    nobs = gvl%no
    call mpi_reduce(  gvl%flg, allflg, nobs, mpi%i8  ,   &
                      mpi_max, 0, mpi%comm, ierr)
    call mpi_reduce(  gvl%inc, allinc, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
    call mpi_reduce(  gvl%b_a, allb_a, nobs, mpi%r8  ,   &
                      mpi_sum, 0, mpi%comm, ierr)
 endif
 if(mpi%myrank.eq.0)then
    write(215) gvl%no
   if(gvl%no.ne.0) then
    gvl%flg(1:gvl%no) = allflg(1:gvl%no)
    gvl%inc(1:gvl%no) = allinc(1:gvl%no)
    gvl%b_a(1:gvl%no) = allb_a(1:gvl%no)
         write (215)                                  &
      gvl%ino(1:gvl%no), gvl%flg(1:gvl%no), gvl%par(1:gvl%no) &
     ,gvl%lon(1:gvl%no), gvl%lat(1:gvl%no)                    &
     ,gvl%dpt(1:gvl%no), gvl%tim(1:gvl%no)                    &
     ,gvl%tms(1:gvl%no), gvl%tme(1:gvl%no)                    &
     ,gvl%val(1:gvl%no), gvl%bac(1:gvl%no)                    &
     ,gvl%err(1:gvl%no), gvl%res(1:gvl%no)                    &
     ,gvl%bia(1:gvl%no), gvl%inc(1:gvl%no)                    &
     ,gvl%b_a(1:gvl%no)
   endif
 endif

 if(mpi%myrank.eq.0)then
    write(215) tra%no
   if(tra%no.ne.0) then
       write (215)                                              &
        tra%dpt                                                 &
       ,tra%ino(1:tra%no), tra%flg(1:tra%no)                    &
       ,tra%loi(1:tra%no), tra%lai(1:tra%no)                    &
       ,tra%lof(1:tra%no), tra%laf(1:tra%no)                    &
       ,tra%lob(tra%nt+1,1:tra%no), tra%lab(tra%nt+1,1:tra%no)  &
       ,tra%rex(1:tra%no), tra%inx(1:tra%no)                    &
       ,tra%rey(1:tra%no), tra%iny(1:tra%no)                    &
       ,tra%loa(1:tra%no), tra%laa(1:tra%no)                    &
       ,tra%erx(1:tra%no), tra%ery(1:tra%no)
   endif
 endif

 if(mpi%myrank.eq.0)then
    write(215) trd%no
   if(trd%no.ne.0) then
       write (215)                                              &
        trd%dpt                                                 &
       ,trd%ino(1:trd%no), trd%flg(1:trd%no)                    &
       ,trd%loi(1:trd%no), trd%lai(1:trd%no)                    &
       ,trd%lof(1:trd%no), trd%laf(1:trd%no)                    &
       ,trd%lob(trd%nt+1,1:trd%no), trd%lab(trd%nt+1,1:trd%no)  &
       ,trd%rex(1:trd%no), trd%inx(1:trd%no)                    &
       ,trd%rey(1:trd%no), trd%iny(1:trd%no)                    &
       ,trd%loa(1:trd%no), trd%laa(1:trd%no)                    &
       ,trd%erx(1:trd%no), trd%ery(1:trd%no)
   endif
 endif


  if(mpi%myrank.eq.0) close (215)
  if(mpi%myrank.eq.0) close (216)
  if(mpi%myrank.eq.0) close (217)
  if(mpi%myrank.eq.0) close (218)
  if(mpi%myrank.eq.0) close (219)

  if(maxno.gt.0) DEALLOCATE ( allflg, allinc, allb_a ) 

end subroutine wrt_dia

subroutine write_inn ( km, fldin, cpr)
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Write outputs and diagnostics communicator                           !
!                                                                      !
! Version 1: S.Dobricic 2010                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str
 use mpi_str
 use netcdf

 implicit none

 include "mpif.h"

 INTEGER(i4)            :: km
 INTEGER(i4)            :: i, j, kk, k
! REAL(r8)                  ::  fldin(grd%im,grd%jm)
 REAL(r8)                  ::  fldin(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,km)
 REAL(r4), allocatable     ::  bff(:)
 REAL(r4), allocatable     ::  bffa(:,:)
 REAL(r4), allocatable     ::  fld(:,:)

 INTEGER                    :: ierr, kproc
 INTEGER                    :: stat, ncid, idvar, idlon, idlat, iddep
 CHARACTER*3                :: cpr

 if(mpi%myrank.eq.0) ALLOCATE ( fld(grd%img,grd%jmg) )

  ALLOCATE ( bff(grd%npsm), bffa(grd%npsm,mpi%nproc) )

  do k = 1, km

 if(mpi%nproc.gt.1)then


     bff(:) = 0.0
  kk = 0
    do j=1,grd%jm
     do i=1,grd%im
      kk = kk + 1
      bff(kk) = fldin(i,j,k)
     enddo
    enddo

  call mpi_gather(bff, grd%npsm, mpi_real, bffa, grd%npsm, mpi_real, 0, mpi%comm, ierr)


  if(mpi%myrank.eq.0)then

       fld(:,:) = 0.0
   do kproc = 1,mpi%nproc
     kk = 0
     do j=grd%aj1(kproc),grd%ajm(kproc)
      do i=grd%ai1(kproc),grd%aim(kproc)
       kk = kk + 1
       fld(i,j) = bffa(kk,kproc)
      enddo
     enddo
    enddo

  endif


 else

    fld(:,:) = fldin(:,:,k)

 endif

 if(mpi%myrank.eq.0) then

 if(k.eq.1)then

  stat = nf90_create('corr_'//cpr//'.nc', nf90_share, ncid)
  stat = nf90_open  ('corr_'//cpr//'.nc', nf90_write, ncid)

  stat = nf90_redef(ncid)

  stat = nf90_def_dim(ncid, 'im', grd%img-2*grd%iex, idlon)
  stat = nf90_def_dim(ncid, 'jm', grd%jmg, idlat)
  stat = nf90_def_dim(ncid, 'km', grd%km , iddep)

  if(km.eq.1)then
   stat = nf90_def_var ( ncid, cpr, nf90_float, (/ idlon, idlat /), idvar)
  else
   stat = nf90_def_var ( ncid, cpr, nf90_float, (/ idlon, idlat, iddep /), idvar)
  endif

  stat = nf90_enddef(ncid)
  if (stat /= nf90_noerr) call netcdf_err(stat)

!  open(101,file='corr_'//cpr//'.dat',form='unformatted')

 endif

  if(km.eq.1)then
   stat = nf90_put_var ( ncid, idvar,                                 &
                         fld(grd%iex+1:grd%img-grd%iex,1:grd%jmg),    &
                         start = (/ 1, 1 /),                          &
                         count = (/ grd%img-2*grd%iex, grd%jmg/) )
  else
   stat = nf90_put_var ( ncid, idvar,                                 &
                         fld(grd%iex+1:grd%img-grd%iex,1:grd%jmg),    &
                         start = (/ 1, 1, k /),                       &
                         count = (/ grd%img-2*grd%iex, grd%jmg, 1/) )
  endif
   if (stat /= nf90_noerr) call netcdf_err(stat)

!  write(101) fld

  if(k.eq.km)then
    close(101)
    stat = nf90_close (ncid)
  endif


 endif

  enddo ! k

  DEALLOCATE ( bff, bffa )

 if(mpi%myrank.eq.0) DEALLOCATE ( fld )

  call mpi_barrier(mpi%comm,ierr)

end subroutine write_inn

