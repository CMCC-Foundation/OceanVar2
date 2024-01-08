subroutine get_obs_trd

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
! Load trajectory observations by surface drifters                     !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use grd_str
 use obs_str
 use mpi_str

 implicit none

  INTEGER(i4)   ::  k, ierr
  INTEGER(i4)   ::  i1, j1, i, j


   trd%no = 0
   trd%nc = 0
   trd%ncs = 0
   trd%ncc = 0


   

    open(511,file='trd_obs.dat',form='unformatted',status='old',err=1111)

    read(511) trd%no, trd%nc, trd%nt, trd%im, trd%jm, trd%km, trd%dpt

! ---
! Allocate memory for observations 

   write(drv%dia,*) 'Number of trajectory observations by surface drifters: ',  trd%no


   if(trd%no.eq.0)then
      close(511)
      return
   endif

   allocate ( trd%ino(trd%no), trd%flg(trd%no), trd%flc(trd%no))
   allocate ( trd%loi(trd%no), trd%lai(trd%no), trd%tim(trd%no), trd%dtm(trd%no))
   allocate ( trd%lof(trd%no), trd%laf(trd%no))
   allocate ( trd%err(trd%no))
   allocate ( trd%lob(trd%nt+1,trd%no), trd%lab(trd%nt+1,trd%no) )
   allocate ( trd%loa(trd%no), trd%laa(trd%no) )
   allocate ( trd%xob(trd%no), trd%xmn(trd%nt+1,trd%no), trd%erx(trd%no) )
   allocate ( trd%yob(trd%no), trd%ymn(trd%nt+1,trd%no), trd%ery(trd%no) )

   allocate ( trd%rex(trd%no), trd%inx(trd%no))
   allocate ( trd%rey(trd%no), trd%iny(trd%no))
   allocate ( trd%xtl(trd%no), trd%ytl(trd%no) )
   allocate ( trd%xtl_ad(trd%no), trd%ytl_ad(trd%no) )
   allocate ( trd%fls(trd%no))
   allocate ( trd%rsx(trd%no))
   allocate ( trd%rsy(trd%no))
   allocate ( trd%isx(trd%no))
   allocate ( trd%isy(trd%no))

   trd%fls(:) = 0.0
   trd%rsx(:) = 0.0
   trd%isx(:) = 0.0
   trd%rsy(:) = 0.0
   trd%isy(:) = 0.0


   read(511)  trd%ino, trd%flg, trd%tim, trd%dtm,      &
              trd%loi, trd%lai, trd%lof, trd%laf,      &
              trd%err, trd%lob, trd%lab,               &
              trd%xob, trd%xmn, trd%erx,               &
              trd%yob, trd%ymn, trd%ery

   close(511)

! ---
! Initialise quality flag
   if(obs%trd.eq.0)then
    do j=1,trd%no
     if(trd%flg(j).ne.0)trd%flg(j)=-1
    enddo
   endif
   trd%flc(:) = trd%flg(:)

    do j=1,trd%no
     trd%rex(j) = trd%xob(j) - trd%xmn(trd%nt+1,j)
     trd%rey(j) = trd%yob(j) - trd%ymn(trd%nt+1,j)
     trd%loa(j) = trd%lob(trd%nt+1,j)
     trd%laa(j) = trd%lab(trd%nt+1,j)
    enddo


! residual check
  do j=1,trd%no
   if(abs(trd%rex(j)).gt.10. .or. abs(trd%rey(j)).gt.10.) trd%flg(j) = 0
  enddo


! ---
! Set quality flags to zero on all processors except the first one

  if(mpi%myrank.gt.0) trd%flg(:) = 0


! ---
! Count good observations
    trd%nc = 0
  do k=1,trd%no
   if(trd%flg(k).eq.1)then
    trd%nc = trd%nc + 1
   endif
  enddo

  trd%flc(:) = trd%flg(:)

  trd%ncc = trd%nc

 write(drv%dia,*) 'Number of good drifter observations: ',  tra%nc

  if(mpi%nproc.gt.1) call mpi_bcast( trd%ncc, 1, mpi%i8, 0, mpi%comm, ierr)

1111 continue

   

end subroutine get_obs_trd

subroutine int_par_trd

!-----------------------------------------------------------------------
!                                                                      !
! Get interpolation parameters for a grid                              !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use grd_str
 use obs_str
 use mpi_str

 implicit none

  INTEGER(i4)   ::  img, jmg, k, km

  if(trd%no.eq.0) return

   if(mpi%myrank.eq.0) then
    img = grd%img
    jmg = grd%jmg
   else
    img = 1
    jmg = 1
   endif

   allocate ( trd%umn(img,jmg), trd%vmn(img,jmg) )
   allocate( trd%uvl(img,jmg), trd%vvl(img,jmg) )
   allocate( trd%uvl_ad(img,jmg), trd%vvl_ad(img,jmg) )
   allocate ( trd%dx(img,jmg), trd%dy(img,jmg) )


  if(mpi%myrank.eq.0) then
    open(511,file='trd_umn.dat',form='unformatted',status='old')
    read(511)  trd%umn
    read(511)  trd%vmn
    close(511)
  endif

    k   = 1
    km  = 1

   if(mpi%nproc.gt.1)then
      call gth_mpi( img, jmg, k, km, grd%dx, trd%dx)
      call gth_mpi( img, jmg, k, km, grd%dy, trd%dy)
   else
      trd%dx(:,:) = grd%dx(:,:)
      trd%dy(:,:) = grd%dy(:,:)
   endif

   trd%lev = 1
   do k=1,grd%km-1
    if( trd%dpt.ge.grd%dep(k) .and. trd%dpt.lt.grd%dep(k+1) ) trd%lev = k
   enddo

  trd%inx(:) = 0.0
  trd%iny(:) = 0.0

end subroutine int_par_trd
