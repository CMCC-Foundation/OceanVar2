subroutine get_obs_tra

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
! Load Argo trajectory observations                                    !
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
  INTEGER(i8)   ::  ncc


   tra%no = 0
   tra%nc = 0
   tra%ncs = 0
   tra%ncc = 0

!  if(mpi%myrank.eq.0) then

   
    open(511,file='tra_obs.dat',form='unformatted',status='old',err=1111)

    read(511) tra%no, tra%nc, tra%nt, tra%im, tra%jm, tra%km, tra%dpt

! ---
! Allocate memory for observations 

   write(drv%dia,*) 'Number of Argo trajectory observations: ',  tra%no,tra%nc


   if(tra%no.eq.0)then
      close(511)
      return
   endif

   allocate ( tra%ino(tra%no), tra%flg(tra%no), tra%flc(tra%no))
   allocate ( tra%loi(tra%no), tra%lai(tra%no), tra%tim(tra%no), tra%dtm(tra%no))
   allocate ( tra%lof(tra%no), tra%laf(tra%no))
   allocate ( tra%err(tra%no))
   allocate ( tra%lob(tra%nt+1,tra%no), tra%lab(tra%nt+1,tra%no) )
   allocate ( tra%loa(tra%no), tra%laa(tra%no) )
   allocate ( tra%xob(tra%no), tra%xmn(tra%nt+1,tra%no), tra%erx(tra%no) )
   allocate ( tra%yob(tra%no), tra%ymn(tra%nt+1,tra%no), tra%ery(tra%no) )

   allocate ( tra%rex(tra%no), tra%inx(tra%no))
   allocate ( tra%rey(tra%no), tra%iny(tra%no))
   allocate ( tra%xtl(tra%no), tra%ytl(tra%no) )
   allocate ( tra%xtl_ad(tra%no), tra%ytl_ad(tra%no) )
   allocate ( tra%fls(tra%no))
   allocate ( tra%rsx(tra%no))
   allocate ( tra%rsy(tra%no))
   allocate ( tra%isx(tra%no))
   allocate ( tra%isy(tra%no))

   tra%fls(:) = 0.0
   tra%rsx(:) = 0.0
   tra%isx(:) = 0.0
   tra%rsy(:) = 0.0
   tra%isy(:) = 0.0


   read(511)  tra%ino, tra%flg, tra%tim, tra%dtm,      &
              tra%loi, tra%lai, tra%lof, tra%laf,      &
              tra%err, tra%lob, tra%lab,               &
              tra%xob, tra%xmn, tra%erx,               &
              tra%yob, tra%ymn, tra%ery

   close(511)

! ---
! Initialise quality flag
   if(obs%tra.eq.0)then
    do j=1,tra%no
     if(tra%flg(j).ne.0)tra%flg(j)=-1
    enddo
   endif
   tra%flc(:) = tra%flg(:)

    do j=1,tra%no
     tra%rex(j) = tra%xob(j) - tra%xmn(tra%nt+1,j)
     tra%rey(j) = tra%yob(j) - tra%ymn(tra%nt+1,j)
     tra%loa(j) = tra%lob(tra%nt+1,j)
     tra%laa(j) = tra%lab(tra%nt+1,j)
    enddo


! residual check
  do j=1,tra%no
   if(abs(tra%rex(j)).gt.10. .or. abs(tra%rey(j)).gt.10.) tra%flg(j) = 0
  enddo

! ---
! Set quality flags to zero on all processors except the first one

  if(mpi%myrank.gt.0) tra%flg(:) = 0

! ---
! Count good observations
    tra%nc = 0
  do k=1,tra%no
   if(tra%flg(k).eq.1)then
    tra%nc = tra%nc + 1
   endif
  enddo

  tra%flc(:) = tra%flg(:)

  tra%ncc = tra%nc

   write(drv%dia,*) 'Number of good Argo trajectory observations: ',  tra%nc


  if(mpi%nproc.gt.1) call mpi_bcast( tra%ncc, 1, mpi%i8, 0, mpi%comm, ierr)

1111 continue



end subroutine get_obs_tra

subroutine int_par_tra

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

  if(tra%no.eq.0) return

   if(mpi%myrank.eq.0) then
    img = grd%img
    jmg = grd%jmg
   else
    img = 1
    jmg = 1
   endif

   allocate ( tra%umn(img,jmg), tra%vmn(img,jmg) )
   allocate( tra%uvl(img,jmg), tra%vvl(img,jmg) )
   allocate( tra%uvl_ad(img,jmg), tra%vvl_ad(img,jmg) )
   allocate ( tra%dx(img,jmg), tra%dy(img,jmg) )

  if(mpi%myrank.eq.0) then
    open(511,file='tra_umn.dat',form='unformatted',status='old')
    read(511)  tra%umn
    read(511)  tra%vmn
    close(511)
  endif


    k   = 1
    km  = 1

   if(mpi%nproc.gt.1)then
      call gth_mpi( img, jmg, k, km, grd%dx, tra%dx)
      call gth_mpi( img, jmg, k, km, grd%dy, tra%dy)
   else
      tra%dx(:,:) = grd%dx(:,:)
      tra%dy(:,:) = grd%dy(:,:)
   endif

   tra%lev = 1
   do k=1,grd%km-1
    if( tra%dpt.ge.grd%dep(k) .and. tra%dpt.lt.grd%dep(k+1) ) tra%lev = k
   enddo

  tra%inx(:) = 0.0
  tra%iny(:) = 0.0


end subroutine int_par_tra
