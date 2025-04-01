!======================================================================
!
! This file is part of Oceanvar.
!
!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version (GPL-3.0-or-later).
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!======================================================================
!-----------------------------------------------------------------------
!                                                                      !
!> Load Argo trajectory observations                                    
!!
!!
!!
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_tra

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : tra,  obs
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k, ierr
   INTEGER(i4)   ::  i1, j1, i, j
   INTEGER(i8)   ::  ncc

   tra%no = 0
   tra%nc = 0
   tra%ncs = 0
   tra%ncc = 0

   OPEN (511,FILE=drv%inpdir//'/tra_obs.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) tra%no, tra%nc, tra%nt, tra%im, tra%jm, tra%km, tra%dpt

   WRITE (drv%dia,*) 'Number of Argo trajectory observations: ',  tra%no,tra%nc

   IF (tra%no.EQ.0) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( tra%ino(tra%no), tra%flg(tra%no), tra%flc(tra%no) )
   ALLOCATE ( tra%loi(tra%no), tra%lai(tra%no), tra%tim(tra%no), tra%dtm(tra%no) )
   ALLOCATE ( tra%lof(tra%no), tra%laf(tra%no) )
   ALLOCATE ( tra%err(tra%no) )
   ALLOCATE ( tra%lob(tra%nt+1,tra%no), tra%lab(tra%nt+1,tra%no) )
   ALLOCATE ( tra%loa(tra%no), tra%laa(tra%no) )
   ALLOCATE ( tra%xob(tra%no), tra%xmn(tra%nt+1,tra%no), tra%erx(tra%no) )
   ALLOCATE ( tra%yob(tra%no), tra%ymn(tra%nt+1,tra%no), tra%ery(tra%no) )
   ALLOCATE ( tra%rex(tra%no), tra%inx(tra%no) )
   ALLOCATE ( tra%rey(tra%no), tra%iny(tra%no) )
   ALLOCATE ( tra%xtl(tra%no), tra%ytl(tra%no) )
   ALLOCATE ( tra%xtl_ad(tra%no), tra%ytl_ad(tra%no) )
   ALLOCATE ( tra%fls(tra%no) )
   ALLOCATE ( tra%rsx(tra%no) )
   ALLOCATE ( tra%rsy(tra%no) )
   ALLOCATE ( tra%isx(tra%no) )
   ALLOCATE ( tra%isy(tra%no), tra%eve(tra%no) )

   tra%fls(:) = 0_i8
   tra%rsx(:) = 0.0_r8
   tra%isx(:) = 0.0_r8
   tra%rsy(:) = 0.0_r8
   tra%isy(:) = 0.0_r8
   tra%eve(:) = 999

   READ (511)  tra%ino, tra%flg, tra%tim, tra%dtm,      &
      tra%loi, tra%lai, tra%lof, tra%laf,      &
      tra%err, tra%lob, tra%lab,               &
      tra%xob, tra%xmn, tra%erx,               &
      tra%yob, tra%ymn, tra%ery

   CLOSE (511)

! ---
! Initialise quality flag
   IF ( obs%tra .EQ. 0 ) THEN
      DO j = 1,tra%no
         IF (tra%flg(j) .NE. 0) tra%flg(j)=-1
      ENDDO
      WRITE (drv%dia,*)'Bad quality flag ',obs%tra
   ENDIF
   tra%flc(:) = tra%flg(:)

   DO j = 1,tra%no
      tra%rex(j) = tra%xob(j) - tra%xmn(tra%nt+1,j)
      tra%rey(j) = tra%yob(j) - tra%ymn(tra%nt+1,j)
      tra%loa(j) = tra%lob(tra%nt+1,j)
      tra%laa(j) = tra%lab(tra%nt+1,j)
   ENDDO


! ---
! Set quality flags to zero on all processors except the first one
   IF ( mpi%myrank .GT. 0 ) tra%flg(:) = 0

! ---
! Count good observations
   tra%nc = 0
   DO k = 1,tra%no
      IF ( tra%flg(k) .EQ. 1 ) THEN
         tra%nc = tra%nc + 1
      ENDIF
   ENDDO

   tra%flc(:) = tra%flg(:)
   WHERE ( tra%flc .EQ. 0 ) tra%eve = 1

   tra%ncc = tra%nc

   WRITE (drv%dia,*) 'Number of good Argo trajectory observations after reading: ',  tra%nc

   IF ( mpi%nproc .GT. 1 ) CALL mpi_bcast( tra%ncc, 1, mpi%i8, 0, mpi%comm, ierr)

1111 CONTINUE


END SUBROUTINE get_obs_tra
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                              
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------
SUBROUTINE int_par_tra

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  img, jmg, k, km

   IF ( tra%no .EQ. 0 ) RETURN

   IF ( mpi%myrank .EQ. 0 ) THEN
      img = grd%img
      jmg = grd%jmg
   ELSE
      img = 1
      jmg = 1
   ENDIF

   ALLOCATE ( tra%umn(img,jmg), tra%vmn(img,jmg) )
   ALLOCATE( tra%uvl(img,jmg), tra%vvl(img,jmg) )
   ALLOCATE( tra%uvl_ad(img,jmg), tra%vvl_ad(img,jmg) )
   ALLOCATE ( tra%dx(img,jmg), tra%dy(img,jmg) )

   IF ( mpi%myrank .EQ. 0 ) THEN
      OPEN (511,FILE=drv%inpdir//'/tra_umn.dat',FORM='unformatted',STATUS='old')
      READ (511)  tra%umn
      READ (511)  tra%vmn
      CLOSE (511)
   ENDIF


   k   = 1
   km  = 1

   IF (mpi%nproc.GT.1) THEN
      CALL gth_mpi( img, jmg, k, km, grd%dx, tra%dx)
      CALL gth_mpi( img, jmg, k, km, grd%dy, tra%dy)
   ELSE
      tra%dx(:,:) = grd%dx(:,:)
      tra%dy(:,:) = grd%dy(:,:)
   ENDIF

   tra%lev = 1
   DO k = 1,grd%km-1
      IF ( tra%dpt .GE. grd%dep(k) .AND. tra%dpt .LT. grd%dep(k+1) ) tra%lev = k
   ENDDO

   tra%inx(:) = 0.0_r8
   tra%iny(:) = 0.0_r8

END SUBROUTINE int_par_tra
