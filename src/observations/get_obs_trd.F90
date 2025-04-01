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
!> Load trajectory observations by surface drifters                     
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic and Vincent Taillandier 2007              !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_trd

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : trd,  obs
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k, ierr
   INTEGER(i4)   ::  i1, j1, i, j

   trd%no = 0
   trd%nc = 0
   trd%ncs = 0
   trd%ncc = 0

   OPEN (511,FILE=drv%inpdir//'/trd_obs.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ(511) trd%no, trd%nc, trd%nt, trd%im, trd%jm, trd%km, trd%dpt

   WRITE (drv%dia,*) 'Number of trajectory observations by surface drifters: ',  trd%no

   IF ( trd%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( trd%ino(trd%no), trd%flg(trd%no), trd%flc(trd%no) )
   ALLOCATE ( trd%loi(trd%no), trd%lai(trd%no), trd%tim(trd%no), trd%dtm(trd%no) )
   ALLOCATE ( trd%lof(trd%no), trd%laf(trd%no) )
   ALLOCATE ( trd%err(trd%no) )
   ALLOCATE ( trd%lob(trd%nt+1,trd%no), trd%lab(trd%nt+1,trd%no) )
   ALLOCATE ( trd%loa(trd%no), trd%laa(trd%no) )
   ALLOCATE ( trd%xob(trd%no), trd%xmn(trd%nt+1,trd%no), trd%erx(trd%no) )
   ALLOCATE ( trd%yob(trd%no), trd%ymn(trd%nt+1,trd%no), trd%ery(trd%no) )
   ALLOCATE ( trd%rex(trd%no), trd%inx(trd%no) )
   ALLOCATE ( trd%rey(trd%no), trd%iny(trd%no) )
   ALLOCATE ( trd%xtl(trd%no), trd%ytl(trd%no) )
   ALLOCATE ( trd%xtl_ad(trd%no), trd%ytl_ad(trd%no) )
   ALLOCATE ( trd%fls(trd%no) )
   ALLOCATE ( trd%rsx(trd%no) )
   ALLOCATE ( trd%rsy(trd%no) )
   ALLOCATE ( trd%isx(trd%no) )
   ALLOCATE ( trd%isy(trd%no), trd%eve(trd%no) )

   trd%fls(:) = 0_i8
   trd%rsx(:) = 0.0_r8
   trd%isx(:) = 0.0_r8
   trd%rsy(:) = 0.0_r8
   trd%isy(:) = 0.0_r8
   trd%eve(:) = 999_i8

   READ(511)  trd%ino, trd%flg, trd%tim, trd%dtm,      &
              trd%loi, trd%lai, trd%lof, trd%laf,      &
              trd%err, trd%lob, trd%lab,               &
              trd%xob, trd%xmn, trd%erx,               &
              trd%yob, trd%ymn, trd%ery
   CLOSE (511)

! ---
! Initialise quality flag
   IF ( obs%trd .EQ. 0 ) THEN
      DO j = 1,trd%no
         IF ( trd%flg(j) .NE. 0 )trd%flg(j)=-1
      ENDDO
      WRITE (drv%dia,*)'Bad quality flag ',obs%trd
   ENDIF
   trd%flc(:) = trd%flg(:)

   DO j = 1,trd%no
      trd%rex(j) = trd%xob(j) - trd%xmn(trd%nt+1,j)
      trd%rey(j) = trd%yob(j) - trd%ymn(trd%nt+1,j)
      trd%loa(j) = trd%lob(trd%nt+1,j)
      trd%laa(j) = trd%lab(trd%nt+1,j)
   ENDDO

! ---
! Set quality flags to zero on all processors except the first one
   IF ( mpi%myrank .GT. 0 ) trd%flg(:) = 0

! ---
! Count good observations
   trd%nc = 0
   DO k = 1,trd%no
      IF ( trd%flg(k) .EQ. 1 ) THEN
         trd%nc = trd%nc + 1
      ENDIF
   ENDDO

   trd%flc(:) = trd%flg(:)
   WHERE ( trd%flc .EQ. 0 ) trd%eve = 1

   trd%ncc = trd%nc

   WRITE (drv%dia,*) 'Number of good drifter observations after reading: ',  trd%nc

   IF ( mpi%nproc .GT. 1 ) CALL mpi_bcast( trd%ncc, 1, mpi%i8, 0, mpi%comm, ierr)

1111 CONTINUE

END SUBROUTINE get_obs_trd
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                              
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_trd

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  img, jmg, k, km

   IF ( trd%no .EQ. 0 ) RETURN

   IF ( mpi%myrank .EQ. 0 ) THEN
      img = grd%img
      jmg = grd%jmg
   ELSE
      img = 1
      jmg = 1
   ENDIF

   ALLOCATE ( trd%umn(img,jmg), trd%vmn(img,jmg) )
   ALLOCATE ( trd%uvl(img,jmg), trd%vvl(img,jmg) )
   ALLOCATE ( trd%uvl_ad(img,jmg), trd%vvl_ad(img,jmg) )
   ALLOCATE ( trd%dx(img,jmg), trd%dy(img,jmg) )


   IF ( mpi%myrank .EQ. 0 ) THEN
      OPEN (511,FILE=drv%inpdir//'/trd_umn.dat',FORM='unformatted',STATUS='old')
      READ (511)  trd%umn
      READ (511)  trd%vmn
      CLOSE (511)
   ENDIF

   k   = 1
   km  = 1

   IF ( mpi%nproc .GT. 1 ) THEN
      CALL gth_mpi( img, jmg, k, km, grd%dx, trd%dx)
      CALL gth_mpi( img, jmg, k, km, grd%dy, trd%dy)
   ELSE
      trd%dx(:,:) = grd%dx(:,:)
      trd%dy(:,:) = grd%dy(:,:)
   ENDIF

   trd%lev = 1
   DO k = 1,grd%km-1
      IF ( trd%dpt .GE. grd%dep(k) .AND. trd%dpt .LT. grd%dep(k+1) ) trd%lev = k
   ENDDO

   trd%inx(:) = 0.0_r8
   trd%iny(:) = 0.0_r8

END SUBROUTINE int_par_trd
