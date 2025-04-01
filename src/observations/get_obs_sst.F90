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
!> Load SST observations                                                
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Jenny Pistoia   2013                                      !
!            Mario Adani 2023 Add events                               !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_sst

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : sst,  obs

   IMPLICIT NONE

   INTEGER       ::  ierr
   INTEGER(i4)   ::  k
   INTEGER(i4)   ::  i1, i, iter
   REAL(r8)      ::  sumt, sumi, timp, dxx, dyy, dsm

   sst%no = 0
   sst%nc = 0


   OPEN (511,FILE=drv%inpdir//'/sst_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) sst%no

   WRITE (drv%dia,*) 'Number of SST observations: ',  sst%no, obs%sst

   IF (sst%no.EQ.0) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( sst%ino(sst%no),sst%par(sst%no), sst%flg(sst%no),sst%flc(sst%no) )
   ALLOCATE ( sst%lon(sst%no), sst%lat(sst%no), sst%dpt(sst%no), sst%tim(sst%no) )
   ALLOCATE ( sst%val(sst%no), sst%bac(sst%no), sst%inc(sst%no) )
   ALLOCATE ( sst%bia(sst%no), sst%err(sst%no), sst%res(sst%no),sst%b_a(sst%no) )
   ALLOCATE ( sst%ib(sst%no),sst%pb(sst%no), sst%jb(sst%no), sst%qb(sst%no) )
   ALLOCATE ( sst%kb(sst%no), sst%rb(sst%no) )
   ALLOCATE ( sst%pq1(sst%no), sst%pq2(sst%no), sst%pq3(sst%no), sst%pq4(sst%no) )
   ALLOCATE ( sst%fls(sst%no),  sst%eve(sst%no) )

! ---
! Initialise quality flag
   sst%flc(:) = 1

   READ (511)                                              &
       sst%ino(1:sst%no), sst%flg(1:sst%no), sst%par(1:sst%no) &
      ,sst%lon(1:sst%no), sst%lat(1:sst%no)                    &
      ,sst%dpt(1:sst%no), sst%tim(1:sst%no)                    &
      ,sst%val(1:sst%no), sst%bac(1:sst%no)                    &
      ,sst%err(1:sst%no), sst%res(1:sst%no)                    &
      ,sst%ib(1:sst%no), sst%jb(1:sst%no), sst%kb(1:sst%no)    &
      ,sst%pb(1:sst%no), sst%qb(1:sst%no), sst%rb(1:sst%no)
   CLOSE (511)

   sst%kb(:) = 1
   sst%eve(:) = 999

! ---
! Initialise quality flag
   IF ( obs%sst .EQ. 0 ) THEN
      sst%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%sst
   ENDIF

! ---
! Count good observations
   sst%nc = 0
   DO k = 1,sst%no
      IF ( sst%flg(k) .EQ. 1 ) THEN
         sst%nc = sst%nc + 1
      ELSE
         sst%inc(k) = 0.0_r8
         sst%b_a(k) = 0.0_r8
         sst%pq1(k) = 0.0_r8
         sst%pq2(k) = 0.0_r8
         sst%pq3(k) = 0.0_r8
         sst%pq4(k) = 0.0_r8
      ENDIF
   ENDDO

   sst%flc(:) = sst%flg(:)
   WHERE ( sst%flc .EQ. 0 ) sst%eve = 1

   WRITE (drv%dia,*) 'Number of good SST observations after interpolation parameters: ',  sst%nc

1111 CONTINUE

END SUBROUTINE get_obs_sst
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                              
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_sst

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k, ierr
   INTEGER(i4)   ::  i1, kk, i, j1, j, idep
   INTEGER(i8)   ::  klev
   REAL(r8)      ::  p1, q1
   REAL(r8)      ::  msk4, div_x, div_y, rmn, dst, dstm
   REAL(r8)      ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
   REAL(r8)      ::  lat_lt_rot, lon_rb_rot

   rmn = 1.e-6

   IF ( sst%no .GT. 0 ) THEN

      sst%flc(:) = sst%flg(:)

      sst%dpt(:) = 0.0_r8

! ---
! Adjust longitudes
      IF ( grd%bwst .GT. 180. ) THEN
         DO k = 1,sst%no
            IF ( sst%lon(k) .LT. 0.0_r8 ) THEN
               sst%lon(k) = sst%lon(k) + 360.
            ENDIF
         ENDDO
      ENDIF

! ---
! Interpolation parameters
      CALL int_obs_hor ( sst%no, sst%lat, sst%lon, sst%flc, sst%eve, sst%ib, sst%jb, sst%pb, sst%qb)

      DO kk = 1,sst%no
         IF ( sst%flc(kk) .EQ. 1 ) THEN
            i1 = sst%ib(kk)
            j1 = sst%jb(kk)
            idep = sst%kb(kk)
            msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
            IF ( msk4 .LT. 1.0) THEN
               sst%flc(kk) = 0
               sst%eve(kk) = 3
            ENDIF
         ENDIF
      ENDDO

! ---
! Horizontal interpolation parameter for each masked grid
      klev = 1
      DO k = 1,sst%no
         IF ( sst%flc(k) .EQ. 1 ) THEN
            CALL int_obs_pq( sst%ib(k), sst%jb(k), klev, sst%pb(k), sst%qb(k),  &
                             sst%pq1(k), sst%pq2(k), sst%pq3(k), sst%pq4(k))
         ENDIF
      ENDDO

! ---
! Count good observations
      sst%nc = 0
      DO k = 1,sst%no
         IF ( sst%flc(k) .EQ. 1 ) THEN
            sst%nc = sst%nc + 1
         ENDIF
      ENDDO

      WRITE (drv%dia,*) 'Number of good SST observations after parameter interpolation: ',  sst%nc

   ENDIF

END SUBROUTINE int_par_sst
