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
!> Load GLIDER observations                                            
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
! Version 2: Paolo Oddo      2015                                      !
! Version 3: Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_gld

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)               ::  k,kk

   gld%no = 0_i8
   gld%nc = 0_i8

   OPEN (511,FILE=drv%inpdir//'/gld_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) gld%no

   WRITE (drv%dia,*) 'Number of glider observations: ',gld%no

   IF ( gld%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( gld%ino(gld%no), gld%flg(gld%no), gld%flc(gld%no), gld%par(gld%no) )
   ALLOCATE ( gld%lon(gld%no), gld%lat(gld%no), gld%dpt(gld%no), gld%tim(gld%no) )
   ALLOCATE ( gld%val(gld%no), gld%bac(gld%no), gld%inc(gld%no) )
   ALLOCATE ( gld%bia(gld%no), gld%err(gld%no) )
   ALLOCATE ( gld%res(gld%no), gld%b_a(gld%no) )
   ALLOCATE ( gld%ib(gld%no), gld%jb(gld%no), gld%kb(gld%no) )
   ALLOCATE ( gld%pb(gld%no), gld%qb(gld%no), gld%rb(gld%no) )
   ALLOCATE ( gld%pq1(gld%no), gld%pq2(gld%no), gld%pq3(gld%no), gld%pq4(gld%no) )
   ALLOCATE ( gld%pq5(gld%no), gld%pq6(gld%no), gld%pq7(gld%no), gld%pq8(gld%no) )
   ALLOCATE ( gld%rss(gld%no), gld%ins(gld%no) )
   ALLOCATE ( gld%fls(gld%no),gld%eve(gld%no) )

   arg%flc(:) = 0_i8
   arg%inc(:) = 0.0_r8
   arg%bia(:) = 0.0_r8
   arg%b_a(:) = 0.0_r8
   arg%pq1(:) = 0.0_r8
   arg%pq2(:) = 0.0_r8
   arg%pq3(:) = 0.0_r8
   arg%pq4(:) = 0.0_r8
   arg%pq5(:) = 0.0_r8
   arg%pq6(:) = 0.0_r8
   arg%pq7(:) = 0.0_r8
   arg%pq8(:) = 0.0_r8
   arg%rss(:) = 0.0_r8
   arg%ins(:) = 0.0_r8
   arg%fls(:) = 0_i8
   arg%eve(:) = 999_i8

   READ (511)                                               &
       gld%ino(1:gld%no), gld%flg(1:gld%no), gld%par(1:gld%no) &
      ,gld%lon(1:gld%no), gld%lat(1:gld%no)                    &
      ,gld%dpt(1:gld%no), gld%tim(1:gld%no)                    &
      ,gld%val(1:gld%no), gld%bac(1:gld%no)                    &
      ,gld%err(1:gld%no), gld%res(1:gld%no)                    &
      ,gld%ib(1:gld%no), gld%jb(1:gld%no), gld%kb(1:gld%no)    &
      ,gld%pb(1:gld%no), gld%qb(1:gld%no), gld%rb(1:gld%no)
   CLOSE (511)

! ---
! Initialise quality flag
   IF ( obs%gld .EQ. 0 ) THEN
      gld%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%gld
   ENDIF

! ---
! Vertical interpolation PARAMETERs
   DO k = 1,gld%no
      IF ( gld%flg(k) .EQ. 1 ) THEN
         gld%kb(k) = grd%km-1
         DO kk = 1,grd%km-1
            IF ( gld%dpt(k) .GE. grd%dep(kk) .AND. gld%dpt(k) .LT. grd%dep(kk+1) ) THEN
               gld%kb(k) = kk
               gld%rb(k) = (gld%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
            ENDIF
         ENDDO
      ENDIF
   ENDDO

! ---
! Count good observations
   gld%nc = 0
   DO k = 1,gld%no
      IF ( gld%flg(k) .EQ. 1 ) THEN
         gld%nc = gld%nc + 1
      ELSE
         gld%bia(k) = 0.0_r8
         gld%res(k) = 0.0_r8
         gld%inc(k) = 0.0_r8
         gld%b_a(k) = 0.0_r8
         gld%pq1(k) = 0.0_r8
         gld%pq2(k) = 0.0_r8
         gld%pq3(k) = 0.0_r8
         gld%pq4(k) = 0.0_r8
         gld%pq5(k) = 0.0_r8
         gld%pq6(k) = 0.0_r8
         gld%pq7(k) = 0.0_r8
         gld%pq8(k) = 0.0_r8
      ENDIF
   ENDDO

   gld%flc(:) = gld%flg(:)
   WHERE ( gld%flc .EQ. 0 ) gld%eve = 1

   WRITE (drv%dia,*) 'Number of good GLIDER observations after reading: ',  gld%nc

1111 CONTINUE

END SUBROUTINE get_obs_gld
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                             
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Paolo Oddo      2015                                      !
! Version 3: Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_gld

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, k
   INTEGER(i8)   ::  klev
   INTEGER(i4)   ::  i1, j1, k1, idep
   REAL(r8)      ::  p1, q1, r1
   REAL(r8)      ::  msk4, div_x, div_y

   IF (gld%no.GT.0) THEN

      gld%flc(:) = gld%flg(:)
! --------------------------------------------
! Adjust longitudes
! --------------------------------------------
      IF ( grd%bwst .GT. 180. ) THEN
         DO k = 1,gld%no
            IF ( gld%lon(k) .LT. 0.0_r8 ) THEN
               gld%lon(k) = gld%lon(k) + 360.
            ENDIF
         ENDDO
      ENDIF

! --------------------------------------------
! Horizontal interpolation PARAMETERs
! --------------------------------------------
      CALL int_obs_hor ( gld%no, gld%lat, gld%lon, gld%flc, gld%eve, gld%ib, gld%jb, gld%pb, gld%qb)

! --------------------------------------------
! Undefine masked for multigrid
! --------------------------------------------
      DO k = 1,gld%no
         IF ( gld%flc(k) .EQ. 1 ) THEN
            i1 = gld%ib(k)
            j1 = gld%jb(k)
            idep = gld%kb(k)+1
            msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
            IF ( msk4 .LT. 1. ) THEN
               gld%flc(k) = 0
               gld%eve(k) = 3
            ENDIF
         ENDIF
      ENDDO

! --------------------------------------------
! Horizontal interpolation PARAMETERs for each masked grid
! --------------------------------------------
      DO k = 1,gld%no
         IF ( gld%flc(k) .EQ. 1 ) THEN

            klev = gld%kb(k)
            CALL int_obs_pq( gld%ib(k), gld%jb(k), klev, gld%pb(k), gld%qb(k),  &
                             gld%pq1(k),gld%pq2(k), gld%pq3(k), gld%pq4(k))
            klev = gld%kb(k) + 1
            CALL int_obs_pq( gld%ib(k), gld%jb(k), klev, gld%pb(k), gld%qb(k),  &
                             gld%pq5(k),gld%pq6(k), gld%pq7(k), gld%pq8(k))
            r1=gld%rb(k)
            gld%pq1(k) = (1.-r1) * gld%pq1(k)
            gld%pq2(k) = (1.-r1) * gld%pq2(k)
            gld%pq3(k) = (1.-r1) * gld%pq3(k)
            gld%pq4(k) = (1.-r1) * gld%pq4(k)
            gld%pq5(k) =     r1  * gld%pq5(k)
            gld%pq6(k) =     r1  * gld%pq6(k)
            gld%pq7(k) =     r1  * gld%pq7(k)
            gld%pq8(k) =     r1  * gld%pq8(k)
         ENDIF
      ENDDO

! ---
! Count good observations
      gld%nc = 0
      DO k = 1,gld%no
         IF ( gld%flc(k) .EQ. 1 ) THEN
            gld%nc = gld%nc + 1
         ENDIF
      ENDDO

      gld%inc(:) = 0.0_r8

   ENDIF

END SUBROUTINE int_par_gld
