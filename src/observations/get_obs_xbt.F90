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
!> Load XBT observations                                                
!!
!!
!!
!                                                                      !
! Version 1  : Srdjan Dobricic 2006                                    !
! Version 1.1: Mario Adani     2024  Add events                        !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_xbt

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : xbt, obs

   IMPLICIT NONE

   INTEGER(i4)   ::  k
   INTEGER(i8)   ::  kk

   xbt%no = 0
   xbt%nc = 0

   OPEN (511,FILE=drv%inpdir//'/xbt_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) xbt%no

   WRITE (drv%dia,*) 'Number of XBT observations: ',  xbt%no, obs%xbt

   IF ( xbt%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( xbt%ino(xbt%no), xbt%flg(xbt%no), xbt%flc(xbt%no), xbt%par(xbt%no) )
   ALLOCATE ( xbt%lon(xbt%no), xbt%lat(xbt%no), xbt%dpt(xbt%no), xbt%tim(xbt%no) )
   ALLOCATE ( xbt%val(xbt%no), xbt%bac(xbt%no), xbt%inc(xbt%no) )
   ALLOCATE ( xbt%bia(xbt%no), xbt%err(xbt%no) )
   ALLOCATE ( xbt%res(xbt%no), xbt%b_a(xbt%no) )
   ALLOCATE ( xbt%ib(xbt%no), xbt%jb(xbt%no), xbt%kb(xbt%no) )
   ALLOCATE ( xbt%pb(xbt%no), xbt%qb(xbt%no), xbt%rb(xbt%no) )
   ALLOCATE ( xbt%pq1(xbt%no), xbt%pq2(xbt%no), xbt%pq3(xbt%no), xbt%pq4(xbt%no) )
   ALLOCATE ( xbt%pq5(xbt%no), xbt%pq6(xbt%no), xbt%pq7(xbt%no), xbt%pq8(xbt%no) )
   ALLOCATE ( xbt%rss(xbt%no), xbt%ins(xbt%no) )
   ALLOCATE ( xbt%fls(xbt%no), xbt%eve(xbt%no) )

   xbt%flc(:) = 0_i8
   xbt%inc(:) = 0.0_r8
   xbt%bia(:) = 0.0_r8
   xbt%b_a(:) = 0.0_r8
   xbt%pq1(:) = 0.0_r8
   xbt%pq2(:) = 0.0_r8
   xbt%pq3(:) = 0.0_r8
   xbt%pq4(:) = 0.0_r8
   xbt%pq5(:) = 0.0_r8
   xbt%pq6(:) = 0.0_r8
   xbt%pq7(:) = 0.0_r8
   xbt%pq8(:) = 0.0_r8
   xbt%rss(:) = 0.0_r8
   xbt%ins(:) = 0.0_r8
   xbt%fls(:) = 0_i8
   xbt%eve(:) = 999_i8


   READ (511)                                               &
       xbt%ino(1:xbt%no), xbt%flg(1:xbt%no), xbt%par(1:xbt%no) &
      ,xbt%lon(1:xbt%no), xbt%lat(1:xbt%no)                    &
      ,xbt%dpt(1:xbt%no), xbt%tim(1:xbt%no)                    &
      ,xbt%val(1:xbt%no), xbt%bac(1:xbt%no)                    &
      ,xbt%err(1:xbt%no), xbt%res(1:xbt%no)                    &
      ,xbt%ib(1:xbt%no), xbt%jb(1:xbt%no), xbt%kb(1:xbt%no)    &
      ,xbt%pb(1:xbt%no), xbt%qb(1:xbt%no), xbt%rb(1:xbt%no)
   CLOSE  (511)

! ---
! Initialise quality flag
   IF ( obs%xbt .EQ. 0 ) THEN
      xbt%flg(:) = -1_i8
      WRITE (drv%dia,*)'Bad quality flag ',obs%xbt
   ENDIF

! ---
! Vertical interpolation parameters
   DO k = 1,xbt%no
      IF ( xbt%flg(k) .EQ. 1 ) THEN
         xbt%kb(k) = grd%km-1
         DO kk = 1,grd%km-1
            IF ( xbt%dpt(k) .GE. grd%dep(kk) .AND. xbt%dpt(k) .LT. grd%dep(kk+1) ) THEN
               xbt%kb(k) = kk
               xbt%rb(k) = (xbt%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
            ENDIF
         ENDDO
      ENDIF
   ENDDO

! Count good observations
   xbt%nc = 0
   DO k = 1,xbt%no
      IF ( xbt%flg(k) .EQ. 1 ) THEN
         xbt%nc = xbt%nc + 1
      ELSE
         xbt%bia(k) = 0._r8
         xbt%res(k) = 0._r8
         xbt%inc(k) = 0._r8
         xbt%b_a(k) = 0._r8
         xbt%pq1(k) = 0._r8
         xbt%pq2(k) = 0._r8
         xbt%pq3(k) = 0._r8
         xbt%pq4(k) = 0._r8
         xbt%pq5(k) = 0._r8
         xbt%pq6(k) = 0._r8
         xbt%pq7(k) = 0._r8
         xbt%pq8(k) = 0._r8
      ENDIF
   ENDDO

   xbt%flc(:) = xbt%flg(:)
   WHERE ( xbt%flc .EQ. 0 ) xbt%eve = 1

   WRITE (drv%dia,*) 'Number of good XBT observations after reading: ',  xbt%nc

1111 CONTINUE

END SUBROUTINE get_obs_xbt
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                              
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_xbt

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k, k1
   INTEGER(i4)   ::  i1, j1, i, idep
   INTEGER(i8)   ::  klev
   REAL(r8)      ::  p1, q1, r1
   REAL(r8)      ::  msk4, div_x, div_y

   IF ( xbt%no .GT. 0 ) THEN

      xbt%flc(:) = xbt%flg(:)

! ---
! Horizontal interpolation PARAMETERs
      CALL int_obs_hor ( xbt%no, xbt%lat, xbt%lon, xbt%flc, xbt%eve, xbt%ib, xbt%jb, xbt%pb, xbt%qb)

! ---
! Undefine masked for multigrid
      DO k = 1,xbt%no
         IF ( xbt%flc(k) .EQ. 1 ) THEN
            i1 = xbt%ib(k)
            j1 = xbt%jb(k)
            idep = xbt%kb(k)+1
            msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
            IF (msk4 .LT. 1.0) THEN
               xbt%flc(k) = 0
               xbt%eve(k) = 3
            ENDIF
         ENDIF
      ENDDO

! ---
! Horizontal interpolation PARAMETERs for each masked grid
      DO k = 1,xbt%no

         IF ( xbt%flc(k) .EQ. 1 ) THEN
            klev = xbt%kb(k)
            CALL int_obs_pq( xbt%ib(k), xbt%jb(k), klev, xbt%pb(k), xbt%qb(k),  &
                             xbt%pq1(k), xbt%pq2(k), xbt%pq3(k), xbt%pq4(k))
            klev = xbt%kb(k) + 1
            CALL int_obs_pq( xbt%ib(k), xbt%jb(k), klev, xbt%pb(k), xbt%qb(k),  &
                            xbt%pq5(k), xbt%pq6(k), xbt%pq7(k), xbt%pq8(k))

            r1=xbt%rb(k)
            xbt%pq1(k) = (1.-r1) * xbt%pq1(k)
            xbt%pq2(k) = (1.-r1) * xbt%pq2(k)
            xbt%pq3(k) = (1.-r1) * xbt%pq3(k)
            xbt%pq4(k) = (1.-r1) * xbt%pq4(k)
            xbt%pq5(k) =     r1  * xbt%pq5(k)
            xbt%pq6(k) =     r1  * xbt%pq6(k)
            xbt%pq7(k) =     r1  * xbt%pq7(k)
            xbt%pq8(k) =     r1  * xbt%pq8(k)
         ENDIF

      ENDDO

! ---
! Count good observations
      xbt%nc = 0
      DO k = 1,xbt%no
         IF ( xbt%flc(k) .EQ. 1 ) THEN
            xbt%nc = xbt%nc + 1
         ENDIF
      ENDDO

      WRITE (drv%dia,*) 'Number of good XBT observations after parameters intrpolation: ',  xbt%nc

      xbt%inc(:) = 0.0_r8

   ENDIF

END SUBROUTINE int_par_xbt
