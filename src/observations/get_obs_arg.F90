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
!> Load ARGO observations                                               
!!
!!
!!
!                                                                      !
! Version 1:   Srdjan Dobricic 2006                                    !
!              Mario Adani     2023                                    !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_arg

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : arg,  obs

   IMPLICIT NONE

   INTEGER(i4)   ::  k, ierr
   INTEGER(i4)   ::  i1, kk, i
   REAL          ::  cnti

   arg%no = 0_i8
   arg%nc = 0_i8


   OPEN (511,FILE=drv%inpdir//'/arg_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) arg%no

   WRITE (drv%dia,*)' --- No of ARGO obs: ',arg%no

   IF ( arg%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
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
   ALLOCATE ( arg%fls(arg%no), arg%eve(arg%no))

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
      arg%ino(1:arg%no), arg%flg(1:arg%no), arg%par(1:arg%no) &
      ,arg%lon(1:arg%no), arg%lat(1:arg%no)                    &
      ,arg%dpt(1:arg%no), arg%tim(1:arg%no)                    &
      ,arg%val(1:arg%no), arg%bac(1:arg%no)                    &
      ,arg%err(1:arg%no), arg%res(1:arg%no)                    &
      ,arg%ib(1:arg%no), arg%jb(1:arg%no), arg%kb(1:arg%no)    &
      ,arg%pb(1:arg%no), arg%qb(1:arg%no), arg%rb(1:arg%no)
   CLOSE (511)

   arg%nc = 0
   DO k=1,arg%no
      IF (arg%flg(k).EQ.1 )arg%nc = arg%nc+1
   ENDDO

   WRITE (drv%dia,*) 'Number of good ARGO observations after reading: ',  arg%nc

! ---
! Initialise quality flag
   IF  ( obs%arg .EQ. 0 ) THEN
      arg%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%arg
   ENDIF

! ---
! Vertical interpolation PARAMETERs
   DO k = 1,arg%no
      IF ( arg%flg(k) .EQ. 1 ) THEN
         arg%kb(k) = grd%km-1
         DO kk = 1,grd%km-1
            IF ( arg%dpt(k) .GE. grd%dep(kk) .AND. arg%dpt(k) .LT. grd%dep(kk+1) ) THEN
               arg%kb(k) = kk
               arg%rb(k) = (arg%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
            ENDIF
         ENDDO
      ENDIF
   ENDDO

! ---
! Count good observations
   arg%nc = 0
   DO k = 1,arg%no
      IF ( arg%flg(k) .EQ. 1 ) THEN
         arg%nc = arg%nc + 1
      ELSE
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
      ENDIF
   ENDDO

   arg%flc(:) = arg%flg(:)
   WHERE ( arg%flc .EQ. 0 ) arg%eve = 1

1111 CONTINUE

END SUBROUTINE get_obs_arg
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                            
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!            Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_arg

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, k,kk, ierr
   INTEGER(i4)   ::  i1, j1, k1, idep
   INTEGER(i8)   ::  klev
   REAL(r8)      ::  p1, q1, r1
   REAL(r8)      ::  msk4, div_x, div_y, rmn

   rmn = 1.e-6

   IF ( arg%no .GT. 0 ) THEN

      arg%flc(:) = arg%flg(:)
! ---
! Adjust longitudes
      IF ( grd%bwst .GT. 180. ) THEN
         DO k = 1,arg%no
            IF ( arg%lon(k) .LT. 0.0_r8 ) THEN
               arg%lon(k) = arg%lon(k) + 360.
            ENDIF
         ENDDO
      ENDIF
! ---
! Horizontal interpolation PARAMETERs
      CALL int_obs_hor ( arg%no, arg%lat, arg%lon, arg%flc, arg%eve, arg%ib, arg%jb, arg%pb, arg%qb)
! ---
! Undefine masked for multigrid
      DO k = 1,arg%no
         IF ( arg%flc(k) .EQ. 1 ) THEN
            i1 = arg%ib(k)
            j1 = arg%jb(k)
            idep = arg%kb(k)+1
            msk4 = grd%msk(i1,j1,idep) + grd%msk(i1+1,j1,idep) + grd%msk(i1,j1+1,idep) + grd%msk(i1+1,j1+1,idep)
            IF ( msk4 .LT. 1. ) THEN
               arg%flc(k) = 0
               arg%eve(k) = 3
            ENDIF
         ENDIF
      ENDDO

! ---
! Horizontal interpolation parameters for each masked grid
      DO k = 1,arg%no
         IF ( arg%flc(k) .EQ. 1 ) THEN

            klev = arg%kb(k)
            CALL int_obs_pq( arg%ib(k), arg%jb(k), klev, arg%pb(k), arg%qb(k),  &
                             arg%pq1(k), arg%pq2(k), arg%pq3(k), arg%pq4(k))
            klev = arg%kb(k) + 1
            CALL int_obs_pq( arg%ib(k), arg%jb(k), klev, arg%pb(k), arg%qb(k),  &
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
         ELSE
            arg%pq1(k) = 0.0_r8
            arg%pq2(k) = 0.0_r8
            arg%pq3(k) = 0.0_r8
            arg%pq4(k) = 0.0_r8
            arg%pq5(k) = 0.0_r8
            arg%pq6(k) = 0.0_r8
            arg%pq7(k) = 0.0_r8
            arg%pq8(k) = 0.0_r8
         ENDIF
      ENDDO

! ---
! Count good observations
      arg%nc = 0
      DO k = 1,arg%no
         IF ( arg%flc(k) .EQ. 1 ) THEN
            arg%nc = arg%nc + 1
         ENDIF
      ENDDO

      WRITE (drv%dia,*) 'Number of good ARGO observations after parametr interpolation: ',  arg%nc

      arg%inc(:) = 0.0_r8

   ENDIF

END SUBROUTINE int_par_arg
