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
!> Load observations of glider velocities                               
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!            Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_gvl

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : gvl,  obs

   IMPLICIT NONE

   INTEGER(i4)   ::  k
   INTEGER(i4)   ::  i1, kk, i
   REAL(r8)      ::  zbo, zbn

   gvl%no = 0_i8
   gvl%nc = 0_i8

   OPEN (511,FILE=drv%inpdir//'/gvl_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ(511) gvl%no

   WRITE (drv%dia,*) 'Number of velocity observations by gliders: ',gvl%no

   IF (gvl%no.EQ.0) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( gvl%ino(gvl%no), gvl%flg(gvl%no), gvl%flc(gvl%no), gvl%par(gvl%no) )
   ALLOCATE ( gvl%lon(gvl%no), gvl%lat(gvl%no), gvl%dpt(gvl%no), gvl%kdp(gvl%no) )
   ALLOCATE ( gvl%tim(gvl%no), gvl%tms(gvl%no), gvl%tme(gvl%no) )
   ALLOCATE ( gvl%val(gvl%no), gvl%bac(gvl%no), gvl%inc(gvl%no) )
   ALLOCATE ( gvl%bia(gvl%no), gvl%err(gvl%no) )
   ALLOCATE ( gvl%res(gvl%no), gvl%b_a(gvl%no) )
   ALLOCATE ( gvl%ib(gvl%no), gvl%jb(gvl%no), gvl%kb(gvl%no) )
   ALLOCATE ( gvl%pb(gvl%no), gvl%qb(gvl%no), gvl%rb(gvl%no) )
   ALLOCATE ( gvl%pq1(gvl%no), gvl%pq2(gvl%no), gvl%pq3(gvl%no), gvl%pq4(gvl%no) )
   ALLOCATE ( gvl%pq5(gvl%no), gvl%pq6(gvl%no), gvl%pq7(gvl%no), gvl%pq8(gvl%no) )
   ALLOCATE ( gvl%nav(gvl%no) )
   ALLOCATE ( gvl%dzr(grd%km,gvl%no) )
   ALLOCATE ( gvl%rss(gvl%no), gvl%ins(gvl%no) )
   ALLOCATE ( gvl%fls(gvl%no),gvl%eve(gvl%no) )

   gvl%rss(:) = 0.0_r8
   gvl%ins(:) = 0.0_r8
   gvl%fls(:) = 0_i8

   gvl%bia(:) = 0.0_r8
   gvl%b_a(:) = 0.0_r8
   gvl%eve(:) = 999_i8

   READ(511)                                               &
       gvl%ino(1:gvl%no), gvl%flg(1:gvl%no), gvl%par(1:gvl%no) &
      ,gvl%lon(1:gvl%no), gvl%lat(1:gvl%no), gvl%dpt(1:gvl%no) &
      ,gvl%tim(1:gvl%no), gvl%tms(1:gvl%no), gvl%tme(1:gvl%no) &
      ,gvl%val(1:gvl%no), gvl%bac(1:gvl%no), gvl%res(1:gvl%no) &
      ,gvl%err(1:gvl%no), gvl%nav(1:gvl%no)                    &
      ,gvl%ib(1:gvl%no), gvl%jb(1:gvl%no)                      &
      ,gvl%pb(1:gvl%no), gvl%qb(1:gvl%no)

   CLOSE (511)

! ---
! Initialise quality flag
   IF ( obs%gvl .EQ. 0 ) THEN
      gvl%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%gvl
   ENDIF


! ---
! Vertical interpolation parameters
   DO k = 1,gvl%no

      IF ( gvl%flg(k) .EQ. 1 ) THEN
         gvl%kb(k) = grd%km-1
         zbo = 0.0_r8
         DO kk = 1,grd%km-1
            zbn = grd%dep(kk)+grd%dz(kk)*0.5
            IF ( gvl%dpt(k) .GT. zbn) THEN
               gvl%dzr(kk,k) = grd%dz(kk)/gvl%dpt(k)
            ELSEIF ( gvl%dpt(k) .GE. zbo) THEN
               gvl%kb(k) = kk
               gvl%dzr(kk,k) = (gvl%dpt(k)-zbo)/gvl%dpt(k)
            ELSE
               gvl%dzr(kk,k) = 0.0_r8
            ENDIF
            zbo=zbn
         ENDDO
      ENDIF

   ENDDO

! ---
! Count good observations
   gvl%nc = 0
   DO k= 1,gvl%no
      IF ( gvl%flg(k) .EQ. 1 ) THEN
         gvl%nc = gvl%nc + 1
      ELSE
         gvl%bia(k) = 0.0_r8
         gvl%res(k) = 0.0_r8
         gvl%inc(k) = 0.0_r8
         gvl%b_a(k) = 0.0_r8
         gvl%pq1(k) = 0.0_r8
         gvl%pq2(k) = 0.0_r8
         gvl%pq3(k) = 0.0_r8
         gvl%pq4(k) = 0.0_r8
         gvl%pq5(k) = 0.0_r8
         gvl%pq6(k) = 0.0_r8
         gvl%pq7(k) = 0.0_r8
         gvl%pq8(k) = 0.0_r8
      ENDIF
   ENDDO

   gvl%flc(:) = gvl%flg(:)
   WHERE ( gvl%flc .EQ. 0 ) gvl%eve = 1

   WRITE (drv%dia,*) 'Number of good glider velocity observations after reading: ',  gvl%nc

1111 CONTINUE

END SUBROUTINE get_obs_gvl
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
SUBROUTINE int_par_gvl

   USE set_knd
   USE grd_str
   USE obs_str
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k
   INTEGER(i4)   ::  i1, j1, k1, i, idep
   REAL(r8)      ::  p1, q1
   REAL(r8)      ::  msk4u, msk4v, r1, div_x, div_y
   LOGICAL       ::  ins

   ins(i,i1) = i.GE.1 .AND. i.LT.i1

   IF (gvl%no.GT.0) THEN

! ---
! Horizontal interpolation PARAMETERs
      DO k = 1,gvl%no
         IF ( gvl%par(k) .EQ. 2 ) THEN
            q1 = (gvl%lat(k) - grd%lat(1,1)-(grd%lat(1,2)-grd%lat(1,1))*0.5) /    &
                 (grd%lat(1,2)-grd%lat(1,1)) + 1.0
         ELSE
            q1 = (gvl%lat(k) - grd%lat(1,1)) / (grd%lat(1,2)-grd%lat(1,1)) + 1.0
         ENDIF
         j1 = INT(q1)
         IF ( gvl%par(k) .EQ. 1 ) THEN
            p1 = (gvl%lon(k) - grd%lon(1,1)-(grd%lon(2,1)-grd%lon(1,1))*0.5) /    &
                 (grd%lon(2,1)-grd%lon(1,1)) + 1.0
         ELSE
            p1 = (gvl%lon(k) - grd%lon(1,1)) / (grd%lon(2,1)-grd%lon(1,1)) + 1.0
         ENDIF
         j1 = INT(q1)
         i1 = INT(p1)
         IF ( ins(j1,grd%jm+grd%jae) .AND. ins(i1,grd%im+grd%iae) ) THEN
            gvl%ib(k) = i1
            gvl%jb(k) = j1
            gvl%pb(k) = (p1-i1)
            gvl%qb(k) = (q1-j1)
         ELSE
            gvl%flc(k) = 0
            gvl%eve(k) = 2
         ENDIF
      ENDDO

! ---
! Undefine masked
      DO k = 1,gvl%no
         IF ( gvl%flc(k) .EQ. 1 ) THEN
            i1 = gvl%ib(k)
            j1 = gvl%jb(k)
            idep = gvl%kb(k)+1
            msk4u = grd%ums(i1,j1  ,idep) + grd%ums(i1+1,j1  ,idep) +      &
                    grd%ums(i1,j1+1,idep) + grd%ums(i1+1,j1+1,idep)
            msk4v = grd%vms(i1,j1  ,idep) + grd%vms(i1+1,j1  ,idep) +      &
                    grd%vms(i1,j1+1,idep) + grd%vms(i1+1,j1+1,idep)
            IF ( msk4u .LT. 4 .OR. msk4v .LT. 4) THEN
               gvl%flc(k) = 0
               gvl%eve(k) = 3
            ENDIF
         ENDIF
      ENDDO

! Horizontal interpolation parameters for each masked grid
      DO k = 1,gvl%no
         IF ( gvl%flc(k) .EQ. 1 .AND. gvl%par(k) .EQ. 1 ) THEN

            i1 = gvl%ib(k)
            p1 = gvl%pb(k)
            j1 = gvl%jb(k)
            q1 = gvl%qb(k)

            k1 = gvl%kb(k)+1
            div_y =  (1.-q1) * MAX(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
            gvl%pq1(k) = grd%ums(i1,j1,k1)                                        &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            gvl%pq2(k) = grd%ums(i1+1,j1,k1)                                      &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
            gvl%pq3(k) = grd%ums(i1,j1+1,k1)                                      &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            gvl%pq4(k) = grd%ums(i1+1,j1+1,k1)                                    &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

         ELSEIF ( gvl%flc(k) .EQ. 1 .AND. gvl%par(k) .EQ. 2 ) THEN

            i1 = gvl%ib(k)
            p1 = gvl%pb(k)
            j1 = gvl%jb(k)
            q1 = gvl%qb(k)

            k1 = gvl%kb(k)+1
            div_y =  (1.-q1) * MAX(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
            gvl%pq1(k) = grd%vms(i1,j1,k1)                                        &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            gvl%pq2(k) = grd%vms(i1+1,j1,k1)                                      &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
            gvl%pq3(k) = grd%vms(i1,j1+1,k1)                                      &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            gvl%pq4(k) = grd%vms(i1+1,j1+1,k1)                                    &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

         ENDIF

      ENDDO


! ---
! Count good observations
      gvl%nc = 0_i8
      DO k = 1,gvl%no
         IF ( gvl%flc(k) .EQ. 1 ) THEN
            gvl%nc = gvl%nc + 1
         ENDIF
      ENDDO

      gvl%inc(:) = 0.0_r8

      WRITE (drv%dia,*) 'Number of good glider velocity observations after parameter interpolation: ',  gvl%nc

   ENDIF

END SUBROUTINE int_par_gvl
