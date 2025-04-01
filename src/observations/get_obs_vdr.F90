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
!> Load observations of drifter velocities                             
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2007                                      !
!          : Mario Adani     2024 Add events                           !
!-----------------------------------------------------------------------
SUBROUTINE get_obs_vdr

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : vdr,  obs

   IMPLICIT NONE

   INTEGER(i4)   ::  k
   INTEGER(i4)   ::  i1, kk, i

   vdr%no = 0_i8
   vdr%nc = 0_i8


   OPEN (511,FILE=drv%inpdir//'/vdr_mis.dat',FORM='unformatted',STATUS='old',ERR=1111)

   READ (511) vdr%no

   WRITE (drv%dia,*) 'Number of velocity observations by gliders: ',vdr%no

   IF ( vdr%no .EQ. 0 ) THEN
      CLOSE (511)
      RETURN
   ENDIF

! ---
! Allocate memory for observations
   ALLOCATE ( vdr%ino(vdr%no), vdr%flg(vdr%no), vdr%flc(vdr%no), vdr%par(vdr%no) )
   ALLOCATE ( vdr%lon(vdr%no), vdr%lat(vdr%no), vdr%dpt(vdr%no), vdr%kdp(vdr%no) )
   ALLOCATE ( vdr%tim(vdr%no), vdr%tms(vdr%no), vdr%tme(vdr%no) )
   ALLOCATE ( vdr%val(vdr%no), vdr%bac(vdr%no), vdr%inc(vdr%no) )
   ALLOCATE ( vdr%bia(vdr%no), vdr%err(vdr%no) )
   ALLOCATE ( vdr%res(vdr%no), vdr%b_a(vdr%no) )
   ALLOCATE ( vdr%ib(vdr%no), vdr%jb(vdr%no), vdr%kb(vdr%no) )
   ALLOCATE ( vdr%pb(vdr%no), vdr%qb(vdr%no), vdr%rb(vdr%no) )
   ALLOCATE ( vdr%pq1(vdr%no), vdr%pq2(vdr%no), vdr%pq3(vdr%no), vdr%pq4(vdr%no) )
   ALLOCATE ( vdr%pq5(vdr%no), vdr%pq6(vdr%no), vdr%pq7(vdr%no), vdr%pq8(vdr%no) )
   ALLOCATE ( vdr%nav(vdr%no) )
   ALLOCATE ( vdr%rss(vdr%no), vdr%ins(vdr%no) )
   ALLOCATE ( vdr%fls(vdr%no), vdr%eve(vdr%no) )

   vdr%rss(:) = 0.0_r8
   vdr%ins(:) = 0.0_r8
   vdr%fls(:) = 0_i8
   vdr%bia(:) = 0.0_r8
   vdr%eve(:) = 999

   READ(511)                                               &
       vdr%ino(1:vdr%no), vdr%flg(1:vdr%no), vdr%par(1:vdr%no) &
      ,vdr%lon(1:vdr%no), vdr%lat(1:vdr%no), vdr%dpt(1:vdr%no) &
      ,vdr%tim(1:vdr%no), vdr%tms(1:vdr%no), vdr%tme(1:vdr%no) &
      ,vdr%val(1:vdr%no), vdr%bac(1:vdr%no)                    &
      ,vdr%err(1:vdr%no), vdr%nav(1:vdr%no)                    &
      ,vdr%ib(1:vdr%no), vdr%jb(1:vdr%no)                      &
      ,vdr%pb(1:vdr%no), vdr%qb(1:vdr%no)
   CLOSE (511)

! ---
! Initialise quality flag
   IF ( obs%vdr .EQ. 0 ) THEN
      vdr%flg(:) = -1
      WRITE (drv%dia,*)'Bad quality flag ',obs%vdr
   ENDIF

! ---
! Vertical interpolation PARAMETERs
   DO k = 1,vdr%no
      IF ( vdr%flg(k) .EQ. 1 ) THEN
         vdr%kb(k) = grd%km-1
         DO kk = 1,grd%km-1
            IF ( vdr%dpt(k).GE.grd%dep(kk) .AND. vdr%dpt(k).LT.grd%dep(kk+1) ) THEN
               vdr%kb(k) = kk
               vdr%rb(k) = (vdr%dpt(k) - grd%dep(kk)) / (grd%dep(kk+1) - grd%dep(kk))
            ENDIF
         ENDDO
      ENDIF
   ENDDO

! calculate residuals
   DO k = 1,vdr%no
      IF ( vdr%nav(k) .GT. 0 ) THEN
         vdr%bac(k) = vdr%bac(k)/vdr%nav(k)
         vdr%res(k) = vdr%val(k)-vdr%bac(k)
      ENDIF
   ENDDO

! ---
! Count good observations
   vdr%nc = 0
   DO k = 1,vdr%no
      IF ( vdr%flg(k) .EQ. 1 ) THEN
         vdr%nc = vdr%nc + 1
      ELSE
         vdr%bia(k) = 0.0_r8
         vdr%res(k) = 0.0_r8
         vdr%inc(k) = 0.0_r8
         vdr%b_a(k) = 0.0_r8
         vdr%pq1(k) = 0.0_r8
         vdr%pq2(k) = 0.0_r8
         vdr%pq3(k) = 0.0_r8
         vdr%pq4(k) = 0.0_r8
         vdr%pq5(k) = 0.0_r8
         vdr%pq6(k) = 0.0_r8
         vdr%pq7(k) = 0.0_r8
         vdr%pq8(k) = 0.0_r8
      ENDIF
   ENDDO

   vdr%flc(:) = vdr%flg(:)
   WHERE ( vdr%flc .EQ. 0 ) vdr%eve = 1

   WRITE (drv%dia,*) 'Number of good drifter velocity observations after reading: ',  vdr%nc

1111 CONTINUE

END SUBROUTINE get_obs_vdr
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                             
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_par_vdr

   USE set_knd
   USE grd_str
   USE obs_str
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)   ::  k
   INTEGER(i4)   ::  i1, j1, k1, i, idep
   REAL(r8)      ::  p1, q1
   REAL(r8)      ::  msk4, r1, div_x, div_y
   LOGICAL       ::  ins

   ins(i,i1) = i.GE.1 .AND. i.LT.i1

   IF ( vdr%no .GT. 0 ) THEN

      vdr%flc(:) = vdr%flg(:)

! ---
! Horizontal interpolation PARAMETERs
      DO k = 1,vdr%no
         IF ( vdr%par(k) .EQ. 2 ) THEN
            q1 = (vdr%lat(k) - grd%lat(1,1)-(grd%lat(1,2)-grd%lat(1,1))*0.5) /    &
                 (grd%lat(1,2)-grd%lat(1,1)) + 1.0
         ELSE
            q1 = (vdr%lat(k) - grd%lat(1,1)) / (grd%lat(1,2)-grd%lat(1,1)) + 1.0
         ENDIF
         j1 = int(q1)
         IF ( vdr%par(k) .EQ. 1 ) THEN
            p1 = (vdr%lon(k) - grd%lon(1,1)-(grd%lon(2,1)-grd%lon(1,1))*0.5) /    &
                 (grd%lon(2,1)-grd%lon(1,1)) + 1.0
         ELSE
            p1 = (vdr%lon(k) - grd%lon(1,1)) / (grd%lon(2,1)-grd%lon(1,1)) + 1.0
         ENDIF
         j1 = int(q1)
         i1 = int(p1)
         IF ( ins(j1,grd%jm+grd%jae) .AND. ins(i1,grd%im+grd%iae) ) THEN
            vdr%ib(k) = i1
            vdr%jb(k) = j1
            vdr%pb(k) = (p1-i1)
            vdr%qb(k) = (q1-j1)
         ELSE
            vdr%flc(k) = 0
            vdr%eve(k) = 2
         ENDIF
      ENDDO

! ---
! Undefine masked
      DO k = 1,vdr%no
         IF ( vdr%flc(k) .EQ. 1 ) THEN
            i1 = vdr%ib(k)
            j1 = vdr%jb(k)
            idep = vdr%kb(k)+1
            IF ( vdr%par(k) .EQ. 1 ) THEN
               msk4 = grd%ums(i1,j1  ,idep) + grd%ums(i1+1,j1  ,idep) +      &
                       grd%ums(i1,j1+1,idep) + grd%ums(i1+1,j1+1,idep)
            ENDIF
            IF ( vdr%par(k) .EQ. 2 ) THEN
               msk4 = grd%vms(i1,j1  ,idep) + grd%vms(i1+1,j1  ,idep) +      &
                      grd%vms(i1,j1+1,idep) + grd%vms(i1+1,j1+1,idep)
            ENDIF
            IF ( msk4 .LT. 2 ) THEN
               vdr%flc(k) = 0
               vdr%eve(k) = 3
            ENDIF
         ENDIF
      ENDDO

! Horizontal interpolation parameters for each masked grid
      DO k = 1,vdr%no
         IF ( vdr%flc(k) .EQ. 1 .AND. vdr%par(k) .EQ. 1 ) THEN

            i1 = vdr%ib(k)
            p1 = vdr%pb(k)
            j1 = vdr%jb(k)
            q1 = vdr%qb(k)
            r1 = vdr%rb(k)

            k1 = vdr%kb(k)
            div_y =  (1.-q1) * MAX(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
            vdr%pq1(k) = grd%ums(i1,j1,k1)                                        &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq2(k) = grd%ums(i1+1,j1,k1)                                      &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
            vdr%pq3(k) = grd%ums(i1,j1+1,k1)                                      &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq4(k) = grd%ums(i1+1,j1+1,k1)                                    &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

            k1 = vdr%kb(k) + 1
            div_y =  (1.-q1) * MAX(grd%ums(i1,j1  ,k1),grd%ums(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%ums(i1  ,j1,k1) + p1 * grd%ums(i1+1,j1,k1)
            vdr%pq5(k) = grd%ums(i1,j1,k1)                                        &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq6(k) = grd%ums(i1+1,j1,k1)                                      &
                         * MAX(grd%ums(i1,j1,k1),grd%ums(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%ums(i1  ,j1+1,k1) + p1 * grd%ums(i1+1,j1+1,k1)
            vdr%pq7(k) = grd%ums(i1,j1+1,k1)                                      &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq8(k) = grd%ums(i1+1,j1+1,k1)                                    &
                         * MAX(grd%ums(i1,j1+1,k1),grd%ums(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

            vdr%pq1(k) = (1.-r1) * vdr%pq1(k)
            vdr%pq2(k) = (1.-r1) * vdr%pq2(k)
            vdr%pq3(k) = (1.-r1) * vdr%pq3(k)
            vdr%pq4(k) = (1.-r1) * vdr%pq4(k)
            vdr%pq5(k) =     r1  * vdr%pq5(k)
            vdr%pq6(k) =     r1  * vdr%pq6(k)
            vdr%pq7(k) =     r1  * vdr%pq7(k)
            vdr%pq8(k) =     r1  * vdr%pq8(k)

         ELSEIF ( vdr%flc(k) .EQ. 1 .AND. vdr%par(k) .EQ. 2 ) THEN

            i1 = vdr%ib(k)
            p1 = vdr%pb(k)
            j1 = vdr%jb(k)
            q1 = vdr%qb(k)
            r1 = vdr%rb(k)

            k1 = vdr%kb(k)
            div_y =  (1.-q1) * MAX(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
            vdr%pq1(k) = grd%vms(i1,j1,k1)                                        &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq2(k) = grd%vms(i1+1,j1,k1)                                      &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
            vdr%pq3(k) = grd%vms(i1,j1+1,k1)                                      &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq4(k) = grd%vms(i1+1,j1+1,k1)                                    &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

            k1 = vdr%kb(k) + 1
            div_y =  (1.-q1) * MAX(grd%vms(i1,j1  ,k1),grd%vms(i1+1,j1  ,k1))     &
                    +    q1  * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))
            div_x =  (1.-p1) * grd%vms(i1  ,j1,k1) + p1 * grd%vms(i1+1,j1,k1)
            vdr%pq5(k) = grd%vms(i1,j1,k1)                                        &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         * (1.-p1) * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq6(k) = grd%vms(i1+1,j1,k1)                                      &
                         * MAX(grd%vms(i1,j1,k1),grd%vms(i1+1,j1,k1))             &
                         *     p1  * (1.-q1)                                      &
                         /( div_x * div_y + 1.e-16 )
            div_x =  (1.-p1) * grd%vms(i1  ,j1+1,k1) + p1 * grd%vms(i1+1,j1+1,k1)
            vdr%pq7(k) = grd%vms(i1,j1+1,k1)                                      &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         * (1.-p1) *     q1                                       &
                         /( div_x * div_y + 1.e-16 )
            vdr%pq8(k) = grd%vms(i1+1,j1+1,k1)                                    &
                         * MAX(grd%vms(i1,j1+1,k1),grd%vms(i1+1,j1+1,k1))         &
                         *     p1  *     q1                                       &
                         /( div_x * div_y + 1.e-16 )

            vdr%pq1(k) = (1.-r1) * vdr%pq1(k)
            vdr%pq2(k) = (1.-r1) * vdr%pq2(k)
            vdr%pq3(k) = (1.-r1) * vdr%pq3(k)
            vdr%pq4(k) = (1.-r1) * vdr%pq4(k)
            vdr%pq5(k) =     r1  * vdr%pq5(k)
            vdr%pq6(k) =     r1  * vdr%pq6(k)
            vdr%pq7(k) =     r1  * vdr%pq7(k)
            vdr%pq8(k) =     r1  * vdr%pq8(k)

         ENDIF

      ENDDO

! ---
! Count good observations
      vdr%nc = 0
      DO k = 1,vdr%no
         IF ( vdr%flc(k) .EQ. 1 ) THEN
            vdr%nc = vdr%nc + 1
         ENDIF
      ENDDO

      vdr%inc(:) = 0.0_r8

      WRITE (drv%dia,*) 'Number of good drifter velocity observations after parameter interpolation: ',  vdr%nc

   ENDIF

END SUBROUTINE int_par_vdr
