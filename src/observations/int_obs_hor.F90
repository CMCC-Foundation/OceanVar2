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
!> Get interpolation parameters for a grid                              
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
! Version 2: Mario Adani     2023                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_obs_hor ( no, olat, olon, flc, eve, ib, jb, pb, qb)

   USE set_knd
   USE drv_str
   USE grd_str

   IMPLICIT NONE


   INTEGER(i8)               ::  no
   REAL(r8)                  ::  olat(no), olon(no), pb(no), qb(no)
   INTEGER(i8)               ::  flc(no), ib(no), jb(no), eve(no)
   INTEGER(i4)               ::  k, ierr, kks, kkp
   INTEGER(i4)               ::  i1, kk, i, j1, j, IF1, jf1
   REAL(r8)                  ::  p1, q1, pf1, qf1
   REAL(r8)                  ::  msk4, div_x, div_y, rmn, dst, dstm
   REAL(r8)                  ::  tga, ang, lat_rot, lon_rot, lat_lb_rot, lon_lb_rot
   REAL(r8)                  ::  lat_lt_rot, lon_rb_rot
   REAL(r8)                  ::  lonv(4),latv(4),beas,bwst
   logical                   :: results,llinvalidcell
   !FUNCTION
   logical                   :: linquad
   REAL(r8),ALLOCATABLE      ::  lon(:,:),lat(:,:)
   REAL(r8),ALLOCATABLE      ::  plon(:),plat(:)


   rmn = 1.e-6

   IF ( grd%prj .EQ. 0 ) THEN
! Equally distant lat-lon grid
      DO kk = 1,no
         q1 = (olat(kk) - grd%lat1_1) / grd%dlt + 1.0
         j1 = INT(q1)
         p1 = (olon(kk) - grd%lon1_1) / grd%dln + 1.0
         i1 = INT(p1)
         IF ( j1 .GE. grd%jgs .AND. j1 .LT. grd%jge+grd%jae .AND. &
              i1 .GE. grd%igs .AND. i1 .LT. grd%ige+grd%iae) THEN
            IF ( flc(kk) .EQ. 1 ) THEN
               ib(kk) = i1-grd%igs+1
               jb(kk) = j1-grd%jgs+1
               pb(kk) = MAX((p1-i1),rmn)
               qb(kk) = MAX((q1-j1),rmn)
            ENDIF
         ELSE
            IF ( flc(kk) .EQ. 1 ) THEN
               eve(kk) = 2
               flc(kk) = 0
            ENDIF
         ENDIF
      ENDDO

   ELSE

      ib(:) = -999
      jb(:) = -999
      pb(:) = -999.9
      qb(:) = -999.9

!-----------------------------------------------------------------------
! Copy grid positions to temporary arrays and normalization 0 to 360
! done in def_grd.F90 to 0 to 360.
!-----------------------------------------------------------------------
      ALLOCATE ( lon (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
      ALLOCATE ( lat (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
      lon(:,:) = grd%lon
      lat(:,:) = grd%lat
      WHERE ( lon(:,:) .LT. 0.0_r8 )
         lon(:,:) = lon(:,:) + 360.0_r8
      END WHERE
      WHERE ( lon(:,:) .GT. 360.0_r8 )
         lon(:,:) = lon(:,:) - 360.0_r8
      END WHERE

!-----------------------------------------------------------------------
! Copy observation positions to temporary arrays and renormalize to 0 to 360.
! It may be not necdssary becaUSE the check is in get_xxx.F90 but for the
! time being we leave it also here.
!-----------------------------------------------------------------------
      ALLOCATE ( plon(no),plat(no) )
      plon(:) = olon(:)
      plat(:) = olat(:)
      WHERE ( plon(:) .LT. 0.0_r8   )
         plon(:) = plon(:) + 360.0_r8
      END WHERE
      WHERE ( plon(:) .GT. 360.0_r8   )
         plon(:) = plon(:) - 360.0_r8
      END WHERE

! Add grid resolution to handle ambiguity of 360
      beas = MAXVAL(MAXVAL(lon,1)) + grd%dln
      bwst = MINVAL(MINVAL(lon,1)) - grd%dln
      DO kk = 1,no
         IF ( plat(kk) .GE. grd%bnrt .OR. plat(kk) .LT. grd%bsth .OR.   &
            plon(kk) .GE. beas .OR. plon(kk) .LT. bwst ) THEN
            flc(kk) = 0
            eve(kk) = 2
         ENDIF
      ENDDO

!------------------------------------------------------------------------
! Master loop for grid search brute force method
!------------------------------------------------------------------------
      k   = 1 ! first level
      kkp = 1
      kks = 0
      obsloop:     DO kk = 1, no
         ! valid obs  ?
         IF ( flc(kk).EQ. 1 )  THEN
            ! same position of observation copy same values and THEN cycle
            IF ( flc(kkp) .EQ. 1 .AND. kks .GT. 0 .AND.          &
               plat(kk) .EQ. plat(kkp) .AND. plon(kk) .EQ. plon(kkp)) THEN
               ib(kk) = ib(kkp)
               jb(kk) = jb(kkp)
               pb(kk) = pb(kkp)
               qb(kk) = qb(kkp)
            ELSE
!------------------------------------------------------------------------
! START Master loop for grid search brute force method
!------------------------------------------------------------------------
               gridloop:   DO j = 1,grd%jm-1
                  DO i = 1,grd%im-1
                     lonv(1)=lon(i  ,j  )
                     lonv(2)=lon(i+1,j  )
                     lonv(3)=lon(i+1,j+1)
                     lonv(4)=lon(i  ,j+1)
                     latv(1)=lat(i  ,j  )
                     latv(2)=lat(i+1,j  )
                     latv(3)=lat(i+1,j+1)
                     latv(4)=lat(i  ,j+1)
                     llinvalidcell = grd%msk(i  ,j  ,k) .EQ. 0.0_r8 .AND. &
                                     grd%msk(i+1,j  ,k) .EQ. 0.0_r8 .AND. &
                                     grd%msk(i+1,j+1,k) .EQ. 0.0_r8 .AND. &
                                     grd%msk(i  ,j+1,k) .EQ. 0.0_r8
!---------------------------------------------------------------------
! Find observations which are on within 1e-6 of a grid point
!---------------------------------------------------------------------
                     IF ( (ABS(lon(i,j) - plon(kk)) < 1e-6 ) .AND. &
                        (ABS(lat(i,j) - plat(kk)) < 1e-6 ) ) THEN
                        IF ( llinvalidcell ) THEN
                           ! Land cycle
                           flc(kk) = 0
                           eve(kk) = 3
                           EXIT gridloop
                        ENDIF
                        ib(kk) = i
                        jb(kk) = j
                        pb(kk) = 1.0_r8
                        qb(kk) = 1.0_r8
                        EXIT gridloop
                     ENDIF
                     results = linquad(plon(kk),plat(kk),lonv,latv)
                     IF ( results ) THEN
                        IF ( llinvalidcell ) THEN
                           ! Land cycle
                           flc(kk) = 0
                           eve(kk) = 3
                           EXIT gridloop
                        ENDIF
                        ib(kk) = i
                        jb(kk) = j
                        tga = ((lon(i,j+1)-lon(i,j))/(lat(i,j+1)-lat(i,j)))
                        ang = ATAN (tga)

                        lon_lb_rot = lon(i,j)*COS(ang) - lat(i,j)*SIN(ang)
                        lat_lb_rot = lon(i,j)*SIN(ang) + lat(i,j)*COS(ang)
                        lon_rb_rot = lon(i+1,j)*COS(ang) - lat(i+1,j)*SIN(ang)
                        lat_lt_rot = lon(i,j+1)*SIN(ang) + lat(i,j+1)*COS(ang)
                        lon_rot = plon(kk) * COS(ang) - plat(kk) * SIN(ang)
                        lat_rot = plon(kk) * SIN(ang) + plat(kk) * COS(ang)
                        qb(kk) = (lat_rot - lat_lb_rot) / (lat_lt_rot - lat_lb_rot)
                        pb(kk) = (lon_rot - lon_lb_rot) / (lon_rb_rot - lon_lb_rot)

                        EXIT gridloop
                     ENDIF

                  ENDDO
               ENDDO  gridloop
!------------------------------------------------------------------------
! END Master loop for grid search brute force method
!------------------------------------------------------------------------

            ENDIF  ! same position of previous obs?
         ENDIF  ! valid flag
         kkp = kk
         kks = kk
      ENDDO  obsloop
!----------------------
! ADANI WARNING!!!!!!!!
! at the border of each subdomain can not find the grid cell
! for now we must procede
!---------------------
      WHERE ( ib(:) < 0 .OR. jb(:) < 0 )
         flc(:) = 0
         eve(:) = 4
      END WHERE
! END ADANI WARNING!!!!!!!!
      DEALLOCATE(lon, lat, plon, plat)
   ENDIF   ! Rotated grid

CONTAINS
   PURE FUNCTION ins(i,j)
      LOGICAL                :: ins
      INTEGER, INTENT(in   ) :: i,j
      ins =  i.GE.1 .AND. i.LT.j
   END FUNCTION ins

END SUBROUTINE int_obs_hor
!-----------------------------------------------------------------------
!                                                                      !
!> Get interpolation parameters for a grid                             
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE int_obs_pq( i1, j1, k1, p1, q1, pq1, pq2, pq3, pq4)

   USE set_knd
   USE grd_str

   IMPLICIT NONE

   INTEGER(i8)   ::  i1, j1, k1
   REAL(r8)      ::  p1, q1
   REAL(r8)      ::  pq1, pq2, pq3, pq4
   REAL(r8)      ::  div_x, div_y

   div_y =  (1.-q1) * MAX(grd%msk(i1,j1  ,k1),grd%msk(i1+1,j1  ,k1))     &
           +    q1  * MAX(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))
   div_x =  (1.-p1) * grd%msk(i1  ,j1,k1) + p1 * grd%msk(i1+1,j1,k1)
   pq1 = grd%msk(i1,j1,k1)                                      &
         * MAX(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))           &
         * (1.-p1) * (1.-q1)                                    &
         /( div_x * div_y + 1.e-16_r8 )
   pq2 = grd%msk(i1+1,j1,k1)                                    &
         * MAX(grd%msk(i1,j1,k1),grd%msk(i1+1,j1,k1))           &
         *     p1  * (1.-q1)                                    &
         /( div_x * div_y + 1.e-16_r8 )
   div_x =  (1.-p1) * grd%msk(i1  ,j1+1,k1) + p1 * grd%msk(i1+1,j1+1,k1)
   pq3 = grd%msk(i1,j1+1,k1)                                    &
         * MAX(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))       &
         * (1.-p1) *     q1                                     &
         /( div_x * div_y + 1.e-16_r8 )
   pq4 = grd%msk(i1+1,j1+1,k1)                                  &
         * MAX(grd%msk(i1,j1+1,k1),grd%msk(i1+1,j1+1,k1))       &
         *     p1  *     q1                                     &
         /( div_x * div_y + 1.e-16_r8 )

END SUBROUTINE int_obs_pq
!!----------------------------------------------------------------------
!>                    ***  function linquad ***
!!
!! ** Purpose : Determine whether a point P(x,y) lies within or on the
!!              boundary of a quadrangle (ABCD) of any shape on a plane.
!!
!! ** Method  : Check if the vectorial products PA x PC, PB x PA,
!!              PC x PD, and PD x PB are all negative.
!
! ** Action  :
!
! History :
!        !  2001-11  (N. Daget, A. Weaver)
!        !  2006-08  (A. Weaver) NEMOVAR migration
!        !  2006-10  (A. Weaver) Cleanup
!!----------------------------------------------------------------------
LOGICAL FUNCTION linquad( px, py, pxv, pyv )

   USE set_knd
   IMPLICIT NONE

   !! * Arguments
   REAL(r8), INTENT(IN) :: px        ! (lon) of the point P(x,y)
   REAL(r8), INTENT(IN) :: py        ! (lat) of the point P(x,y)
   REAL(r8), DIMENSION(4), INTENT(IN) :: &
   & pxv,  &                  ! (lon, lat) of the surrounding cell
   & pyv

   !! * Local declarations
   REAL(r8) :: zst1
   REAL(r8) :: zst2
   REAL(r8) :: zst3
   REAL(r8) :: zst4

   !-----------------------------------------------------------------------
   ! Test to see IF the point is within the cell
   !-----------------------------------------------------------------------
   linquad = .FALSE.
   zst1 =   ( px - pxv(1) ) * ( py - pyv(4) ) &
   &   - ( py - pyv(1) ) * ( px - pxv(4) )
   IF ( zst1 <= 0.0_r8 ) THEN
      zst2 =   ( px - pxv(4) ) * ( py - pyv(3) ) &
      &   - ( py - pyv(4) ) * ( px - pxv(3) )
      IF ( zst2 <= 0.0_r8 ) THEN
         zst3 =   ( px - pxv(3) ) * ( py - pyv(2) ) &
         &   - ( py - pyv(3) ) * ( px - pxv(2) )
         IF ( zst3 <= 0.0_r8) THEN
            zst4 =   ( px - pxv(2) ) * ( py - pyv(1) ) &
            &   - ( py - pyv(2) ) * ( px - pxv(1) )
            IF ( zst4 <= 0.0_r8 ) linquad = .TRUE.
         ENDIF
      ENDIF
   ENDIF

END FUNCTION linquad

