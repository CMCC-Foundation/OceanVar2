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
!> glider velocities quality control                                   
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab from namelist only
!! - Climatology qc abs[y-h(xclim)] > threshold,
!! - Too close tp the coast,
!! - If profile is only subsurface ,
!                                                                      !
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_gvl

   USE set_knd
   USE obs_str, ONLY : gvl, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)              :: kk, k, j, i
   REAL(r8),ALLOCATABLE :: treshold(:)
   REAL(r8)                 :: weightsum,pq1,pq2,pq3,pq4,val
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- GLIDER velocity observations:     ',gvl%nc

!Initialize
   ALLOCATE ( gvl%bgerr(gvl%no) )
   gvl%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( gvl%flc .EQ. 1 .AND. ABS(gvl%res) .GT. qck%res_vel )
         gvl%flc = 0
         gvl%eve = 5
      ENDWHERE
      CALL reset_gvl
      WRITE (drv%dia,*)' ---- after residual check:             ',gvl%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
! with constant background from NAMELIST
      ALLOCATE ( treshold(gvl%no) )
      ALLOCATE ( flc(gvl%no) )
      flc = gvl%flc
      gvl%bgerr(:) = qck%bgr_vel(1)
      treshold (:) = qck%bgr_vel(2)
      CALL qc_bgerr(gvl%no,gvl%res,gvl%err,gvl%bgerr,treshold,gvl%flc)
      CALL reset_gvl
      WHERE ( ( flc .EQ. 1 ) .AND. ( gvl%flc .EQ. 0 ) ) gvl%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',gvl%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Background error computed from EOFs not implemented yet.              '
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! IF greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Climatological check for glider velocity not implemented yet.               '
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%gvl ) THEN
      DO kk = 1,gvl%no
         IF ( gvl%flc(kk) .EQ. 1 ) THEN
            i = gvl%ib(kk)
            j = gvl%jb(kk)
            weightsum = gvl%pq1(kk) + gvl%pq2(kk) + gvl%pq3(kk) + gvl%pq4(kk)
            pq1 = gvl%pq1(kk) / weightsum
            pq2 = gvl%pq2(kk) / weightsum
            pq3 = gvl%pq3(kk) / weightsum
            pq4 = gvl%pq4(kk) / weightsum
            val = pq1 * coastrej%distc(i  ,j  ) +       &
                  pq2 * coastrej%distc(i+1,j  ) +       &
                  pq3 * coastrej%distc(i  ,j+1) +       &
                  pq4 * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_gvl*1000._r8 ) THEN
               gvl%flc(kk) = 0
               gvl%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_gvl
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',gvl%nc
   ENDIF

! 5 ) Perform vertical check
!!  - If min depth with good flag > MAX depth wih not good flag
!!  then get rid of the profile
   IF ( qck%vert ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Vertical check for glider velocity not implemented yet.               '
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_gvl
!-----------------------------------------------------------------------
!                                                                      !
!> Reset observation
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_gvl

   USE obs_str, ONLY : gvl

   IMPLICIT NONE

   WHERE (gvl%flc .EQ. 0)
      gvl%bia(:) = 0.
      gvl%res(:) = 0.
      gvl%inc(:) = 0.
      gvl%b_a(:) = 0.
      gvl%pq1(:) = 0.
      gvl%pq2(:) = 0.
      gvl%pq3(:) = 0.
      gvl%pq4(:) = 0.
      gvl%pq5(:) = 0.
      gvl%pq6(:) = 0.
      gvl%pq7(:) = 0.
      gvl%pq8(:) = 0.
   END WHERE
   gvl%nc = SUM(gvl%flc)

END SUBROUTINE reset_gvl
