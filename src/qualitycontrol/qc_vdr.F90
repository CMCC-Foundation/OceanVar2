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
!> drifter velocities quality control                                 
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist
!! - Too close tp the coast,
!!
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_vdr

   USE set_knd
   USE obs_str, ONLY : vdr, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)              :: kk,i,j
   REAL(r8),ALLOCATABLE     :: treshold(:)
   REAL(r8)                 :: val,pq1,pq2,pq3,pq4,weightsum
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- Drifter velocity  observations:   ',vdr%nc

!Initialize
   ALLOCATE ( vdr%bgerr(vdr%no) )
   vdr%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( vdr%flc .EQ. 1 .AND. ABS(vdr%res) .GT. qck%res_vel )
         vdr%flc(:) = 0
         vdr%eve(:) = 5
      ENDWHERE
      CALL reset_vdr
      WRITE (drv%dia,*)' ---- after residual check:             ',vdr%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
      ALLOCATE ( treshold(vdr%no) )
      ALLOCATE ( flc(vdr%no) )
      flc = vdr%flc
      vdr%bgerr(:) = qck%bgr_vel(1)
      treshold (:) = qck%bgr_vel(2)
      CALL qc_bgerr(vdr%no,vdr%res,vdr%err,vdr%bgerr,treshold,vdr%flc)
      CALL reset_vdr
      WHERE ( ( flc .EQ. 1 ) .AND. ( vdr%flc .EQ. 0 ) ) vdr%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',vdr%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Background error computed from EOFs not implemented yet               '
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! IF greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Climatological check for  drifter velocity not implemented yet.       '
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%vdr ) THEN
      DO kk = 1,vdr%no
         IF ( vdr%flc(kk) .EQ. 1 ) THEN
            i = vdr%ib(kk)
            j = vdr%jb(kk)
            weightsum = vdr%pq1(kk) + vdr%pq2(kk) + vdr%pq3(kk) + vdr%pq4(kk)
            pq1 = vdr%pq1(kk) / weightsum
            pq2 = vdr%pq2(kk) / weightsum
            pq3 = vdr%pq3(kk) / weightsum
            pq4 = vdr%pq4(kk) / weightsum
            val = pq1 * coastrej%distc(i  ,j  ) +       &
                  pq2 * coastrej%distc(i+1,j  ) +       &
                  pq3 * coastrej%distc(i  ,j+1) +       &
                  pq4 * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_vdr*1000._r8 ) THEN
               vdr%flc(kk) = 0
               vdr%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_vdr
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',vdr%nc
   ENDIF

! 5 ) Perform vertical check
!!  - if min depth with good flag > max depth wih not good flag
!!  then get rid of the profile
   IF ( qck%vert ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Vertical check for  drifter velocity not implemented yet.             '
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_vdr
!-----------------------------------------------------------------------
!> Reset observation        
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_vdr

   USE obs_str, ONLY : vdr

   IMPLICIT NONE

   WHERE (vdr%flc .EQ. 0)
      vdr%bia(:) = 0.
      vdr%res(:) = 0.
      vdr%inc(:) = 0.
      vdr%b_a(:) = 0.
      vdr%pq1(:) = 0.
      vdr%pq2(:) = 0.
      vdr%pq3(:) = 0.
      vdr%pq4(:) = 0.
      vdr%pq5(:) = 0.
      vdr%pq6(:) = 0.
      vdr%pq7(:) = 0.
      vdr%pq8(:) = 0.
   END WHERE

   vdr%nc = SUM(vdr%flc)

END SUBROUTINE reset_vdr

