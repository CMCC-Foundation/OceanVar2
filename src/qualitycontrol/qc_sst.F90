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
!> SST quality control                                                  
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist or computed from EOFs
!! - Too close tp the coast,
!!
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_sst

   USE set_knd
   USE obs_str, ONLY : sst, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)          :: kk, k, j, i
   REAL(r8),ALLOCATABLE :: treshold(:)
   REAL(r8)             :: val
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- SST observations:                ',sst%nc

!Initialize
   ALLOCATE ( sst%bgerr(sst%no) )
   sst%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( sst%flc .EQ. 1 .AND. ABS(sst%res) .GT. qck%res_sst)
         sst%flc = 0
         sst%eve = 5
      ENDWHERE
      CALL reset_sst
      WRITE (drv%dia,*)' ---- after residual check:             ',sst%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
! with constant background from namelist
      ALLOCATE ( treshold(sst%no) )
      ALLOCATE ( flc(sst%no) )
      flc = sst%flc
      sst%bgerr(:) = qck%bgr_sst(1)
      treshold (:) = qck%bgr_sst(2)
      CALL qc_bgerr(sst%no,sst%res,sst%err,sst%bgerr,treshold,sst%flc)
      CALL reset_sst
      WHERE ( ( flc .EQ. 1 ) .AND. ( sst%flc .EQ. 0 ) ) sst%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',sst%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Background error computed from EOFs not implementd yet.               '
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! IF greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      DO kk = 1,sst%no
         IF ( sst%flc(kk).EQ.1 ) THEN
            i = sst%ib(kk)
            j = sst%jb(kk)
            k = sst%kb(kk)
            ! for safe interpolation
            IF ( ANY(qck%climtem(i:i+1,j:j+1,k:k+1) .LE. 0) ) cycle
            val = sst%pq1(kk) * qck%climtem(i  ,j ,k ) +       &
                  sst%pq2(kk) * qck%climtem(i+1,j ,k ) +       &
                  sst%pq3(kk) * qck%climtem(i  ,j+1,k) +       &
                  sst%pq4(kk) * qck%climtem(i+1,j+1,k)
            IF ( ABS(val-sst%val(kk)) .GT. qck%clm_lim(1) ) THEN
               sst%flc(kk) = 0
               sst%eve(kk) = 8
            ENDIF
         ENDIF
      ENDDO
      CALL reset_sst
      WRITE (drv%dia,*)' ---- after climatological check:       ',sst%nc
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%sst ) THEN
      DO kk = 1,sst%no
         IF ( sst%flc(kk) .EQ. 1) THEN
            i = sst%ib(kk)
            j = sst%jb(kk)
            val = sst%pq1(kk) * coastrej%distc(i  ,j  ) +       &
                  sst%pq2(kk) * coastrej%distc(i+1,j  ) +       &
                  sst%pq3(kk) * coastrej%distc(i  ,j+1) +       &
                  sst%pq4(kk) * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_sst*1000._r8 ) THEN
               sst%flc(kk) = 0
               sst%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_sst
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',sst%nc
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_sst
!-----------------------------------------------------------------------
!> Reset observation        
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_sst

   USE obs_str, ONLY : sst

   IMPLICIT NONE

   WHERE ( sst%flc .EQ. 0 )
      sst%inc(:) = 0.
      sst%b_a(:) = 0.
      sst%pq1(:) = 0.
      sst%pq2(:) = 0.
      sst%pq3(:) = 0.
      sst%pq4(:) = 0.
   END WHERE

   sst%nc = SUM(sst%flc)

END SUBROUTINE reset_sst
