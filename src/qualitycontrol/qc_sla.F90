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
!> SLA quality control                                                 
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist or computed from EOFs
!! - Too close tp the coast,
!! - Obs in a too shallow area ,
!! - Obs too close to equator,
!!
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_sla

   USE set_knd
   USE obs_str, ONLY : sla, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)  :: kk,i,j
   REAL(r8),ALLOCATABLE :: treshold(:)
   REAL(r8)     :: val
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- SLA  observations:                ',sla%nc

!Initialize
   ALLOCATE ( sla%bgerr(sla%no) )
   sla%bgerr(:) = -999

! 0 ) Remove observations close to equator if dynamic height operator is active
   IF ( ANY(drv%bal(:) .EQ. 1) ) THEN
      WHERE ( ABS(sla%lat) .LE. 2._r8 )
         sla%flc = 0
         sla%eve = 10
      ENDWHERE
      CALL reset_sla
      WRITE (drv%dia,*)' ---- close to equator check:           ',sla%nc
   ENDIF

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( sla%flc .EQ. 1 .AND. ABS(sla%res) .GT. qck%res_sla )
         sla%flc = 0
         sla%eve = 5
      ENDWHERE
      CALL reset_sla
      WRITE (drv%dia,*)' ---- after residual check:             ',sla%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr ) THEN
! with constant background from namelist
      ALLOCATE ( treshold(sla%no) )
      ALLOCATE ( flc(sla%no) )
      flc = sla%flc
      sla%bgerr(:) = qck%bgr_sla(1)
      treshold (:) = qck%bgr_sla(2)
      CALL qc_bgerr(sla%no,sla%res,sla%err,sla%bgerr,treshold,sla%flc)
      CALL reset_sla
      WHERE ( (flc .EQ. 1) .AND. (sla%flc .EQ. 0 ) ) sla%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',sla%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
      IF ( drv%bmd(drv%ktr) .EQ. 1 ) THEN
         WRITE (drv%dia,*)' ------------------------------------- WARNING ------------------------------------------------'
         WRITE (drv%dia,*)' Background quality check for based on EOFs for SLA in case of barotropic model not implemented'
      ELSE
! with computed background from eofs
         IF ( .NOT. ALLOCATED( treshold ) ) ALLOCATE ( treshold(sla%no)  )
         ALLOCATE ( flc(sla%no) )
         flc = sla%flc
         CALL int_bgerr_sla
         treshold (:) = qck%bgr_sla(2)
         CALL qc_bgerr(sla%no,sla%res,sla%err,sla%bgerr,treshold,sla%flc)
         CALL reset_sla
         WHERE ( ( flc .EQ. 1 ) .AND. ( sla%flc .EQ. 0 ) ) sla%eve(:) =  7
         DEALLOCATE ( flc )
         WRITE (drv%dia,*)' ---- after EOFs background check:      ',sla%nc
      ENDIF
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! IF greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Climatological check for SLA not implemented yet.               '
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%sla ) THEN
      DO kk = 1,sla%no
         IF ( sla%flc(kk) .EQ. 1) THEN
            i = sla%ib(kk)
            j = sla%jb(kk)
            val = sla%pq1(kk) * coastrej%distc(i  ,j  ) +       &
                  sla%pq2(kk) * coastrej%distc(i+1,j  ) +       &
                  sla%pq3(kk) * coastrej%distc(i  ,j+1) +       &
                  sla%pq4(kk) * coastrej%distc(i+1,j+1)
            IF ( val .LT. coastrej%km_sla*1000._r8 ) THEN
               sla%flc(kk) = 0
               sla%eve(kk) = 9
            ENDIF
         ENDIF
      ENDDO
      CALL reset_sla
      WRITE (drv%dia,*)' ---- after coastal rejection:          ',sla%nc
   ENDIF

! 5 ) Remove observations shallower than level of no motion
   WHERE ( sla%flc .EQ. 1 .AND. sla%dpt < sla%dep )
      sla%flc(:) = 0
      sla%eve(:) = 12
   ENDWHERE
   CALL reset_sla
   WRITE (drv%dia,*)' ---- after level of no motion:        ',sla%nc

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_sla
!-----------------------------------------------------------------------
!                                                                      !
!> Interpolation background error
!                                                                      !
! Version 1: Mario Adani 2023                                          !             
!-----------------------------------------------------------------------
SUBROUTINE int_bgerr_sla

   USE set_knd
   USE obs_str, ONLY : sla, qck

   IMPLICIT NONE

   INTEGER(i4)   :: i, j, k

   DO k = 1, sla%no
      IF ( sla%flc(k) .EQ. 1 ) THEN
         i = sla%ib(k)
         j = sla%jb(k)
         sla%bgerr(k) = sla%pq1(k) * qck%eta(i  ,j  ) +       &
                        sla%pq2(k) * qck%eta(i+1,j  ) +       &
                        sla%pq3(k) * qck%eta(i  ,j+1) +       &
                        sla%pq4(k) * qck%eta(i+1,j+1)
      ENDIF
   ENDDO

END SUBROUTINE int_bgerr_sla
!-----------------------------------------------------------------------
!> Reset observation        
!                                                                      !
! Version 1: Mario Adani 2023                                        !             
!-----------------------------------------------------------------------
SUBROUTINE reset_sla

   USE obs_str, ONLY : sla

   IMPLICIT NONE

   WHERE ( sla%flc .EQ. 0 )
      sla%inc(:) = 0.
      sla%b_a(:) = 0.
      sla%pq1(:) = 0.
      sla%pq2(:) = 0.
      sla%pq3(:) = 0.
      sla%pq4(:) = 0.
   END WHERE

   sla%nc = SUM(sla%flc)

END SUBROUTINE reset_sla
