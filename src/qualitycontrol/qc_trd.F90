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
!> Surface drifter trajectory quality control                          
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist 
!!
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_trd

   USE set_knd
   USE obs_str, ONLY : trd, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)  :: k
   REAL(r8),ALLOCATABLE :: treshold(:)
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- Surf. drif. traj. observations:   ',trd%nc

!Initialize
   ALLOCATE ( trd%bgerr(trd%no) )
   trd%bgerr(:) = -999
!
! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( trd%flc .EQ. 1 .AND. ABS(trd%rex) .GT. qck%res_dis .OR. ABS(trd%rey) .GT. qck%res_dis)
         trd%flc(:) = 0
         trd%eve(:) = 5
      ENDWHERE
      CALL reset_trd
      WRITE (drv%dia,*)' ---- after residual check:             ',trd%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr )  THEN
! with constant background from NAMELIST
      ALLOCATE ( treshold(trd%no) )
      ALLOCATE ( flc(trd%no) )
      flc = trd%flc
      trd%bgerr(:) = qck%bgr_dis(1)
      treshold (:) = qck%bgr_dis(2)
      CALL qc_bgerr(trd%no,trd%rex,trd%err,trd%bgerr,treshold,trd%flc)
      CALL qc_bgerr(trd%no,trd%rey,trd%err,trd%bgerr,treshold,trd%flc)
      CALL reset_trd
      WHERE ( ( flc .EQ. 1 ) .AND. ( trd%flc .EQ. 0 ) ) trd%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',trd%nc
   ENDIF
   IF ( qck%eofbgr ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING ----------------------------------------- '
      WRITE (drv%dia,*)' Background error computed from EOFs not implemented yet.                '
   ENDIF

! 3 ) Perform climatological quality check
!! - evaluate abs(clim-obs)
!! IF greater then a user-defined factor the obs is rejected.
   IF ( qck%clm ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Climatological check for drifter trajectory not implemented yet.      '
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%trd ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Coastal rejection for drifter trajectory not implemented yet.         '
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_trd
!-----------------------------------------------------------------------
!> Count observation        
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_trd

   USE obs_str, ONLY : trd

   IMPLICIT NONE

   trd%nc = SUM(trd%flc)

END SUBROUTINE reset_trd

