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
!> ARGO trajectory quality control                                    
!!
!! - Misfit qc      abs[y-h(xb)] > threshold,
!! - background qc  [y-h(xb)]**2/( sigmab**2 + sigmao**2) > threshold ,
!! sigmab can be from namelist 
!!
!                                                                      !
! Version 1: Andrea Storto 2022                                        !
!            Mario Adani   2023                                        !
!-----------------------------------------------------------------------
SUBROUTINE qc_tra

   USE set_knd
   USE obs_str, ONLY : tra, qck, coastrej
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)  :: k
   REAL(r8),ALLOCATABLE :: treshold(:)
   INTEGER(i8),ALLOCATABLE  :: flc(:)

   WRITE (drv%dia,*)' -------------------------------------- '
   WRITE (drv%dia,*)' ---- ARGO trajectory observations:     ',tra%nc

!Initialize
   ALLOCATE ( tra%bgerr(tra%no) )
   tra%bgerr(:) = -999

! 1 ) Remove obseravtaions with large residuals
! residual check
   IF ( qck%res ) THEN
      WHERE ( tra%flc .EQ. 1 .AND. ABS(tra%rex) .GT. qck%res_dis .OR. ABS(tra%rey) .GT. qck%res_dis )
         tra%flc(:) = 0
         tra%eve(:) = 5
      ENDWHERE
      CALL reset_tra
      WRITE (drv%dia,*)' ---- after residual check:             ',tra%nc
   ENDIF

! 2 ) Perform background quality check and reject observations accordingly
!!  - evaluate [y-h(xb)]**2 > alfa * ( sigmab**2 + sigmao**2)
!!  where sigmab and sigmao are background and observational
!!  errors, respectively, and alfa is a user-defined factor.
!!  when the relation is satisfied, the obs is rejected.
   IF ( qck%conbgr )  THEN
! with constant background from namelist
      ALLOCATE ( treshold(tra%no) )
      ALLOCATE ( flc(tra%no) )
      flc = tra%flc
      tra%bgerr(:) = qck%bgr_dis(1)
      treshold (:) = qck%bgr_dis(2)
      CALL qc_bgerr(tra%no,tra%rex,tra%err,tra%bgerr,treshold,tra%flc)
      CALL qc_bgerr(tra%no,tra%rey,tra%err,tra%bgerr,treshold,tra%flc)
      CALL reset_tra
      WHERE ( ( flc .EQ. 1 ) .AND. ( tra%flc .EQ. 0 ) ) tra%eve(:) =  6
      DEALLOCATE ( flc )
      WRITE (drv%dia,*)' ---- after constant background check:  ',tra%nc
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
      WRITE (drv%dia,*)' Climatological check for ARGO trajectory not implemented yet.         '
   ENDIF

! 4 ) Perform Coastal Rejection
   IF ( coastrej%tra ) THEN
      WRITE (drv%dia,*)' --------------------- WARNING --------------------------------------- '
      WRITE (drv%dia,*)' Coastal rejection for ARGO trajectory not implemented yet.            '
   ENDIF

   IF ( ALLOCATED(treshold) )  DEALLOCATE ( treshold )

END SUBROUTINE qc_tra
!-----------------------------------------------------------------------
!> Count observation        
!                                                                      !
! Version 1: Mario Adani 2023                                          ! 
!-----------------------------------------------------------------------
SUBROUTINE reset_tra

   USE obs_str, ONLY : tra

   IMPLICIT NONE

   tra%nc = SUM(tra%flc)

END SUBROUTINE reset_tra

