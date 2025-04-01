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
!> Apply thinning                                                     
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thinning

   USE obs_str
   USE drv_str

   IMPLICIT NONE

   IF ( thin%arg .AND. obs%arg .NE. 0 .AND. arg%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to arg...'
      CALL thin_arg
   ENDIF

   IF ( thin%gld .AND. obs%gld .NE. 0 .AND. gld%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to gld...'
      CALL thin_gld
   ENDIF

   IF ( thin%gvl .AND. obs%gvl .NE. 0 .AND. gvl%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to gvl...'
      CALL thin_gvl
   ENDIF

   IF ( thin%sla .AND. obs%sla .NE. 0 .AND. sla%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to sla...'
      CALL thin_sla
   ENDIF

   IF ( thin%sst .AND. obs%sst .NE. 0 .AND. sst%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to sst...'
      CALL thin_sst
   ENDIF

   IF ( thin%tra .AND. obs%tra .NE. 0 .AND. tra%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to tra...'
      CALL thin_tra
   ENDIF

   IF ( thin%trd .AND. obs%trd .NE. 0 .AND. trd%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to trd...'
      CALL thin_trd
   ENDIF

   IF ( thin%vdr .AND. obs%vdr .NE. 0 .AND. vdr%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to vdr...'
      CALL thin_vdr
   ENDIF

   IF ( thin%xbt .AND. obs%xbt .NE. 0 .AND. xbt%no .GT. 0 ) THEN
      WRITE (drv%dia,*) ' ...to xbt...'
      CALL thin_xbt
   ENDIF

END SUBROUTINE thinning
