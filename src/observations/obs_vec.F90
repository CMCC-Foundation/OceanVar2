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
! Create the observational vector                                      !
!                                                                      !
! Version 1: Srdjan Dobricic  2006                                     !
! Version 2: Francesco Carere 2024 Reproducibility                     !
!-----------------------------------------------------------------------
SUBROUTINE obs_vec

   USE set_knd
   USE drv_str
   USE obs_str
#ifdef REPRO
   USE rpr_str
#endif

   IMPLICIT NONE

   INTEGER(i4)    ::  k, i

! -------
! Define observational vector
   WRITE (drv%dia,*) ' ---- Defining the Observational Vector'
   WRITE (drv%dia,*) ' ---- number of good observations (sla%nc): ', sla%nc
   WRITE (drv%dia,*) ' ---- number of good observations (arg%nc): ', arg%nc
   WRITE (drv%dia,*) ' ---- number of good observations (xbt%nc): ', xbt%nc
   WRITE (drv%dia,*) ' ---- number of good observations (gld%nc): ', gld%nc
   WRITE (drv%dia,*) ' ---- number of good observations (tra%nc): ', tra%nc
   WRITE (drv%dia,*) ' ---- number of good observations (trd%nc): ', trd%nc
   WRITE (drv%dia,*) ' ---- number of good observations (vdr%nc): ', vdr%nc
   WRITE (drv%dia,*) ' ---- number of good observations (gvk%nc): ', gvl%nc
   WRITE (drv%dia,*) ' ---- number of good observations (sst%nc): ', sst%nc
   obs%no = sla%nc + arg%nc + xbt%nc + gld%nc + 2 * tra%nc + 2 * trd%nc  &
          + vdr%nc + gvl%nc + sst%nc

#ifdef REPRO
   rpr%obs_nog = sla%no + arg%no + xbt%no + gld%no + 2 * tra%no + 2 * trd%no  &
               + vdr%no + gvl%no + sst%no
#endif

   WRITE (drv%dia,*) ' ---- Total number of good observations: ', obs%no
   WRITE (drv%dia,*) ' -------------------------------------- '
   ALLOCATE ( obs%inc(obs%no), obs%amo(obs%no), obs%res(obs%no) )
   ALLOCATE ( obs%err(obs%no), obs%gra(obs%no) )
   IF ( huberqc%any ) ALLOCATE ( obs%ahub(obs%no,2), obs%ahub2(obs%no,2) )

   k = 0

! SLA observations
   DO i = 1,sla%no
      IF ( sla%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = sla%res(i)
         obs%err(k) = sla%err(i)
      ENDIF
   ENDDO

! ARGO observations
   DO i = 1,arg%no
      IF ( arg%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = arg%res(i)
         obs%err(k) = arg%err(i)
      ENDIF
   ENDDO

! XBT observations
   DO i = 1,xbt%no
      IF ( xbt%flc(i) .EQ. 1 ) THEN
         k =  k+1
         obs%res(k) = xbt%res(i)
         obs%err(k) = xbt%err(i)
      ENDIF
   ENDDO

! Glider observations
   DO i = 1,gld%no
      IF ( gld%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = gld%res(i)
         obs%err(k) = gld%err(i)
      ENDIF
   ENDDO

! Argo trajectory observations
   DO i = 1,tra%no
      IF ( tra%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = tra%rex(i)
         obs%err(k) = tra%erx(i)
      ENDIF
   ENDDO
   DO i = 1,tra%no
      IF ( tra%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = tra%rey(i)
         obs%err(k) = tra%ery(i)
      ENDIF
   ENDDO

! Trajectory observations of surface drIFters
   DO i = 1,trd%no
      IF ( trd%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = trd%rex(i)
         obs%err(k) = trd%erx(i)
      ENDIF
   ENDDO
   DO i = 1,trd%no
      IF ( trd%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = trd%rey(i)
         obs%err(k) = trd%ery(i)
      ENDIF
   ENDDO

! Observations of drIFter velocities
   DO i = 1,vdr%no
      IF ( vdr%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = vdr%res(i)
         obs%err(k) = vdr%err(i)
      ENDIF
   ENDDO

! Observations of glider velocities
   DO i = 1,gvl%no
      IF ( gvl%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = gvl%res(i)
         obs%err(k) = gvl%err(i)
      ENDIF
   ENDDO

! SST observations
   DO i = 1,sst%no
      IF ( sst%flc(i) .EQ. 1 ) THEN
         k = k+1
         obs%res(k) = sst%res(i)
         obs%err(k) = sst%err(i)
      ENDIF
   ENDDO

END SUBROUTINE obs_vec
