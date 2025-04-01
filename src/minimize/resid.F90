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
!> Calculate analysis - observation                                     
!!
!! It computes the analysis observation differences weighted for the 
!! observational error
!!
!                                                                      !
! Version 1: Srdjan Dobricic  2006                                     !
!            Francesco Carere 2024 : Reproducibility                   !
!-----------------------------------------------------------------------
SUBROUTINE resid

   USE set_knd
   USE obs_str
#ifdef REPRO
   USE mpi_str
   USE rpr_str
#endif

   IMPLICIT NONE

   INTEGER(i4)   :: i, k
#ifdef REPRO
   INTEGER       :: ii
#endif

#ifdef REPRO
   ii = 1
   IF ( mpi%nproc .GT. 1 ) THEN
      IF ( ALLOCATED(rpr%obs_amo) ) DEALLOCATE ( rpr%obs_amo )
      ALLOCATE( rpr%obs_amo(rpr%obs_nog) )
      rpr%obs_amo(:) = 0.0_r8
   ENDIF
#endif

   k = 0
! ---
! Satellite observations of SLA
   DO i = 1,sla%no
      IF ( sla%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = sla%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! ARGO observations
   DO i = 1,arg%no
      IF ( arg%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = arg%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! XBT observations
   DO i = 1,xbt%no
      IF ( xbt%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = xbt%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Glider observations
   DO i = 1,gld%no
      IF ( gld%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = gld%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Observations of Argo float positions
   DO i = 1,tra%no
      IF ( tra%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = tra%inx(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO
   DO i = 1,tra%no
      IF ( tra%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = tra%iny(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Observations of positions of surface drifters
   DO i = 1,trd%no
      IF ( trd%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = trd%inx(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO
   DO i = 1,trd%no
      IF ( trd%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = trd%iny(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Observations of velocity by drIFters
   DO i = 1,vdr%no
      IF ( vdr%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = vdr%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Observations of velocity by gliders
   DO i = 1,gvl%no
      IF ( gvl%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = gvl%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

! ---
! Observations of SST
   DO i = 1,sst%no
      IF ( sst%flc(i) .EQ. 1 ) THEN
         k = k + 1
         obs%inc(k) = sst%inc(i)
         obs%amo(k) = ( obs%inc(k) - obs%res(k) ) / obs%err(k)
#ifdef REPRO
         IF ( mpi%nproc .GT. 1 ) THEN
            rpr%obs_amo(ii) = obs%amo(k)
         ENDIF
      ENDIF
      ii= ii + 1
#else
      ENDIF
#endif
   ENDDO

END SUBROUTINE resid
