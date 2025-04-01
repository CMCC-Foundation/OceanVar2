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
!> Apply  thinning to SLA                                              
!!
!! if more observations lay in a grid cell the one closest 
!! to the analysis time is retained
!!                                                                     !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_sla

   USE set_knd
   USE drv_str
   USE obs_str, ONLY : thin, sla
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER(i4)   ::  k,kk
   REAL(r8)      ::  cnti
   REAL(r8)      ::  time
   REAL(r8)      ::  dist_obs
   REAL(r8)      ::  dxx, dyy
   INTEGER       :: ierr
   INTEGER(i8),ALLOCATABLE  :: allflag(:)

   GOTO 999

   ! Time from second to fraction of day
   time = thin%tim/86400._r8

   ALLOCATE ( allflag(sla%no) )
   !Just in case a profile starts in a region and ends in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  sla%flc,  allflag, sla%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( sla%nc .EQ. 0 ) return
   ELSE
      allflag = sla%flc
   ENDIF

   ! Thinning
   DO k = 1,sla%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1,sla%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (sla%lon(k)/cnti-sla%lon(kk)) *COS(sla%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (sla%lat(k)/cnti-sla%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (        sla%ksat(k)               .EQ. sla%ksat(kk))  .AND. &
                 (          dist_obs                .LE. thin%spc    )  .AND. &
                 (ABS(sla%tim(k)/cnti-sla%tim(kk))  .LE. time        )  .AND. &
                 (          allflag(kk)             .EQ. 1           ) ) THEN
               sla%lon(k)   = sla%lon(k)   + sla%lon(kk)
               sla%lat(k)   = sla%lat(k)   + sla%lat(kk)
               sla%tim(k)   = sla%tim(k)   + sla%tim(kk)
               sla%val(k)   = sla%val(k)   + sla%val(kk)
               sla%bac(k)   = sla%bac(k)   + sla%bac(kk)
               sla%res(k)   = sla%res(k)   + sla%res(kk)
               sla%dpt(k)   = sla%dpt(k)   + sla%dpt(kk)
               sla%bgerr(k) = sla%bgerr(k) + sla%bgerr(kk)
               sla%err(k)   = sla%err(k)   + sla%err(kk)
               sla%dtm(k)   = sla%dtm(k)   + sla%dtm(kk)
               sla%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         sla%eve(k)   = 998
         sla%lon(k)   = sla%lon(k)  /cnti
         sla%lat(k)   = sla%lat(k)  /cnti
         sla%tim(k)   = sla%tim(k)  /cnti
         sla%val(k)   = sla%val(k)  /cnti
         sla%bac(k)   = sla%bac(k)  /cnti
         sla%res(k)   = sla%res(k)  /cnti
         sla%dpt(k)   = sla%dpt(k)  /cnti
         sla%bgerr(k) = sla%bgerr(k)/cnti
         sla%err(k)   = sla%err(k)  /cnti
         sla%dtm(k)   = sla%dtm(k)  /cnti
      ENDIF
   ENDDO

   WHERE ( allflag .LT. sla%flc ) sla%flc = allflag
   CALL reset_sla
   WRITE (drv%dia,*) 'Number of good SLA observations after thinning: ',  sla%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   WRITE (drv%dia,*) 'WARNING: Currently SLA thinning does not take into accout Satellite differences'
   WRITE (drv%dia,*) 'It averges all the data within the time and space windows                     '
   WRITE (drv%dia,*) 'Maybe in future it would be usefull to take into account different satellites  '
   ! Subsample based on obs error
   WRITE (drv%dia,*) 'Subsample based obs error'

   DO k = 1,sla%no-1
      IF ( sla%flc(k) .EQ. 1 ) THEN
         DO kk = k+1,sla%no
            ! IF in the same grid point
            IF ( sla%jb(k)   .EQ. sla%jb(kk)       .AND. &
                 sla%ib(k)   .EQ. sla%ib(kk)       .AND. &
                 sla%flc(kk) .EQ. 1 )              THEN
               ! choose  the one with the smaller error
               IF ( sla%err(kk) .GT. sla%err(k) ) THEN
                  sla%eve(kk) = 13
                  sla%flc(kk) = 0
               ELSE
                  sla%eve(k) = 13
                  sla%flc(k) = 0
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDDO

   CALL reset_sla
   WRITE (drv%dia,*) 'Number of good SLA observations after thinning: ',  sla%nc

END SUBROUTINE thin_sla
