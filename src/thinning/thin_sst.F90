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
!> Apply  thinning to SST                                               
!!
!! if more observations lay in a grid cell the one closest 
!! to the analysis time is retained
!!                                                                     !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_sst

   USE set_knd
   USE drv_str
   USE obs_str, ONLY : thin, sst
   USE mpi_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INCLUDE "mpif.h"

   INTEGER(i4)   ::  k,kk
   REAL(r8)      ::  cnti
   REAL(r8)      ::  time
   REAL(r8)      ::  dist_obs
   REAL(r8)      ::  dxx, dyy
   INTEGER                  :: ierr
   INTEGER(i8),ALLOCATABLE  :: allflag(:)

   GOTO 999
   WRITE (drv%dia,*) 'WARNING: Currently SST thinning does not take into accout Satellite differences'
   WRITE (drv%dia,*) 'It averges all the data within the time and space windows                     '
   WRITE (drv%dia,*) 'Maybe in future it would be USEfull to take into account different satellites  '

   ! Time from second to fraction of day
   time = thin%tim/86400._r8

   ALLOCATE ( allflag(sst%no) )
   !Just in case a profile starts in a region and ENDs in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  sst%flc,  allflag, sst%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF (sst%nc .EQ. 0 ) return
   ELSE
      allflag = sst%flc
   ENDIF

   ! Thinning
   DO k = 1,sst%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1,sst%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (sst%lon(k)/cnti-sst%lon(kk)) *COS(sst%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (sst%lat(k)/cnti-sst%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (        dist_obs                .LE. thin%spc    )  .AND. &
                 (ABS(sst%tim(k)/cnti-sst%tim(kk)).LE. time        )  .AND. &
                 (          allflag(kk)           .EQ. 1           ) ) THEN
               sst%lon(k)   = sst%lon(k)   + sst%lon(kk)
               sst%lat(k)   = sst%lat(k)   + sst%lat(kk)
               sst%tim(k)   = sst%tim(k)   + sst%tim(kk)
               sst%val(k)   = sst%val(k)   + sst%val(kk)
               sst%bac(k)   = sst%bac(k)   + sst%bac(kk)
               sst%res(k)   = sst%res(k)   + sst%res(kk)
               sst%bgerr(k) = sst%bgerr(k) + sst%bgerr(kk)
               sst%err(k)   = sst%err(k)   + sst%err(kk)
               sst%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         sst%eve(k)   = 998
         sst%lon(k)   = sst%lon(k)  /cnti
         sst%lat(k)   = sst%lat(k)  /cnti
         sst%tim(k)   = sst%tim(k)  /cnti
         sst%val(k)   = sst%val(k)  /cnti
         sst%bac(k)   = sst%bac(k)  /cnti
         sst%res(k)   = sst%res(k)  /cnti
         sst%bgerr(k) = sst%bgerr(k)/cnti
         sst%err(k)   = sst%err(k)  /cnti
      ENDIF
   ENDDO

   WHERE ( allflag .LT. sst%flc ) sst%flc = allflag
   CALL reset_sst
   WRITE (drv%dia,*) 'Number of good SST observations after thinning: ',  sst%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   WRITE (drv%dia,*) 'WARNING: Currently SST thinning does not take into accout Satellite differences'
   WRITE (drv%dia,*) 'It averges all the data within the time and space windows                     '
   WRITE (drv%dia,*) 'Maybe in future it would be usefull to take into account different satellites  '

   ! Subsample based on obs error
   WRITE (drv%dia,*) 'Subsample based obs error'

   DO k = 1,sst%no-1
      IF ( sst%flc(k) .EQ. 1 ) THEN
         DO kk = k+1,sst%no
            ! if in the same grid point
            IF ( sst%jb(k)   .EQ. sst%jb(kk)       .AND. &
                 sst%ib(k)   .EQ. sst%ib(kk)       .AND. &
                 sst%flc(kk) .EQ. 1 )              THEN
               ! choose  the one with the smaller error
               IF ( sst%err(kk) .GT. sst%err(k) ) THEN
                  sst%eve(kk) = 13
                  sst%flc(kk) = 0
               ELSE
                  sst%eve(k) = 13
                  sst%flc(k) = 0
               ENDIF
            ENDIF
         ENDDO
      ENDIF
   ENDDO

   CALL reset_sst
   WRITE (drv%dia,*) 'Number of good SST observations after thinning: ',  sst%nc

END SUBROUTINE thin_sst
