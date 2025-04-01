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
!> Apply  thinning to DRIFTER velocity                                  
!!
!! if more observations lay in a grid cell they are averaged
!!                                                                     !
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_vdr

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : thin, vdr
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

   ALLOCATE ( allflag(vdr%no) )
   !Just in case a profile starts in a region and ends in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  vdr%flc,  allflag, vdr%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( vdr%nc .EQ. 0 ) return
   ELSE
      allflag = vdr%flc
   ENDIF

   ! Thinning
   DO k = 1,vdr%no-1
      IF (allflag(k).EQ.1) THEN
         cnti = 1.
         DO kk = k+1,vdr%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (vdr%lon(k)/cnti-vdr%lon(kk)) *COS(vdr%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (vdr%lat(k)/cnti-vdr%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (          vdr%par(k)             .EQ. vdr%par(kk))  .AND. &
                 (          vdr%kb(k)              .EQ. vdr%kb(kk) )  .AND. &
                 (          dist_obs               .LE. thin%spc   )  .AND. &
                 (ABS(vdr%tim(k)/cnti-vdr%tim(kk)) .LE. time       )  .AND. &
                 (          allflag(kk)            .EQ. 1          ) ) THEN
               vdr%lon(k)   = vdr%lon(k)   + vdr%lon(kk)
               vdr%lat(k)   = vdr%lat(k)   + vdr%lat(kk)
               vdr%tim(k)   = vdr%tim(k)   + vdr%tim(kk)
               vdr%val(k)   = vdr%val(k)   + vdr%val(kk)
               vdr%bac(k)   = vdr%bac(k)   + vdr%bac(kk)
               vdr%res(k)   = vdr%res(k)   + vdr%res(kk)
               vdr%dpt(k)   = vdr%dpt(k)   + vdr%dpt(kk)
               vdr%bgerr(k) = vdr%bgerr(k) + vdr%bgerr(kk)
               vdr%err(k)   = vdr%err(k)   + vdr%err(kk)
               vdr%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         vdr%eve(k)   = 998
         vdr%lon(k)   = vdr%lon(k)  /cnti
         vdr%lat(k)   = vdr%lat(k)  /cnti
         vdr%tim(k)   = vdr%tim(k)  /cnti
         vdr%val(k)   = vdr%val(k)  /cnti
         vdr%bac(k)   = vdr%bac(k)  /cnti
         vdr%res(k)   = vdr%res(k)  /cnti
         vdr%dpt(k)   = vdr%dpt(k)  /cnti
         vdr%bgerr(k) = vdr%bgerr(k)/cnti
         vdr%err(k)   = vdr%err(k)  /cnti
         vdr%rb(k)    =  (vdr%dpt(k) - grd%dep(vdr%kb(k))) / (grd%dep(vdr%kb(k)+1) - grd%dep(vdr%kb(k)))
      ENDIF
   ENDDO

   WHERE ( allflag .LT. vdr%flc ) vdr%flc = allflag
   CALL reset_vdr
   WRITE (drv%dia,*) 'Number of good DrIF. vel. obs.  after thinning: ',  vdr%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   ! first subsampling if more profile are in the same grid point
   ! choose the one closest to the analysis time.
   ! it assumes all obs of the same profile have same time
   WRITE (drv%dia,*) 'Subsample based obs analysis time difference'
   DO k = 1,vdr%no-1
      time = ABS(vdr%tim(k)-drv%zanjul1950)
      IF ( vdr%flc(k) .EQ. 1 ) THEN
         WHERE ( vdr%par(k) .EQ. vdr%par                     .AND. &  ! Same parameter
                 vdr%jb(k)  .EQ. vdr%jb                      .AND. &  ! Same grid point x
                 vdr%ib(k)  .EQ. vdr%ib                      .AND. &  ! Same grid point y
                 vdr%ino(k) .NE. vdr%ino                     .AND. &  ! Not Same profile
                 time       .LT. ABS(vdr%tim-drv%zanjul1950) .AND. &  ! Time of k obs is closer to analysis time
                 vdr%flc    .EQ. 1 )                                  ! Good flag
            vdr%eve = 13
            vdr%flc = 0
         END WHERE
      ENDIF
   ENDDO

   ! ... then thinning in vertical for the same profile
   DO k = 1,vdr%no-1
      cnti = 1.
      IF ( vdr%flc(k) .EQ. 1 ) THEN
         DO kk = k+1,vdr%no
            IF ( vdr%par(k) .EQ. vdr%par(kk)      .AND. &      ! Same parameter
                 vdr%kb(k)  .EQ. vdr%kb(kk)       .AND. &      ! Same.LE.el
                 vdr%jb(k)  .EQ. vdr%jb(kk)       .AND. &      ! Same grid point x
                 vdr%ib(k)  .EQ. vdr%ib(kk)       .AND. &      ! Same grid point y
                 vdr%tim(k) .EQ. vdr%tim(kk)      .AND. &      ! Same time
                 vdr%ino(k) .EQ. vdr%ino(kk)      .AND. &      ! Same profile
                 vdr%flc(kk).EQ. 1 )              THEN         ! Good flag
               vdr%val(k)   = vdr%val(k)   + vdr%val(kk)        ! Obs value
               vdr%bac(k)   = vdr%bac(k)   + vdr%bac(kk)        ! Bkg value
               vdr%res(k)   = vdr%res(k)   + vdr%res(kk)        ! Misfits
               vdr%dpt(k)   = vdr%dpt(k)   + vdr%dpt(kk)        ! Depth
               vdr%bgerr(k) = vdr%bgerr(k) + vdr%bgerr(kk)      ! Bkg error
               vdr%err(k)   = vdr%err(k)   + vdr%err(kk)        ! Obs error
               vdr%flc(kk)  = 0
               vdr%eve(kk)  = 997
               cnti         = cnti + 1.
            ENDIF  ! Same obs
         ENDDO   ! kk
         vdr%eve(k)   = 998
         vdr%val(k)   = vdr%val(k)  /cnti                      ! Obs value
         vdr%bac(k)   = vdr%bac(k)  /cnti                      ! Bkg value
         vdr%res(k)   = vdr%res(k)  /cnti                      ! Misfits
         vdr%dpt(k)   = vdr%dpt(k)  /cnti                      ! Depth
         vdr%bgerr(k) = vdr%bgerr(k)/cnti                      ! Bkg error
         vdr%err(k)   = vdr%err(k)  /cnti                      ! Obs error
         vdr%rb(k)    =  (vdr%dpt(k) - grd%dep(vdr%kb(k))) / (grd%dep(vdr%kb(k)+1) - grd%dep(vdr%kb(k)))
      ENDIF
   ENDDO

   CALL reset_vdr
   WRITE (drv%dia,*) 'Number of good Drifter vel. obs.  after thinning: ',  vdr%nc

END SUBROUTINE thin_vdr

