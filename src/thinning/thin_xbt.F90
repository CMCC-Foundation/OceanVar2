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
!> Apply  thinning to XBT                                               
!!
!! if more observations lay in a grid cell they are averaged
!!                                                                     !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_xbt

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : thin, xbt
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
   ! Time from second to fraction of day
   time = thin%tim/86400._r8

   ALLOCATE ( allflag(xbt%no) )
   !Just in case a profile starts in a region and ENDs in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  xbt%flc,  allflag, xbt%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( xbt%nc .EQ. 0 ) return
   ELSE
      allflag = xbt%flc
   ENDIF

   ! Thinning
   DO k = 1,xbt%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1,xbt%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (xbt%lon(k)/cnti-xbt%lon(kk)) *COS(xbt%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (xbt%lat(k)/cnti-xbt%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (          xbt%kb(k)              .EQ. xbt%kb(kk) )  .AND. &
                 (          dist_obs               .LE. thin%spc   )  .AND. &
                 (ABS(xbt%tim(k)/cnti-xbt%tim(kk)) .LE. time       )  .AND. &
                 (          allflag(kk)            .EQ. 1          ) ) THEN
               xbt%lon(k)   = xbt%lon(k)   + xbt%lon(kk)
               xbt%lat(k)   = xbt%lat(k)   + xbt%lat(kk)
               xbt%tim(k)   = xbt%tim(k)   + xbt%tim(kk)
               xbt%val(k)   = xbt%val(k)   + xbt%val(kk)
               xbt%bac(k)   = xbt%bac(k)   + xbt%bac(kk)
               xbt%res(k)   = xbt%res(k)   + xbt%res(kk)
               xbt%dpt(k)   = xbt%dpt(k)   + xbt%dpt(kk)
               xbt%bgerr(k) = xbt%bgerr(k) + xbt%bgerr(kk)
               xbt%err(k)   = xbt%err(k)   + xbt%err(kk)
               xbt%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         xbt%eve(k)   = 998
         xbt%lon(k)   = xbt%lon(k)  /cnti
         xbt%lat(k)   = xbt%lat(k)  /cnti
         xbt%tim(k)   = xbt%tim(k)  /cnti
         xbt%val(k)   = xbt%val(k)  /cnti
         xbt%bac(k)   = xbt%bac(k)  /cnti
         xbt%res(k)   = xbt%res(k)  /cnti
         xbt%dpt(k)   = xbt%dpt(k)  /cnti
         xbt%bgerr(k) = xbt%bgerr(k)/cnti
         xbt%err(k)   = xbt%err(k)  /cnti
         xbt%rb(k)    =  (xbt%dpt(k) - grd%dep(xbt%kb(k))) / (grd%dep(xbt%kb(k)+1) - grd%dep(xbt%kb(k)))
      ENDIF
   ENDDO

   WHERE ( allflag .LT. xbt%flc ) xbt%flc = allflag
   CALL reset_xbt
   WRITE (drv%dia,*) 'Number of good XBT observations after thinning: ',  xbt%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   ! first subsampling if more profile are in the same grid point
   ! choose the one closest to the analysis time.
   ! it assumes all obs of the same profile have same time
   WRITE (drv%dia,*) 'Subsample based obs analysis time difference'
   DO k = 1,xbt%no-1
      time = ABS(xbt%tim(k)-drv%zanjul1950)
      IF ( xbt%flc(k) .EQ. 1 ) THEN
         WHERE ( xbt%jb(k)  .EQ. xbt%jb                      .AND. &  ! Same grid point x
                 xbt%ib(k)  .EQ. xbt%ib                      .AND. &  ! Same grid point y
                 xbt%ino(k) .NE. xbt%ino                     .AND. &  ! not Same profile
                 time       .LT. ABS(xbt%tim-drv%zanjul1950) .AND. &  ! Time of k obs is closer to analysis time
                 xbt%flc    .EQ. 1 )                                  ! Good flag
            xbt%eve = 13
            xbt%flc = 0
         END WHERE
      ENDIF
   ENDDO

   ! ... then thinning in vertical for the same profile
   DO k = 1,xbt%no-1
      cnti = 1.
      IF ( xbt%flc(k) .EQ. 1 ) THEN
         DO kk = k+1,xbt%no
            IF ( xbt%kb(k)  .EQ. xbt%kb(kk)       .AND. &      ! Same level
                 xbt%jb(k)  .EQ. xbt%jb(kk)       .AND. &      ! Same grid point x
                 xbt%ib(k)  .EQ. xbt%ib(kk)       .AND. &      ! Same grid point y
                 xbt%tim(k) .EQ. xbt%tim(kk)      .AND. &      ! Same time
                 xbt%ino(k) .EQ. xbt%ino(kk)      .AND. &      ! Same profile
                 xbt%flc(kk).EQ. 1 )              THEN         ! Good flag
               xbt%val(k)   = xbt%val(k)   + xbt%val(kk)        ! Obs value
               xbt%bac(k)   = xbt%bac(k)   + xbt%bac(kk)        ! Bkg value
               xbt%res(k)   = xbt%res(k)   + xbt%res(kk)        ! Misfits
               xbt%dpt(k)   = xbt%dpt(k)   + xbt%dpt(kk)        ! Depth
               xbt%bgerr(k) = xbt%bgerr(k) + xbt%bgerr(kk)      ! Bkg error
               xbt%err(k)   = xbt%err(k)   + xbt%err(kk)        ! Obs error
               xbt%flc(kk)  = 0
               xbt%eve(kk)  = 997
               cnti         = cnti + 1.
            ENDIF  ! Same obs
         ENDDO   ! kk
         xbt%eve(k)   = 998
         xbt%val(k)   = xbt%val(k)  /cnti                      ! Obs value
         xbt%bac(k)   = xbt%bac(k)  /cnti                      ! Bkg value
         xbt%res(k)   = xbt%res(k)  /cnti                      ! Misfits
         xbt%dpt(k)   = xbt%dpt(k)  /cnti                      ! Depth
         xbt%bgerr(k) = xbt%bgerr(k)/cnti                      ! Bkg error
         xbt%err(k)   = xbt%err(k)  /cnti                      ! Obs error
         xbt%rb(k)    =  (xbt%dpt(k) - grd%dep(xbt%kb(k))) / (grd%dep(xbt%kb(k)+1) - grd%dep(xbt%kb(k)))
      ENDIF
   ENDDO

   CALL reset_xbt
   WRITE (drv%dia,*) 'Number of good XBT observations after thinning:',  xbt%nc

END SUBROUTINE thin_xbt
