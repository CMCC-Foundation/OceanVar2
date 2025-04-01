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
!> Apply  thinning to GLIDER velocity                                   !
!!
!! if more observations lay in a grid cell they are averaged
!!                                                                     !
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_gvl

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : thin, gvl
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

   ALLOCATE ( allflag(gvl%no) )
   !Just in case a profile starts in a region and ends in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  gvl%flc,  allflag, gvl%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( gvl%nc .EQ. 0 ) return
   ELSE
      allflag = gvl%flc
   ENDIF

   ! Thinning
   DO k = 1,gvl%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1,gvl%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (gvl%lon(k)/cnti-gvl%lon(kk)) *COS(gvl%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (gvl%lat(k)/cnti-gvl%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (          gvl%par(k)             .EQ. gvl%par(kk))  .AND. &
               (            gvl%kb(k)              .EQ. gvl%kb(kk) )  .AND. &
               (            dist_obs               .LE. thin%spc   )  .AND. &
               (ABS(gvl%tim(k)/cnti-gvl%tim(kk))   .LE. time       )  .AND. &
               (          allflag(kk)              .EQ. 1          ) ) THEN
               gvl%lon(k)   = gvl%lon(k)   + gvl%lon(kk)
               gvl%lat(k)   = gvl%lat(k)   + gvl%lat(kk)
               gvl%tim(k)   = gvl%tim(k)   + gvl%tim(kk)
               gvl%val(k)   = gvl%val(k)   + gvl%val(kk)
               gvl%bac(k)   = gvl%bac(k)   + gvl%bac(kk)
               gvl%res(k)   = gvl%res(k)   + gvl%res(kk)
               gvl%dpt(k)   = gvl%dpt(k)   + gvl%dpt(kk)
               gvl%bgerr(k) = gvl%bgerr(k) + gvl%bgerr(kk)
               gvl%err(k)   = gvl%err(k)   + gvl%err(kk)
               gvl%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         gvl%eve(k)   = 998
         gvl%lon(k)   = gvl%lon(k)  /cnti
         gvl%lat(k)   = gvl%lat(k)  /cnti
         gvl%tim(k)   = gvl%tim(k)  /cnti
         gvl%val(k)   = gvl%val(k)  /cnti
         gvl%bac(k)   = gvl%bac(k)  /cnti
         gvl%res(k)   = gvl%res(k)  /cnti
         gvl%dpt(k)   = gvl%dpt(k)  /cnti
         gvl%bgerr(k) = gvl%bgerr(k)/cnti
         gvl%err(k)   = gvl%err(k)  /cnti
         gvl%rb(k)    =  (gvl%dpt(k) - grd%dep(gvl%kb(k))) / (grd%dep(gvl%kb(k)+1) - grd%dep(gvl%kb(k)))
      ENDIF
   ENDDO

   WHERE ( allflag .LT. gvl%flc ) gvl%flc = allflag
   CALL reset_gvl
   WRITE (drv%dia,*) 'Number of good GLIDER vel. obs. after thinning: ',  gvl%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   ! first subsampling if more profile are in the same grid point
   ! choose the one closest to the analysis time.
   ! it assumes all obs of the same profile have same time
   WRITE (drv%dia,*) 'Subsample based obs analysis time difference'
   DO k = 1,gvl%no-1
      time = ABS(gvl%tim(k)-drv%zanjul1950)
      IF ( gvl%flc(k) .EQ. 1 ) THEN
         WHERE ( gvl%par(k) .EQ. gvl%par                     .AND. &  ! Same parameter
                 gvl%jb(k)  .EQ. gvl%jb                      .AND. &  ! Same grid point x
                 gvl%ib(k)  .EQ. gvl%ib                      .AND. &  ! Same grid point y
                 gvl%ino(k) .NE. gvl%ino                     .AND. &  ! not Same profile
                 time       .LT. ABS(gvl%tim-drv%zanjul1950) .AND. &  ! Time of k obs is closer to analysis time
                 gvl%flc    .EQ. 1 )                                  ! Good flag
            gvl%eve = 13
            gvl%flc = 0
         END WHERE
      ENDIF
   ENDDO

   ! ... then thinning in vertical for the same profile
   DO k = 1,gvl%no-1
      cnti = 1.
      IF (gvl%flc(k).EQ.1) THEN
         DO kk = k+1, gvl%no
            IF ( gvl%par(k) .EQ. gvl%par(kk)      .AND. &      ! Same parameter
                 gvl%kb(k)  .EQ. gvl%kb(kk)       .AND. &      ! Same level
                 gvl%jb(k)  .EQ. gvl%jb(kk)       .AND. &      ! Same grid point x
                 gvl%ib(k)  .EQ. gvl%ib(kk)       .AND. &      ! Same grid point y
                 gvl%tim(k) .EQ. gvl%tim(kk)      .AND. &      ! Same time
                 gvl%ino(k) .EQ. gvl%ino(kk)      .AND. &      ! Same profile
                 gvl%flc(kk).EQ. 1 )              THEN         ! Good flag
               gvl%val(k)   = gvl%val(k)   + gvl%val(kk)        ! Obs value
               gvl%bac(k)   = gvl%bac(k)   + gvl%bac(kk)        ! Bkg value
               gvl%res(k)   = gvl%res(k)   + gvl%res(kk)        ! Misfits
               gvl%dpt(k)   = gvl%dpt(k)   + gvl%dpt(kk)        ! Depth
               gvl%bgerr(k) = gvl%bgerr(k) + gvl%bgerr(kk)      ! Bkg error
               gvl%err(k)   = gvl%err(k)   + gvl%err(kk)        ! Obs error
               gvl%flc(kk)  = 0
               gvl%eve(kk)  = 997
               cnti         = cnti + 1.
            ENDIF  ! Same obs
         ENDDO   ! kk
         gvl%eve(k)   = 998
         gvl%val(k)   = gvl%val(k)  /cnti                      ! Obs value
         gvl%bac(k)   = gvl%bac(k)  /cnti                      ! Bkg value
         gvl%res(k)   = gvl%res(k)  /cnti                      ! Misfits
         gvl%dpt(k)   = gvl%dpt(k)  /cnti                      ! Depth
         gvl%bgerr(k) = gvl%bgerr(k)/cnti                      ! Bkg error
         gvl%err(k)   = gvl%err(k)  /cnti                      ! Obs error
         gvl%rb(k)    =  (gvl%dpt(k) - grd%dep(gvl%kb(k))) / (grd%dep(gvl%kb(k)+1) - grd%dep(gvl%kb(k)))
      ENDIF
   ENDDO

   CALL reset_gvl
   WRITE (drv%dia,*) 'Number of good GLIDER observations after thinning:',  gvl%nc

END SUBROUTINE thin_gvl

