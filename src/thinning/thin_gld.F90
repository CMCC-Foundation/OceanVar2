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
!> Apply  thinning to GLIDER                                            
!!
!! if more observations lay in a grid cell they are averaged
!!                                                                     !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_gld

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : thin, gld
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

   ALLOCATE ( allflag(gld%no) )
   !Just in case a prof.LE.starts in a region and ends in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  gld%flc,  allflag, gld%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( gld%nc .EQ. 0 ) return
   ELSE
      allflag = gld%flc
   ENDIF

   ! Thinning
   DO k = 1,gld%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1, gld%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (gld%lon(k)/cnti-gld%lon(kk)) *COS(gld%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (gld%lat(k)/cnti-gld%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (          gld%par(k)             .EQ. gld%par(kk))  .AND. &
               (            gld%kb(k)              .EQ. gld%kb(kk) )  .AND. &
               (            dist_obs               .LE. thin%spc   )  .AND. &
               (ABS(gld%tim(k)/cnti-gld%tim(kk))   .LE. time       )  .AND. &
               (            allflag(kk)            .EQ. 1          ) ) THEN
               gld%lon(k)   = gld%lon(k)   + gld%lon(kk)
               gld%lat(k)   = gld%lat(k)   + gld%lat(kk)
               gld%tim(k)   = gld%tim(k)   + gld%tim(kk)
               gld%val(k)   = gld%val(k)   + gld%val(kk)
               gld%bac(k)   = gld%bac(k)   + gld%bac(kk)
               gld%res(k)   = gld%res(k)   + gld%res(kk)
               gld%dpt(k)   = gld%dpt(k)   + gld%dpt(kk)
               gld%bgerr(k) = gld%bgerr(k) + gld%bgerr(kk)
               gld%err(k)   = gld%err(k)   + gld%err(kk)
               gld%eve(kk) = 997
               allflag(kk) = 0
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         gld%eve(k)   = 998
         gld%lon(k)   = gld%lon(k)  /cnti
         gld%lat(k)   = gld%lat(k)  /cnti
         gld%tim(k)   = gld%tim(k)  /cnti
         gld%val(k)   = gld%val(k)  /cnti
         gld%bac(k)   = gld%bac(k)  /cnti
         gld%res(k)   = gld%res(k)  /cnti
         gld%dpt(k)   = gld%dpt(k)  /cnti
         gld%bgerr(k) = gld%bgerr(k)/cnti
         gld%err(k)   = gld%err(k)  /cnti
         gld%rb(k)    =  (gld%dpt(k) - grd%dep(gld%kb(k))) / (grd%dep(gld%kb(k)+1) - grd%dep(gld%kb(k)))
      ENDIF
   ENDDO

   WHERE ( allflag .LT. gld%flc ) gld%flc = allflag
   CALL reset_gld
   WRITE (drv%dia,*) 'Number of good GLIDER observations after thinning:',  gld%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   ! first subsampling if more profile are in the same grid point
   ! choose the one closest to the analysis time.
   ! it assumes all obs of the same profile have same time
   WRITE (drv%dia,*) 'Subsample based obs analysis time difference'
   DO k = 1,gld%no-1
      time = ABS(gld%tim(k)-drv%zanjul1950)
      IF ( gld%flc(k) .EQ. 1 ) THEN
         WHERE ( gld%par(k) .EQ. gld%par                     .AND. &  ! Same parameter
                 gld%jb(k)  .EQ. gld%jb                      .AND. &  ! Same grid point x
                 gld%ib(k)  .EQ. gld%ib                      .AND. &  ! Same grid point y
                 gld%ino(k) .NE. gld%ino                     .AND. &  ! not Same profile
                 time       .LT. ABS(gld%tim-drv%zanjul1950) .AND. &  ! Time of k obs is closer to analysis time
                 gld%flc    .EQ. 1 )                                  ! Good flag
            gld%eve = 13
            gld%flc = 0
         END WHERE
      ENDIF
   ENDDO

   ! ... THEN thinning in vertical for the same profile
   DO k = 1,gld%no-1
      cnti = 1.
      IF ( gld%flc(k) .EQ. 1 ) THEN
         DO kk = k+1,gld%no
            IF ( gld%par(k) .EQ. gld%par(kk)      .AND. &      ! Same parameter
                 gld%kb(k)  .EQ. gld%kb(kk)       .AND. &      ! Same level
                 gld%jb(k)  .EQ. gld%jb(kk)       .AND. &      ! Same grid point x
                 gld%ib(k)  .EQ. gld%ib(kk)       .AND. &      ! Same grid point y
                 gld%tim(k) .EQ. gld%tim(kk)      .AND. &      ! Same time
                 gld%ino(k) .EQ. gld%ino(kk)      .AND. &      ! Same profile
                 gld%flc(kk).EQ. 1 )              THEN         ! Good flag
               gld%val(k)   = gld%val(k)   + gld%val(kk)        ! Obs value
               gld%bac(k)   = gld%bac(k)   + gld%bac(kk)        ! Bkg value
               gld%res(k)   = gld%res(k)   + gld%res(kk)        ! Misfits
               gld%dpt(k)   = gld%dpt(k)   + gld%dpt(kk)        ! Depth
               gld%bgerr(k) = gld%bgerr(k) + gld%bgerr(kk)      ! Bkg error
               gld%err(k)   = gld%err(k)   + gld%err(kk)        ! Obs error
               gld%flc(kk)  = 0
               gld%eve(kk)  = 997
               cnti         = cnti + 1.
            ENDIF  ! Same obs
         ENDDO   ! kk
         gld%eve(k)   = 998
         gld%val(k)   = gld%val(k)  /cnti                      ! Obs value
         gld%bac(k)   = gld%bac(k)  /cnti                      ! Bkg value
         gld%res(k)   = gld%res(k)  /cnti                      ! Misfits
         gld%dpt(k)   = gld%dpt(k)  /cnti                      ! Depth
         gld%bgerr(k) = gld%bgerr(k)/cnti                      ! Bkg error
         gld%err(k)   = gld%err(k)  /cnti                      ! Obs error
         gld%rb(k)    =  (gld%dpt(k) - grd%dep(gld%kb(k))) / (grd%dep(gld%kb(k)+1) - grd%dep(gld%kb(k)))
      ENDIF
   ENDDO

   CALL reset_gld
   WRITE (drv%dia,*) 'Number of good GLIDER observations after thinning:',  gld%nc

END SUBROUTINE thin_gld
