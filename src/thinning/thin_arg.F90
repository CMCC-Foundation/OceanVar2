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
!> Apply  thinning to ARGO                                            
!!
!! if more observations lay in a grid cell they are averaged
!!                                                                     !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE thin_arg

   USE set_knd
   USE drv_str
   USE grd_str
   USE obs_str, ONLY : thin, arg
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

   ALLOCATE ( allflag(arg%no) )
   !Just in case a profile starts in a region and ends in another one
   IF ( mpi%nproc .GT. 1 ) THEN
      CALL MPI_ALLREDUCE(  arg%flc,  allflag, arg%no, mpi%i8  ,   &
         MPI_MAX,  mpi%comm, ierr)
      IF ( arg%nc .EQ. 0 ) return
   ELSE
      allflag = arg%flc
   ENDIF

   ! Thinning
   DO k = 1,arg%no-1
      IF ( allflag(k) .EQ. 1 ) THEN
         cnti = 1.
         DO kk = k+1,arg%no
            ! compute horizontal distance between obs
            dxx = phy%re*phy%d2r * (arg%lon(k)/cnti-arg%lon(kk)) *COS(arg%lat(k)/cnti*phy%d2r)
            dyy = phy%re*phy%d2r * (arg%lat(k)/cnti-arg%lat(kk))
            dist_obs = DSQRT(dxx**2+dyy**2)
            IF ( (          arg%par(k)           .EQ. arg%par(kk))  .AND. &
               (          arg%kb(k)              .EQ. arg%kb(kk) )  .AND. &
               (          dist_obs               .LE. thin%spc   )  .AND. &
               (ABS(arg%tim(k)/cnti-arg%tim(kk)) .LE. time       )  .AND. &
               (          allflag(kk)            .EQ. 1          ) ) THEN
               arg%lon(k)   = arg%lon(k)   + arg%lon(kk)
               arg%lat(k)   = arg%lat(k)   + arg%lat(kk)
               arg%tim(k)   = arg%tim(k)   + arg%tim(kk)
               arg%val(k)   = arg%val(k)   + arg%val(kk)
               arg%bac(k)   = arg%bac(k)   + arg%bac(kk)
               arg%res(k)   = arg%res(k)   + arg%res(kk)
               arg%dpt(k)   = arg%dpt(k)   + arg%dpt(kk)
               arg%bgerr(k) = arg%bgerr(k) + arg%bgerr(kk)
               arg%err(k)   = arg%err(k)   + arg%err(kk)
               allflag(kk) = 0
               arg%eve(kk) = 997
               cnti = cnti + 1.
            ENDIF ! Same obs
         ENDDO    ! kk
         arg%eve(k)   = 998
         arg%lon(k)   = arg%lon(k)  /cnti
         arg%lat(k)   = arg%lat(k)  /cnti
         arg%tim(k)   = arg%tim(k)  /cnti
         arg%val(k)   = arg%val(k)  /cnti
         arg%bac(k)   = arg%bac(k)  /cnti
         arg%res(k)   = arg%res(k)  /cnti
         arg%dpt(k)   = arg%dpt(k)  /cnti
         arg%bgerr(k) = arg%bgerr(k) + arg%bgerr(kk)
         arg%err(k)   = arg%err(k)  /cnti
         arg%rb(k)    =  (arg%dpt(k) - grd%dep(arg%kb(k))) / (grd%dep(arg%kb(k)+1) - grd%dep(arg%kb(k)))
      ENDIF
   ENDDO

   WHERE ( allflag .LT. arg%flc ) arg%flc = allflag
   CALL reset_arg
   WRITE (drv%dia,*) 'Number of good ARGO observations after thinning:',  arg%nc

   DEALLOCATE ( allflag )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

999 CONTINUE

   ! first subsampling if more profile are in the same grid point
   ! choose the one closest to the analysis time.
   ! it assumes all obs of the same profile have same time
   WRITE (drv%dia,*) 'Subsample based obs analysis time difference'
   DO k = 1,arg%no-1
      time = ABS(arg%tim(k)-drv%zanjul1950)
      IF ( arg%flc(k) .EQ. 1 ) THEN
         WHERE ( arg%par(k) .EQ. arg%par                     .AND. &  ! Same parameter
             arg%jb(k)      .EQ. arg%jb                      .AND. &  ! Same grid point x
             arg%ib(k)      .EQ. arg%ib                      .AND. &  ! Same grid point y
             arg%ino(k)     .NE. arg%ino                     .AND. &  ! not Same profile
             time           .LT. ABS(arg%tim-drv%zanjul1950) .AND. &  ! Time of k obs is closer to analysis time
             arg%flc        .EQ. 1 )                                  ! Good flag
             arg%eve = 13
             arg%flc = 0
         END WHERE
      ENDIF
   ENDDO

   ! ... then thinning in vertical for the same profile
   DO k = 1,arg%no-1
      cnti = 1.
      IF ( arg%flc(k) .EQ. 1 ) THEN
         DO kk = k+1, arg%no
            IF ( arg%par(k) .EQ. arg%par(kk)      .AND. &      ! Same parameter
                 arg%kb(k)  .EQ. arg%kb(kk)       .AND. &      ! Same leveel
                 arg%jb(k)  .EQ. arg%jb(kk)       .AND. &      ! Same grid point x
                 arg%ib(k)  .EQ. arg%ib(kk)       .AND. &      ! Same grid point y
                 arg%tim(k) .EQ. arg%tim(kk)      .AND. &      ! Same time
                 arg%ino(k) .EQ. arg%ino(kk)      .AND. &      ! Same profile
                 arg%flc(kk).EQ. 1 )              THEN         ! Good flag

               arg%val(k)   = arg%val(k)   + arg%val(kk)        ! Obs value
               arg%bac(k)   = arg%bac(k)   + arg%bac(kk)        ! Bkg value
               arg%res(k)   = arg%res(k)   + arg%res(kk)        ! Misfits
               arg%dpt(k)   = arg%dpt(k)   + arg%dpt(kk)        ! Depth
               arg%bgerr(k) = arg%bgerr(k) + arg%bgerr(kk)      ! Background error
               arg%err(k)   = arg%err(k)   + arg%err(kk)        ! Obs error
               arg%flc(kk)  = 0
               arg%eve(kk)  = 997
               cnti         = cnti + 1.
            ENDIF  ! Same obs
         ENDDO   ! kk
         arg%eve(k)   = 998
         arg%val(k)   = arg%val(k)  /cnti                      ! Obs value
         arg%bac(k)   = arg%bac(k)  /cnti                      ! Bkg value
         arg%res(k)   = arg%res(k)  /cnti                      ! Misfits
         arg%dpt(k)   = arg%dpt(k)  /cnti                      ! Depth
         arg%bgerr(k) = arg%bgerr(k)/cnti                      ! Bkg error
         arg%err(k)   = arg%err(k)  /cnti                      ! Obs error
         arg%rb(k)    =  (arg%dpt(k) - grd%dep(arg%kb(k))) / (grd%dep(arg%kb(k)+1) - grd%dep(arg%kb(k)))
      ENDIF
   ENDDO

   CALL reset_arg
   WRITE (drv%dia,*) 'Number of good ARGO observations after thinning:',  arg%nc

END SUBROUTINE thin_arg
