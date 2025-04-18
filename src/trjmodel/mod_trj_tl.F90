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
!>  Trajectory model                                     
!                                                                      !
! Version 1: V. Taillandier 2007                                       !
!                                                                      !
!-----------------------------------------------------------------------
SUBROUTINE mod_trj_tl( jpi,jpj,umod,vmod,e1u,e2v, flc, fls,         &
   jpt,jpn,pimod,pjmod,ptime,utan,vtan,xtan,ytan )
!!----------------------------------------------------------------------
!! Arguments
!! =========
   IMPLICIT NONE
   INTEGER*4 :: jpi,jpj
   INTEGER*8 :: jpt,jpn
   INTEGER*8 :: flc(jpn), fls(jpn)
   REAL*8, DIMENSION(jpt+1,jpn) :: pimod,pjmod
   REAL*8, DIMENSION(jpn) :: ptime
   REAL*8, DIMENSION(jpn) :: xtan,ytan
   REAL*8, DIMENSION(jpi,jpj) :: umod,vmod,e1u,e2v
   REAL*8, DIMENSION(jpi,jpj) :: utan,vtan
!!----------------------------------------------------------------------
!! Local declarations
!! ==================
   INTEGER :: jn,jt
   REAL*8 :: zpin,zpjn,zpitan,zpjtan,zrdt,zutan,zvtan
!!----------------------------------------------------------------------

!*... Trajectories computation
   DO jn = 1, jpn
      if(flc(jn).eq.1 .or. fls(jn).eq.1)then
         zpitan = 0.d0
         zpjtan = 0.d0
         zrdt = ptime(jn)*3600.d0/(jpt)
         DO jt = 1, jpt
            zpin = pimod(jt,jn)
            zpjn = pjmod(jt,jn)
            CALL floitptan( jpi,jpj,umod,vmod,utan,vtan,e1u,e2v,    &
               zpin,zpjn,zpitan,zpjtan,zutan,zvtan )
            zutan  = zpitan + zrdt * zutan
            zvtan  = zpjtan + zrdt * zvtan
            zpitan = zutan
            zpjtan = zvtan
         END DO
         xtan(jn) = zpitan
         ytan(jn) = zpjtan
      endif
   END DO

   RETURN
END SUBROUTINE mod_trj_tl

!!======================================================================

SUBROUTINE floitptan( jpi,jpj,umod,vmod,utan,vtan,e1u,e2v,     &
   pifl,pjfl,pifltl,pjfltl,pufltl,pvfltl )
!!----------------------------------------------------------------------
!! Arguments
!! =========
   IMPLICIT NONE
   INTEGER*4 :: jpi,jpj
   REAL*8 :: pifl,pjfl,pifltl,pjfltl,pufltl,pvfltl
   REAL*8, DIMENSION(jpi,jpj) :: umod,vmod,utan,vtan,e1u,e2v
!!----------------------------------------------------------------------
!! Local declarations
!! ==================
   INTEGER :: iil,ijl,jind1,jind2
   INTEGER, DIMENSION(2) :: iid,ijd
   REAL*8, DIMENSION(2) :: zlagx,zlagy,zlagxtl,zlagytl
   REAL*8, DIMENSION(2,2) :: zuv,zuvtl
!!----------------------------------------------------------------------

!! 1. Interpolation of the zonal velocity
!! ======================================

!*... Neighbooring points (background)
   iil = INT(pifl-.5d0)
   ijl = INT(pjfl   )
   DO jind1 = 1, 2
      iid(jind1) = iil + jind1 - 1
      ijd(jind1) = ijl + jind1 - 1
   END DO

!*... Lagrange coefficients (background)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         IF( jind1.NE.jind2 ) THEN
            zlagx(jind1) = ( pifl - ((iid(jind2))+.5d0) ) / ( iid(jind1)-iid(jind2) )
            zlagy(jind1) = ( pjfl -  (ijd(jind2))     ) / ( ijd(jind1)-ijd(jind2) )
         ENDIF
      END DO
   END DO

!*... Value of the zonal velocity (background)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         zuv(jind1,jind2) = umod(iid(jind1),ijd(jind2)) / e1u(iid(jind1),ijd(jind2))
      END DO
   END DO

!*... Lagrange coefficients (tangent)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         IF( jind1.NE.jind2 ) THEN
            zlagxtl(jind1) = pifltl / ( iid(jind1)-iid(jind2) )
            zlagytl(jind1) = pjfltl / ( ijd(jind1)-ijd(jind2) )
         ENDIF
      END DO
   END DO

!*... Value of the zonal velocity (tangent)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         zuvtl(jind1,jind2) = utan(iid(jind1),ijd(jind2)) / e1u(iid(jind1),ijd(jind2))
      END DO
   END DO

!*... Interpolation of the zonal velocity
   pufltl = 0.e0
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         pufltl = pufltl                                              &
            + zuv(jind1,jind2)   * zlagxtl(jind1) * zlagy(jind2)      &
            + zuv(jind1,jind2)   * zlagx(jind1)   * zlagytl(jind2)    &
            + zuvtl(jind1,jind2) * zlagx(jind1)   * zlagy(jind2)
      END DO
   END DO

!! 2. Interpolation of the meridian velocity
!! =========================================

!*... Neighbooring points (background)
   iil = INT(pifl   )
   ijl = INT(pjfl-.5)
   DO jind1 = 1, 2
      iid(jind1) = iil + jind1 - 1
      ijd(jind1) = ijl + jind1 - 1
   END DO

!*... Lagrange coefficients (background)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         IF( jind1.NE.jind2 ) THEN
            zlagx(jind1) = ( pifl -  (iid(jind2))     ) / ( iid(jind1)-iid(jind2) )
            zlagy(jind1) = ( pjfl - ((ijd(jind2))+.5d0) ) / ( ijd(jind1)-ijd(jind2) )
         ENDIF
      END DO
   END DO

!*... Value of the meridian velocity (background)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         zuv(jind1,jind2) = vmod(iid(jind1),ijd(jind2)) / e2v(iid(jind1),ijd(jind2))
      END DO
   END DO

!*... Lagrange coefficients (tangent)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         IF( jind1.NE.jind2 ) THEN
            zlagxtl(jind1) = pifltl / ( iid(jind1)-iid(jind2) )
            zlagytl(jind1) = pjfltl / ( ijd(jind1)-ijd(jind2) )
         ENDIF
      END DO
   END DO

!*... Value of the meridian velocity (tangent)
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         zuvtl(jind1,jind2) = vtan(iid(jind1),ijd(jind2)) / e2v(iid(jind1),ijd(jind2))
      END DO
   END DO

!*... Interpolation of the meridian velocity
   pvfltl = 0.d0
   DO jind1 = 1, 2
      DO jind2 = 1, 2
         pvfltl = pvfltl                                             &
            + zuv(jind1,jind2)   * zlagxtl(jind1) * zlagy(jind2)     &
            + zuv(jind1,jind2)   * zlagx(jind1)   * zlagytl(jind2)   &
            + zuvtl(jind1,jind2) * zlagx(jind1)   * zlagy(jind2)
      END DO
   END DO

END SUBROUTINE floitptan
