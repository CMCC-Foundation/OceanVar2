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
!> Call drifter trajectory model                                        
!!
!!
!!
!                                                                      !
! Version 1: Vincent Taillandier, Srdjan Dobricic 2007                 !
!                                                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_trd

   USE set_knd
   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, img, jmg, k, km

   IF ( trd%ncc .GT. 0 .OR. trd%ncs .GT. 0 ) THEN

      IF ( mpi%myrank .EQ. 0 ) THEN
         img = grd%img
         jmg = grd%jmg
      ELSE
         img = 1
         jmg = 1
      ENDIF

      k   = trd%lev
      km  = grd%km

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL gth_mpi( img, jmg, k, km, grd%uvl, trd%uvl)
         CALL gth_mpi( img, jmg, k, km, grd%vvl, trd%vvl)
      ELSE
         trd%uvl(:,:) = grd%uvl(:,:,k)
         trd%vvl(:,:) = grd%vvl(:,:,k)
      ENDIF

      IF ( mpi%myrank .EQ. 0 ) THEN

         CALL mod_trj_tl( trd%im,trd%jm,trd%umn,trd%vmn,trd%dx,trd%dy,trd%flc, trd%fls, &
            trd%nt,trd%no,trd%xmn,trd%ymn,trd%dtm,                        &
            trd%uvl,trd%vvl,trd%xtl,trd%ytl )

         DO k = 1,trd%no
            IF ( trd%flc(k) .EQ. 1 .OR. trd%fls(k) .EQ. 1 ) THEN
               trd%inx(k) = trd%xtl(k)
               trd%iny(k) = trd%ytl(k)
               i = INT(trd%xmn(trd%nt+1,k)+trd%xtl(k))
               j = INT(trd%ymn(trd%nt+1,k)+trd%ytl(k))
            ENDIF
         ENDDO

      ENDIF

   ENDIF

END SUBROUTINE obs_trd
