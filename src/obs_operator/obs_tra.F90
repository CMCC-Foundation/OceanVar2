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
!> Call Argo trajectory model                                          
!!
!!
!!
!                                                                      !
! Version 1: Vincent Taillandier, Srdjan Dobricic 2007                 !
!                                                                      !
!-----------------------------------------------------------------------
SUBROUTINE obs_tra

   USE set_knd
   USE grd_str
   USE obs_str
   USE mpi_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, img, jmg, k, km

   IF ( tra%ncc .GT. 0 .OR. tra%ncs .GT. 0) THEN

      IF ( mpi%myrank .EQ. 0 ) THEN
         img = grd%img
         jmg = grd%jmg
      ELSE
         img = 1
         jmg = 1
      ENDIF

      k   = tra%lev
      km  = grd%km

      IF ( mpi%nproc .GT. 1 ) THEN
         CALL gth_mpi( img, jmg, k, km, grd%uvl, tra%uvl)
         CALL gth_mpi( img, jmg, k, km, grd%vvl, tra%vvl)
      ELSE
         tra%uvl(:,:) = grd%uvl(:,:,k)
         tra%vvl(:,:) = grd%vvl(:,:,k)
      ENDIF

      IF ( mpi%myrank .EQ. 0 ) THEN
         CALL mod_trj_tl( tra%im,tra%jm,tra%umn,tra%vmn,tra%dx,tra%dy,tra%flc, tra%fls, &
                          tra%nt,tra%no,tra%xmn,tra%ymn,tra%dtm,                        &
                          tra%uvl,tra%vvl,tra%xtl,tra%ytl )
         DO k = 1,tra%no
            IF ( tra%flc(k) .EQ. 1 .OR. tra%fls(k) .EQ. 1 ) THEN
               tra%inx(k) = tra%xtl(k)
               tra%iny(k) = tra%ytl(k)
               i = INT(tra%xmn(tra%nt+1,k)+tra%xtl(k))
               j = INT(tra%ymn(tra%nt+1,k)+tra%ytl(k))
            ENDIF
         ENDDO
      ENDIF

   ENDIF

END SUBROUTINE obs_tra
