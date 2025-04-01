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
!> Define filter constants, EOFs, etc.                               
!!
!! If activated from namelist it initialize recursive filter.
!! If activated from namelist it initialize diffusive filter.
!! It reads empirical orthogonal function for the model error background 
!! covariance matrix.
!!
!                                                                      !
! Version 1: Srdjan Dobricic              2006                         !
! Version 2: Srdjan Dobricic and R.Farina 2013                         !
! Version 3: Mario Adani                  2023                         !
!-----------------------------------------------------------------------
SUBROUTINE def_cov

   USE set_knd
   USE drv_str
   USE grd_str
   USE eof_str

   IMPLICIT NONE

! ---
! Horizontal Filter
   IF ( drv%filter .EQ. 1) THEN
!Recursive filter
      CALL ini_rcfl

   ELSEIF ( drv%filter .EQ. 2) THEN
!Diffusive filter
      CALL ini_dflt

   ELSEIF ( drv%filter .EQ. 3) THEN
! No filter
      WRITE (drv%dia,*)' Running with no filter '

   ELSE
! Error
      WRITE (drv%dia,*)' -------------------------------------------------------- '
      WRITE (drv%dia,*)' ERROR!!! Possible choice for filter 1,2,3. Value chosen: ',drv%filter
      WRITE (drv%dia,*)' -------------------------------------------------------- '
      CALL abort
   ENDIF

! ---
! Vertical EOFs
   ros%kmt = grd%km * 2 + 1

   CALL rdeofs

   ALLOCATE ( grd%ro( grd%im, grd%jm, ros%neof) )
   ALLOCATE ( grd%ro_ad( grd%im, grd%jm, ros%neof) )

END SUBROUTINE def_cov
