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
!> Read Normalization weight of the diffusife filter                   
!!
!! It reads normalization weights for the diffusive filter that can 
!! be different from the analytical solution close to the coasts
!!
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE rdwgh

   USE set_knd
   USE grd_str
   USE dfl_str
   USE drv_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4) :: stat,ncid,idvar
   INTEGER(i4) :: im,jm,km
   INTEGER(i4) :: start(3),count(3)

   stat = NF90_OPEN(drv%inpdir//'/'//dfl%wgh_flname, NF90_NOWRITE, ncid)

   stat = NF90_INQ_DIMID (ncid, 'im', idvar)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = im)
   stat = NF90_INQ_DIMID (ncid, 'jm', idvar)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = jm)
   stat = NF90_INQ_DIMID (ncid, 'km', idvar)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = km)

   IF ( ( im .NE. grd%img) .OR. (jm .NE. grd%jmg) .OR. (km .NE. grd%km) ) THEN
      WRITE (drv%dia,*),'Dimensions of weight grid differs'
      STOP
   ENDIF

   ALLOCATE(dfl%wgh(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))

   start(1) = grd%igs-grd%ias
   start(2) = grd%jgs-grd%jas
   start(3) = 1
   count(1) = grd%im + grd%ias + grd%iae
   count(2) = grd%jm + grd%jas + grd%jae
   count(3) = grd%km

   stat = NF90_INQ_VARID (ncid, 'coef', idvar)
   stat = NF90_GET_VAR (ncid,idvar,dfl%wgh,start(1:3),count(1:3))
   stat = NF90_CLOSE(ncid)

END SUBROUTINE rdwgh
