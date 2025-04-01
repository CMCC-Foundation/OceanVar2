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
!----------------------------------------------------------------------------------
!                                                                      
!> Read Temperature Salinity field  to compute expantion/contraction coefficients.
!!
!!
!!
!                                                                                 !
! Version 1: Mario Adani 2023                                                     !
!----------------------------------------------------------------------------------
SUBROUTINE rdeos

   USE set_knd
   USE drv_str
   USE grd_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                           :: stat, ncid, idvar
   INTEGER(i4)                           :: img, jmg, km
   REAL(r4), ALLOCATABLE                 :: x3(:,:,:)
   INTEGER(i4)                           :: start(3), count(3)
   INTEGER(i4)                           :: i,j,k

! Open file
   stat = NF90_OPEN(drv%inpdir//'/'//drv%eosflname, NF90_NOWRITE, ncid)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
! Get dimensions
   stat = NF90_INQ_DIMID (ncid, 'im', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = img)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'jm', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = jmg)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'km', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = km)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

   IF ( ( grd%img .NE. img ) .OR. &
        ( grd%jmg .NE. jmg ) .OR. &
        ( grd%km  .NE. km ) ) THEN
      WRITE (drv%dia,*)' ------------------------------------------ '
      WRITE (drv%dia,*)' Wrong dimension of the climatological grid '
      WRITE (drv%dia,*)' ------------------------------------------ '
      CALL abort
   ENDIF

   start(1) = grd%igs-grd%ias
   start(2) = grd%jgs-grd%jas
   start(3) = 1
   count(1) = grd%im + grd%ias + grd%iae
   count(2) = grd%jm + grd%jas + grd%jae
   count(3) = grd%km

   ALLOCATE ( x3         (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%temb(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%salb(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )

   stat = NF90_INQ_VARID (ncid, 'vosaline', idvar)
   stat = NF90_GET_VAR (ncid,idvar,x3,start,count)
   grd%salb(:,:,:) = DBLE(x3(:,:,:))

   stat = NF90_INQ_VARID (ncid, 'votemper', idvar)
   stat = NF90_GET_VAR (ncid,idvar,x3,start,count)
   grd%temb(:,:,:) = DBLE(x3(:,:,:))

   stat = NF90_CLOSE(ncid)

   DEALLOCATE ( x3 )

   DO k = 1,grd%km
      DO j = 1-grd%jas,grd%jm+grd%jae
         DO i = 1-grd%ias,grd%im+grd%iae
            grd%temb(i,j,k) = grd%temb(i,j,k) * grd%msk(i,j,k)
            grd%salb(i,j,k) = grd%salb(i,j,k) * grd%msk(i,j,k)
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE rdeos
