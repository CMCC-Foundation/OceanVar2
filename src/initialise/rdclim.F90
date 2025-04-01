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
!> Read climatology                                                    
!!
!! It reads the climatology for quality control of observations
!!
!                                                                      !
! Version 1: Mario Adani 2024                                          !
!-----------------------------------------------------------------------
subroutine rdclim

   use set_knd
   use grd_str
   use obs_str,      only : qck
   use drv_str
   use netcdf
   use datetime_module, only : datetime, date2num, num2date


   implicit none


   INTEGER(i4)                           :: stat, ncid, idvar
   INTEGER(i4)                           :: img, jmg, km
   REAL(r4), ALLOCATABLE                 :: x3(:,:,:)
   INTEGER(i4)                           :: start(4), count(4)
   TYPE(datetime)                        :: daterun

   daterun = num2date(drv%zanjul1950 + date2num(datetime(1950,1,01,0)))

! Open file
   stat = nf90_open(drv%inpdir//'/'//qck%flname, NF90_NOWRITE, ncid)
   if (stat /= nf90_noerr) call netcdf_err(stat)
! Get dimensions
   stat = nf90_inq_dimid (ncid, 'im', idvar)
   if (stat /= nf90_noerr) call netcdf_err(stat)
   stat = nf90_inquire_dimension (ncid, idvar, len = img)
   if (stat /= nf90_noerr) call netcdf_err(stat)
   stat = nf90_inq_dimid (ncid, 'jm', idvar)
   if (stat /= nf90_noerr) call netcdf_err(stat)
   stat = nf90_inquire_dimension (ncid, idvar, len = jmg)
   if (stat /= nf90_noerr) call netcdf_err(stat)
   stat = nf90_inq_dimid (ncid, 'km', idvar)
   if (stat /= nf90_noerr) call netcdf_err(stat)
   stat = nf90_inquire_dimension (ncid, idvar, len = km)
   if (stat /= nf90_noerr) call netcdf_err(stat)

   if( ( grd%img .ne. img ) .OR. &
      ( grd%jmg .ne. jmg ) .OR. &
      ( grd%km  .ne. km ) ) then
      write(drv%dia,*)' ------------------------------------------ '
      write(drv%dia,*)' Wrong dimension of the climatological grid '
      write(drv%dia,*)' ------------------------------------------ '
      call abort
   endif

   start(1) = grd%igs-grd%ias
   start(2) = grd%jgs-grd%jas
   start(3) = 1
   start(4) = daterun%getMonth()
   count(1) = grd%im + grd%ias + grd%iae
   count(2) = grd%jm + grd%jas + grd%jae
   count(3) = grd%km
   count(4) = 1

   ALLOCATE ( x3         (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
   ALLOCATE ( qck%climtem(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
   ALLOCATE ( qck%climsal(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))

   stat = nf90_inq_varid (ncid, 'vosaline', idvar)
   stat = nf90_get_var (ncid,idvar,x3,start,count)
   qck%climsal(:,:,:) = dble(x3(:,:,:))

   stat = nf90_inq_varid (ncid, 'votemper', idvar)
   stat = nf90_get_var (ncid,idvar,x3,start,count)
   qck%climtem(:,:,:) = dble(x3(:,:,:))

   stat = nf90_close(ncid)

   DEALLOCATE ( x3 )

end subroutine rdclim
