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
!> Write outputs and diagnostics for SLA
!!
!!
!!
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE wrt_sla

   USE set_knd
   USE obs_str, ONLY : sla
   USE mpi_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i8), ALLOCATABLE      ::  allflg(:), alleve(:)
   REAL(r8),    ALLOCATABLE      ::  allinc(:)
   INTEGER(i4)                   ::  stat, ncid, dimid, ierr
   INTEGER(i4)                   ::  varid(15) ! number of field to be saved
   INTEGER(i4)                   ::  nobs

   INCLUDE "mpif.h"

   ALLOCATE ( allflg(sla%no), allinc(sla%no),  alleve(sla%no) )

   alleve(:)  = 0_i8
   allflg(:)  = 0_i8
   allinc(:)  = 0.0_r8
   nobs = sla%no
   CALL MPI_REDUCE( sla%inc, allinc, nobs, mpi%r8  ,   &
                    MPI_SUM, 0, mpi%comm, ierr)
   CALL MPI_REDUCE( sla%flc, allflg, nobs, mpi%i8  ,   &
                    MPI_MAX, 0, mpi%comm, ierr)
   CALL MPI_REDUCE( sla%eve, alleve, nobs, mpi%i8  ,   &
                    MPI_MAX, 0, mpi%comm, ierr)

   IF ( mpi%myrank .EQ. 0 ) THEN

      sla%flg(1:sla%no) = allflg(1:sla%no)
      sla%inc(1:sla%no) = allinc(1:sla%no)
      sla%eve(1:sla%no) = alleve(1:sla%no)

      ! Create file
      stat = NF90_CREATE('obsstat_sla.nc', NF90_SHARE, ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_OPEN  ('obsstat_sla.nc', NF90_WRITE, ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      stat = NF90_REDEF(ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      ! Define dimension
      stat = NF90_DEF_DIM(ncid, 'nobs',   nobs,   dimid )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      ! Define variables
      stat = NF90_DEF_VAR(ncid, 'lon',    NF90_FLOAT, (/dimid/), varid(1) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(1), 'units', 'degree_east')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(1), 'long_name', 'longitude')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(1), 'standard_name', 'longitude')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'lat',    NF90_FLOAT, (/dimid/), varid(2) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(2), 'units', 'degree_north')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(2), 'long_name', 'latitude')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(2), 'standard_name', 'latitude')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'time',   NF90_FLOAT, (/dimid/), varid(3) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(3), 'units', 'days since 1950-01-01')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(3), 'long_name', 'time')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(3), 'standard_name', 'time')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(3), 'calendar', 'standard')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'track',   NF90_INT, (/dimid/), varid(4) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(4), 'long_name', 'track number')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(4), 'standard_name', 'track number')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'satnm',   NF90_INT,  (/dimid/), varid(5) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(5), 'long_name', 'satellite number')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(5), 'standard_name', 'satellite number')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'flag',    NF90_INT, (/dimid/),  varid(6) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(6), 'long_name', 'quality flag')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(6), 'standard_name', 'quality flag')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'eve',    NF90_INT, (/dimid/),  varid(7) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(7), 'long_name', 'event')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(7), 'standard_name', 'event')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'value',  NF90_FLOAT, (/dimid/),  varid(8) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(8), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(8), 'long_name', 'observational value')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(8), 'standard_name', 'observational value')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'bkg  ',  NF90_FLOAT, (/dimid/),  varid(9) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(9), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(9), 'long_name', 'background value')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(9), 'standard_name', 'background value')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'error',  NF90_FLOAT, (/dimid/),  varid(10) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(10), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(10), 'long_name', 'observational error')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(10), 'standard_name', 'observational error')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'res',    NF90_FLOAT, (/dimid/), varid(11) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(11), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(11), 'long_name', 'observation-analysis')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(11), 'standard_name', 'observation-analysis')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'bias',   NF90_FLOAT, (/dimid/), varid(12) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(12), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(12), 'long_name', 'bias')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(12), 'standard_name', 'bias')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'incr'   ,NF90_FLOAT, (/dimid/), varid(13) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(13), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(13), 'long_name', 'increment')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(13), 'standard_name', 'increment')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'dpt',    NF90_FLOAT, (/dimid/), varid(14) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(14), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(14), 'long_name', 'water column depth')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(14), 'standard_name', 'water column depth')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_DEF_VAR(ncid, 'mdt',    NF90_FLOAT, (/dimid/), varid(15) )
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(15), 'units', 'meter')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(15), 'long_name', 'mean dynamic topography')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_ATT(ncid, varid(15), 'standard_name', 'mean dynamic topography')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)


      ! Close definition
      stat = NF90_PUT_ATT( ncid, NF90_GLOBAL, 'Product',  'OceanVar Variational Analysis')
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_ENDDEF(ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      ! Fill in variables
      stat = NF90_PUT_VAR(ncid, varid(1), sla%lon)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(2), sla%lat)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(3), sla%tim)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(4), sla%ino)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(5), sla%ksat)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(6), sla%flg)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(7), sla%eve)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(8), sla%val)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(9), sla%bac)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(10),sla%err)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(11),sla%res)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(12),sla%bia)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(13),sla%inc)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(14),sla%dpt)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_PUT_VAR(ncid, varid(15),sla%dtm)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

      ! Close file
      stat = NF90_CLOSE(ncid)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

   ENDIF ! rank

   DEALLOCATE ( allflg, allinc,  alleve )
   CALL mpi_barrier(mpi%comm,ierr)

END SUBROUTINE wrt_sla
