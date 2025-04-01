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
!> Initialize time of analysis                                         
!!
!! It initializes the time variable to weight the observational error
!! dependent from the time distance between the analysis and the 
!! the observation retrival.
!! Reference date in julian day start from January, 1st 1950
!!
!                                                                      !
! Version 1: Mario Adani 2023                                          !
!-----------------------------------------------------------------------
SUBROUTINE ini_time

   USE set_knd
   USE drv_str
   USE datetime_module, ONLY : datetime, date2num

   IMPLICIT NONE

   REAL(r8)          :: zjul1950
   CHARACTER*8       :: cdate
   CHARACTER*4       :: cyear
   CHARACTER*2       :: cmonth
   CHARACTER*2       :: cday
   CHARACTER*2       :: chour
   INTEGER(i4)       :: year, month, day, hour

   WRITE (cdate,"(I8)")drv%sdat
   WRITE (chour,"(I2)")drv%shou

   cyear  = cdate(1:4)
   cmonth = cdate(5:6)
   cday   = cdate(7:8)

   READ (cyear,"(I4)" )year
   READ (cmonth,"(I2)")month
   READ (cday,"(I2)"  )day
   READ (chour,"(I2)" )hour

!Reference date
   zjul1950   = date2num(datetime(1950,1,01,0))
   drv%zanjul1950 = date2num(datetime(year,month,day,hour)) - zjul1950

END SUBROUTINE ini_time

