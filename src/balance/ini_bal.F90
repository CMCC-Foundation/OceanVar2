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
!                                                                      
!>   Inizialization of simplified balance operator based on dynamic height
!!
!! It defines reference level based on sla_dep value read from namelist
!! and layer width.
!!
!                                                                    
! Version 1: Andrea Storto 2021                                       
!            Mario  Adani  2024                                        
!-----------------------------------------------------------------------
SUBROUTINE ini_bal

   USE set_knd
   USE bal_str
   USE obs_str, ONLY : sla
   USE grd_str
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)    :: k,i,j
   REAL(r8)       :: zdz

!sla%dep becomes level of no motion > 1000 !!!
   bal%lnm = sla%dep

   ALLOCATE ( bal%dhdz(grd%km) )

   zdz = 0._r8
   bal%dhdz(:) = grd%dz(:)

   DO k=1,grd%km
      zdz = zdz + grd%dz(k)
      IF ( zdz .LT. bal%lnm ) THEN
         bal%nlevs = k
      ENDIF
   ENDDO

   bal%nlevs   = bal%nlevs+1
   bal%dhdz(bal%nlevs) = SUM( grd%dz(1:bal%nlevs) ) - bal%lnm

   WRITE(drv%dia,*) ' Level of no motion: Level/Depth: ',bal%nlevs,'/',bal%lnm
   IF ( bal%lnm < 1000._r8) WRITE(drv%dia,*) 'WARINIG: No motion level is less than 1000 meters'

   CALL FLUSH(drv%dia)

END SUBROUTINE ini_bal
