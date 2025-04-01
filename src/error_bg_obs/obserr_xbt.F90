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
!> XBT observational error                                              
!!
!! It applys the observational error. 
!! It can be spatially variable if obserr%ts_ver_dep_err is true in 
!! namelist. It consists of a vertical profile for temperature and
!! salinity obserr%tem multiplied for a horizontal dependend
!! factor tem_fc.
!! If obserr%ts_ver_dep_err is false error is spatially constant and it
!! is defined in the namelist tem_con_err.
!! If time dependent error (tim_dep_err) is true in the namelist it 
!! increase the observational error far from the analysis time
!                                                                      !
! Version 1:  Andrea Storto 2022                                       !
!             Mario Adani   2024                                       !
!-----------------------------------------------------------------------
SUBROUTINE obserr_xbt

   USE set_knd
   USE obs_str, ONLY : xbt, obserr
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k
   REAL(r8)      :: rdep, tdist, cf
   INTEGER(i4)   :: ndepths, idep
!FUNCTION
   REAL(r8)      :: tim_dep_errfact

   IF (obserr%ts_ver_dep_err) THEN
      rdep    = obserr%depth(2)-obserr%depth(1)
      ndepths = size(obserr%depth,1)
   ENDIF

   DO k = 1,xbt%no
      IF ( xbt%flc(k) .EQ. 1 ) THEN

         IF (obserr%ts_ver_dep_err) THEN ! Vertical dependent error
            idep=nint(xbt%dpt(k)/rdep)+1
            IF ( idep .GT. ndepths ) idep=ndepths
            i = xbt%ib(k)
            j = xbt%jb(k)

            cf  = xbt%pq1(k) * obserr%tem_fc(i  ,j  ) +       &
                  xbt%pq2(k) * obserr%tem_fc(i+1,j  ) +       &
                  xbt%pq3(k) * obserr%tem_fc(i  ,j+1) +       &
                  xbt%pq4(k) * obserr%tem_fc(i+1,j+1) +       &
                  xbt%pq5(k) * obserr%tem_fc(i  ,j  ) +       &
                  xbt%pq6(k) * obserr%tem_fc(i+1,j  ) +       &
                  xbt%pq7(k) * obserr%tem_fc(i  ,j+1) +       &
                  xbt%pq8(k) * obserr%tem_fc(i+1,j+1)
            xbt%err(k) = obserr%tem(idep) * cf
         ELSE   ! constant error
            xbt%err(k) = obserr%tem_con_err
         ENDIF 

! Time dependent factor
         IF ( obserr%tim_dep_err ) THEN
            tdist = xbt%tim(k) - drv%zanjul1950
            xbt%err(k)= tim_dep_errfact(tdist) * xbt%err(k)
         ENDIF

      ENDIF   ! flc
   ENDDO

END SUBROUTINE obserr_xbt
