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
!> SLA observational error                                            
!!
!! It applys the observational error. 
!! It can be satellite dependent if obserr%sla_sat_dep_err is true. 
!! It can be spatially variable if obserr%sla_hor_dep_err is true in 
!! namelist. 
!! If time dependent error (tim_dep_err) is true in the namelist it 
!! increase the observational error far from the analysis time
!! Error values are defined in the namelist.
!! Errors are additive.
!                                                                      !
! Version 1:  Andrea Storto 2022                                       !
!             Mario Adani   2024                                       !
!-----------------------------------------------------------------------
SUBROUTINE obserr_sla

   USE set_knd
   USE obs_str, ONLY : sla, obserr
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)   ::  i, j, k
   REAL(r8)      :: tdist, err
!FUNCTION
   REAL(r8)      :: tim_dep_errfact

!Satellite dependent error
   IF ( obserr%sla_sat_dep_err ) THEN
      DO i = 1,obserr%sla_sat_nu
         WHERE (sla%ksat .EQ. i ) sla%err = obserr%sla_sat_err(i)
      ENDDO
   ELSE
      sla%err(:) = obserr%sla_con_err
   ENDIF

!Horizontaly dependent error
   IF ( obserr%sla_hor_dep_err ) THEN
      DO k = 1,sla%no
         IF ( sla%flc(k) .EQ. 1) THEN
            i = sla%ib(k)
            j = sla%jb(k)
            err        = sla%pq1(k) * obserr%stde(i  ,j  ) +       &
                         sla%pq2(k) * obserr%stde(i+1,j  ) +       &
                         sla%pq3(k) * obserr%stde(i  ,j+1) +       &
                         sla%pq4(k) * obserr%stde(i+1,j+1)
            err        = obserr%sla_hor_dep_ecf * MAX(MIN(err,0.5_r8),0.005_r8)
            sla%err(k) = sla%err(k) + err
         ENDIF
      ENDDO
   ENDIF

! Time dependent factor
   IF ( obserr%tim_dep_err ) THEN
      DO k = 1,sla%no
         tdist = sla%tim(k) - drv%zanjul1950
         sla%err(k)= tim_dep_errfact(tdist) * sla%err(k)
      ENDDO
   ENDIF

END SUBROUTINE obserr_sla

