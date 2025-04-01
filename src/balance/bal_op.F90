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
!>  Balance operator based on dynamic height
!!
!! It computes density anomalies from temperature and salinity
!! using tangent linear version of equation of state.
!! Eventually, it computes the dynamic height integrating density
!! anomalies from top to the reference level:
!! sla_dep variable in namelist 
!!
!                                                                      !
! Version 1: Andrea Storto 2021                                        !
!            Mario  Adani  2024                                        !
!-----------------------------------------------------------------------
SUBROUTINE bal_op

   USE set_knd
   USE grd_str
   USE drv_str
   USE bal_str
   USE cns_str, ONLY : phy

   IMPLICIT NONE

   INTEGER(i4)             :: k,j,i
   REAL(r8),ALLOCATABLE    :: rhtl(:,:)
!Function
   REAL(r8)                :: rho_unescotl

! Initilize
   ALLOCATE ( rhtl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )

!-----------------------------
!Compute density
   IF ( drv%nneos(drv%ktr) .EQ. 1 ) THEN
      grd%dns(:,:,:) = (-grd%alpha * grd%tem(:,:,:)  +                    &
                         grd%beta  * grd%sal(:,:,:) ) *  grd%msk(:,:,:)
   ELSEIF ( drv%nneos(drv%ktr) .EQ. 2 ) THEN
      grd%dns(:,:,:) = (-grd%alpha3d(:,:,:) * grd%tem(:,:,:)   +          &
                         grd%beta3d(:,:,:)  * grd%sal(:,:,:) ) *  grd%msk(:,:,:)
   ELSEIF ( drv%nneos(drv%ktr) .EQ. 3 ) THEN
      DO k = 1,grd%km
         DO j = 1-grd%jas,grd%jm+grd%jae
            DO i = 1-grd%ias,grd%im+grd%iae
               grd%dns(i,j,k) = rho_unescotl(grd%salb(i,j,k),grd%temb(i,j,k), &
                                              grd%sal(i,j,k),grd%tem(i,j,k) ) &
                                              * grd%msk(i,j,k)
            ENDDO
         ENDDO
      ENDDO
   ELSE
      write(drv%dia,*)' -----------------------------------------------'
      write(drv%dia,*)' bal_op: Unsupported  equation of state option. '
      write(drv%dia,*)' Please choose drv%nneos from 1 to 3.           '
      write(drv%dia,*)' -----------------------------------------------'
      CALL abort
   ENDIF

!-----------------------------
!Compute dynamic height
   rhtl(:,:) = 0._r8
   IF( .NOT. drv%ssh_unbalanced(drv%ktr) ) grd%eta = 0._r8
   DO k = bal%nlevs,1,-1
      rhtl(:,:) = rhtl(:,:)  + grd%dns(:,:,k)*bal%dhdz(k)
   ENDDO
   grd%eta(:,:) = grd%eta(:,:) - rhtl(:,:) / phy%rho0

   DEALLOCATE ( rhtl )

END SUBROUTINE bal_op
