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
!> Simplified balance operator adjoint check 
!                                                                     
! It checks if the linear tangent and the adjoint are consistent.
! It can only be activated without domain decomposition.
! It is used for debugging purposes only.
!
! Version 1: M.Adani    2023                                           !
!-----------------------------------------------------------------------
SUBROUTINE adjck_balop

   USE set_knd
   USE grd_str
   USE drv_str

   IMPLICIT NONE

   REAL(r8),ALLOCATABLE ::  tem(:,:,:)      ! Temperature increment
   REAL(r8),ALLOCATABLE ::  sal(:,:,:)      ! Salinity increment
   REAL(r8),ALLOCATABLE ::  eta(:,:)        ! Sea level increment
   REAL(r8),ALLOCATABLE ::  tem_ad(:,:,:)   ! Temperature adjoint
   REAL(r8),ALLOCATABLE ::  sal_ad(:,:,:)   ! Salinity adjoint
   REAL(r8),ALLOCATABLE ::  eta_ad(:,:)     ! Sea.LE.el adjoint
   REAL(r8),ALLOCATABLE ::  eta_1d(:)
   REAL(r8),ALLOCATABLE ::  temsal_1d(:)
   REAL(r8),ALLOCATABLE ::  eta_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  temsal_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  dummy(:,:,:,:)
   REAL(r8)             ::  a,b

!Allocate Variables
   ALLOCATE ( dummy(grd%img,grd%jmg,grd%km,2)        )
   ALLOCATE ( temsal_1d(grd%img*grd%jmg*grd%km*2)    )
   ALLOCATE ( temsal_ad_1d(grd%img*grd%jmg*grd%km*2) )
   ALLOCATE ( eta_1d(grd%img*grd%jmg)                )
   ALLOCATE ( eta_ad_1d(grd%img*grd%jmg)             )
   ALLOCATE ( tem(grd%img,grd%jmg,grd%km)            )
   ALLOCATE ( sal(grd%img,grd%jmg,grd%km)            )
   ALLOCATE ( eta(grd%img,grd%jmg)                   )
   ALLOCATE ( tem_ad(grd%img,grd%jmg,grd%km)         )
   ALLOCATE ( sal_ad(grd%img,grd%jmg,grd%km)         )
   ALLOCATE ( eta_ad(grd%img,grd%jmg)                )

!Store Variables
   tem(   :,:,:) = grd%tem(:,:,:)
   sal(   :,:,:) = grd%sal(:,:,:)
   eta(   :,:  ) = grd%eta(:,:)
   tem_ad(:,:,:) = grd%tem_ad(:,:,:)
   sal_ad(:,:,:) = grd%sal_ad(:,:,:)
   eta_ad(:,:  ) = grd%eta_ad(:,:)

!Initialize Variable
   grd%tem   (:,:,:) = 0.0_r8
   grd%sal   (:,:,:) = 0.0_r8
   grd%eta     (:,:) = 0.0_r8
   grd%tem_ad(:,:,:) = 0.0_r8
   grd%sal_ad(:,:,:) = 0.0_r8
   grd%eta_ad  (:,:) = 0.0_r8

!Perturbation
   CALL RANDOM_NUMBER(grd%tem(:,:,:))
   CALL RANDOM_NUMBER(grd%sal(:,:,:))
   CALL RANDOM_NUMBER(grd%eta_ad(:,:))
   grd%tem(:,:,:)  = grd%tem (:,:,:) *grd%msk(:,:,:)
   grd%sal(:,:,:)  = grd%sal (:,:,:) *grd%msk(:,:,:)
   grd%eta_ad(:,:) = grd%eta_ad(:,:) *grd%msk(:,:,1)

!Getting grd%eta
   CALL bal_op

!Getting grd%tem_ad, grd%sal_ad
   CALL bal_op_ad

!Reshaping
   eta_ad_1d(:)    = RESHAPE(grd%eta_ad,(/grd%img*grd%jmg/))
   eta_1d   (:)    = RESHAPE(grd%eta,(/grd%img*grd%jmg/))
   dummy(:,:,:,1)  = grd%tem(:,:,:)
   dummy(:,:,:,2)  = grd%sal(:,:,:)
   temsal_1d  (:)  = RESHAPE(dummy,(/grd%img*grd%jmg*grd%km*2/))
   dummy(:,:,:,1)  = grd%tem_ad(:,:,:)
   dummy(:,:,:,2)  = grd%sal_ad(:,:,:)
   temsal_ad_1d(:) = RESHAPE(dummy,(/grd%img*grd%jmg*grd%km*2/))

   a = DOT_PRODUCT(temsal_1d,temsal_ad_1d)
   b = DOT_PRODUCT(eta_1d,eta_ad_1d)

   WRITE (drv%dia,*) '-------------------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK SIMPLIFIED BALANCE OPERATOR:      '
   WRITE (drv%dia,*) ' <grd%temsal,grd%temsal_ad>                    = ',a
   WRITE (drv%dia,*) ' <grd%eta,grd%eta_ad>                          = ',b
   WRITE (drv%dia,*) ' Absolute and relative difference              = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '-------------------------------------------------'
   CALL FLUSH(drv%dia)

!Restore Variables
   grd%tem(   :,:,:) = tem(:,:,:)
   grd%sal(   :,:,:) = sal(:,:,:)
   grd%eta(   :,:  ) = eta(:,:)
   grd%tem_ad(:,:,:) = tem_ad(:,:,:)
   grd%sal_ad(:,:,:) = sal_ad(:,:,:)
   grd%eta_ad(:,:  ) = eta_ad(:,:)

!DeAllocate Variables
   DEALLOCATE ( tem )
   DEALLOCATE ( sal )
   DEALLOCATE ( eta )
   DEALLOCATE ( tem_ad )
   DEALLOCATE ( sal_ad )
   DEALLOCATE ( eta_ad )
   DEALLOCATE ( dummy, temsal_1d, temsal_ad_1d, eta_1d, eta_ad_1d )

END SUBROUTINE adjck_balop
