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
!>  Buoyancy adjoint check                                            
!!
!! It checks if the linear tangent and the adjoint are consistent.
!! It can only be activated without domain decomposition.
!! It is used for debugging purposes only.
!!
!                                                                      !
! Version 1: M.Adani    2023                                           !
!-----------------------------------------------------------------------
SUBROUTINE adjck_byg

   USE set_knd
   USE grd_str
   USE drv_str

   IMPLICIT NONE

   REAL(r8),ALLOCATABLE ::  by(:,:)         ! Buoyancy
   REAL(r8),ALLOCATABLE ::  bx(:,:)         ! Buoyancy
   REAL(r8),ALLOCATABLE ::  b_y(:,:,:)      ! Buoyancy
   REAL(r8),ALLOCATABLE ::  b_x(:,:,:)      ! Buoyancy
   REAL(r8),ALLOCATABLE ::  tem(:,:,:)      ! Temperature increment
   REAL(r8),ALLOCATABLE ::  sal(:,:,:)      ! Salinity increment
   REAL(r8),ALLOCATABLE ::  tem_ad(:,:,:)   ! Temperature adjoint
   REAL(r8),ALLOCATABLE ::  sal_ad(:,:,:)   ! Salinity adjoint
   REAL(r8),ALLOCATABLE ::  temsal_1d(:)
   REAL(r8),ALLOCATABLE ::  temsal_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  bxby_1d(:)
   REAL(r8),ALLOCATABLE ::  bxby_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  dummy3(:,:,:)
   REAL(r8),ALLOCATABLE ::  dummy4(:,:,:,:)
   REAL(r8)             ::  a,b

!Allocate Variables
   ALLOCATE ( dummy4(grd%img,grd%jmg,grd%km,2)       )
   ALLOCATE ( dummy3(grd%img,grd%jmg,2)              )
   ALLOCATE ( temsal_1d(grd%img*grd%jmg*grd%km*2)    )
   ALLOCATE ( temsal_ad_1d(grd%img*grd%jmg*grd%km*2) )
   ALLOCATE ( bxby_1d(grd%img*grd%jmg*2)             )
   ALLOCATE ( bxby_ad_1d(grd%img*grd%jmg*2)          )
   ALLOCATE ( tem(grd%img,grd%jmg,grd%km)            )
   ALLOCATE ( sal(grd%img,grd%jmg,grd%km)            )
   ALLOCATE ( tem_ad(grd%img,grd%jmg,grd%km)         )
   ALLOCATE ( sal_ad(grd%img,grd%jmg,grd%km)         )
   ALLOCATE ( bx(grd%img,grd%jmg)                    )
   ALLOCATE ( by(grd%img,grd%jmg)                    )
   ALLOCATE ( b_x(grd%img,grd%jmg,grd%km)            )
   ALLOCATE ( b_y(grd%img,grd%jmg,grd%km)            )

!Store Variables
   b_x(   :,:,:) = grd%b_x(:,:,:)
   b_y(   :,:,:) = grd%b_y(:,:,:)
   tem(   :,:,:) = grd%tem(:,:,:)
   sal(   :,:,:) = grd%sal(:,:,:)
   tem_ad(:,:,:) = grd%tem_ad(:,:,:)
   sal_ad(:,:,:) = grd%sal_ad(:,:,:)
   bx      (:,:) = grd%bx(:,:)
   by      (:,:) = grd%by(:,:)

   grd%b_x   (:,:,:) = 0.0_r8
   grd%b_y   (:,:,:) = 0.0_r8
   grd%tem   (:,:,:) = 0.0_r8
   grd%sal   (:,:,:) = 0.0_r8
   grd%tem_ad(:,:,:) = 0.0_r8
   grd%sal_ad(:,:,:) = 0.0_r8
   grd%dns(   :,:,:) = 0.0_r8
   grd%bx      (:,:) = 0.0_r8
   grd%by      (:,:) = 0.0_r8
   dummy4  (:,:,:,:) = 0.0_r8
   dummy3    (:,:,:) = 0.0_r8

   CALL RANDOM_NUMBER(grd%tem(:,:,:))
   CALL RANDOM_NUMBER(grd%sal(:,:,:))

   grd%tem(:,:,:)  = grd%tem (:,:,:) * grd%msk(:,:,:)
   grd%sal(:,:,:)  = grd%sal (:,:,:) * grd%msk(:,:,:)
   dummy4(:,:,:,1) = grd%tem (:,:,:)
   dummy4(:,:,:,2) = grd%sal (:,:,:)
   temsal_1d(:)  = RESHAPE(dummy4,(/grd%img*grd%jmg*grd%km*2/))

   CALL get_byg
   dummy3(:,:,1) = grd%bx(:,:)
   dummy3(:,:,2) = grd%by(:,:)
   bxby_1d   (:)  = RESHAPE(dummy3,(/grd%img*grd%jmg*2/))

   grd%b_x   (:,:,:) = 0.0_r8
   grd%b_y   (:,:,:) = 0.0_r8
   grd%tem   (:,:,:) = 0.0_r8
   grd%sal   (:,:,:) = 0.0_r8
   grd%dns   (:,:,:) = 0.0_r8
   grd%bx      (:,:) = 0.0_r8
   grd%by      (:,:) = 0.0_r8
   dummy4  (:,:,:,:) = 0.0_r8
   dummy3    (:,:,:) = 0.0_r8

   CALL RANDOM_NUMBER(grd%bx(:,:))
   CALL RANDOM_NUMBER(grd%by(:,:))
   grd%bx(:,:)   = grd%bx(:,:) * grd%msk(:,:,1)
   grd%by(:,:)   = grd%by(:,:) * grd%msk(:,:,1)
   dummy3(:,:,1) = grd%bx(:,:)
   dummy3(:,:,2) = grd%by(:,:)
   bxby_ad_1d(:)    = RESHAPE(dummy3,(/grd%img*grd%jmg*2/))

   b = DOT_PRODUCT(bxby_1d,bxby_ad_1d)

   CALL get_byg_ad

   dummy4(:,:,:,1)  = grd%tem_ad(:,:,:)
   dummy4(:,:,:,2)  = grd%sal_ad(:,:,:)
   temsal_ad_1d(:) = RESHAPE(dummy4,(/grd%img*grd%jmg*grd%km*2/))

   a = DOT_PRODUCT(temsal_1d,temsal_ad_1d)

   WRITE (drv%dia,*) '-------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK BUOYANCY MODEL        '
   WRITE (drv%dia,*) ' <grd%temsal,grd%temsal_ad>        = ',a
   WRITE (drv%dia,*) ' <grd%bxby,grd%bxby_ad>            = ',b
   WRITE (drv%dia,*) ' Absolute and relative dIFference  = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '-------------------------------------'
   CALL FLUSH(drv%dia)

!ReStore Variables
   grd%tem(   :,:,:) = tem(:,:,:)
   grd%sal(   :,:,:) = sal(:,:,:)
   grd%tem_ad(:,:,:) = tem_ad(:,:,:)
   grd%sal_ad(:,:,:) = sal_ad(:,:,:)
   grd%bx    (:,:  ) = bx(:,:)
   grd%by    (:,:  ) = by(:,:)
   grd%b_x   (:,:,:) = b_x(:,:,:)
   grd%b_y   (:,:,:) = b_y(:,:,:)

!DeAllocate Variables
   DEALLOCATE ( tem )
   DEALLOCATE ( sal )
   DEALLOCATE ( b_x )
   DEALLOCATE ( b_y )
   DEALLOCATE ( tem_ad )
   DEALLOCATE ( sal_ad )
   DEALLOCATE ( by )
   DEALLOCATE ( bx )
   DEALLOCATE ( dummy4, dummy3, temsal_1d, temsal_ad_1d, bxby_1d, bxby_ad_1d )

END SUBROUTINE adjck_byg
