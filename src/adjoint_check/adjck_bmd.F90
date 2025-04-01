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
!>  Barotropic Model adjoint check                     
!!
!! It checks if the linear tangent and the adjoint are consistent.
!! It can only be activated without domain decomposition.
!! It is used for debugging purposes only.
!!
!                                                                      !
! Version 1: M.Adani    2023                                           !
!-----------------------------------------------------------------------
SUBROUTINE adjck_bmd

   USE set_knd
   USE grd_str
   USE drv_str

   IMPLICIT NONE

   REAL(r8),ALLOCATABLE ::  eta(:,:)        ! Sea.LE.el increment
   REAL(r8),ALLOCATABLE ::  bx(:,:)         ! Buoyancy
   REAL(r8),ALLOCATABLE ::  by(:,:)         ! Buoyancy
   REAL(r8),ALLOCATABLE ::  eta_ad(:,:)     ! Sea.LE.el adjoint
   REAL(r8),ALLOCATABLE ::  eta_1d(:)
   REAL(r8),ALLOCATABLE ::  eta_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  bxby_1d(:)
   REAL(r8),ALLOCATABLE ::  bxby_ad_1d(:)
   REAL(r8),ALLOCATABLE ::  dummy3(:,:,:)
   REAL(r8)             ::  a,b

!Allocate Variables
   ALLOCATE ( eta_1d(grd%img*grd%jmg) )
   ALLOCATE ( eta_ad_1d(grd%img*grd%jmg) )
   ALLOCATE ( dummy3(grd%img,grd%jmg,2) )
   ALLOCATE ( bxby_1d(grd%img*grd%jmg*2) )
   ALLOCATE ( bxby_ad_1d(grd%img*grd%jmg*2) )
   ALLOCATE ( bx(grd%img,grd%jmg) )
   ALLOCATE ( by(grd%img,grd%jmg) )
   ALLOCATE ( eta(grd%img,grd%jmg) )
   ALLOCATE ( eta_ad(grd%img,grd%jmg) )

!Store Variables
   eta(   :,:  ) = grd%eta(:,:)
   bx (   :,:  ) = grd%bx(:,:)
   by (   :,:  ) = grd%by(:,:)
   eta_ad(:,:  ) = grd%eta_ad(:,:)

   grd%bx    (:,:) = 0.0_r8
   grd%by    (:,:) = 0.0_r8
   grd%eta   (:,:) = 0.0_r8
   grd%eta_ad(:,:) = 0.0_r8
   dummy3  (:,:,:) = 0.0_r8

   CALL RANDOM_NUMBER(grd%bx(:,:))
   CALL RANDOM_NUMBER(grd%by(:,:))
   grd%bx(:,:)   = grd%bx (:,:) *grd%msk(:,:,1)
   grd%by(:,:)   = grd%by (:,:) *grd%msk(:,:,1)
   dummy3(:,:,1) = grd%bx(:,:)
   dummy3(:,:,2) = grd%by(:,:)
   bxby_1d(:)    = RESHAPE(dummy3,(/grd%img*grd%jmg*2/))

   CALL bar_mod

   eta_1d      (:) = RESHAPE(grd%eta,(/grd%img*grd%jmg/))
   grd%bx  (:,:  ) = 0.0_r8
   grd%by  (:,:  ) = 0.0_r8
   grd%eta   (:,:) = 0.0_r8
   grd%eta_ad(:,:) = 0.0_r8
   dummy3  (:,:,:) = 0.0_r8

   CALL RANDOM_NUMBER(grd%eta_ad(:,:))
   grd%eta_ad(:,:) = grd%eta_ad(:,:) *grd%msk(:,:,1)
   eta_ad_1d(:)    = RESHAPE(grd%eta_ad,(/grd%img*grd%jmg/))

   b = DOT_PRODUCT(eta_1d,eta_ad_1d)

   CALL bar_mod_ad

   dummy3(:,:,1)  = grd%bx(:,:)
   dummy3(:,:,2)  = grd%by(:,:)
   bxby_ad_1d(:) = RESHAPE(dummy3,(/grd%img*grd%jmg*2/))

   a = DOT_PRODUCT(bxby_1d,bxby_ad_1d)

   WRITE (drv%dia,*) '---------------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK BAROTROPIC MODEL              '
   WRITE (drv%dia,*) ' <grd%bxby,grd%bxby_ad>                    = ',a
   WRITE (drv%dia,*) ' <grd%eta,grd%eta_ad>                      = ',b
   WRITE (drv%dia,*) ' Absolute and relative difference          = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '---------------------------------------------'
   CALL FLUSH(drv%dia)

!ReStore Variables
   grd%eta   (:,:) = eta(:,:)
   grd%eta_ad(:,:) = eta_ad(:,:)
   grd%bx    (:,:) = bx(:,:)
   grd%by    (:,:) = by(:,:)

!DeAllocate Variables
   DEALLOCATE ( eta )
   DEALLOCATE ( eta_ad )
   DEALLOCATE ( bx, by )
   DEALLOCATE ( dummy3, bxby_1d, bxby_ad_1d, eta_1d, eta_ad_1d )

END SUBROUTINE adjck_bmd
