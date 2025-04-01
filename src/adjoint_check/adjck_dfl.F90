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
!> Diffusive Filter adjoint check                                     
!!
!! It checks if the linear tangent and the adjoint are consistent.
!! It can only be activated without domain decomposition.
!! It is used for debugging purposes only.
!!
!                                                                      !
! Version 1: M.Adani    2023                                           !
!-----------------------------------------------------------------------
SUBROUTINE adjck_dfl

   USE set_knd
   USE grd_str
   USE dfl_str
   USE drv_str

   IMPLICIT NONE

   REAL(r8),POINTER :: tem(:,:,:)
   REAL(r8),POINTER :: sal(:,:,:)
   REAL(r8),POINTER :: eta(:,:)
   REAL(r8),POINTER :: tem_ad(:,:,:)
   REAL(r8),POINTER :: sal_ad(:,:,:)
   REAL(r8),POINTER :: eta_ad(:,:)
   REAL(r8),POINTER :: tem1d_in(:)
   REAL(r8),POINTER :: sal1d_in(:)
   REAL(r8),POINTER :: eta1d_in(:)
   REAL(r8),POINTER :: tem1d_ou(:)
   REAL(r8),POINTER :: sal1d_ou(:)
   REAL(r8),POINTER :: eta1d_ou(:)
   REAL(r8),POINTER :: tem1d_ad_in(:)
   REAL(r8),POINTER :: sal1d_ad_in(:)
   REAL(r8),POINTER :: eta1d_ad_in(:)
   REAL(r8),POINTER :: tem1d_ad_ou(:)
   REAL(r8),POINTER :: sal1d_ad_ou(:)
   REAL(r8),POINTER :: eta1d_ad_ou(:)
   REAL(r8)         :: a,b

   ALLOCATE ( tem(grd%img,grd%jmg,grd%km),         &
              sal(grd%img,grd%jmg,grd%km),         &
              eta(grd%img,grd%jmg)         )
   ALLOCATE ( tem_ad(grd%img,grd%jmg,grd%km),      &
              sal_ad(grd%img,grd%jmg,grd%km),      &
              eta_ad(grd%img,grd%jmg)      )
   ALLOCATE ( tem1d_in(grd%img*grd%jmg*grd%km),    &
              sal1d_in(grd%img*grd%jmg*grd%km),    &
              eta1d_in(grd%img*grd%jmg)    )
   ALLOCATE ( tem1d_ad_in(grd%img*grd%jmg*grd%km), &
              sal1d_ad_in(grd%img*grd%jmg*grd%km), &
              eta1d_ad_in(grd%img*grd%jmg) )
   ALLOCATE ( tem1d_ou(grd%img*grd%jmg*grd%km),    &
              sal1d_ou(grd%img*grd%jmg*grd%km),    &
              eta1d_ou(grd%img*grd%jmg)    )
   ALLOCATE (tem1d_ad_ou(grd%img*grd%jmg*grd%km),  &
             sal1d_ad_ou(grd%img*grd%jmg*grd%km),  &
             eta1d_ad_ou(grd%img*grd%jmg)  )

!Store Variables
   tem(:,:,:) = grd%tem(:,:,:)
   sal(:,:,:) = grd%sal(:,:,:)
   eta(:,:) = grd%eta(:,:)
   tem_ad(:,:,:) = grd%tem_ad(:,:,:)
   sal_ad(:,:,:) = grd%sal_ad(:,:,:)
   eta_ad(:,:) = grd%eta_ad(:,:)

! Random perturbation
   CALL RANDOM_NUMBER(grd%tem(:,:,:))
   CALL RANDOM_NUMBER(grd%sal(:,:,:))
   CALL RANDOM_NUMBER(grd%eta(:,:))

! Mask
   grd%tem(:,:,:) = grd%tem(:,:,:) * grd%msk(:,:,:)
   grd%sal(:,:,:) = grd%sal(:,:,:) * grd%msk(:,:,:)
   grd%eta(:,:)   = grd%eta(:,:)   * grd%msk(:,:,1)

! Save input
   tem1d_in =  RESHAPE(grd%tem,(/grd%img*grd%jmg*grd%km/))
   sal1d_in =  RESHAPE(grd%sal,(/grd%img*grd%jmg*grd%km/))
   eta1d_in =  RESHAPE(grd%eta,(/grd%img*grd%jmg/))

! start dIFfusion
   CALL diffusive_filter

! save output farward
   tem1d_ou =  RESHAPE(grd%tem,(/grd%img*grd%jmg*grd%km/))
   sal1d_ou =  RESHAPE(grd%sal,(/grd%img*grd%jmg*grd%km/))
   eta1d_ou =  RESHAPE(grd%eta,(/grd%img*grd%jmg/))

!--- adjoint part ---

! Random perturbation
   CALL RANDOM_NUMBER(grd%tem_ad(:,:,:))
   CALL RANDOM_NUMBER(grd%sal_ad(:,:,:))
   CALL RANDOM_NUMBER(grd%eta_ad(:,:))

! Mask
   grd%tem_ad(:,:,:) = grd%tem_ad(:,:,:) * grd%msk(:,:,:)
   grd%sal_ad(:,:,:) = grd%sal_ad(:,:,:) * grd%msk(:,:,:)
   grd%eta_ad(:,:)   = grd%eta_ad(:,:)   * grd%msk(:,:,1)

! Save input
   tem1d_ad_in =  RESHAPE(grd%tem_ad,(/grd%img*grd%jmg*grd%km/))
   sal1d_ad_in =  RESHAPE(grd%sal_ad,(/grd%img*grd%jmg*grd%km/))
   eta1d_ad_in =  RESHAPE(grd%eta_ad,(/grd%img*grd%jmg/))

! ... adjoint
   CALL diffusive_filter_ad

! save output backward
   tem1d_ad_ou =  RESHAPE(grd%tem_ad,(/grd%img*grd%jmg*grd%km/))
   sal1d_ad_ou =  RESHAPE(grd%sal_ad,(/grd%img*grd%jmg*grd%km/))
   eta1d_ad_ou =  RESHAPE(grd%eta_ad,(/grd%img*grd%jmg/))

   a = DOT_PRODUCT(tem1d_in,tem1d_ad_ou)
   b = DOT_PRODUCT(tem1d_ad_in,tem1d_ou)

   WRITE (drv%dia,*) '--------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK DIFFUSIVE FILTER:      '
   WRITE (drv%dia,*) ' TEMPERATURE:                         '
   WRITE (drv%dia,*) ' <fld1d_in,fld1d_ad_out>            = ',a
   WRITE (drv%dia,*) ' <fld1d_ad_in,fld1d_out>            = ',b
   WRITE (drv%dia,*) ' Absolute and relative difference   = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '--------------------------------------'

   a = DOT_PRODUCT(sal1d_in,sal1d_ad_ou)
   b = DOT_PRODUCT(sal1d_ad_in,sal1d_ou)

   WRITE (drv%dia,*) '--------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK DIFFUSIVE FILTER:      '
   WRITE (drv%dia,*) ' SALINITY:                            '
   WRITE (drv%dia,*) ' <fld1d_in,fld1d_ad_out>            = ',a
   WRITE (drv%dia,*) ' <fld1d_ad_in,fld1d_out>            = ',b
   WRITE (drv%dia,*) ' Absolute and relative difference   = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '--------------------------------------'

   a = DOT_PRODUCT(eta1d_in,eta1d_ad_ou)
   b = DOT_PRODUCT(eta1d_ad_in,eta1d_ou)

   WRITE (drv%dia,*) '--------------------------------------'
   WRITE (drv%dia,*) ' ADJOINT CHECK DIFFUSIVE FILTER:      '
   WRITE (drv%dia,*) ' SLA:                                 '
   WRITE (drv%dia,*) ' <fld1d_in,fld1d_ad_out>            = ',a
   WRITE (drv%dia,*) ' <fld1d_ad_in,fld1d_out>            = ',b
   WRITE (drv%dia,*) ' Absolute and relative difference   = ',a-b,' / ',(a-b)/b
   WRITE (drv%dia,*) '--------------------------------------'

! ReStore Variables
   grd%tem(   :,:,:) = tem(:,:,:)
   grd%sal(   :,:,:) = sal(:,:,:)
   grd%eta(   :,:  ) = eta(:,:)
   grd%tem_ad(:,:,:) = tem_ad(:,:,:)
   grd%sal_ad(:,:,:) = sal_ad(:,:,:)
   grd%eta_ad(:,:  ) = eta_ad(:,:)

!DeAllocate Variables
   DEALLOCATE ( tem, sal, eta  )
   DEALLOCATE ( tem_ad, sal_ad, eta_ad  )
   DEALLOCATE ( tem1d_in, sal1d_in, eta1d_in )
   DEALLOCATE ( tem1d_ou, sal1d_ou, eta1d_ou )
   DEALLOCATE ( tem1d_ad_in, sal1d_ad_in, eta1d_ad_in )
   DEALLOCATE ( tem1d_ad_ou, sal1d_ad_ou, eta1d_ad_ou )

END SUBROUTINE adjck_dfl
