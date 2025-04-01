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
!> Apply horizontal filter                                              
!                                                                      !
! Version 1: S.Dobricic               2006                             !
! Version 2: S.Dobricic               2007                             !
! Version 3: S.Dobricic and R. Farina 2013                             !
!     Symmetric calculation in presence of coastal boundaries          !
!     eta_ad, tem_ad, and sal_ad are here temporary arrays             !
!-----------------------------------------------------------------------
SUBROUTINE recursive_filter

   USE set_knd
   USE grd_str
   USE cns_str
   USE drv_str

   IMPLICIT NONE

   INTEGER(i4)    :: k, ione,   iter, niter
   LOGICAL        :: eta_condition

   ione = 1

   eta_condition = ( (drv%bmd(drv%ktr)+drv%bal(drv%ktr) .EQ. 0)                     .OR.  &
                    ((drv%bal(drv%ktr) .EQ. 1) .AND. (drv%ssh_unbalanced(drv%ktr))) )

!---
   IF ( drv%mask(drv%ktr) .GT. 1 ) THEN
      niter = 2
   ELSE
      niter = 1
   ENDIF
!---

   DO iter=1,niter
! ---
! Load temporary arrays
      IF ( drv%mask(drv%ktr) .GT. 1 ) THEN
         IF ( eta_condition ) THEN
            grd%eta_ad(:,:)   = grd%eta(:,:  )
         ENDIF
         grd%tem_ad(:,:,:)    = grd%tem(:,:,:)
         grd%sal_ad(:,:,:)    = grd%sal(:,:,:)
      ENDIF

! ---
! x direction
      IF ( eta_condition ) THEN
         CALL rcfl_x( grd%im, grd%jm, ione, grd%eta(1:grd%im,1:grd%jm), grd%alx(1,1,iter), grd%btx(1,1,iter),      &
                      grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
      ENDIF
      CALL rcfl_x( grd%im, grd%jm, grd%km, grd%tem(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter), &
                   grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
      CALL rcfl_x( grd%im, grd%jm, grd%km, grd%sal(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter), &
                   grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
! ---------------------------
! Scale by the scaling factor
      IF ( eta_condition ) THEN
         grd%eta(1:grd%im,1:grd%jm) = grd%eta(1:grd%im,1:grd%jm)   * grd%scx(1:grd%im,1:grd%jm,iter)
      ENDIF
      DO k = 1,grd%km
         grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
         grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
      ENDDO
! ---
! y direction
      IF ( eta_condition ) THEN
         CALL rcfl_y( grd%im, grd%jm, ione, grd%eta(1:grd%im,1:grd%jm), grd%aly(1,1,iter), grd%bty(1,1,iter),  &
                      grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
      ENDIF
      CALL rcfl_y( grd%im, grd%jm, grd%km, grd%tem(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter), &
                   grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
      CALL rcfl_y( grd%im, grd%jm, grd%km, grd%sal(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter),   &
                   grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
! ---
! Scale by the scaling factor
      IF ( eta_condition ) THEN
         grd%eta(1:grd%im,1:grd%jm) = grd%eta(1:grd%im,1:grd%jm)   * grd%scy(1:grd%im,1:grd%jm,iter)
      ENDIF
      DO k=1,grd%km
         grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
         grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
      ENDDO

! ---
! Transpose calculation in the presense of coastal boundaries
      IF ( drv%mask(drv%ktr) .GT. 1 ) THEN

! --------------------------------
! y direction
         IF ( eta_condition ) THEN
            CALL rcfl_y( grd%im, grd%jm,ione, grd%eta_ad(1:grd%im,1:grd%jm), grd%aly(1,1,iter), grd%bty(1,1,iter),    &
                         grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
         ENDIF
         CALL rcfl_y( grd%im, grd%jm, grd%km, grd%tem_ad(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter), &
                      grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
         CALL rcfl_y( grd%im, grd%jm, grd%km, grd%sal_ad(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter),  &
                      grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
! ---
! Scale by the scaling factor
         IF ( eta_condition ) THEN
            grd%eta_ad(1:grd%im,1:grd%jm) = grd%eta_ad(1:grd%im,1:grd%jm)   * grd%scy(1:grd%im,1:grd%jm,iter)
         ENDIF
         DO k=1,grd%km
            grd%tem_ad(1:grd%im,1:grd%jm,k) = grd%tem_ad(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
            grd%sal_ad(1:grd%im,1:grd%jm,k) = grd%sal_ad(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
         ENDDO
! ---
! x direction
         IF ( eta_condition ) THEN
            CALL rcfl_x( grd%im, grd%jm,ione, grd%eta_ad(1:grd%im,1:grd%jm), grd%alx(1,1,iter), grd%btx(1,1,iter),              &
                         grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
         ENDIF
         CALL rcfl_x( grd%im, grd%jm, grd%km, grd%tem_ad(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter),  &
                      grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
         CALL rcfl_x( grd%im, grd%jm, grd%km, grd%sal_ad(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter),  &
                      grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
! ---------------------------
! Scale by the scaling factor
         IF (eta_condition) THEN
            grd%eta_ad(1:grd%im,1:grd%jm) = grd%eta_ad(1:grd%im,1:grd%jm)   * grd%scx(1:grd%im,1:grd%jm,iter)
         ENDIF
         DO k=1,grd%km
            grd%tem_ad(1:grd%im,1:grd%jm,k) = grd%tem_ad(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
            grd%sal_ad(1:grd%im,1:grd%jm,k) = grd%sal_ad(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
         ENDDO
! ---
! Average
         IF (eta_condition) THEN
            grd%eta(1:grd%im,1:grd%jm)  = (grd%eta(1:grd%im,1:grd%jm  ) + grd%eta_ad(1:grd%im,1:grd%jm  ) ) * 0.5_r8
         ENDIF
         grd%tem(1:grd%im,1:grd%jm,:)   = (grd%tem(1:grd%im,1:grd%jm,:) + grd%tem_ad(1:grd%im,1:grd%jm,:) ) * 0.5_r8
         grd%sal(1:grd%im,1:grd%jm,:)   = (grd%sal(1:grd%im,1:grd%jm,:) + grd%sal_ad(1:grd%im,1:grd%jm,:) ) * 0.5_r8

      ENDIF 
! ---
! Mask
      IF ( eta_condition ) THEN
         grd%eta(1:grd%im,1:grd%jm)= grd%eta(1:grd%im,1:grd%jm)   * grd%msk(1:grd%im,1:grd%jm,1)
      ENDIF
      grd%tem(1:grd%im,1:grd%jm,:) = grd%tem(1:grd%im,1:grd%jm,:) * grd%msk(1:grd%im,1:grd%jm,:)
      grd%sal(1:grd%im,1:grd%jm,:) = grd%sal(1:grd%im,1:grd%jm,:) * grd%msk(1:grd%im,1:grd%jm,:)

   ENDDO  ! iter

END SUBROUTINE recursive_filter
