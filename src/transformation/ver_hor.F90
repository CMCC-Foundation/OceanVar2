subroutine ver_hor

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                          !
!                                                                          !
!    OceanVar is free software: you can redistribute it and/or modify.     !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    OceanVar is distributed in the hope that it will be useful,           !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.       !
!                                                                          !
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!                                                                      !
! Apply horizontal filter                                              !
!                                                                      !
! Version 1: S.Dobricic               2006                             !
! Version 2: S.Dobricic               2007                             !
! Version 3: S.Dobricic and R. Farina 2013                             !
!     Symmetric calculation in presence of coastal boundaries          !
!     eta_ad, tem_ad, and sal_ad are here temporary arrays             !
!-----------------------------------------------------------------------


 use set_knd
 use grd_str
 use eof_str
 use cns_str
 use drv_str
 use mpi_str

 implicit none

 INTEGER(i4)    :: k, ione, ierr, i, j, iter, niter
 REAL(r8)       :: time1,time2
  

 ione = 1


! ---
! Vertical EOFs           
    call veof

!---
    if(drv%mask(drv%ktr) .gt. 1)then
      niter = 2
    else
      niter = 1
    endif
!---

     do iter=1,niter

! ---
! Load temporary arrays
    if(drv%mask(drv%ktr) .gt. 1) then
       if(drv%bmd(drv%ktr) .ne. 1)      &
          grd%eta_ad(:,:  )    = grd%eta(:,:  )
          grd%tem_ad(:,:,:)    = grd%tem(:,:,:)
          grd%sal_ad(:,:,:)    = grd%sal(:,:,:)
    endif


! ---
! x direction
     if(drv%bmd(drv%ktr) .ne. 1)       &
      call rcfl_x( grd%im, grd%jm, ione, grd%eta(1:grd%im,1:grd%jm), grd%alx(1,1,iter), grd%btx(1,1,iter),      &
                           grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))  
      call rcfl_x( grd%im, grd%jm, grd%km, grd%tem(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter), &
                          grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
      call rcfl_x( grd%im, grd%jm, grd%km, grd%sal(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter), &
                          grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
    
! ---------------------------
! Scale by the scaling factor
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta(1:grd%im,1:grd%jm)   = grd%eta(1:grd%im,1:grd%jm)   * grd%scx(1:grd%im,1:grd%jm,iter)
   do k=1,grd%km
    grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
    grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
   enddo
! ---
! y direction
     if(drv%bmd(drv%ktr) .ne. 1)        &
       call rcfl_y( grd%im, grd%jm, ione, grd%eta(1:grd%im,1:grd%jm), grd%aly(1,1,iter), grd%bty(1,1,iter),  &
                         grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
     call rcfl_y( grd%im, grd%jm, grd%km, grd%tem(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter), &
                           grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
     call rcfl_y( grd%im, grd%jm, grd%km, grd%sal(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter),   &
                           grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))

! ---
! Scale by the scaling factor
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta(1:grd%im,1:grd%jm)   = grd%eta(1:grd%im,1:grd%jm)   * grd%scy(1:grd%im,1:grd%jm,iter)
   do k=1,grd%km
    grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
    grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
   enddo

! ---
! Transpose calculation in the presense of coastal boundaries
 if(drv%mask(drv%ktr) .gt. 1) then

! --------------------------------
! y direction
       if(drv%bmd(drv%ktr) .ne. 1)         &
       call rcfl_y( grd%im, grd%jm,ione, grd%eta_ad(1:grd%im,1:grd%jm), grd%aly(1,1,iter), grd%bty(1,1,iter),    &
                          grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter)) 
       call rcfl_y( grd%im, grd%jm, grd%km, grd%tem_ad(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter), &
                          grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))
       call rcfl_y( grd%im, grd%jm, grd%km, grd%sal_ad(1:grd%im,1:grd%jm,1:grd%km), grd%aly(1,1,iter), grd%bty(1,1,iter),  &
                          grd%gmy(1,1,iter), grd%dly(1,1,iter),grd%mat_bc_y(1,1,1,iter))

! ---
! Scale by the scaling factor
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta_ad(1:grd%im,1:grd%jm)   = grd%eta_ad(1:grd%im,1:grd%jm)   * grd%scy(1:grd%im,1:grd%jm,iter)
   do k=1,grd%km
    grd%tem_ad(1:grd%im,1:grd%jm,k) = grd%tem_ad(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
    grd%sal_ad(1:grd%im,1:grd%jm,k) = grd%sal_ad(1:grd%im,1:grd%jm,k) * grd%scy(1:grd%im,1:grd%jm,iter)
   enddo
! ---
! x direction
  if(drv%bmd(drv%ktr) .ne. 1)       &
    call rcfl_x( grd%im, grd%jm,ione, grd%eta_ad(1:grd%im,1:grd%jm), grd%alx(1,1,iter), grd%btx(1,1,iter),              &
                           grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
    call rcfl_x( grd%im, grd%jm, grd%km, grd%tem_ad(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter),  &
                          grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter))
    call rcfl_x( grd%im, grd%jm, grd%km, grd%sal_ad(1:grd%im,1:grd%jm,1:grd%km), grd%alx(1,1,iter), grd%btx(1,1,iter),  &
                           grd%gmx(1,1,iter), grd%dlx(1,1,iter),grd%mat_bc_x(1,1,1,iter)) 

! ---------------------------
! Scale by the scaling factor
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta_ad(1:grd%im,1:grd%jm)   = grd%eta_ad(1:grd%im,1:grd%jm)   * grd%scx(1:grd%im,1:grd%jm,iter)
   do k=1,grd%km
    grd%tem_ad(1:grd%im,1:grd%jm,k) = grd%tem_ad(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
    grd%sal_ad(1:grd%im,1:grd%jm,k) = grd%sal_ad(1:grd%im,1:grd%jm,k) * grd%scx(1:grd%im,1:grd%jm,iter)
   enddo

! ---
! Average
    if(drv%bmd(drv%ktr) .ne. 1)         &
        grd%eta(1:grd%im,1:grd%jm  )   = (grd%eta(1:grd%im,1:grd%jm  ) + grd%eta_ad(1:grd%im,1:grd%jm  ) ) * 0.5
        grd%tem(1:grd%im,1:grd%jm,:)   = (grd%tem(1:grd%im,1:grd%jm,:) + grd%tem_ad(1:grd%im,1:grd%jm,:) ) * 0.5
        grd%sal(1:grd%im,1:grd%jm,:)   = (grd%sal(1:grd%im,1:grd%jm,:) + grd%sal_ad(1:grd%im,1:grd%jm,:) ) * 0.5

 endif

! ---
! Mask
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta(1:grd%im,1:grd%jm)   = grd%eta(1:grd%im,1:grd%jm)   * grd%msk(1:grd%im,1:grd%jm,1)
    grd%tem(1:grd%im,1:grd%jm,:) = grd%tem(1:grd%im,1:grd%jm,:) * grd%msk(1:grd%im,1:grd%jm,:)
    grd%sal(1:grd%im,1:grd%jm,:) = grd%sal(1:grd%im,1:grd%jm,:) * grd%msk(1:grd%im,1:grd%jm,:)

   enddo 

! ---
! Horizontal localization
    if(drv%bmd(drv%ktr) .ne. 1)         &
       grd%eta(1:grd%im,1:grd%jm)   = grd%eta(1:grd%im,1:grd%jm)   * grd%loc(1:grd%im,1:grd%jm)
   do k=1,grd%km
    grd%tem(1:grd%im,1:grd%jm,k) = grd%tem(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
    grd%sal(1:grd%im,1:grd%jm,k) = grd%sal(1:grd%im,1:grd%jm,k) * grd%loc(1:grd%im,1:grd%jm)
   enddo

1000 continue


! ---
! Bouyancy force
       call get_byg
! ---
! Barotropic model
    if(drv%bmd(drv%ktr) .eq. 1) then
       call bar_mod
    endif

! ---
! Barotropic model
       call get_vel

! ---
! Divergence damping
!    if(drv%dda(drv%ktr) .eq. 1) then
    if( drv%dda(drv%ktr).eq.1 .or. (drv%ddi(drv%ktr).eq.1 .and. drv%ktr.eq.drv%ntr) ) then
        call div_dmp
    endif


end subroutine ver_hor
