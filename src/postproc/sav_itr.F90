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
!> Save the result on the coarse grid                                   
!!
!!
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE sav_itr

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str
   USE eof_str
   USE ctl_str
   USE bmd_str
   USE bal_str
   USE mpi_str

   IMPLICIT NONE

! ---
! Save grid dimensions
   drv%im = grd%im
   drv%jm = grd%jm

   ALLOCATE ( drv%ro(drv%im,drv%jm,ros%neof), drv%ro_ad(drv%im,drv%jm,ros%neof) )
   ALLOCATE ( drv%msk(drv%im,drv%jm) )

! ---
! Save eigenvalues
   drv%ro   (:,:,:) = grd%ro   (:,:,:)
   drv%ro_ad(:,:,:) = grd%ro_ad(:,:,:)
   drv%msk  (:,:)   = grd%msr  (:,:,1)

! ---
! Deallocate everithing related to the old grid

! Grid structure
   DEALLOCATE ( grd%reg )
   DEALLOCATE ( grd%msk )
   DEALLOCATE ( grd%hgt )
   DEALLOCATE ( grd%f )
   DEALLOCATE ( grd%tem, grd%sal )
   DEALLOCATE ( grd%uvl, grd%vvl )
   DEALLOCATE ( grd%uvl_ad, grd%vvl_ad )
   DEALLOCATE ( grd%b_x, grd%b_y )
   DEALLOCATE ( grd%dns )
   DEALLOCATE ( grd%bx, grd%by )
   DEALLOCATE ( grd%eta )
   IF ( ANY(drv%nneos .GE. 2) ) THEN
      DEALLOCATE ( grd%temb, grd%salb )
   ENDIF
   DEALLOCATE ( grd%tem_ad, grd%sal_ad )
   DEALLOCATE ( grd%eta_ad )
   DEALLOCATE ( grd%lon, grd%lat, grd%dep )
   DEALLOCATE ( grd%dx, grd%dy, grd%dz )
   DEALLOCATE ( grd%dxdy )
   IF ( drv%filter .EQ. 1 ) THEN
      DEALLOCATE ( grd%alx )
      DEALLOCATE ( grd%aly )
      DEALLOCATE ( grd%btx )
      DEALLOCATE ( grd%bty )
      DEALLOCATE ( grd%gmx )
      DEALLOCATE ( grd%gmy )
      DEALLOCATE ( grd%dlx )
      DEALLOCATE ( grd%dly )
      DEALLOCATE ( grd%scx )
      DEALLOCATE ( grd%scy )
      DEALLOCATE ( grd%msr )
   ENDIF
   DEALLOCATE ( grd%mxd )
   DEALLOCATE ( grd%lcl )

! Observational vector
   DEALLOCATE ( obs%inc, obs%amo, obs%res )
   DEALLOCATE ( obs%err, obs%gra )
! Covariances structure
   DEALLOCATE ( grd%ro)
   DEALLOCATE ( grd%ro_ad)
   DEALLOCATE ( ros%evc, ros%eva )
   DEALLOCATE ( ros%cor )
! Control structure
   DEALLOCATE( ctl%nbd, ctl%iwa )
   DEALLOCATE( ctl%x_c, ctl%g_c )
   DEALLOCATE( ctl%l_c, ctl%u_c )
   IF ( mpi%myrank .EQ. 0 .AND. mpi%flg_min .EQ. 1 ) THEN
      DEALLOCATE( ctl_glob%nbd, ctl_glob%iwa )
      DEALLOCATE( ctl_glob%x_c, ctl_glob%g_c )
      DEALLOCATE( ctl_glob%l_c, ctl_glob%u_c )
      DEALLOCATE( ctl_glob%ws, ctl_glob%wy )
      DEALLOCATE( ctl_glob%z_c, ctl_glob%r_c, ctl_glob%d_c, ctl_glob%t_c )
   ENDIF
   DEALLOCATE( ctl%wa, ctl%sg, ctl%sgo, ctl%yg, ctl%ygo )
   DEALLOCATE( ctl%ws, ctl%wy )
   DEALLOCATE( ctl%sy, ctl%ss, ctl%yy )
   DEALLOCATE( ctl%wt, ctl%wn, ctl%snd )
   DEALLOCATE( ctl%z_c, ctl%r_c, ctl%d_c, ctl%t_c )
! Barotropic model
   IF ( drv%bmd(drv%ktr) .EQ. 1 ) THEN
      DEALLOCATE ( bmd%itr )
      DEALLOCATE ( bmd%mst, bmd%msu, bmd%msv )
      DEALLOCATE ( bmd%hgt, bmd%hgu, bmd%hgv )
      DEALLOCATE ( bmd%dxu, bmd%dxv )
      DEALLOCATE ( bmd%dyu, bmd%dyv )
      DEALLOCATE ( bmd%a1, bmd%a2, bmd%a3 )
      DEALLOCATE ( bmd%a4, bmd%a0 )
      DEALLOCATE ( bmd%bxby, bmd%rgh )
      DEALLOCATE ( bmd%etb, bmd%ub, bmd%vb )
      DEALLOCATE ( bmd%un, bmd%vn )
      DEALLOCATE ( bmd%eta, bmd%ua, bmd%va )
      DEALLOCATE ( bmd%etm, bmd%um, bmd%vm )
      DEALLOCATE ( bmd%div, bmd%cu, bmd%cv )
      DEALLOCATE ( bmd%dux, bmd%duy )
      DEALLOCATE ( bmd%dvx, bmd%dvy )
      DEALLOCATE ( bmd%etx, bmd%ety )
   ENDIF
! Simplified balanced operator
   IF (drv%bal(drv%ktr) .EQ. 1 ) THEN
      DEALLOCATE( bal%dhdz )
   ENDIF
!ADANI New features not implemented yet on multiple grids !!

END SUBROUTINE sav_itr
