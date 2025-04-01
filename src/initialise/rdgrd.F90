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
!> Read grid                                                           
!!
!! It define mpi grid variables and 
!! readsa all the variables realated to the grid.
!!
!                                                                      !
! Version 1: Srjan Dobricic 2006                                       !
! Version 2: Mario Adani    2024                                       !
!-----------------------------------------------------------------------
SUBROUTINE rdgrd

   USE set_knd
   USE drv_str
   USE grd_str
   USE bal_str
   USE bmd_str
   USE obs_str, ONLY : obserr,coastrej
   USE eof_str
   USE netcdf

   IMPLICIT NONE

   INTEGER(i4)                :: stat, ncid, idvar,k, ierr
   INTEGER(i4)                :: start(3), count(3)
   REAL(r4), ALLOCATABLE      :: x3(:,:,:), x2(:,:), x1(:)
   INTEGER(i4)                :: km
   CHARACTER*1                :: cgrd

   WRITE (drv%dia,*)' -----In rdgrd grd%grd_mod= ',grd%grd_mod
   WRITE (cgrd,'(i1)') grd%grd_mod

! Vertical levels ------------------------------------------------
   stat = NF90_OPEN(drv%inpdir//'/'//'grid'//cgrd//'.nc', NF90_NOWRITE, ncid)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

! Get dimensions
   stat = NF90_INQ_DIMID (ncid, 'im', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = grd%img)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'jm', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = grd%jmg)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQ_DIMID (ncid, 'km', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = grd%km)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)

   stat = NF90_INQ_DIMID (ncid, 'iex', idvar)
   IF (stat /= NF90_NOERR) THEN
      grd%iex = 0
   ELSE
      stat = NF90_INQUIRE_DIMENSION (ncid, idvar, len = grd%iex)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   ENDIF

   WRITE (drv%dia,*)' ---- Grid dims are:'
   WRITE (drv%dia,*)' ---- jpidta = ',grd%img,'jpjdta=',grd%jmg,'jpkdta=',grd%km
   WRITE (drv%dia,*)' '
   WRITE (drv%dia,*)' May want to check with NEMO for consistency'
   WRITE (drv%dia,*)'-----------------------------------------------------------------'

   CALL def_mpid

!  Allocate grid arrays
   ALLOCATE ( grd%hgt(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%ums(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%vms(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%f(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%tem(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%sal(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%uvl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%vvl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%eta(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%tem_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%sal_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%uvl_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%vvl_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%eta_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%dx  (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%dy  (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%dxdy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%lon (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%lat (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%b_x(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) ); grd%b_x=0.0_r8
   ALLOCATE ( grd%b_y(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) ); grd%b_y=0.0_r8
   ALLOCATE ( grd%dns(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( grd%bx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%by(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   IF ( coastrej%any ) THEN
      ALLOCATE ( coastrej%distc(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ENDIF
   IF ( drv%filter .EQ. 2 ) THEN
      ALLOCATE ( grd%distc3d(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ENDIF
   ALLOCATE ( grd%bx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%by(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%temb(grd%im,grd%jm,grd%km), grd%salb(grd%im,grd%jm,grd%km) )
   ALLOCATE ( grd%dep(grd%km) )
   ALLOCATE ( grd%dz(grd%km) )
   IF ( drv%filter .EQ. 1 ) THEN
      ALLOCATE ( grd%alx(grd%im,grd%jm,2) )
      ALLOCATE ( grd%aly(grd%jm,grd%im,2) )
      ALLOCATE ( grd%btx(grd%im,grd%jm,2) )
      ALLOCATE ( grd%bty(grd%jm,grd%im,2) )
      ALLOCATE ( grd%gmx(grd%im,grd%jm,2) )
      ALLOCATE ( grd%gmy(grd%jm,grd%im,2) )
      ALLOCATE ( grd%dlx(grd%im,grd%jm,2) ); grd%dlx=0.0_r8
      ALLOCATE ( grd%dly(grd%jm,grd%im,2) ); grd%dly=0.0_r8
      ALLOCATE ( grd%mat_bc_x(9,grd%im,grd%jm,2) )
      ALLOCATE ( grd%mat_bc_y(9,grd%jm,grd%im,2) )
      ALLOCATE ( grd%tmx(grd%jm*grd%km) )
      ALLOCATE ( grd%tmy(grd%im*grd%km) )
      ALLOCATE ( grd%scx(grd%im,grd%jm,2) )
      ALLOCATE ( grd%scy(grd%im,grd%jm,2) )
      ALLOCATE ( grd%msr(1-3:grd%im+3,1-3:grd%jm+3,grd%km) )
   ENDIF
   ALLOCATE ( grd%mxd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( grd%lcl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) ) 
   ALLOCATE ( grd%loc(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( x3(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km) )
   ALLOCATE ( x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
   ALLOCATE ( x1(grd%km) )

   stat = NF90_INQ_VARID (ncid, 'dz', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x1)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   grd%dz(:) = DBLE(x1(:))
   stat = NF90_INQ_VARID (ncid, 'dep', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid, idvar, x1)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   grd%dep(:) = DBLE(x1(:))

   start(1) = grd%igs-grd%ias
   start(2) = grd%jgs-grd%jas
   start(3) = 1
   count(1) = grd%im + grd%ias + grd%iae
   count(2) = grd%jm + grd%jas + grd%jae
   count(3) = grd%km

   stat = NF90_INQ_VARID (ncid, 'lon', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x2,start(1:2),count(1:2))
   grd%lon(:,:) = DBLE(x2(:,:))

   stat = NF90_INQ_VARID (ncid, 'lat', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x2,start(1:2),count(1:2))
   grd%lat(:,:) = DBLE(x2(:,:))

   stat = NF90_INQ_VARID (ncid, 'dx', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x2,start(1:2),count(1:2))
   grd%dx(:,:) = DBLE(x2(:,:))

   stat = NF90_INQ_VARID (ncid, 'dy', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x2,start(1:2),count(1:2))
   grd%dy(:,:) = DBLE(x2(:,:))

   stat = NF90_INQ_VARID (ncid, 'topo', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x2,start(1:2),count(1:2))
   grd%hgt(:,:) = DBLE(x2(:,:))

   stat = NF90_INQ_VARID (ncid, 'tmsk', idvar)
   IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
   stat = NF90_GET_VAR (ncid,idvar,x3,start,count)
   grd%msk(:,:,:) = DBLE(x3(:,:,:))

   ros%already_read = .False.
   stat = NF90_INQ_VARID (ncid, 'regs', idvar)
   IF (stat == NF90_NOERR) THEN
      ALLOCATE ( grd%reg(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
      stat = NF90_GET_VAR (ncid, idvar, x2,start(1:2),count(1:2))
      grd%reg(:,:) = INT(x2(:,:))
      ros%already_read = .True.
      WHERE( grd%reg(:,:) .EQ. 0 ) grd%reg(:,:) = 1
   ENDIF

   IF ( coastrej%any )  THEN
      stat = NF90_INQ_VARID (ncid, 'dcoast2d', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, x2,start(1:2),count(1:2))
      coastrej%distc(:,:) = DBLE(x2(:,:))
   ENDIF

   IF ( drv%filter .EQ. 2 )  THEN
      stat = NF90_INQ_VARID (ncid, 'dcoast', idvar)
      IF (stat /= NF90_NOERR) CALL netcdf_err(stat)
      stat = NF90_GET_VAR (ncid, idvar, x3,start(1:3),count(1:3))
      grd%distc3d(:,:,:) = DBLE(x3(:,:,:))
   ENDIF

   DEALLOCATE ( x3, x2, x1  )

   stat = NF90_CLOSE(ncid)

! ----------------------------------------------------------------

END SUBROUTINE rdgrd
