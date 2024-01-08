subroutine rdgrd

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2006 Srdjan Dobricic, CMCC, Bologna                         !
!                                                                          !
!    This file is part of OceanVar.                                        !
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
!    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.     !
!                                                                          !
!--------------------------------------------------------------------------- 

!-----------------------------------------------------------------------
!                                                                      !
! READ parameters of the MFS_16_72 grid                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use set_knd
  use drv_str
  use grd_str
  use netcdf

  implicit none

  INTEGER(i4)                :: stat, ncid, idvar,k, ierr
  INTEGER(i4)                :: start(3), count(3)
  real(r4), ALLOCATABLE      :: x3(:,:,:), x2(:,:), x1(:)

  character*1                    :: cgrd

    write(drv%dia,*)' -----In rdgrd grd%grd_mod= ',grd%grd_mod
    write(cgrd,'(i1)') grd%grd_mod

! Vertical levels ------------------------------------------------
    stat = nf90_open('grid'//cgrd//'.nc', NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)


! Get dimensions 
      stat = nf90_inq_dimid (ncid, 'im', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%img)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'jm', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%jmg)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'km', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%km)
    if (stat /= nf90_noerr) call netcdf_err(stat)

      stat = nf90_inq_dimid (ncid, 'iex', idvar)
    if (stat /= nf90_noerr) then
      grd%iex = 0
    else
      stat = nf90_inquire_dimension (ncid, idvar, len = grd%iex)
    if (stat /= nf90_noerr) call netcdf_err(stat)
    endif

    write(drv%dia,*)' ---- Grid dims are:'
    write(drv%dia,*)' ---- jpidta = ',grd%img,'jpjdta=',grd%jmg,'jpkdta=',grd%km
    write(drv%dia,*)' '
    write(drv%dia,*)' May want to check with NEMO for consistency'
    write(drv%dia,*)'-----------------------------------------------------------------'

      call def_mpid

!  Allocate grid arrays
     ALLOCATE ( grd%reg(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%hgt(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%ums(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%vms(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%f(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     ALLOCATE ( grd%tem(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%sal(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%uvl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%vvl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%eta(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     ALLOCATE ( grd%tem_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%sal_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%uvl_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%vvl_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%eta_ad(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     ALLOCATE ( grd%dx  (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%dy  (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%dxdy(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     ALLOCATE ( grd%lon (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%lat (1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))

     ALLOCATE ( grd%b_x(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%b_y(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%dns(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%bx(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( grd%by(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae) )
     ALLOCATE ( grd%temb(grd%im,grd%jm,grd%km), grd%salb(grd%im,grd%jm,grd%km))
     ALLOCATE ( grd%etab(grd%im,grd%jm))
     ALLOCATE ( grd%dep(grd%km))
     ALLOCATE ( grd%dz(grd%km))
     ALLOCATE ( grd%alx(grd%im,grd%jm,2) )
     ALLOCATE ( grd%aly(grd%jm,grd%im,2) )
     ALLOCATE ( grd%btx(grd%im,grd%jm,2) )
     ALLOCATE ( grd%bty(grd%jm,grd%im,2) )
     ALLOCATE ( grd%gmx(grd%im,grd%jm,2) )
     ALLOCATE ( grd%gmy(grd%jm,grd%im,2) )
     ALLOCATE ( grd%dlx(grd%im,grd%jm,2) )
     ALLOCATE ( grd%dly(grd%jm,grd%im,2) )
     ALLOCATE ( grd%mat_bc_x(9,grd%im,grd%jm,2) )
     ALLOCATE ( grd%mat_bc_y(9,grd%jm,grd%im,2) )
     ALLOCATE ( grd%tmx(grd%jm*grd%km) )
     ALLOCATE ( grd%tmy(grd%im*grd%km) )


     ALLOCATE ( grd%scx(grd%im,grd%jm,2) )
     ALLOCATE ( grd%scy(grd%im,grd%jm,2) )
     ALLOCATE ( grd%msr(1-3:grd%im+3,1-3:grd%jm+3,grd%km) )
     ALLOCATE ( grd%imx(grd%jm*grd%km,2), grd%jmx(grd%im*grd%km,2))
     ALLOCATE ( grd%istp(grd%im,grd%jm), grd%jstp(grd%im,grd%jm))
     ALLOCATE ( grd%inx(grd%im,grd%jm,grd%km,2), grd%jnx(grd%im,grd%jm,grd%km,2))
     ALLOCATE ( grd%fct(grd%im,grd%jm,grd%km) )

     ALLOCATE ( grd%mxd(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))
     ALLOCATE ( grd%lcl(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km))
     ALLOCATE ( grd%loc(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae))


     ALLOCATE ( x3(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,grd%km)) 
     ALLOCATE ( x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae)) 
     ALLOCATE ( x1(grd%km) ) 

      stat = nf90_inq_varid (ncid, 'dz', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x1)
    if (stat /= nf90_noerr) call netcdf_err(stat)
!ADANI    grd%dz(:) = x1(:)
    grd%dz(:) = dble(x1(:))
     stat = nf90_inq_varid (ncid, 'dep', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
     stat = nf90_get_var (ncid, idvar, x1)
    if (stat /= nf90_noerr) call netcdf_err(stat)
!ADANI    grd%dep(:) = x1(:)
    grd%dep(:) = dble(x1(:))

      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      start(3) = 1
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae
      count(3) = grd%km

      stat = nf90_inq_varid (ncid, 'lon', idvar)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
!ADANI    grd%lon(:,:) = x2(:,:)
    grd%lon(:,:) = dble(x2(:,:))

      stat = nf90_inq_varid (ncid, 'lat', idvar)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
!ADANI    grd%lat(:,:) = x2(:,:)
    grd%lat(:,:) = dble(x2(:,:))

      stat = nf90_inq_varid (ncid, 'dx', idvar)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
!ADANI    grd%dx(:,:) = x2(:,:)
    grd%dx(:,:) = dble(x2(:,:))

      stat = nf90_inq_varid (ncid, 'dy', idvar)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
!ADANI    grd%dy(:,:) = x2(:,:)
    grd%dy(:,:) = dble(x2(:,:))

      stat = nf90_inq_varid (ncid, 'topo', idvar)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
!ADANI    grd%hgt(:,:) = x2(:,:)
    grd%hgt(:,:) = dble(x2(:,:))

      stat = nf90_inq_varid (ncid, 'tmsk', idvar)
      stat = nf90_get_var (ncid,idvar,x3,start,count)
!ADANI    grd%msk(:,:,:) = x3(:,:,:)
    grd%msk(:,:,:) = dble(x3(:,:,:))

      stat = nf90_inq_varid (ncid, 'regs', idvar)
      stat = nf90_get_var (ncid, idvar, x2,start,count)
    grd%reg(:,:) = int(x2(:,:))

    DEALLOCATE ( x3, x2, x1 )

    stat = nf90_close(ncid)

! ----------------------------------------------------------------

end subroutine rdgrd


