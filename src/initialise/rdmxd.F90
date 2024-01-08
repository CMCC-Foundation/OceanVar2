subroutine rdmxd

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
! READ parameters of the MFS_16_72 grid                                !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
! This routine will have effect only if compiled with netcdf library.  !
!-----------------------------------------------------------------------

  use set_knd
  use drv_str
  use grd_str
  use eof_str
  use mpi_str
  use netcdf

  implicit none

  INTEGER(i4)                :: stat, ncid, idvar,k, ierr, i, j
  INTEGER(i4)                :: start(3), count(3)
  real(r4), ALLOCATABLE      :: x3(:,:,:), x2(:,:), x1(:)


 if( drv%kts.eq.1 .and. drv%nts.eq.2 .and. ros%mld.eq.1 ) then

    stat = nf90_open('mldn.nc', NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)

     ALLOCATE ( x2(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae)) 

      start(1) = grd%igs-grd%ias
      start(2) = grd%jgs-grd%jas
      start(3) = 1
      count(1) = grd%im + grd%ias + grd%iae
      count(2) = grd%jm + grd%jas + grd%jae
      count(3) = 1 !grd%km

      stat = nf90_inq_varid (ncid, 'mldd', idvar)
       if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x2,start,count)
      if (stat /= nf90_noerr) call netcdf_err(stat)
      grd%mxd(:,:) =  x2(:,:)

    DEALLOCATE ( x2 )

    stat = nf90_close(ncid)



   do j = 1-grd%jas,grd%jm+grd%jae
    do i = 1-grd%ias,grd%im+grd%iae

       grd%mxd(i,j) = max(2,min(grd%mxd(i,j),grd%km-1))

      do k=1,grd%mxd(i,j)-1
       grd%lcl(i,j,k) = 1.0
      enddo
         k=grd%mxd(i,j)
       grd%lcl(i,j,k) = 0.5
      do k=grd%mxd(i,j)+1,grd%km
       grd%lcl(i,j,k) = 0.0
      enddo

    enddo
   enddo

 else

      grd%lcl(:,:,:) = 1.0
      write(*,*)'Mixlayer loc=1'

 endif

end subroutine rdmxd


