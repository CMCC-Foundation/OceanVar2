subroutine rdeofs

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
  use eof_str
  use grd_str
  use mpi_str
  use netcdf

  implicit none

  INTEGER(i4)                    :: stat, ncid, idvar, neofs, nlevs
  INTEGER(i4)                    :: nec, k, nrg
  INTEGER(i4)                :: start(4), count(4)
  real(r4), ALLOCATABLE      :: x4(:,:,:,:), x3(:,:,:), x2(:,:) 
  INTEGER                    :: i, j, ierr, kt, ke
  INTEGER                    :: tmp


    stat = nf90_open('eofs.nc', NF90_NOWRITE, ncid)
    if (stat /= nf90_noerr) call netcdf_err(stat)

! Get dimensions 
      stat = nf90_inq_dimid (ncid, 'nreg', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = ros%nreg)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'nlev', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = nlevs)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_dimid (ncid, 'neof', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inquire_dimension (ncid, idvar, len = neofs)
    if (stat /= nf90_noerr) call netcdf_err(stat)

    write(drv%dia,*)' ---- Eof dimensions are: ',ros%nreg, ros%kmt, neofs
    write(drv%dia,*)' ---- Uses ',ros%neof,' eofs.'

    if(ros%neof .gt. neofs) then
      write(drv%dia,*)'Error: Requires more Eofs than available in the input file.'
      stop
    endif
    if(ros%kmt .ne. nlevs) then
!      write(drv%dia,*)'Error: Vertical dimension different than in the input file.'
!      stop
    endif

! Check the size of nreg. If it is equal to the size of the global grid split the 
! matrix on processors. Otherwise read the whole matrix on each processor.

  if( ros%nreg .eq. grd%img*grd%jmg ) then
      
    ros%nreg = grd%im * grd%jm
    write(drv%dia,*)' ---- NREG equal to ', ros%nreg

     ALLOCATE ( ros%evc( ros%nreg, ros%kmt, neofs), ros%eva( ros%nreg, neofs) )
     ALLOCATE ( ros%cor( ros%nreg, ros%neof, neofs) )

!!     ALLOCATE ( x2( grd%im, neofs), x3( grd%im, ros%kmt, neofs) )
     ALLOCATE ( x3( grd%im, grd%jm, neofs), x4( grd%im, grd%jm, ros%kmt, neofs) )


         k = 0
       do j=1,grd%jm
        do i=1,grd%im
         k = k + 1
         grd%reg(i,j) = k
        enddo
       enddo

      start(1) = grd%igs
      start(2) = grd%jgs
      start(3) = 1
! Oddo
      start(4) = 0

      count(1) = grd%im
      count(2) = grd%jm
      count(3) = neofs
! Oddo
      count(4) = 0

      stat = nf90_inq_varid (ncid, 'eva', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)


      stat = nf90_get_var (ncid,idvar,x3,start(1:3),count(1:3))
    if (stat /= nf90_noerr) call netcdf_err(stat)
      do ke=1,neofs
         k = 0
       do j=1,grd%jm
        do i=1,grd%im
         k = k + 1
         ros%eva(k,ke)  = x3(i,j,ke) 
        enddo
       enddo
      enddo

      start(1) = grd%igs
      start(2) = grd%jgs
      start(3) = 1
      start(4) = 1
      count(1) = grd%im
      count(2) = grd%jm
      count(3) = ros%kmt
      count(4) = neofs
      stat = nf90_inq_varid (ncid, 'evc', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,x4,start,count)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      do ke=1,neofs
      do kt=1,ros%kmt
         k = 0
       do j=1,grd%jm
        do i=1,grd%im
         k = k + 1
         ros%evc(k,kt,ke)  = x4(i,j,kt,ke)
        enddo
       enddo
      enddo
      enddo

  

!!     DEALLOCATE ( x2, x3 )
     DEALLOCATE ( x3, x4 )

  else
!  Allocate eof arrays
     ALLOCATE ( ros%evc( ros%nreg, ros%kmt, neofs), ros%eva( ros%nreg, neofs) )
     ALLOCATE ( ros%cor( ros%nreg, ros%neof, neofs) )

      stat = nf90_inq_varid (ncid, 'eva', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,ros%eva)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_inq_varid (ncid, 'evc', idvar)
    if (stat /= nf90_noerr) call netcdf_err(stat)
      stat = nf90_get_var (ncid,idvar,ros%evc)
    if (stat /= nf90_noerr) call netcdf_err(stat)

  endif

    stat = nf90_close(ncid)

    do nec=1,ros%neof
     do k=1,ros%neof
      do nrg=1,ros%nreg
       if(k.eq.nec)then
        ros%cor(nrg,k,nec) = ros%eva(nrg,nec)
        else
        ros%cor(nrg,k,nec) = 0.0
        endif
      enddo
     enddo
    enddo


!!
!Test for salinity in the case of the SST assimilation
   if( drv%kts.eq.1 .and. drv%nts.eq.2 ) then
    do nec=1,ros%neof
     do k=1,grd%km
      do nrg=1,ros%nreg
!        ros%evc( nrg, k+1+grd%km, nec)  = 0.0
      enddo
     enddo
    enddo
   endif
!!


end subroutine rdeofs


