subroutine def_grd

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
! Define the grid                                                      !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use mpi_str

  implicit none

  include 'mpif.h'

  INTEGER(I4)    :: i, j, k , ierr, ich
  REAL(r8)       :: dxdyl, dxdyp
!ADANI FROM OCEANVAR
  REAL(r8) :: dxmin,dxmax,dymin,dymax
!ADANI FROM OCEANVAR


! ---
! Define grid 

 grd%grd_mod  = drv%grid (drv%ktr)
 grd%read_grd = drv%read_grd (drv%ktr)

!Read grid definition
      write(drv%dia,*)' ----'
      write(drv%dia,*)' ---- Reading GRIDS : '
      call rdgrd

!ADANI FROM OCEANVAR
    dxmin = minval(grd%dx,mask=(int(grd%msk(:,:,1)).eq.1))
    dxmax = maxval(grd%dx,mask=(int(grd%msk(:,:,1)).eq.1))
    dymin = minval(grd%dy,mask=(int(grd%msk(:,:,1)).eq.1))
    dymax = maxval(grd%dy,mask=(int(grd%msk(:,:,1)).eq.1))
    where(grd%dx.lt.dxmin) grd%dx=dxmin
    where(grd%dy.lt.dymin) grd%dy=dymin
    where(grd%dx.gt.dxmax) grd%dx=dxmax
    where(grd%dy.gt.dymax) grd%dy=dymax
!ADANI FROM OCEANVAR

     grd%dxdy(:,:) =  grd%dy(:,:) * grd%dx(:,:)

     grd%dlt = (grd%lat(1,2) - grd%lat(1,1))
     grd%dln = (grd%lon(2,1) - grd%lon(1,1))

     if(mpi%ir.eq.1      ) grd%msk(1,:,:)      = 0.0
     if(mpi%jr.eq.1      ) grd%msk(:,1,:)      = 0.0
     if(mpi%ir.eq.mpi%irm) grd%msk(grd%im,:,:) = 0.0
     if(mpi%jr.eq.mpi%jrm) grd%msk(:,grd%jm,:) = 0.0


! ---
      grd%npsa = grd%img*grd%jmg
      do j=1,grd%jm
         do i=1,grd%im
!ADANI       dxdyl = sum(grd%dxdy(1:grd%im,1:grd%jm))
             dxdyl = dxdyl + grd%dxdy(i,j)
          enddo
      enddo

 if(mpi%nproc.gt.1)then
  call mpi_reduce( dxdyl, dxdyp, 1, mpi%r8, mpi_sum, 0, mpi%comm, ierr)
  if(mpi%myrank.eq.0) dxdyl = dxdyp
  call mpi_bcast( dxdyl, 1, mpi%r8, 0, mpi%comm, ierr)
 endif
      grd%adxdy = dxdyl/dble(grd%npsa)
      grd%adxdy = dsqrt(grd%adxdy)
      grd%dxdy(:,:) = dsqrt(grd%dxdy(:,:))
! ---

      grd%nps = grd%im*grd%jm

     do k=1,grd%km
        grd%ums(1:grd%im+grd%iae-1,1:grd%jm+grd%jae,k) =                      &
          grd%msk(1:grd%im+grd%iae-1,1:grd%jm+grd%jae,k) * grd%msk(2:grd%im+grd%iae,1:grd%jm+grd%jae,k)
        grd%vms(1:grd%im+grd%iae,1:grd%jm+grd%jae-1,k) =                      &
          grd%msk(1:grd%im+grd%iae,1:grd%jm+grd%jae-1,k) * grd%msk(1:grd%im+grd%iae,2:grd%jm+grd%jae,k)
     enddo

 ! Define grid for horizontal covariances
      if( drv%mask(drv%ktr).eq.1)then
          grd%msr(:,:,:) = 1.0
      else if( drv%mask(drv%ktr).eq.2)then
         do k=1,grd%km
          grd%msr(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,k) =    &
              grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,1)
        enddo
      else if( drv%mask(drv%ktr).eq.3)then
          grd%msr(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,:) =    &
              grd%msk(1-grd%ias:grd%im+grd%iae,1-grd%jas:grd%jm+grd%jae,:)
      else
          write(drv%dia,*)'Wrong mask for horizontal covariances ',  &
                          drv%mask(drv%ktr)
         stop
      endif
 
      if( drv%mask(drv%ktr).eq.2 .or. drv%mask(drv%ktr).eq.3 )then
       call exa_mpi ( 1_i4, 2_i4, grd%im-1, 1-2, grd%im+2,   &
                            2_i4, grd%jm-1, 1-2, grd%jm+2,   &
                            1-3,grd%im+3,1-3,grd%jm+3, grd%km, grd%msr)
       call exa_mpi ( 1_i4, 3_i4, grd%im-2, 1-3, grd%im+3,   &
                            3_i4, grd%jm-2, 1-3, grd%jm+3,   &
                            1-3,grd%im+3,1-3,grd%jm+3, grd%km, grd%msr)
      endif


    do j=1-grd%jas,grd%jm+grd%jae      ! 1,grd%jm
     do i=1-grd%ias,grd%im+grd%iae      ! 1,grd%im
       grd%f(i,j) = 0.00014584_r8 *                &
                    sin(grd%lat(i,j)*3.141592654_r8/180._r8)
      enddo
     enddo


         ich = 0
       do j=1,grd%jm+grd%jae
        do i=1,grd%im+grd%iae-1
         if(grd%lon(i+1,j)-grd%lon(i,j) .lt. -300.) ich = 1
        enddo
       enddo
     if(ich.eq.1)then
       do j=1-grd%jas,grd%jm+grd%jae
        do i=1-grd%ias,grd%im+grd%iae
         if(grd%lon(i,j).lt.0.0) grd%lon(i,j) = grd%lon(i,j) + 360.
        enddo
       enddo
     endif

      grd%bsth =  1.e20
      grd%bwst =  1.e20
      grd%bnrt = -1.e20
      grd%beas = -1.e20
    do j=1,grd%jm+grd%jae
     do i=1,grd%im+grd%iae
        grd%bsth = min(grd%bsth,grd%lat(i,j))
        grd%bnrt = max(grd%bnrt,grd%lat(i,j))
        grd%bwst = min(grd%bwst,grd%lon(i,j))
        grd%beas = max(grd%beas,grd%lon(i,j))
      enddo
     enddo


end subroutine def_grd
