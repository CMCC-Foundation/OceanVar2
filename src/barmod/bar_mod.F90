subroutine bar_mod

!---------------------------------------------------------------------------
!                                                                          !
!    Copyright 2007 Srdjan Dobricic, CMCC, Bologna                         !
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
! Barotropic model                                                     !
!                                                                      !
! Version 1: S.Dobricic 2007                                           !
!-----------------------------------------------------------------------


  use set_knd
  use drv_str
  use grd_str
  use bmd_str
  use mpi_str

  implicit none
  INTEGER(i4)    :: i, j, kstp

 if(mpi%nproc.gt.1)then

  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, 1_i4, grd%jm+1,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%bx )
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, 1_i4, grd%im+1,                               &
                -1_i4, 1_i4, grd%jm+1,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, grd%by )

 endif

!---
! Bouyancy forcing
    bmd%bxby(:,:)=0._r8
    bmd%rgh(:,:)=0._r8
    do j=2-grd%jas,grd%jm-1+grd%jae                 ! 2,grd%jm-1
     do i=2-grd%ias,grd%im-1+grd%iae                ! 2,grd%im-1
        bmd%bxby(i,j) = - (bmd%dt**2)*((grd%bx(i+1,j)-grd%bx(i,j))/grd%dx(i,j) +               &
                                       (grd%by(i,j+1)-grd%by(i,j))/grd%dy(i,j)) * bmd%mst(i,j)
        bmd%rgh(i,j) = bmd%bxby(i,j) * bmd%alp1**2
       enddo
      enddo

       bmd%etb(:,:) = 0.0_r8
        bmd%ub(:,:) = 0.0_r8
        bmd%vb(:,:) = 0.0_r8
       bmd%etn(:,:) = 0.0_r8
        bmd%un(:,:) = 0.0_r8
        bmd%vn(:,:) = 0.0_r8
       bmd%eta(:,:) = 0.0_r8
        bmd%ua(:,:) = 0.0_r8
        bmd%va(:,:) = 0.0_r8

          bmd%etx(:,:) = 0.0_r8
          bmd%ety(:,:) = 0.0_r8
          bmd%div(:,:) = 0.0_r8
           bmd%cu(:,:) = 0.0_r8
           bmd%cv(:,:) = 0.0_r8
          bmd%dux(:,:) = 0.0_r8
          bmd%dvx(:,:) = 0.0_r8
          bmd%duy(:,:) = 0.0_r8
          bmd%dvy(:,:) = 0.0_r8

         bmd%etm(:,:) = 0.0_r8
          bmd%um(:,:) = 0.0_r8
          bmd%vm(:,:) = 0.0_r8


!---
! Time loop

     do kstp = 1,bmd%nstps


       do j=2-grd%jas,grd%jm-1+grd%jae                ! 2,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae                ! 2,grd%im-1
          bmd%div(i,j) = ((bmd%ub(i+1,j)-bmd%ub(i,j))/grd%dx(i,j) + (bmd%vb(i,j+1)-bmd%vb(i,j))/grd%dy(i,j)) * bmd%mst(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm-1+grd%jae    ! 2,grd%jm-1
        do i=2-grd%ias,grd%im             ! 2,grd%im
          bmd%etx(i,j) = bmd%alp2*bmd%dt*bmd%g*bmd%hgu(i,j)*(bmd%etb(i,j)-bmd%etb(i-1,j  ))/bmd%dxu(i,j) * bmd%msu(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm              ! 2,grd%jm
        do i=2-grd%ias,grd%im-1+grd%iae   ! 2,grd%im-1
          bmd%ety(i,j) = bmd%alp2*bmd%dt*bmd%g*bmd%hgv(i,j)*(bmd%etb(i,j)-bmd%etb(i  ,j-1))/bmd%dyv(i,j) * bmd%msv(i,j)
        enddo
       enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, 1_i4, grd%jm+1,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%etx)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, 1_i4, grd%im+1,                               &
                -1_i4, 1_i4, grd%jm+1,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ety)
 endif

       do j=2-grd%jas,grd%jm-1+grd%jae        ! 2,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae       ! 2,grd%im-1
         bmd%rgh(i,j) = bmd%bxby(i,j) - bmd%etb(i,j) + bmd%dt*bmd%div(i,j)        &
                   - bmd%alp1*bmd%dt*(bmd%etx(i+1,j)-bmd%etx(i,j))/grd%dx(i,j)    &
                   - bmd%alp1*bmd%dt*(bmd%ety(i,j+1)-bmd%ety(i,j))/grd%dy(i,j)
        enddo
       enddo

!---
! Calculate sea level
       call invrt( kstp )


       do j=1,grd%jm-1+grd%jae              ! 1,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae     ! 2,grd%im-1
         bmd%cu(i,j) = -((bmd%vn(i,j  )*bmd%dxv(i,j  ) + bmd%vn(i-1,j  )*bmd%dxv(i-1,j  ))*grd%f(i,j  ) +     &
                         (bmd%vn(i,j+1)*bmd%dxv(i,j+1) + bmd%vn(i-1,j+1)*bmd%dxv(i-1,j+1))*grd%f(i,j+1))      &
                        * 0.25_r8 / bmd%dxu(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm-1+grd%jae    ! 2,grd%jm-1
        do i=1,grd%im-1+grd%iae           ! 1,grd%im-1
         bmd%cv(i,j) =  ((bmd%un(i  ,j)*bmd%dyu(i  ,j) + bmd%un(i  ,j-1)*bmd%dyu(i  ,j-1))*grd%f(i  ,j) +      &
                         (bmd%un(i+1,j)*bmd%dyu(i+1,j) + bmd%un(i+1,j-1)*bmd%dyu(i+1,j-1))*grd%f(i+1,j))       &
                        * 0.25_r8 / bmd%dyv(i,j)
        enddo
       enddo
       do j=1,grd%jm                              ! 1,grd%jm
        do i=2-grd%ias,grd%im                     ! 2,grd%im
          bmd%dux(i,j) = (bmd%ub(i,j) - bmd%ub(i-1,j))/bmd%dxu(i,j)
          bmd%dvx(i,j) = (bmd%vb(i,j) - bmd%vb(i-1,j))/bmd%dxu(i,j)
        enddo
       enddo

       do j=2-grd%jas,grd%jm                      ! 2,grd%jm
        do i=1,grd%im                             ! 1,grd%im
          bmd%duy(i,j) = (bmd%ub(i,j) - bmd%ub(i,j-1))/bmd%dyv(i,j)
          bmd%dvy(i,j) = (bmd%vb(i,j) - bmd%vb(i,j-1))/bmd%dyv(i,j)
        enddo
       enddo

! MPI

 if(mpi%nproc.gt.1)then
  call exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, grd%im, 0_i4,                               &
                -1_i4, 1_i4, grd%jm+1,                              &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, grd%im, 0_i4,                               &
                -1_i4, 1_i4, grd%jm+1,                              &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
 endif

!---
! Calculate new velocity
       do j=1,grd%jm-1+grd%jae                    ! 1,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae                   ! 2,grd%im
          bmd%ua(i,j) = bmd%ub(i,j) - bmd%dt*(bmd%cu(i,j) +                                           &
                        bmd%alp1*bmd%g*bmd%hgu(i,j)*(bmd%eta(i,j)-bmd%eta(i-1,j  ))/bmd%dxu(i,j) +    &
                        grd%bx(i,j))*bmd%msu(i,j)  -                                                  &
                        bmd%etx(i,j) +                                                                &
                        bmd%df1 * ((bmd%dux(i+1,j)-bmd%dux(i,j))/bmd%dxu(i,j) +                       &
                                   (bmd%duy(i,j+1)-bmd%duy(i,j))/bmd%dyu(i,j)) * bmd%msu(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm-1+grd%jae                    ! 2,grd%jm
        do i=1,grd%im-1+grd%iae                   ! 1,grd%im-1
          bmd%va(i,j) = bmd%vb(i,j) - bmd%dt*(bmd%cv(i,j) +                                           &
                        bmd%alp1*bmd%g*bmd%hgv(i,j)*(bmd%eta(i,j)-bmd%eta(i  ,j-1))/bmd%dyv(i,j) +    &
                        grd%by(i,j))*bmd%msv(i,j) -                                                   &
                        bmd%ety(i,j) +                                                                &
                        bmd%df1 * ((bmd%dvx(i+1,j)-bmd%dvx(i,j))/bmd%dxv(i,j) +                       &
                                   (bmd%dvy(i,j+1)-bmd%dvy(i,j))/bmd%dyv(i,j)) * bmd%msv(i,j)
        enddo
       enddo

 if(mpi%nproc.gt.1)then
  call exo_mpi( 1_i4, 1_i4,                                       &
                 1_i4, grd%im, 0_i4,                               &
                 1_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%ua)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 1_i4, grd%im, 0_i4,                               &
                 1_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%va)
 endif

       do j=1,grd%jm                          ! 1,grd%jm
        do i=2-grd%ias,grd%im                 ! 2,grd%im
          bmd%dux(i,j) = (bmd%ua(i,j) - bmd%ua(i-1,j))/bmd%dxu(i,j)
          bmd%dvx(i,j) = (bmd%va(i,j) - bmd%va(i-1,j))/bmd%dxv(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm                  ! 2,grd%jm
        do i=1,grd%im                         ! 1,grd%im
          bmd%duy(i,j) = (bmd%ua(i,j) - bmd%ua(i,j-1))/bmd%dyu(i,j)
          bmd%dvy(i,j) = (bmd%va(i,j) - bmd%va(i,j-1))/bmd%dyv(i,j)
        enddo
       enddo


 if(mpi%nproc.gt.1)then
  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dux)
  call exo_mpi( 1_i4, 1_i4,                                       &
                -1_i4, 1_i4, grd%im+1,                               &
                 0_i4, grd%jm, 0_i4,                               &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvx)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, grd%im, 0_i4,                               &
                -1_i4, 1_i4, grd%jm+1,                              &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%duy)
  call exo_mpi( 1_i4, 1_i4,                                       &
                 0_i4, grd%im, 0_i4,                               &
                -1_i4, 1_i4, grd%jm+1,                              &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae, 1_i4, bmd%dvy)
 endif


       do j=2-grd%jas,grd%jm-1+grd%jae        ! 2,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae       ! 2,grd%im-1
          bmd%ua(i,j) = bmd%ua(i,j) + bmd%df2 * ((bmd%dux(i+1,j)-bmd%dux(i,j))/bmd%dxu(i,j) +      &
                                                 (bmd%duy(i,j+1)-bmd%duy(i,j))/bmd%dyu(i,j))*bmd%msu(i,j)
        enddo
       enddo
       do j=2-grd%jas,grd%jm-1+grd%jae        ! 2,grd%jm-1
        do i=2-grd%ias,grd%im-1+grd%iae       ! 2,grd%im-1
          bmd%va(i,j) = bmd%va(i,j) + bmd%df2 * ((bmd%dvx(i+1,j)-bmd%dvx(i,j))/bmd%dxv(i,j) +      &
                                                 (bmd%dvy(i,j+1)-bmd%dvy(i,j))/bmd%dyv(i,j))*bmd%msv(i,j)
        enddo
       enddo


 if(mpi%nproc.gt.1)then
  call exa_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%eta)
  call eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%ua )
  call eao_mpi( 1_i4, 1_i4, grd%im, 0_i4, grd%im+1, 1_i4, grd%jm, 0_i4, grd%jm+1,     &
                1-grd%ias, grd%im+grd%iae, 1-grd%jas, grd%jm+grd%jae , 1_i4, bmd%va )
 endif

!---
! Asselin filter
           bmd%un(:,:) =  bmd%un(:,:) + ( bmd%ub(:,:) + bmd%ua(:,:) - 2.0_r8*bmd%un(:,:) ) *0.05_r8
           bmd%vn(:,:) =  bmd%vn(:,:) + ( bmd%vb(:,:) + bmd%va(:,:) - 2.0_r8*bmd%vn(:,:) ) *0.05_r8

          bmd%etb(:,:) = bmd%eta(:,:)
           bmd%ub(:,:) =  bmd%un(:,:)
           bmd%vb(:,:) =  bmd%vn(:,:)

           bmd%un(:,:) =  bmd%ua(:,:)
           bmd%vn(:,:) =  bmd%va(:,:)

          bmd%etm(:,:) = bmd%etm(:,:) + bmd%etb(:,:)
           bmd%um(:,:) =  bmd%um(:,:) +  bmd%ub(:,:)
           bmd%vm(:,:) =  bmd%vm(:,:) +  bmd%vb(:,:)


!---
! Temporal average
       if( mod(kstp,bmd%nstpa).eq.0 ) then

          bmd%etm(:,:) = bmd%etm(:,:)/bmd%nstpa * bmd%mst(:,:)
          bmd%um(:,:) =  bmd%um(:,:)/bmd%nstpa
          bmd%vm(:,:) =  bmd%vm(:,:)/bmd%nstpa

          grd%eta(:,:) = bmd%etm(:,:)

          bmd%etm(:,:) = 0.0_r8
          bmd%um(:,:) = 0.0_r8
          bmd%vm(:,:) = 0.0_r8

       endif   !Temporal average
     enddo   ! kstp

end subroutine bar_mod
