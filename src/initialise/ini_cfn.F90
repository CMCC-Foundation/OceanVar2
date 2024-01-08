subroutine ini_cfn

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
! Initialise the minimisation                                          !
!                                                                      !
! Version 1: S.Dobricic 2006                                           !
!-----------------------------------------------------------------------

 use set_knd
 use drv_str
 use obs_str
 use grd_str
 use eof_str
 use ctl_str

 implicit none

 INTEGER(i4)  :: i


  ctl%task = 'START'

  ctl%iprint = 1
  ctl%iprint = 1 !99

! Oddo change the tollerance factor
  !ctl%factr=1.0d+7
  ctl%factr=1.0d+6

!  ctl%pgtol=1.0d-2


  if( drv%ktr.eq.1 .or. drv%ratio(drv%ktr).ne.1.0 ) then

! ---
! Allocate memory for optimization arrays

    ctl%n = grd%nps * ros%neof 

    write(drv%dia,*) ' ---- Size of the control vector: ',ctl%n
    ALLOCATE( ctl%nbd(ctl%n), ctl%iwa(3*ctl%n))
    ALLOCATE( ctl%x_c(ctl%n), ctl%g_c(ctl%n))
    ALLOCATE( ctl%l_c(ctl%n), ctl%u_c(ctl%n))
    ALLOCATE( ctl%wa(8*ctl%m), ctl%sg(ctl%m), ctl%sgo(ctl%m), ctl%yg(ctl%m), ctl%ygo(ctl%m))
    ALLOCATE( ctl%ws(ctl%n,ctl%m), ctl%wy(ctl%n,ctl%m))
    ALLOCATE( ctl%sy(ctl%m,ctl%m), ctl%ss(ctl%m,ctl%m), ctl%yy(ctl%m,ctl%m))
    ALLOCATE( ctl%wt(ctl%m,ctl%m), ctl%wn(2*ctl%m,2*ctl%m), ctl%snd(2*ctl%m,2*ctl%m))
    ALLOCATE( ctl%z_c(ctl%n), ctl%r_c(ctl%n), ctl%d_c(ctl%n), ctl%t_c(ctl%n))

    ctl%nbd(:) = 0
    ctl%iwa(:) = 0
    ctl%x_c(:) = 0.0
    ctl%g_c(:) = 0.0
    ctl%l_c(:) = 0.0
    ctl%u_c(:) = 0.0
    ctl%wa(:) = 0.0
    ctl%sg(:) = 0.0
    ctl%sgo(:) = 0.0
    ctl%yg(:) = 0.0
    ctl%ygo(:) = 0.0
    ctl%ws(:,:) = 0.0
    ctl%wy(:,:) = 0.0
    ctl%sy(:,:) = 0.0
    ctl%ss(:,:) = 0.0
    ctl%yy(:,:) = 0.0
    ctl%wt(:,:) = 0.0
    ctl%wn(:,:) = 0.0
    ctl%snd(:,:) = 0.0
    ctl%z_c(:) = 0.0
    ctl%r_c(:) = 0.0
    ctl%d_c(:) = 0.0
    ctl%t_c(:) = 0.0

! ---
! Initialise the arrays
      do i=1,ctl%n,2
         ctl%nbd(i)=0
         ctl%l_c(i)=-1.0d3
         ctl%u_c(i)=1.0d3
      enddo

      do i=2,ctl%n,2
         ctl%nbd(i)=0
         ctl%l_c(i)=-1.0d3
         ctl%u_c(i)=1.0d3
      enddo

      do i=1,ctl%n
         ctl%x_c(i)=0.0d0
      enddo

  endif

end subroutine ini_cfn
