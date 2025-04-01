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
!> Initialize the minimisation                                        
!!
!! It initializes the variables used in L-BFGS-B algorithm
!!
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!-----------------------------------------------------------------------
SUBROUTINE ini_cfn

   USE set_knd
   USE drv_str
   USE obs_str
   USE grd_str
   USE eof_str
   USE ctl_str

   IMPLICIT NONE

   INTEGER(i4)  :: i

   ctl%task = 'START'
   ctl%iprint = 1
   ctl%iprint = 1 
   ctl%factr=1.0d+6

   IF ( drv%ktr .EQ. 1 .OR. drv%ratio(drv%ktr) .NE. 1.0 ) THEN

! ---
! Allocate memory for optimization arrays
      ctl%n = grd%nps * ros%neof
      WRITE (drv%dia,*) ' ---- Size of the control vector: ',ctl%n
      ALLOCATE( ctl%nbd(ctl%n), ctl%iwa(3*ctl%n))
      ALLOCATE( ctl%x_c(ctl%n), ctl%g_c(ctl%n))
      ALLOCATE( ctl%l_c(ctl%n), ctl%u_c(ctl%n))
      ALLOCATE( ctl%wa(8*ctl%m), ctl%sg(ctl%m), ctl%sgo(ctl%m), ctl%yg(ctl%m), ctl%ygo(ctl%m) )
      ALLOCATE( ctl%ws(ctl%n,ctl%m), ctl%wy(ctl%n,ctl%m) )
      ALLOCATE( ctl%sy(ctl%m,ctl%m), ctl%ss(ctl%m,ctl%m), ctl%yy(ctl%m,ctl%m) )
      ALLOCATE( ctl%wt(ctl%m,ctl%m), ctl%wn(2*ctl%m,2*ctl%m), ctl%snd(2*ctl%m,2*ctl%m) )
      ALLOCATE( ctl%z_c(ctl%n), ctl%r_c(ctl%n), ctl%d_c(ctl%n), ctl%t_c(ctl%n) )

      ctl%nbd(:)   = 0
      ctl%iwa(:)   = 0
      ctl%x_c(:)   = 0.0_r8
      ctl%g_c(:)   = 0.0_r8
      ctl%l_c(:)   = 0.0_r8
      ctl%u_c(:)   = 0.0_r8
      ctl%wa(:)    = 0.0_r8
      ctl%sg(:)    = 0.0_r8
      ctl%sgo(:)   = 0.0_r8
      ctl%yg(:)    = 0.0_r8
      ctl%ygo(:)   = 0.0_r8
      ctl%ws(:,:)  = 0.0_r8
      ctl%wy(:,:)  = 0.0_r8
      ctl%sy(:,:)  = 0.0_r8
      ctl%ss(:,:)  = 0.0_r8
      ctl%yy(:,:)  = 0.0_r8
      ctl%wt(:,:)  = 0.0_r8
      ctl%wn(:,:)  = 0.0_r8
      ctl%snd(:,:) = 0.0_r8
      ctl%z_c(:)   = 0.0_r8
      ctl%r_c(:)   = 0.0_r8
      ctl%d_c(:)   = 0.0_r8
      ctl%t_c(:)   = 0.0_r8

! ---
! Initialise the arrays
      DO i = 1,ctl%n,2
         ctl%l_c(i) = -1.0d3
         ctl%u_c(i) =  1.0d3
      ENDDO

      DO i = 2,ctl%n,2
         ctl%l_c(i) = -1.0d3
         ctl%u_c(i) =  1.0d3
      ENDDO

      DO i = 1,ctl%n
         ctl%x_c(i) = 0.0_r8
      ENDDO
#ifdef REPRO
      CALL get_indexes
#endif

   ENDIF

END SUBROUTINE ini_cfn
