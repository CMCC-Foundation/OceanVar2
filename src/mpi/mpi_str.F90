MODULE mpi_str

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
! Structure of mpi                                                     !
!                                                                      !
! Version 1: S.Dobricic 2009                                           !
!-----------------------------------------------------------------------

 use set_knd

implicit none

public

   TYPE mpi_t

        INTEGER              ::  comm         ! Communicator
        INTEGER              ::  nproc        ! Number of processes
        INTEGER              ::  myrank       ! Rank of the process
        INTEGER              ::  top          ! Top processor
        INTEGER              ::  bot          ! Bottom processor
        INTEGER              ::  lft          ! Left processor
        INTEGER              ::  rgh          ! Right processor
        INTEGER              ::  r4           ! size of real*4 variables
        INTEGER              ::  r8           ! size of real*8 variables
        INTEGER              ::  i4           ! size of integer*4 variables
        INTEGER              ::  i8           ! size of integer*8 variables
        INTEGER              ::  irm          ! number of tiles in x direction
        INTEGER              ::  jrm          ! number of tiles in y direction
        INTEGER              ::  ir           ! position of the tile in x direction
        INTEGER              ::  jr           ! position of the tile in y direction
        INTEGER              ::  thx          ! Mult. factor for recursive filter threads (thi=nproc*thx) 
        INTEGER              ::  thy          ! Mult. factor for recursive filter threads (thj=nproc*thy) 
        INTEGER, POINTER     ::  thi(:)       ! Number of threads in i direction (recursive filter)
        INTEGER, POINTER     ::  thj(:)       ! Number of threads in j direction (recursive filter)
        INTEGER, POINTER     ::  psi(:)       ! position of the subtile in i direction
        INTEGER, POINTER     ::  psj(:)       ! position of the subtile in j direction
        INTEGER, POINTER     ::  topr(:)      ! Top processor for recursive filter
        INTEGER, POINTER     ::  botr(:)      ! Bottom processor for recursive filter
        INTEGER, POINTER     ::  lftr(:)      ! Left processor for recursive filter
        INTEGER, POINTER     ::  rghr(:)      ! Right processor for recursive filter

   END TYPE mpi_t

   TYPE (mpi_t)              :: mpi

END MODULE mpi_str
