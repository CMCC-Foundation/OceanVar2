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
!> Structure of EOFs 
!                                                                      !
! Version 1: Srdjan Dobricic 2024                                      !
!-----------------------------------------------------------------------
MODULE eof_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Variable related to EOFs
   TYPE eof_t

      CHARACTER(len=:),ALLOCATABLE ::  flname       !< EOF filename
      LOGICAL                      ::  already_read !< Flag if EOFS are already read in grd file
      INTEGER(i4)                  ::  neof         !< No. of EOFs
      INTEGER(i4)                  ::  nreg         !< No. of regions
      INTEGER(i4)                  ::  kmt          !< No. of levels of EOFs
      REAL(r8),    POINTER         ::  evcr(:,:,:)  !< Eigenvectors on regions
      REAL(r8),    POINTER         ::  evar(:,:)    !< Eigenvalues on regions
      REAL(r8),    POINTER         ::  corr(:,:,:)  !< Corelations on regions
#ifdef opt_huge_memory
      REAL(r8),    POINTER         ::  evc(:,:,:,:) !< Eigenvectors
      REAL(r8),    POINTER         ::  eva(:,:,:)   !< Eigenvalues
      REAL(r8),    POINTER         ::  cor(:,:,:,:) !< Corelation matrix
#else
      REAL(r8),    POINTER         ::  evc(:,:,:)   !< Eigenvectors
      REAL(r8),    POINTER         ::  eva(:,:)     !< Eigenvalues
      REAL(r8),    POINTER         ::  cor(:,:,:)   !< Corelation matrix
#endif
      INTEGER(i4)                  ::  mld          !< Indicator to use the mixed layer depth information

   END TYPE eof_t

   TYPE (eof_t)                 :: ros              !< initialize derived type

END MODULE eof_str
