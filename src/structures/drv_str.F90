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
!> Structure for the driver of the outer loop 
!                                                                      !
! Version 1: Srdjan Dobricic 2006                                      !
!            Mario Adani     2024                                      !
!-----------------------------------------------------------------------
MODULE drv_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

!---
!> Variables related to the driver
   TYPE drv_t

      CHARACTER(len=:),ALLOCATABLE ::  flag              !< Flag for the analysis
      CHARACTER(len=:),ALLOCATABLE ::  eosflname         !< Expantion/contraction coefficent filename
      CHARACTER(len=:),ALLOCATABLE ::  inpdir            !< Input directory
      INTEGER(I4)                  ::  sdat              !< Starting date of the forecast
      INTEGER(I4)                  ::  shou              !< Starting hour of the forecast
      INTEGER(I4)                  ::  lhou              !< Length of the forecast
      REAL(r8)                     ::  zanjul1950        !< Analysis time from 1950/01/01
      INTEGER(i4)                  ::  dia               !< No. of diagnostic output file
      INTEGER(i4)                  ::  nts               !< No. of outer iterations - SST assimilation
      INTEGER(i4)                  ::  kts               !< Outer iteration - SST assimilation
      INTEGER(i4)                  ::  ntr               !< No. of outer iterations - multigrid
      INTEGER(i4)                  ::  ktr               !< Outer iteration- multigrid
      INTEGER(i4)                  ::  im                !< Dimension of the coarse grid
      INTEGER(i4)                  ::  jm                !< Dimension of the coarse grid
      INTEGER(i4)                  ::  filter            !< Filter type: 1-recursive, 2-dIFfusive, 3-None
      INTEGER(i4), POINTER         ::  nneos(:)          !< Expantion coefficient
      LOGICAL    , POINTER         ::  ssh_unbalanced(:) !< Use EOFs component
      INTEGER(i4), POINTER         ::  grid(:)           !< grid number for the current iterration
      REAL(r8),    POINTER         ::  ratco(:)          !< Ratio between model grid and the current grid
      REAL(r8),    POINTER         ::  ratio(:)          !< Ratio between successive grids
      INTEGER(i4), POINTER         ::  mask(:)           !< Mask USEd for horizontal covariances
      INTEGER(i4), POINTER         ::  bmd(:)            !< 1 - run barotropic model, else - do not run
      INTEGER(i4), POINTER         ::  bal(:)            !< 1 - run simplified balance operator (dynamic height), else - do not run
      INTEGER(i4), POINTER         ::  dda(:)            !< 1 - divergence damping in analysis, else no filter
      INTEGER(i4), POINTER         ::  ddi(:)            !< 1 - divergence damping in initialisation, else no filter
      INTEGER(i4), POINTER         ::  cntm(:)           !< Maximum number of cost function evaluations
      REAL(r8)                     ::  f_ci              !< Inital cost function
      REAL(r8),    POINTER         ::  ro(:,:,:)         !< Vector v
      REAL(r8)                     ::  f_c               !< Cost function
      REAL(r8),    POINTER         ::  ro_ad(:,:,:)      !< Observational part of the cost function gradient
      REAL(r8),    POINTER         ::  msk(:,:)          !< Mask of the old grid

   END TYPE drv_t

   TYPE (drv_t)                    :: drv                !< initialize derived type

END MODULE drv_str
