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
!> Structure of diffusive filter 
!                                                                      !
! Version 1: Mario Adani 2024                                          !
!-----------------------------------------------------------------------
MODULE dfl_str


   USE set_knd

   IMPLICIT NONE

   PUBLIC
! ---
!>  Diffusive filter variables
   TYPE dfl_t

      INTEGER(i4)                   :: nt          !< number of step the filter is applyed
      LOGICAL                       :: rd_corr     !< activate read of spatial varing radius
      CHARACTER(len=:),ALLOCATABLE  :: crl_flname  !< filename of spatial varing radius
      REAL(r8)                      :: rx          !< constant ratdius x
      REAL(r8)                      :: ry          !< constant ratdius ry
      LOGICAL                       :: use_bc      !< activate boundary condition on land
      CHARACTER(len=:),ALLOCATABLE  :: bc_type     !< type of boundary condition applyed on land
      LOGICAL                       :: rd_wgh      !< activate read of the weighting factors
      CHARACTER(len=:),ALLOCATABLE  :: wgh_flname  !< filename of weighting factors
      LOGICAL                       :: use_cst     !< activate reduction of correlation radius close to coast
      REAL(r8)                      :: cst_dst     !< distance from coast
      REAL(r8)                      :: level       !< level to be filtered
      REAL(r8),POINTER              :: wgh(:,:,:)  !< Numerical weight for diffusion filter
      REAL(r8),POINTER              :: Lx(:,:,:)   !< Lower Triangular in x
      REAL(r8),POINTER              :: Ux(:,:,:)   !< Upper Triangular in x
      REAL(r8),POINTER              :: Ax(:,:,:)   !< Subdiagonal coefficients in x
      REAL(r8),POINTER              :: Ly(:,:,:)   !< Lower Triangular in y
      REAL(r8),POINTER              :: Uy(:,:,:)   !< Upper Triangular in y
      REAL(r8),POINTER              :: Ay(:,:,:)   !< Subdiagonal coefficients in y
      REAL(r8),POINTER              :: rx3d(:,:,:) !< 3D correlation radius x
      REAL(r8),POINTER              :: ry3d(:,:,:) !< 3D correlation radius y
      REAL(r8),POINTER              :: kx(:,:,:)   !< Diffusion coefficients x
      REAL(r8),POINTER              :: ky(:,:,:)   !< Diffusion coefficients y

   END TYPE dfl_t

   TYPE (dfl_t)              :: dfl                !< initialize derived type

END MODULE dfl_str
