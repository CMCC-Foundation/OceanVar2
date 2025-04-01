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
!> Structure of the grid 
!                                                                      !
! Version 1: Srdjan Dobricic   2006                                    !
!            Mario Adani       2023                                    !
!            Francesco Carere  2023                                    !
!-----------------------------------------------------------------------
MODULE grd_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Variables related to the grid
   TYPE grid_t

      INTEGER(i4)              ::  grd_mod           !< Grid model
      INTEGER(i4)              ::  prj               !< Projection (0 for lat-lon, 1 for general)
      INTEGER(i4)              ::  img               !< No. points in x direction on global grid
      INTEGER(i4)              ::  jmg               !< No. points in y direction on global grid
      INTEGER(i4)              ::  igs               !< Starting point in x direction on global grid
      INTEGER(i4)              ::  ige               !< End point in x direction on global grid
      INTEGER(i4)              ::  jgs               !< Starting point in y direction on global grid
      INTEGER(i4)              ::  jge               !< End point in y direction on global grid
      INTEGER(i4)              ::  ias               !< Additional starting point in x direction on inner tiles
      INTEGER(i4)              ::  iae               !< Additional END point in x direction on inner tiles
      INTEGER(i4)              ::  jas               !< Additional starting point in y direction on inner tiles
      INTEGER(i4)              ::  jae               !< Additional END point in y direction on inner tiles
      INTEGER(I4)              ::  is2               !< Starting point jump in i direction (invrt.F90)
      INTEGER(I4)              ::  is3               !< Starting point jump in i direction (invrt.F90)
      INTEGER(I4)              ::  js2               !< Starting point jump in j direction (invrt.F90)
      INTEGER(I4)              ::  js3               !< Starting point jump in j direction (invrt.F90)
      INTEGER(i4)              ::  im                !< No. points in x direction
      INTEGER(i4)              ::  jm                !< No. points in y direction
      INTEGER(i4)              ::  km                !< No. points in z direction
      INTEGER(i4)              ::  npsa              !< No. of ocean points on global grid
      INTEGER(i4)              ::  nps               !< No. of ocean points
      INTEGER(i4)              ::  iex               !< No. of extended points for the global grid

      INTEGER(I4), POINTER     ::  irs(:,:)          !< Starting point on MPI subtiles in x direction
      INTEGER(I4), POINTER     ::  ire(:,:)          !< Ending point on MPI subtiles in x direction
      INTEGER(I4), POINTER     ::  imr(:,:)          !< Number of points on MPI subtiles in x direction
      INTEGER(I4), POINTER     ::  jrs(:,:)          !< Starting point on MPI subtiles in y direction
      INTEGER(I4), POINTER     ::  jre(:,:)          !< Ending point on MPI subtiles in y direction
      INTEGER(I4), POINTER     ::  jmr(:,:)          !< Number of points on MPI subtiles in y direction

      INTEGER(I4), POINTER     ::  igsp(:)           !< Starting point on MPI tile in x direction (writing the output)
      INTEGER(I4), POINTER     ::  igep(:)           !< Ending point on MPI tile in x direction (writing the output)
      INTEGER(I4), POINTER     ::  jgsp(:)           !< Starting point on MPI tile in y direction (writing the output)
      INTEGER(I4), POINTER     ::  jgep(:)           !< Ending point on MPI tile in y direction (writing the output)
      INTEGER(I4)              ::  npsm              !< Maximum number of points on a tile

      REAL(r8)                 ::  dln               !< Resolution in the x direction
      REAL(r8)                 ::  dlt               !< Resolution in the y direction
      REAL(r8)                 ::  bsth              !< Latitude of the most southern point
      REAL(r8)                 ::  bnrt              !< Latitude of the most northern point
      REAL(r8)                 ::  bwst              !< Longitude of the most western point
      REAL(r8)                 ::  beas              !< Longitude of the most eastern point

      REAL(r8),    POINTER     ::  ro(:,:,:)         !< Reduced order control vector
      INTEGER(i4), POINTER     ::  reg(:,:)          !< Mask for EOF regions
      REAL(r8),    POINTER     ::  msk(:,:,:)        !< Sea-Land mask for scalar points
      REAL(r8),    POINTER     ::  ums(:,:,:)        !< Sea-Land mask for u points
      REAL(r8),    POINTER     ::  vms(:,:,:)        !< Sea-Land mask for v points
      REAL(r8),    POINTER     ::  hgt(:,:)          !< Topography
      REAL(r8),    POINTER     ::    f(:,:)          !< Coriolis term

      REAL(r8),    POINTER     ::  tem(:,:,:)        !< Temperature increment
      REAL(r8),    POINTER     ::  sal(:,:,:)        !< Salinity increment
      REAL(r8),    POINTER     ::  uvl(:,:,:)        !< u componnet of velocity increment
      REAL(r8),    POINTER     ::  vvl(:,:,:)        !< v componnet of velocity increment
      REAL(r8),    POINTER     ::  eta(:,:)          !< Sea level increment

      REAL(r8),    POINTER     ::  tes(:,:,:)        !< Saved temperature increment
      REAL(r8),    POINTER     ::  sas(:,:,:)        !< Saved salinity increment
      REAL(r8),    POINTER     ::  uvs(:,:,:)        !< Saved u componnet of velocity increment
      REAL(r8),    POINTER     ::  vvs(:,:,:)        !< Saved v componnet of velocity increment
      REAL(r8),    POINTER     ::  ets(:,:)          !< Saved sea level increment

      REAL(r8),    POINTER     ::  temb(:,:,:)       !< Temperature background
      REAL(r8),    POINTER     ::  salb(:,:,:)       !< Salinity background

      REAL(r8),    POINTER     ::  mdt(:,:)          !< Mean dynamic topography

      REAL(r8),    POINTER     ::  ro_ad(:,:,:)      !< Reduced order control vector adjoint
      REAL(r8),    POINTER     ::  tem_ad(:,:,:)     !< Temperature adjoint
      REAL(r8),    POINTER     ::  sal_ad(:,:,:)     !< Salinity adjoint
      REAL(r8),    POINTER     ::  uvl_ad(:,:,:)     !< u componnet of velocity adjoint
      REAL(r8),    POINTER     ::  vvl_ad(:,:,:)     !< v componnet of velocity adjoint
      REAL(r8),    POINTER     ::  eta_ad(:,:)       !< Sea level adjoint

      REAL(r8),    POINTER     ::  dns(:,:,:)        !< density
      REAL(r8),    POINTER     ::  b_x(:,:,:)        !< bouyancy force
      REAL(r8),    POINTER     ::  b_y(:,:,:)        !< bouyancy force
      REAL(r8),    POINTER     ::  bx(:,:)           !< bouyancy force integral
      REAL(r8),    POINTER     ::  by(:,:)           !< bouyancy force integral

      REAL(r8),    POINTER     ::  lon(:,:)          !< Longitude
      REAL(r8),    POINTER     ::  lat(:,:)          !< Latitude
      REAL(r8),    POINTER     ::  dep(:)            !< Depth

      REAL(r8)                 ::  lon1_1            !< Longitude at global gridpoint (1,1)
      REAL(r8)                 ::  lat1_1            !< Latitude at global gridpoint (1,1)

      REAL(r8),    POINTER     ::  loc(:,:)          !< Horizontal localisation


      REAL(r8),    POINTER     ::  dx(:,:)           !< dx
      REAL(r8),    POINTER     ::  dy(:,:)           !< dy
      REAL(r8),    POINTER     ::  dz(:)             !< dz
      REAL(r8),    POINTER     ::  dxdy(:,:)         !< dx*dy
      REAL(r8)                 ::  adxdy             !< Mean dx*dy


      REAL(r8),    POINTER     ::  alx(:,:,:)        !< Coefficient of the recursive filter
      REAL(r8),    POINTER     ::  aly(:,:,:)        !< Coefficient of the recursive filter
      REAL(r8),    POINTER     ::  btx(:,:,:)        !< Coefficient of the recursive filter
      REAL(r8),    POINTER     ::  bty(:,:,:)        !< Coefficient of the recursive filter
      REAL(r8),    POINTER     ::  gmx(:,:,:)        !< Coefficient  of the recursive filter
      REAL(r8),    POINTER     ::  gmy(:,:,:)        !< Coefficient  of the recursive filter
      REAL(r8),    POINTER     ::  dlx(:,:,:)        !< Coefficient  of the recursive filter
      REAL(r8),    POINTER     ::  dly(:,:,:)        !< Coefficient  of the recursive filter
      REAL(r8),    POINTER     ::  mat_bc_x(:,:,:,:) !< Boundary  Matrix in x direction  of the recursive filter
      REAL(r8),    POINTER     ::  mat_bc_y(:,:,:,:) !< Boundary  Matrix in y direction  of the recursive filter
      REAL(r8),    POINTER     ::  scx(:,:,:)        !< Scaling factor for x direction
      REAL(r8),    POINTER     ::  scy(:,:,:)        !< Scaling factor for y direction
      REAL(r8),    ALLOCATABLE ::  msr(:,:,:)        !< Sea-land mask USEd in the recursive filter
      INTEGER(i4), POINTER     ::  tmx(:)            !< Land or sea rows
      INTEGER(i4), POINTER     ::  tmy(:)            !< Land or sea columns

      INTEGER(i4), POINTER     ::  mxd(:,:)          !< Level at the bottom of the mixed layer
      REAL(r8),    POINTER     ::  lcl(:,:,:)        !< Vertical localization FUNCTIONs

      REAL(r8)                 ::  alpha             !< Thermal expansion coefficient
      REAL(r8)                 ::  beta              !< Haline  contraction coefficient
      REAL(r8),    POINTER     ::  alpha3d(:,:,:)    !< Thermal expansion coefficient
      REAL(r8),    POINTER     ::  beta3d(:,:,:)     !< Haline  contraction coefficient

      REAL(r8),    POINTER     ::  lon_rad(:,:)      !< Longitude in radians
      REAL(r8),    POINTER     ::  lat_rad(:,:)      !< Latitude  in radians
      REAL(r8),    POINTER     ::  dlon(:,:)         !< Grid spacing in radians
      REAL(r8),    POINTER     ::  dlat(:,:)         !< Grid spacing in radians
      REAL(r8),    POINTER     ::  distc3d(:,:,:)    !< distance from coast

   END TYPE grid_t

   TYPE (grid_t)                 :: grd              !< initialize derived type

END MODULE grd_str
