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
!> Structure for observational vectorss
!!
!! Data structure for the observational vectors, huber norm,
!! observational error, thinning, quality control, coastal rejection
!!
!                                                                      !
! Version 1: Srdjan Dobricic  2006                                     !
!            Adani Mario      2024                                     !
!-----------------------------------------------------------------------
MODULE obs_str

   USE set_knd

   IMPLICIT NONE

   PUBLIC

! ---
!> Observational vector in the cost function
   TYPE obs_t

      INTEGER(i8)              ::  no         !< Number of local observations
      INTEGER(i8)              ::  k          !< Observation index
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  amo(:)     !< Analysis - observation
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  gra(:)     !< Observational gradient
      INTEGER(i8)              ::  sla        !< Flag for assimilation of SLA
      INTEGER(i8)              ::  arg        !< Flag for assimilation of ARGO floats
      INTEGER(i8)              ::  xbt        !< Flag for assimilation of XBTs
      INTEGER(i8)              ::  gld        !< Flag for assimilation of gliders
      INTEGER(i8)              ::  tra        !< Flag for assimilation of Argo trajectories
      INTEGER(i8)              ::  trd        !< Flag for assimilation of drIFtres trajectories
      INTEGER(i8)              ::  vdr        !< Flag for assimilation of velocities by drIFters
      INTEGER(i8)              ::  gvl        !< Flag for assimilation of velocities by gliders
      INTEGER(i8)              ::  sst        !< Flag for assimilation of SST from satellite
      REAL(r8),    POINTER     ::  ahub(:,:)  !< Coefficient of Huber Norm Distribution
      REAL(r8),    POINTER     ::  ahub2(:,:) !< Coefficient of Huber Norm Distribution

   END TYPE obs_t

   TYPE (obs_t)                :: obs         !< Inizialization of the observational vector

! ---
!> Observational vector for SLA
   TYPE sla_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      REAL(r8)                 ::  dsm        !< Maximum distance between obs sla point to be considered of the same track
      REAL(r8)                 ::  minobspt   !< Minimum number of observation per track

      LOGICAL                  ::  unbias     !< T - Remove global bias
      LOGICAL                  ::  bias_at    !< Remove bias along track
      INTEGER(i8)              ::  kdp        !< Model level corresponding to dep
      INTEGER(i8), POINTER     ::  ino(:)     !< Track number
      INTEGER(i8), POINTER     ::  ksat(:)    !< Satellite ID ( following the same sequence in the NAMELIST)
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      REAL(r8),    POINTER     ::  lon(:)     !< Longitute
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      ! j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  dpt(:)     !< Maximum depth of surrounding points
      REAL(r8)   , POINTER     ::  dtm(:)     !< Mean of the dynamic topography
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE sla_t

   TYPE (sla_t)                 :: sla        !< Inizialization of sla observational vector

! ---
!> Observational vector for ARGO floats
   TYPE arg_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      REAL(r8)                 ::  dsm        !< Max distance allowed between sla obs to be considered of the same track
      INTEGER(i8)              ::  kdp        !< Model level corresponding to dep
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1-temperature, 2-salinity)
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      REAL(r8),    POINTER     ::  lon(:)     !< Longitute
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  dpt(:)     !< Depth
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      INTEGER(i8), POINTER     ::  kb(:)      !< k index of the nearest point below
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq5(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq6(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq7(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq8(:)     !< Interpolation parameter for masked grids
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE arg_t

   TYPE (arg_t)                 :: arg        !< Inizialization of argo observational vector

! ---
!> Observational vector for XBT profiles
   TYPE xbt_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      INTEGER(i8)              ::  kdp        !< Model level corresponding to dep
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1-temperature)
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      REAL(r8),    POINTER     ::  lon(:)     !< Longitute
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  dpt(:)     !< Depth
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      INTEGER(i8), POINTER     ::  kb(:)      !< k index of the nearest point below
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq5(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq6(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq7(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq8(:)     !< Interpolation parameter for masked grids
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE xbt_t

   TYPE (xbt_t)                 :: xbt        !< Inizialization of xbt observational vector

! ---
!> Observational vector for gliders
   TYPE gld_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      INTEGER(i8)              ::  kdp        !< Model level corresponding to dep
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1-temperature, 2-salinity)
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      REAL(r8),    POINTER     ::  lon(:)     !< Longitute
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  dpt(:)     !< Depth
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      INTEGER(i8), POINTER     ::  kb(:)      !< k index of the nearest point below
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq5(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq6(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq7(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq8(:)     !< Interpolation parameter for masked grids
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE gld_t

   TYPE (gld_t)                :: gld         !< Inizialization of the glider observational vector

! ---
!> Observational vector for Argo trajectories
   TYPE tra_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      INTEGER(i8)              ::  ncc        !< Global information on number of good observations
      INTEGER(i8)              ::  ncs        !< Saved global information on number of good observations
      INTEGER(i8)              ::  nt         !< Number of time steps
      REAL(r8)                 ::  dpt        !< Depth of observations
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      REAL(r8),    POINTER     ::  dtm(:)     !< Time spent under surface
      REAL(r8),    POINTER     ::  loi(:)     !< Initial observed longitude
      REAL(r8),    POINTER     ::  lai(:)     !< Initial observed latitude
      REAL(r8),    POINTER     ::  lof(:)     !< Final observed longitude
      REAL(r8),    POINTER     ::  laf(:)     !< Final observed latitude
      REAL(r8),    POINTER     ::  lob(:,:)   !< Longitudes of simulated positions
      REAL(r8),    POINTER     ::  lab(:,:)   !< Latitudes of simulated positions
      REAL(r8),    POINTER     ::  loa(:)     !< Longitudes of analysed positions
      REAL(r8),    POINTER     ::  laa(:)     !< Latitudes of analysed positions
      REAL(r8),    POINTER     ::  xob(:)     !< Observed value in x direction
      REAL(r8),    POINTER     ::  erx(:)     !< Observational error in x direction
      REAL(r8),    POINTER     ::  rex(:)     !< Misfit in x direction
      REAL(r8),    POINTER     ::  rsx(:)     !< Saved misfit in x direction
      REAL(r8),    POINTER     ::  inx(:)     !< Increments in x direction
      REAL(r8),    POINTER     ::  isx(:)     !< Saved increments in x direction
      REAL(r8),    POINTER     ::  yob(:)     !< Observed value in y direction
      REAL(r8),    POINTER     ::  ery(:)     !< Observational error in y direction
      REAL(r8),    POINTER     ::  rey(:)     !< Misfit in y direction
      REAL(r8),    POINTER     ::  rsy(:)     !< Saved misfit in y direction
      REAL(r8),    POINTER     ::  iny(:)     !< Increments in y direction
      REAL(r8),    POINTER     ::  isy(:)     !< Saved increments in y direction
      REAL(r8),    POINTER     ::  err(:)     !< Observational error in meters
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8)   , POINTER     ::  xmn(:,:)   !< X coordinate of mean trajectory position
      REAL(r8)   , POINTER     ::  ymn(:,:)   !< Y coordinate of mean trajectory position
      REAL(r8)   , POINTER     ::  tim(:)     !< Time of the duration of the trajectory
      REAL(r8)   , POINTER     ::  xtl(:)     !< Delta X of trajectory correction
      REAL(r8)   , POINTER     ::  ytl(:)     !< Delta Y of trajectory correction
      REAL(r8)   , POINTER     ::  xtl_ad(:)  !< Delta X of trajectory correction (adjoint)
      REAL(r8)   , POINTER     ::  ytl_ad(:)  !< Delta Y of trajectory correction (adjoint)
      INTEGER(i8), POINTER     ::  i1(:,:)    !< i index for interpolation between grids
      INTEGER(i8), POINTER     ::  j1(:,:)    !< j index for interpolation between grids
      INTEGER(I4)              ::  im         !< I dimension of the grid
      INTEGER(I4)              ::  jm         !< J dimension of the grid
      INTEGER(I4)              ::  km         !< K dimension of the grid
      INTEGER(I4)              ::  lev        !< level on the grid for the trajectory
      REAL(r8)   , POINTER     ::  umn(:,:)   !< U component of background velocity
      REAL(r8)   , POINTER     ::  vmn(:,:)   !< V component of background velocity
      REAL(r8)   , POINTER     ::  uvl(:,:)   !< Delta u on trajectory grid
      REAL(r8)   , POINTER     ::  vvl(:,:)   !< Delta v on trajectory grid
      REAL(r8)   , POINTER     ::  uvl_ad(:,:)!< Delta u on trajectory grid (adjoint)
      REAL(r8)   , POINTER     ::  vvl_ad(:,:)!< Delta v on trajectory grid (adjoint)
      REAL(r8)   , POINTER     ::   dx(:,:)   !< Delta x on trajectory grid
      REAL(r8)   , POINTER     ::   dy(:,:)   !< Delta y on trajectory grid
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE tra_t

   TYPE (tra_t)                :: tra         !< Inizialization of the argo trajectory observational vector

! ---
!> Observational vector for trajectories of surface drifters
   TYPE trd_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      INTEGER(i8)              ::  ncc        !< Global information on number of good observations
      INTEGER(i8)              ::  ncs        !< Saved global information on number of good observations
      INTEGER(i8)              ::  nt         !< Number of time steps
      REAL(r8)                 ::  dpt        !< Depth of observations
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      REAL(r8),    POINTER     ::  dtm(:)     !< Time spent under surface
      REAL(r8),    POINTER     ::  loi(:)     !< Initial observed longitude
      REAL(r8),    POINTER     ::  lai(:)     !< Initial observed latitude
      REAL(r8),    POINTER     ::  lof(:)     !< Final observed longitude
      REAL(r8),    POINTER     ::  laf(:)     !< Final observed latitude
      REAL(r8),    POINTER     ::  lob(:,:)   !< Longitudes of simulated positions
      REAL(r8),    POINTER     ::  lab(:,:)   !< Latitudes of simulated positions
      REAL(r8),    POINTER     ::  loa(:)     !< Longitudes of analysed positions
      REAL(r8),    POINTER     ::  laa(:)     !< Latitudes of analysed positions
      REAL(r8),    POINTER     ::  xob(:)     !< Observed value in x direction
      REAL(r8),    POINTER     ::  erx(:)     !< Observational error in x direction
      REAL(r8),    POINTER     ::  rex(:)     !< Misfit in x direction
      REAL(r8),    POINTER     ::  rsx(:)     !< Saved misfit in x direction
      REAL(r8),    POINTER     ::  inx(:)     !< Increments in x direction
      REAL(r8),    POINTER     ::  isx(:)     !< Saved increments in x direction
      REAL(r8),    POINTER     ::  yob(:)     !< Observed value in y direction
      REAL(r8),    POINTER     ::  ery(:)     !< Observational error in y direction
      REAL(r8),    POINTER     ::  rey(:)     !< Misfit in y direction
      REAL(r8),    POINTER     ::  rsy(:)     !< Saved misfit in y direction
      REAL(r8),    POINTER     ::  iny(:)     !< Increments in y direction
      REAL(r8),    POINTER     ::  isy(:)     !< Saved increments in y direction
      REAL(r8),    POINTER     ::  err(:)     !< Observational error in meters
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8)   , POINTER     ::  xmn(:,:)   !< X coordinate of mean trajectory position
      REAL(r8)   , POINTER     ::  ymn(:,:)   !< Y coordinate of mean trajectory position
      REAL(r8)   , POINTER     ::  tim(:)     !< Time of the duration of the trajectory
      REAL(r8)   , POINTER     ::  xtl(:)     !< Delta X of trajectory correction
      REAL(r8)   , POINTER     ::  ytl(:)     !< Delta Y of trajectory correction
      REAL(r8)   , POINTER     ::  xtl_ad(:)  !< Delta X of trajectory correction (adjoint)
      REAL(r8)   , POINTER     ::  ytl_ad(:)  !< Delta Y of trajectory correction (adjoint)
      INTEGER(i8), POINTER     ::  i1(:,:)    !< i index for interpolation between grids
      INTEGER(i8), POINTER     ::  j1(:,:)    !< j index for interpolation between grids
      INTEGER(I4)              ::  im         !< I dimension of the grid
      INTEGER(I4)              ::  jm         !< J dimension of the grid
      INTEGER(I4)              ::  km         !< K dimension of the grid
      INTEGER(I4)              ::  lev        !< level on the grid for the trajectory
      REAL(r4)   , POINTER     ::  umn(:,:)   !< U component of background velocity
      REAL(r4)   , POINTER     ::  vmn(:,:)   !< V component of background velocity
      REAL(r4)   , POINTER     ::  uvl(:,:)   !< Delta u on trajectory grid
      REAL(r4)   , POINTER     ::  vvl(:,:)   !< Delta v on trajectory grid
      REAL(r4)   , POINTER     ::  uvl_ad(:,:)!< Delta u on trajectory grid (adjoint)
      REAL(r4)   , POINTER     ::  vvl_ad(:,:)!< Delta v on trajectory grid (adjoint)
      REAL(r4)   , POINTER     ::  lon(:,:)   !< Longitudes of grid points
      REAL(r4)   , POINTER     ::  lat(:,:)   !< Latitudes of grid points
      REAL(r4)   , POINTER     ::   dx(:,:)   !< Delta x on trajectory grid
      REAL(r4)   , POINTER     ::   dy(:,:)   !< Delta y on trajectory grid
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE trd_t

   TYPE (trd_t)                :: trd         !< Inizialization of drifters trajectory observational vector

! ---
!> Observational vector for velocity from drifters
   TYPE vdr_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1 - u component, 2 - v component)
      REAL(r8),    POINTER     ::  lon(:)     !< Longitude
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  dpt(:)     !< Depth
      INTEGER(i8), POINTER     ::  kdp(:)     !< Model level corresponding to dep
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  tms(:)     !< Starting time for averaging
      REAL(r8),    POINTER     ::  tme(:)     !< Final time for averaging
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      INTEGER(i8), POINTER     ::  nav(:)     !< Number of time steps for averaging
      INTEGER(i8), POINTER     ::  kb(:)      !< k index of the nearest point below
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq5(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq6(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq7(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq8(:)     !< Interpolation parameter for masked grids
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE vdr_t

   TYPE (vdr_t)                :: vdr         !< Inizialization of drifter velocity observational vector

! ---
!> Observational vector for velocity from gliders
   TYPE gvl_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      INTEGER(i8), POINTER     ::  ino(:)     !< Profile number
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1 - u component, 2 - v component)
      REAL(r8),    POINTER     ::  lon(:)     !< Longitude
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  dpt(:)     !< Depth
      REAL(r8),    POINTER     ::  dzr(:,:)   !< Relative thickness
      INTEGER(i8), POINTER     ::  kdp(:)     !< Model level corresponding to dep
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  tms(:)     !< Starting time for averaging
      REAL(r8),    POINTER     ::  tme(:)     !< Final time for averaging
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      REAL(r8),    POINTER     ::  ins(:)     !< Saved increments
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  rss(:)     !< Saved misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      INTEGER(i8), POINTER     ::  nav(:)     !< Number of time steps for averaging
      INTEGER(i8), POINTER     ::  kb(:)      !< k index of the nearest point below
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq5(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq6(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq7(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq8(:)     !< Interpolation parameter for masked grids
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE gvl_t

   TYPE (gvl_t)                :: gvl         !< Inizialization of glider velocity observational vector

! ---
!> Observational vector for SST
   TYPE sst_t

      INTEGER(i8)              ::  no         !< Number of all observations
      INTEGER(i8)              ::  nc         !< Number of good observations
      INTEGER(i8)              ::  ns         !< Saved number of good observations
      REAL(r8)                 ::  dep        !< Minimum depth for observations
      INTEGER(i8)              ::  kdp        !< Model level corresponding to dep
      INTEGER(i8), POINTER     ::  ino(:)     !< Instrument
      INTEGER(i8), POINTER     ::  flg(:)     !< Quality flag
      INTEGER(i8), POINTER     ::  flc(:)     !< Temporary flag for multigrid
      INTEGER(i8), POINTER     ::  fls(:)     !< Saved temporary flag for multigrid
      REAL(r8),    POINTER     ::  lon(:)     !< Longitute
      REAL(r8),    POINTER     ::  lat(:)     !< Latitude
      REAL(r8),    POINTER     ::  tim(:)     !< Time
      REAL(r8),    POINTER     ::  val(:)     !< Observed value
      REAL(r8),    POINTER     ::  bac(:)     !< Background value
      REAL(r8),    POINTER     ::  inc(:)     !< Increments
      INTEGER(i8), POINTER     ::  par(:)     !< Parameter flag (1-temperature)
      REAL(r8),    POINTER     ::  bia(:)     !< Bias
      REAL(r8),    POINTER     ::  err(:)     !< Observational error
      REAL(r8),    POINTER     ::  bgerr(:)   !< Background error
      REAL(r8),    POINTER     ::  res(:)     !< Misfit
      REAL(r8),    POINTER     ::  b_a(:)     !< Background - analyses
      INTEGER(i8), POINTER     ::  ib(:)      !< i index of the nearest west point
      REAL(r8)   , POINTER     ::  pb(:)      !< distance from the nearest west point
      REAL(r8)   , POINTER     ::  rb(:)      !< distance from the nearest point below
      INTEGER(i8), POINTER     ::  jb(:)      !< j index of the nearest south point
      REAL(r8)   , POINTER     ::  qb(:)      !< distance from the nearest south point
      REAL(r8)   , POINTER     ::  pq1(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq2(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq3(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  pq4(:)     !< Interpolation parameter for masked grids
      REAL(r8)   , POINTER     ::  dpt(:)     !< Maximum depth of surrounding points
      REAL(r8)   , POINTER     ::  dtm(:)     !< Mean of the dynamic topography
      INTEGER(i8), POINTER     ::  kb(:)      !< Mean of the dynamic topography
      INTEGER(i8), POINTER     ::  eve(:)     !< Events

   END TYPE sst_t

   TYPE (sst_t)                :: sst         !< Inizialization of the SST observational vector

! ---
!> Error handling
   TYPE obserr_t

      LOGICAL                     :: rd_err_ff      !< Use error READ in observational file
      LOGICAL                     :: tim_dep_err    !< Apply temporal depENDency to the error
      REAL(r8)                    :: ztime_weigth   !< Temporal decay (in days) for obs error
      LOGICAL                     :: sla_sat_dep_err!< Apply           satellite  depENDency to the error
      LOGICAL                     :: sla_hor_dep_err!< Apply  horizontal spatial  depENDency to the error
      LOGICAL                     :: ts_ver_dep_err !< Apply  insitu vertical     depENDency to the error
      REAL(r8)                    :: sla_hor_dep_ecf!< Corrective factor for horizontal depENDent error
      INTEGER(i8)                 :: sla_sat_nu     !< Number of satellites
      CHARACTER(len=:),ALLOCATABLE:: sla_flname     !< File name of horizontal depENDent error for SLA
      CHARACTER(len=:),ALLOCATABLE:: ts_flname      !< File name of vertical depENDent error for TS
      CHARACTER*13,POINTER        :: sla_sat_name(:)!< Name   of satellites
      REAL(r8),    POINTER        :: sla_sat_err (:)!< Error associated to each satellite
      REAL(r8)                    :: sla_con_err    !< SLA         constant error
      REAL(r8)                    :: tem_con_err    !< Temperature constant error
      REAL(r8)                    :: sal_con_err    !< Salinity    constant error
      REAL(r8),    POINTER        :: stde(:,:)      !< SLA horizontal representation error
      REAL(r8),    POINTER        :: tem_fc(:,:)    !< Temperature horizontal correction factor
      REAL(r8),    POINTER        :: sal_fc(:,:)    !< Salinity    horizontal correction factor
      REAL(r8),    POINTER        :: tem(:)         !< Temperature vertical depENDent error
      REAL(r8),    POINTER        :: sal(:)         !< Salinity    vertical depENDent error
      REAL(r8),    POINTER        :: depth(:)       !< Depth of T/S error

   END TYPE obserr_t

   TYPE (obserr_t)                :: obserr         !< Inizialization of observational error vector

! ---
!> Huber Norm Quality control
   TYPE huberqc_t

      LOGICAL                     :: any    !< ANY of the follwing
      LOGICAL                     :: arg    !< Activate it for ARGO
      LOGICAL                     :: gld    !< Activate it for GLIDER
      LOGICAL                     :: gvl    !< Activate it for GLIDER VELOCITY
      LOGICAL                     :: sla    !< Activate it for SLA
      LOGICAL                     :: sst    !< Activate it for SST
      LOGICAL                     :: tra    !< Activate it for ARGO TRAJECTORY
      LOGICAL                     :: trd    !< Activate it for Suface DrIFter Trajectory
      LOGICAL                     :: vdr    !< Activate it for DrIFter Velocity
      LOGICAL                     :: xbt    !< Activate it for XBT
      LOGICAL                     :: asymm  !< Asymmetry of Huber Norm distribution
      LOGICAL                     :: L05    !< Activate L05 Norm
      CHARACTER(len=:),ALLOCATABLE:: flname !< File name for coefficients
      INTEGER                     :: iter   !< Iteration at which huber start

   END TYPE huberqc_t

   TYPE (huberqc_t)               :: huberqc !< Inizialization of the huber norm data structure

! ---
!> Quality check
   TYPE qck_t

      LOGICAL                     ::  res            !< Absolute misfits quality control
      LOGICAL                     ::  conbgr         !< Background/Observation ratio quality control (Constant Value)
      LOGICAL                     ::  eofbgr         !< Background/Observation ratio quality control (EOF)
      LOGICAL                     ::  vert           !< Reject measurements IF some observations above within the same profile were rejected
      LOGICAL                     ::  clm            !< Climatological   quality control
      REAL(r8)                    ::  res_sla        !< Max SLA          residual
      REAL(r8)                    ::  res_tem        !< Max Temperature  residual
      REAL(r8)                    ::  res_sal        !< Max Salinity     residual
      REAL(r8)                    ::  res_vel        !< Max Velocity     residual
      REAL(r8)                    ::  res_sst        !< Max SST          residual
      REAL(r8)                    ::  res_dis        !< Max distance     residual
      REAL(r8)                    ::  bgr_sla(2)     !< Background and Treshold SLA
      REAL(r8)                    ::  bgr_tem(2)     !< Background and Treshold Temperature
      REAL(r8)                    ::  bgr_sal(2)     !< Background and Treshold Salinity
      REAL(r8)                    ::  bgr_vel(2)     !< Background and Treshold Velocity
      REAL(r8)                    ::  bgr_sst(2)     !< Background and Treshold SST
      REAL(r8)                    ::  bgr_dis(2)     !< Background and Treshold distance
      REAL(r8),POINTER            ::  tem(:,:,:)     !< Background depENDent from i,j,k
      REAL(r8),POINTER            ::  sal(:,:,:)     !< Background depENDent from i,j,k
      REAL(r8),POINTER            ::  eta(:,:)       !< Background depENDent from i,j
      REAL(r8),POINTER            ::  climtem(:,:,:) !< Climatology for T
      REAL(r8),POINTER            ::  climsal(:,:,:) !< Climatology for S
      REAL(r8)                    ::  clm_lim(2)     !< Climatological Treshold ABS(clim-obs) for Temperature clm_lim(1) and Salimity clm_lim(2)
      CHARACTER(len=:),ALLOCATABLE::  flname         !< Climatology file name


   END TYPE qck_t

   TYPE(qck_t)                    :: qck    !< Inizialization of the quality check data structure

! ---
!> Huber Norm Quality control
   TYPE coastrej_t

      LOGICAL                  :: any       !< ANY of the follwing
      LOGICAL                  :: arg       !< Activate it for ARGO
      LOGICAL                  :: gld       !< Activate it for GLIDER
      LOGICAL                  :: gvl       !< Activate it for GLIDER VELOCITY
      LOGICAL                  :: sla       !< Activate it for SLA
      LOGICAL                  :: sst       !< Activate it for SST
      LOGICAL                  :: tra       !< Activate it for ARGO TRAJECTORY
      LOGICAL                  :: trd       !< Activate it for Suface DrIFter Trajectory
      LOGICAL                  :: vdr       !< Activate it for DrIFter Velocity
      LOGICAL                  :: xbt       !< Activate it for XBT
      REAL(r8)                 :: km_arg    !< km form the coast for ARGO
      REAL(r8)                 :: km_gld    !< km form the coast for GLIDER
      REAL(r8)                 :: km_gvl    !< km form the coast for GLIDER VELOCITY
      REAL(r8)                 :: km_sla    !< km form the coast for SLA
      REAL(r8)                 :: km_sst    !< km form the coast for SST
      REAL(r8)                 :: km_tra    !< km form the coast for ARGO TRAJECTORY
      REAL(r8)                 :: km_trd    !< km form the coast for Suface DrIFter Trajectory
      REAL(r8)                 :: km_vdr    !< km form the coast for DrIFter Velocity
      REAL(r8)                 :: km_xbt    !< km form the coast for XBT
      REAL(r8),POINTER         :: distc(:,:)!< km form the coast

   END TYPE coastrej_t

   TYPE (coastrej_t)           :: coastrej     !< Inizialization of the coastal rejection data structure

! ---
!> Thinning of observations
   TYPE thin_t

      LOGICAL                  :: any       !< ANY of the follwing
      LOGICAL                  :: arg       !< Activate it for ARGO
      LOGICAL                  :: gld       !< Activate it for GLIDER
      LOGICAL                  :: gvl       !< Activate it for GLIDER VELOCITY
      LOGICAL                  :: sla       !< Activate it for SLA
      LOGICAL                  :: sst       !< Activate it for SST
      LOGICAL                  :: tra       !< Activate it for ARGO TRAJECTORY
      LOGICAL                  :: trd       !< Activate it for Suface DrIFter Trajectory
      LOGICAL                  :: vdr       !< Activate it for DrIFter Velocity
      LOGICAL                  :: xbt       !< Activate it for XBT
      REAL(r8)                 :: tim       !< Time  winDOw for thinning
      REAL(r8)                 :: spc       !< Space winDOw for thinning

   END TYPE thin_t

   TYPE (thin_t)                :: thin     !< Inizialization of the thinning data structure

END MODULE obs_str
