#!/bin/bash
gfortran -O3 -fopenmp  \
 set_knd.F90 grd_str.F90 \
 read_nml.F90 read_grid.F90  read_corrad.F90 \
 writenc.F90  mk_weights.F90  \
-L/Users/marioadani/miniforge3/envs/MyEnv/lib \
-I/Users/marioadani/miniforge3/envs/MyEnv/include \
-lnetcdff -lnetcdf -o mk_weights.x

