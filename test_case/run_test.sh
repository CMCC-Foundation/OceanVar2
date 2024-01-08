#!/bin/bash
 exe_dir=${PWD}
# repository
# ------------
 var_dir='../bin/'
 var_nml_dir='../bin/'

# files
 var_exe='var_3d'
 var_nml='var_3d_nml'
 eofs_file="eofs.nc"
 grid_file='grid1.nc'

# Link the  namelist
 [ -L ${var_nml} ] && rm ${var_nml}
 ln -s ${var_nml_dir}${var_nml} var_3d_nml

# Link the  exe
 [ -L  ${var_exe} ] && rm ${var_exe}
 ln -s ${var_dir}${var_exe} var_3d

# run
    mpirun -np 1 var_3d 
