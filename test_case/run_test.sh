#!/bin/bash
 exe_dir=${PWD}
# repository
# ------------
 var_dir='../bin/'
 var_nml_dir='../bin/'

# files
 var_exe='OceanVar'
 var_nml='OceanVar_nml'
 eofs_file="eofs.nc"
 grid_file='grid1.nc'

# Link the  exe
 [ -L  ${var_exe} ] && rm ${var_exe}
 ln -s ${var_dir}${var_exe} OceanVar

# run
    yyyymmdd=20181001
    mpi_irm=2
    mpi_jrm=2
    np=$(( mpi_jrm * mpi_irm))
   sed -e "s/@YYYYMMDD@/${yyyymmdd}/g"   \
       -e "s/@MPII@/${mpi_irm}/g" \
       -e "s/@MPIJ@/${mpi_jrm}/g" \
       ${var_nml_dir}${var_nml} > ${var_nml}

    mpirun -np $np ./OceanVar
