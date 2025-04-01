#!/bin/bash

#BSUB -P 0603
#BSUB -q @Q@_medium
#BSUB -J @JOBNAME@
#BSUB -o %Jout
#BSUB -e @SCRIPTDIR@/../err/%Jerr
#BSUB -n @NPROC@
#BSUB -M 10G
@ADD_DEPENDENCY@


module purge
module load --auto impi-2021.6.0/2021.6.0 intel-2021.6.0/impi-2021.6.0/netcdf-cxx-threadsafe intel-2021.6.0/impi-2021.6.0/netcdf-fortran-threadsafe
#module load oneapi-2022.1.0/compiler-rt/2022.1.0 intel-2021.6.0/2021.6.0 impi-2021.6.0/2021.6.0 intel-2021.6.0/cdo-threadsafe intel-2021.6.0/impi-2021.6.0/netcdf-c-threadsafe
#intel-2021.6.0/impi-2021.6.0/parallel-netcdf  intel-2021.6.0/curl intel-2021.6.0/impi-2021.6.0/hdf5-threadsafe
#intel-2021.6.0/impi-2021.6.0/xios intel-2021.6.0/magics-threadsafe intel-2021.6.0/eccodes-threadsafe #intel-2021.6.0/cmake oneapi-2022.1.0/tbb/2021.6.0 oneapi-2022.1.0/mkl/2022.1.0

# ------------
# files
 OceanVar_exe="var_3d"
# OceanVar_nml='OceanVar_nml_seq'
 eofs_file="eofs.nc"
 grid_file='grid1.nc'


 errdir="@SCRIPTDIR@/../err"
 mkdir -p $errdir

echo "Running from directory $PWD"

# Link the  exe
 ln -vfs "@EXE@" $OceanVar_exe


# run
    mpirun -np @NPROC@ ./$OceanVar_exe

### Move profiling
#tempidx=""
#mv "profile_@NPROCX@_@NPROCY@.txt" "profile_@NPROCX@_@NPROCY@_$tempidx.txt" 2>dev/null
