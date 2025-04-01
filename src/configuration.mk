#!======================================================================
#!
#! This file is part of Oceanvar.
#!
#!  Copyright (C) 2025 OceanVar System Team ( oceanvar@cmcc.it )
#!
#! This program is free software: you can redistribute it and/or modify
#! it under the terms of the GNU General Public License as published by
#! the Free Software Foundation, either version 3 of the License, or
#! any later version (GPL-3.0-or-later).
#!
#! This program is distributed in the hope that it will be useful,
#! but WITHOUT ANY WARRANTY; without even the implied warranty of
#! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#! GNU General Public License for more details.
#!
#! You should have received a copy of the GNU General Public License
#! along with this program. If not, see <https://www.gnu.org/licenses/>.
#!======================================================================
#In this file, the user may specify variables needed for compiling
#Compiling should take a few minute


#Name of executable
EXE=OceanVar

#Set reproducibility to "FALSE" or "TRUE"
REPRO=FALSE

#Specify compiler (Fortran 90 and Fortran 77)
F90=mpiifort
F77=mpiifort

#Specify preprocessor options
P_P= 

#Specify compiler and linker options
F_O= -c -O2 -fpp 
F_L= $(F90) -O2 -fpp 

#Define external libraries (likely only NetCDF/NetCDF-fortran libraries. Checking the NetCDF command "nc-config --all" or "nf-config --all" command in your terminal you may find the paths which should be included (-I) and linked (-L) )
EXTINC=-Ipath_to_netcdf/include
EXTLIB=-Lpath_to_netcdf/lib

#Optinally, you may also use a predefined architecture (Zeus/Juno supercomputer, mac+homebrew, ipg compilter, ibm compiler)
#--WARNING--, uncommenting/setting this variable will override previously set variables in this document (except for 'EXE' and 'REPRO')
ARCHITECTURE=ZEUS_JUNO
