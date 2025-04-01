#!/bin/bash
OceanVar_exe="var_3d"
ln -vfs "@EXE@" $OceanVar_exe

mpirun -np @NPROC@ ./$OceanVar_exe\
   1>"@JOBNAME@out"\
   2>"@SCRIPTDIR@/../err/@JOBNAME@err"
