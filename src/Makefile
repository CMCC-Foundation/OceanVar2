SHELL = /bin/sh

############################################################################
#
#        Makefile for OceanVar 
#
############################################################################
#
#    Copyright 2006 Srdjan Dobricic, CMCC, Bologna
#
#    This file is part of OceanVar.
#
#    OceanVar is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OceanVar is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OceanVar.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

# IBM
F90 = mpiifort  #mpxlf90
F77 = mpiifort  -double-size 64 #mpxlf
P_P = #-WF,-Dopt_huge_memory
F_O = -c -O3 # -convert big_endian #-qreport -qmaxmem=-1 # -p # -C #-pg #-C #-p #-qsmp=omp -qreport #-p # -qnosmp -C
F_L = $(F90) -O3 # -qmaxmem=-1 # -p #-pg #-qsmp=omp #-p #-qnosmp

# PGI
#F90 = mpif90
#F77 = mpif90
#P_P = -Dopt_netcdf
#P_C = cpp -traditional-cpp
#F_O = -c -O2 #-Minfo -Mprof=func 
#F_L = $(F90) -O2 #-Mprof=func  


# IFORT
F90 = mpiifort
F77 = mpiifort  -double-size 64
P_P = -Dopt_netcdf
P_C = cpp -traditional-cpp
F_O = -c -O3 #-check bounds -check pointers #-Minfo -Mprof=func
F_L = $(F90) -O3 #-check bounds -check pointers #-Mprof=func

# MY-MACHINE
F90 = mpif90
F77 = mpif77 
P_P = 
P_C = cpp -traditional-cpp
#F_O = -c -O3 -std=legacy #-check bounds -check pointers #-Minfo -Mprof=func
#F_L = $(F90) -O3 #-check bounds -check pointers #-Mprof=func
F_O = -c -O3 -std=legacy #-check bounds -check pointers #-Minfo -Mprof=func
F_L = $(F90) -O3 #-check bounds -check pointers #-Mprof=func
L_O =

RM      = rm -f
MV      = mv -f

# my machine
EXTINC =  -I/usr/local/include 
EXTLIB =  -L/usr/local/lib -lnetcdff -lnetcdf

# ekman
#EXTINC =  -I/srv/lib/netcdf-last/include
#EXTLIB =  -L/srv/lib/netcdf-last/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lcurl

# athena
#EXTINC = -I/users/home/opt/netcdf/netcdf-3.6.3/include
#EXTLIB = -L/users/home/opt/netcdf/netcdf-3.6.3/lib -lnetcdf

EXEC = ../bin/var_3d

RM      = rm -f

KNDSTR  =  \
	set_knd.o
OBJSTR  =  \
	drv_str.o\
	cns_str.o\
	obs_str.o\
	grd_str.o\
	bmd_str.o\
	eof_str.o\
	mpi_str.o\
	ctl_str.o
OBJS    =  \
	routines.o\
	def_mpid.o\
	exc_mpi.o\
	def_nml.o\
	def_grd.o\
	sav_itr.o\
	ini_itr.o\
	ini_bmd.o\
	rdgrd.o\
	rdeofs.o\
	rdmxd.o\
	def_loc.o\
	netcdf_err.o\
	get_obs.o\
	get_obs_sla.o\
	get_obs_arg.o\
	get_obs_xbt.o\
	get_obs_gld.o\
	get_obs_vdr.o\
	get_obs_gvl.o\
	get_obs_tra.o\
	get_obs_trd.o\
	get_obs_sst.o\
	int_par.o\
	int_obs_hor.o\
	obs_vec.o\
	def_cov.o\
	def_coef.o\
	def_coef_bc.o\
	ini_cfn.o\
	ini_nrm.o\
	min_cfn.o\
	sav_msf.o\
	mod_msf.o\
	mod_inc.o\
	costf.o\
	cnv_ctv.o\
	ver_hor.o\
	get_byg.o\
	rcfl_1.o\
	rcfl_2.o\
	rcfl_x.o\
	rcfl_y.o\
	veof.o\
	bar_mod.o\
	get_vel.o\
	div_dmp.o\
	bar_mod_ad.o\
	invrt.o\
	invrt_ad.o\
	obsop.o\
	obs_sla.o\
	obs_arg.o\
	obs_xbt.o\
	obs_gld.o\
	obs_vdr.o\
	obs_gvl.o\
	obs_tra.o\
	obs_trd.o\
	obs_sst.o\
	resid.o\
	res_inc.o\
	obsop_ad.o\
	obs_sla_ad.o\
	obs_arg_ad.o\
	obs_xbt_ad.o\
	obs_gld_ad.o\
	obs_vdr_ad.o\
	obs_gvl_ad.o\
	obs_tra_ad.o\
	obs_trd_ad.o\
	obs_sst_ad.o\
	veof_ad.o\
	ver_hor_ad.o\
	rcfl_1_ad.o\
	rcfl_2_ad.o\
	rcfl_x_ad.o\
	rcfl_y_ad.o\
	get_vel_ad.o\
	div_dmp_ad.o\
	get_byg_ad.o\
	mod_trj_tl.o\
	mod_trj_ad.o\
	cnv_ctv_ad.o\
	cnv_inn.o\
	wrt_dia.o\
        oceanvar.o

.SUFFIXES: .F90 .F


all:  $(EXEC)
	@echo $(EXEC) is compiled

$(EXEC) :	$(KNDSTR) $(OBJSTR) $(OBJS) 
	$(F_L) $(L_X) -o $(EXEC) $(OBJSTR) $(OBJS) $(EXTLIB)

.DEFAULTS:

.F90.o :
	$(F90) $(P_P) $(F_O) $*.F90 $(EXTINC)

.F.o :
	$(F77) $(P_P) $(F_O) $*.F $(EXTINC)

clean:
	/bin/rm -f *.o *.mod $(EXEC)
