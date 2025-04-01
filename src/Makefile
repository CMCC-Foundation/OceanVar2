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
SHELL = /bin/sh

############################################################################
#
#        Makefile for OceanVar
#
############################################################################

#Debug automatically set to zero
#0: no debug; [1-3]: different debug checks
DEBUG_CHKS=0

#--------------------  USER OPTIONS in "configuration.mk" -------------------------

#Read user input
-include configuration.mk

#Parse user input: reproducibility
ifndef REPRO
   REPRO=TRUE
endif
ifeq ($(REPRO),TRUE)
   KEYREPRO = -DREPRO

   #Parse intel compiler
   ifeq ($(F90), mpiifort)
      CO_REP= -fp-model consistent -fimf-use-svml=true
   endif
   ifeq ($(ARCH), ZEUS_JUNO)
      CO_REP= -fp-model consistent -fimf-use-svml=true
   endif
endif

#Set intel precompiler key for traceback call with intel compiler
ifeq ($(F90), mpiifort)
   P_P:=$(P_P) -DINTEL
else ifeq ($(ARCH), ZEUS_JUNO)
   P_P:=$(P_P) -DINTEL
endif

#-----------------------  COMPILER     -------------------------------

#Define default architecture
ifndef ARCHITECTURE
   ifndef F90
      ARCHITECTURE=MAC_PC
   else
      ARCHITECTURE='User defined'
   endif
endif

# IBM
ifeq ($(ARCHITECTURE),IBM)
   F90 = mpiifort  
   F77 = mpiifort  -double-size 64 
   P_P = 
   F_O = -c -O3 
   F_L = $(F90) -O3 
endif

# PGI
ifeq ($(ARCHITECTURE),PGI)
   F90 = mpif90
   F77 = mpif77
   P_P = -Dopt_netcdf
   F_O = -c -O2 
   F_L = $(F90) -O2
endif

# IFORT/ZEUS_JUNO
ifeq ($(ARCHITECTURE),ZEUS_JUNO)
   #KEYTRACEBACKQQ = -DINTEL
   F90 = mpiifort
   F77 = mpiifort -double-size 64
   P_P = -Dopt_netcdf -DINTEL
   F_O = -c -O3 -fpp $(CO_REP) -fp-model precise $(CO_OPTS) 
   F_L = $(F90) -O3 -fpp $(CO_REP) -fp-model precise $(CO_OPTS) 
endif

# MY-MACHINE/Mac_PersonalComputer
ifeq ($(ARCHITECTURE),MAC_PC)
   F90 = mpif90
   F77 = mpif77
   P_P =
   F_O = -c -O0 -std=legacy  -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow 
   F_L = $(F90) -O0 -cpp -finit-local-zero  -fcheck=all -ffpe-trap=invalid,zero,overflow 
endif

RM= rm -f
MV= mv -f


#-----------------------  COMPILER DEBUG FOR ZEUS_JUNO ARCHITECTURE  -----------------------------
### N.B. the debug versions (DEBUG_CHKS=) 1,2,3 are complementary (different checks), 4 turns on all
ifeq ($(F90),mpiifort) #Debugging currently done with interl compiler
   ifeq ($(DEBUG_CHKS),1)
      F_O        :=    $(F_O) -c -O0 -fpp -check bounds,contiguous,uninit,shape,stack,udio_iostat,pointers -traceback 
      F_L        :=    $(F90) $(F_L) -O0 -check bounds,contiguous,uninit,shape,stack,udio_iostat,pointers -traceback 
   else ifeq ($(DEBUG_CHKS),2)
      F_O        :=    $(F_O) -c -O0 -fpp  -traceback -check bound -ftrapuv 
      F_L        :=    $(F90) $(F_L) -O0  -traceback -check bound -ftrapuv 
   else ifeq ($(DEBUG_CHKS),3)
      F_O        :=    $(F_O) -c -O0 -fpp  -traceback -init=array,minus_huge 
      F_L        :=    $(F90) $(F_L) -O0  -traceback -init=array,minus_huge 
   else ifeq ($(DEBUG_CHKS),4) #Turn on all previous debug checks
      F_O1        =    $(F_O) -c -O0 -fpp  -check bounds,contiguous,uninit,shape,stack,udio_iostat,pointers -traceback 
      F_O2        =    $(F_O1) -check bound -ftrapuv 
      F_O        :=    $(F_O2) -init=array,minus_huge 
      F_L1        =    $(F90) -O0  -check bounds,contiguous,uninit,shape,stack,udio_iostat,pointers -traceback 
      F_L2        =    $(F_L1) -check bound -ftrapuv 
      F_L        =    $(F_L2) -init=array,minus_huge 
   endif
endif

F_O:=$(F_O)$(KEYREPRO)
F_L:=$(F_L)$(KEYREPRO)

#---------------------	NETCDF	-----------------------------

ifeq ($(ARCHITECTURE),ZEUS_JUNO)
   EXTINC =  -I${NETCDF}/include
   EXTLIB =  -L${NETCDF}/lib -lnetcdff -lnetcdf
else ifeq ($(ARCHITECTURE),MAC_PC)
# homebrew
#   EXTINC = -I/opt/homebrew/Cellar/netcdf/4.9.2_1/include/ -I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include/
#   EXTLIB = -L/opt/homebrew/Cellar/netcdf/4.9.2_1/lib -L/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib -lnetcdff -lnetcdf
# my old machine
#       EXTINC =  -I/usr/local/include
#       EXTLIB =  -L/usr/local/lib -lnetcdff -lnetcdf
# my new machine
   EXTINC =  -I/Users/marioadani/miniforge3/envs/MyEnv/include
   EXTLIB =  -L/Users/marioadani/miniforge3/envs/MyEnv/lib -lnetcdff -lnetcdf
endif

ifndef EXE
   EXE = OceanVar
endif
EXEC = ../bin/$(EXE)


#----------------------------------------------------------------------
KNDSTR  =  \
	set_knd.o
OBJSTR  =  \
	drv_str.o\
	cns_str.o\
	obs_str.o\
	grd_str.o\
	bmd_str.o\
	bal_str.o\
	eof_str.o\
	mpi_str.o\
	ctl_str.o\
	dfl_str.o\
	adjck_str.o\
	rpr_str.o
OBJS    =  \
        datetime_module.o\
	routines.o\
	def_mpid.o\
	exc_mpi.o\
	def_nml.o\
	def_grd.o\
	sav_itr.o\
        ini_time.o\
        ini_dflt.o\
        ini_huber.o\
	ini_itr.o\
	ini_bmd.o\
	ini_bal.o\
	seqmin.o\
	rdgrd.o\
	rdeofs.o\
	rdclim.o\
	rdwgh.o\
	rdcrl.o\
	rdeos.o\
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
	obserrors.o\
        obserr_arg.o\
        obserr_gld.o\
        obserr_gvl.o\
        obserr_sla.o\
        obserr_sst.o\
        obserr_tra.o\
        obserr_trd.o\
        obserr_vdr.o\
        obserr_xbt.o\
        huber_costf.o\
        get_bgerr.o\
        qualitycheck.o\
        qc_arg.o\
        qc_gld.o\
        qc_gvl.o\
        qc_sla.o\
        qc_sst.o\
        qc_tra.o\
        qc_trd.o\
        qc_vdr.o\
        qc_xbt.o\
        thinning.o\
        thin_arg.o\
        thin_gld.o\
        thin_gvl.o\
        thin_sla.o\
        thin_sst.o\
        thin_tra.o\
        thin_trd.o\
        thin_vdr.o\
        thin_xbt.o\
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
        ini_rcfl.o\
	recursive_filter.o\
	recursive_filter_ad.o\
	diffusive_filter.o\
	diffusive_filter_ad.o\
	rcfl_2.o\
	rcfl_x.o\
	rcfl_y.o\
	veof.o\
	bal_op.o\
	bal_op_ad.o\
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
	wrt_sla.o\
	wrt_arg.o\
	wrt_xbt.o\
	wrt_gld.o\
	wrt_sst.o\
	oceantools.o\
        adjck_balop.o\
        adjck_dfl.o\
        adjck_bmd.o\
        adjck_byg.o\
        oceanvar.o

.SUFFIXES: .F90 .F

.PHONY: help
help:
	@echo "-----------------------------------"
	@echo "The following targets are available:"
	@echo
	@echo "  help        Display this output"
	@echo "  clean       Clean up the \".o\", \".mod and binary/executable file"
	@echo "  all         Compile all of OceanVar"
	@echo "  \"file\".o  Compile a single file"
	@echo
	@echo "  Makefile options should be set in the \"configuration.mk\" file in \"/src\""
	@echo "-----------------------------------"

.PHONY: all
all:  $(EXEC)
	@echo $(EXEC) is compiled
	$(info ARCHITECTURE is $(ARCHITECTURE))
	$(info PREPROCESSOR is $(P_P))
	$(info REPRODUCIBILTY COMPILER OPTIONS ON? = $(REPRO))
	$(info PREPROCESSOR KEYS REPRODUCIBILITY  = $(KEYREPRO))
	$(info COMPILER OPTIONS REPRODUCIBILITY = $(CO_REP))

$(EXEC) :	$(KNDSTR) $(OBJSTR) $(OBJS)
	$(F_L) -o $(EXEC) $(OBJSTR) $(OBJS) $(EXTLIB)


.F90.o :
	$(F90) $(P_P) $(F_O) $*.F90 $(EXTINC)

.F.o :
	$(F77) $(P_P) $(F_O) $*.F $(EXTINC)

clean:
	/bin/rm -f *.o *.mod $(EXEC)
