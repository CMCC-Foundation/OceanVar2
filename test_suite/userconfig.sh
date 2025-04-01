#!/usr/bin/env bash
### CHOOSE ARCHITECTURE, may be 'local' (default) or 'JUNO' for Juno supercomputer
ARCH="local"  # Set to "JUNO" if testing on CMCC Juno supercomputer
#### In case local architecture used, make sure to have the prerequisites below installed
#### If 'JUNO' the environment modules will be loaded automatically (last checked: 15/09/2024)

#PREREQUISITES
#'local': (only tested on MacOS Ventura 13.*))
# Bash >=4.x
# A conda environment with python>=3.11.5
# CDO (climate data operators) installed in the Shell
# Graphviz
# python libraries: matplotlib, pandas


#### Test names #####
#####################################
FULL_CONFIG_SUBDIR="example" #Name of directory containing configuration
TEST_NAME="first_example" #Reference name testing runs (e.g. used for directory naming of tests)
TESTINGSETTINGS_NAME="tests_settings.txt" #Must be in the "full_config/$FULL_CONFIG_SUBDIR" directory, example in "full_config/example/test_settings.txt"
NAMELIST_TEMPLATE_NAME="OceanVar_nml_template"  #Must be in the "full_config/$FULL_CONFIG_SUBDIR" directory, example in "full_config/example/results_settings.txt"
RUNTEST_TEMPLATE_NAME="runfile_local.sh"  #Must be in the "full_config/FULL_CONFIG_SUBDIR" directory, example in "full_config/example/runfile_mac.sh"



#### GIVE INFO ABOUT RESULTS TO BE CALCULATED #####
###################################################
RESULTS_NAME="$TEST_NAME" #Reference name results directory
PLOTS_NAME="$TEST_NAME" #Reference name plot directory
RESULTSSETTINGS_NAME="results_settings.txt" #Must be in the "full_config/$FULL_CONFIGSUBDIR/" directory, example in "full_config/example/results_settings.txt"
CONDA_ENV="@YOUR_CONDA_ENVIRONMENT@" # Specify your conda environment which supports '*'-list unpacking (e.g. python version >=3.11.5)



#Output files to check for bitwise reproducibility, usually you would not change these
FILES_OUT=("corr_eta.nc" "corr_sal.nc" "corr_tem.nc" "corr_uvl.nc" "corr_vvl.nc" "iterate.dat" "sla_stat.txt" "arg_stat.txt" "obs.dat")

### EXTRA OPTIONS
TYPE_RESULTS="max" # This can be changed to "mean", "sum", etc., for different types of results (computing max, mean, sum of L1 differences over grid)
#TODO:
#Also calculate mean over grid
#PRECISION='64' #Precision of "cdo" (climate data operators) operations on NetCDF "output files to check"/"FILES_OUT", current default is 64
#CALCULATIONS='FAST' #Slow/fast calculations of cdo. If slow then the intermediate NetCDF files are saved (memory-heavy and slower).
#  current default is slow calculations
