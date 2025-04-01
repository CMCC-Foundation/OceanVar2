#!/usr/bin/env bash
#This file is part of the OceanVar testing suite
#File for scheduling jobs for automated testing of OceanVar
#Created 21 may 2024 12:39 by F. Carere (francesco.carere@cmcc.it)

set -u
set -e


######################
### Checking input ###
######################
function usage() {
   echo "-------"
   echo "Please pass the following arguments:"
   echo "1. name/reference of test which will be run"
   echo "2. name of setttings for test in a \"full_config/\" subdir"
   echo "3. name of namelist template in a \"full_config\" subdir"
   echo "4. name of executescript template in a \"full_config\" subdir"
   echo "5. name of subdirectory in \"full_config\" directory in which these files are"
   echo "6. List of output files which should be checked for bitwise equality"
   echo "7. Architecture (\"JUNO\" supercomputer or other)"
   echo "8. Amount of rerun for timing"
   #echo ""
   #echo "The following files in a \"full_config/\" subdir  may be edited:"
   #echo a file for the variables in the namelist that need to be changed across runs"
   #echo the same file for the MPI grid-decomposition that need to be changed across runs"
   #echo "\"OceanVar_nml_template\" for setting variables in the namelist that are constant across runs"
   #echo "\"OceanVar_nml_template\" for indicating which variables need to changed across runs"
   #echo Similar to"\"run_test_template.sh\" for the file passed to LSF using bsub"
   #SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
   #echo ""
   #full_configsdir=$(realpath $SCRIPT_DIR"/..")
   #echo "the \"full_configs\" dir should be located in $full_configsdir"
   exit 1
}

###############################
nrargs_ncsry=8
if [ "$#" -ne "$nrargs_ncsry" ]; then
   echo "Error passing arguments in $0" >&2
   echo "$nrargs_ncsry arguments should be passed, you passed $#" >&2
   echo "-------" >&2
   echo "Arguments passed:" >&2
   for((i=0;i<$#;i++)); do
      echo "arg $i: $(eval echo \$$i)" >&2
   done
   usage >&2
fi



#################################
########### SET DIRS ############
#################################
#Directory of current script which is executed
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#Directory to put tests in
TEST_DIR=$(realpath "$SCRIPT_DIR/../tests")

#################################
######### PARSE INPUT ###########
#################################
#directories
USR_DIR=$(realpath "$SCRIPT_DIR/../full_config/$5")
TESTNR=$1
OUTT_DIR="$TEST_DIR/$1"
#files and others
INPUTFILE="$USR_DIR/$2"
NMLFILE="$USR_DIR/$3"
RUNTEST="$USR_DIR/$4"
mapfile -d " " OUTPUTFILES <<< $6
ARCH="$7"
RERUN=$8



#Diagnostics script
metadatafile="$OUTT_DIR/metadata_$1.txt"

#####################################
######## READ INPUT FROM FILE ######
#####################################
echo ""
echo "-----------------"
echo "START READING INPUT"
echo "-----------------"
source "$SCRIPT_DIR/utils/Read_input.sh" "Start" $metadatafile $1
#Reset scriptdir, may be changed in the sourced script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#################################
##### DONE READING INPUT ########
#################################




#################################
####### OUTPUT: DO BSUBS ########
#################################

##Small script for debugging
#for ((j=0;j<${#STAT[@]};j++)); do
   #echo "------------"
   #echo "j = $j"
   #echo "STAT[$j] = ${STAT[$j]}"
   #echo "CAP[$j] = ${CAP[$j]}"
   #echo "NML_ID[$j] = ${NML_ID[$j]}"
   #for ((i=0;i<${#NML_NM[@]};i++)); do
   #   echo "NML_NM[$i] = ${NML_NM[i]}"
   #   curgrp=$( awk -v idt="${STAT[$j]}"  ' { print $idt }' <<< "${NML_NM[$j]}" )
   #   echo "curgrp = $curgrp"
   #done
#done

echo ""
echo "------------------------"
echo "START RUNNING JOBS"
echo "------------------------"
echo ""

source "$SCRIPT_DIR/utils/runtestjob.sh"
