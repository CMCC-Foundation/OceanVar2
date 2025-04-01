#!/usr/bin/env bash
#This file is part of the OceanVar testing suite: this file is the start of a complete test:
#from bsubbing all the tests, to computing the outputs, to plotting the outputs.
#All is automatic and should be relatively flexible, although some hard-coding is present
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a high-level overview of the testing suite
set -e
set -u


##########################################
#########    ADD USER INPUT   ############
source userconfig.sh
##########################################























########### DO NOT CHANGE THE FOLLOWING CODE##############
######## IF YOU DO NOT KNOW WHAT YOU ARE DOING  ###########
if [ ${BASH_VERSINFO[0]} -le 3  ]; then
   echo "##############################################"
   echo "WARNING: The current bash version is ${BASH_VERSINFO[@]}"
   echo "Therefore this script may not run properly"
   echo "This script \($0\) should be run with bash version >=4"
   echo "##############################################"
fi


RERUN=1
#Names (PLEASE DO NOT REMOVE THE RANDOM PART OF THE NAMES. Necessary for LSF job dependency (option #BSUB -w)

TESTID=$RANDOM
TEST_NAME="$TEST_NAME"_"$TESTID"
TESTS_DIR="$SCRIPT_DIR/tests/$TEST_NAME"
while [ ! -d "$TESTS_DIR" ]; do
   TESTID=$RANDOM
   TEST_NAME="$TEST_NAME"_"$TESTID"
   TESTS_DIR="$SCRIPT_DIR/tests/$TEST_NAME"
done
RESULTS_NAME="$RESULTS_NAME"_"$TESTID"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BSUB_DIR="$SCRIPT_DIR/.scripts/utils/bsub"
RESULTS_DIR="$SCRIPT_DIR/results/$RESULTS_NAME"

echo "-----------------------------------"
echo "Starting test suite with id $TESTID"
#Checks if tests already done

mkdir -p "$TESTS_DIR"


#Set substitute values tests
SUB0="$SCRIPT_DIR"
SUB1="$TEST_NAME"
SUB2="$TESTINGSETTINGS_NAME"
SUB3="$NAMELIST_TEMPLATE_NAME"
SUB4="$RUNTEST_TEMPLATE_NAME"
SUB5="$FULL_CONFIG_SUBDIR"
SUB6="\"${FILES_OUT[*]}\""
SUB7="$ARCH"
SUB8="$RERUN"
SUBPOSTFIX="$TESTID"

##############################
#Set above user given input to the part of the testing that initializes up to and including bsub of all tests
bsub_start_file="$TESTS_DIR/rerun_start.sh"


echo "Preparing directories"
OUTDIR="$SCRIPT_DIR/out"
ERRDIR="$SCRIPT_DIR/err"
mkdir -p "$OUTDIR"
mkdir -p "$ERRDIR"



echo "-----------------------------------"
if [ "$ARCH" == "JUNO" ]; then       #Run tests on Juno supercomputer
   echo "submitting all testjobs (silently)"
   sed -e "s,@PATH@,$SUB0,g" "$BSUB_DIR/bsub_start.sh" | \
      sed -e 's,$1,'"$SUB1,g" | \
      sed -e "s,\$2,$SUB2,g" | \
      sed -e "s,\$3,$SUB3,g" | \
      sed -e "s,\$4,$SUB4,g"  | \
      sed -e "s,\$5,$SUB5,g" | \
      sed -e "s,\$6,$SUB6,g" | \
      sed -e "s,\$7,$SUB7,g" | \
      sed -e "s,\$8,$SUB8,g" | \
      sed -e "s,@POSTFIX@,$SUBPOSTFIX,g" > $bsub_start_file
      bsub < $bsub_start_file
      perl -i -pe '/.*BSUB.*-K.*/d' $bsub_start_file
else
   echo "Starting local tests"
   echo "$SCRIPT_DIR/.scripts/StartTestOceanVar.sh $SUB1 $SUB2 $SUB3 $SUB4 $SUB5 $SUB6 \"$ARCH\" $SUB8" > $bsub_start_file
   bash $bsub_start_file 1>"$OUTDIR/Startout$SUBPOSTFIX" 2>"$ERRDIR/Starterr$SUBPOSTFIX"
#else                                 #Choose architecture to run test on
#   echo "Please choose an available architecture"
#   exit 1
fi
echo "-----------------------------------"
echo "Done submitting jobs"
#echo "bsub $bsub_start_file"




##############################
#Set above user given input to the part of the testing that calculates and plots results of the tests
#if [ -d "$RESULTS_DIR" ]; then
#   echo "ERROR: Directory with results $RESULTS_DIR already exists." >&2
#   exit 1
#fi
mkdir -p "$RESULTS_DIR"

bsub_results_file="$RESULTS_DIR/rerun_results.sh"

DEPENDENCY="#BSUB -w done($TEST_NAME*)"

#Set substitute values results
SUB0="$SCRIPT_DIR"
SUB1="$TEST_NAME"
SUB2="$RESULTS_NAME"
SUB3="$PLOTS_NAME"
SUB4="$RESULTSSETTINGS_NAME"
SUB5="$FULL_CONFIG_SUBDIR"
SUB6="$TYPE_RESULTS"
SUB7="\"${FILES_OUT[*]}\""
SUB8="\"$CONDA_ENV\""
SUB9="\"$ARCH\""
SUB10=$RERUN
SUBDEP="$DEPENDENCY"
SUBPOSTFIX="$TESTID"

echo "-----------------------------------"
if [ "$ARCH" == "JUNO" ]; then #Run results on Juno supercomputer
   echo "Submitting job which calculates and plots results"
   echo "NOTE: this job starts only if all testjobs succesfully return (\"done\" state), otherwise infinitely pending"
   sed -e "s,@PATH@,$SUB0,g" "$BSUB_DIR/bsub_results.sh" | \
      sed -e "s,\$1,$SUB1,g" | \
      sed -e "s,\$2,$SUB2,g" | \
      sed -e "s,\$3,$SUB3,g" | \
      sed -e "s,\$4,$SUB4,g" | \
      sed -e "s,\$5,$SUB5,g" | \
      sed -e "s,\$6,$SUB6,g" | \
      sed -e "s,\$7,$SUB7,g" | \
      sed -e "s,\$8,$SUB8,g" | \
      sed -e "s,\$9,$SUB9,g" | \
      sed -e "s,\$(10),$SUB10,g" | \
      sed -e "s,@ADD_DEPENDENCY@,$SUBDEP,g" | \
      sed -e "s,@POSTFIX@,$SUBPOSTFIX,g" > $bsub_results_file
   bsub < $bsub_results_file
   perl -i -pe "/.*BSUB -w.*/d" "$bsub_results_file"
else #if [ "$ARCH" == "LOCAL_PC" ]; then
   echo "Starting to calculate results"
   echo "$SCRIPT_DIR""/.scripts/ResultsAfterStart.sh $SUB1 $SUB2 $SUB3 $SUB4 $SUB5 $SUB6 $SUB7 $SUB8 \"$ARCH\" $SUB10 ">"$bsub_results_file"
   bash $bsub_results_file 1>"$OUTDIR/Resultsout$SUBPOSTFIX" 2>"$ERRDIR/Resultserr$SUBPOSTFIX"
fi
