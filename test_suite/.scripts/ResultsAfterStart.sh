#!/usr/bin/env bash
#This file is part of the OceanVar testing suite
#File for creating diagnostics of tested jobs
#Created 29 may 2024 11:01 by F. Carere (francesco.carere@cmcc.it)

set -u
set -e


###############################


function usage() {
   echo "-------"
   echo "Please pass the following arguments:"
   echo "1. name/reference of test which have been run (in tests/)"
   echo "2. name/reference of test diagnostics (out to results/\$2)"
   echo "3. name/reference of plots (out to plots/\$3)"
   echo "4. name of file in a \"full_config/\" subdir containing what should be calculated/plotted"
   echo "5. name of subdirectory in \"full_config\" directory in which the above file is"
   echo "6. If the tests computed should compute the max or the sum over L1 grid differences"
   echo "7. List of output files which should be checked for bitwise equality"
   echo "8. Conda environment with python version >3.12"
   echo "9. Architecture (\"JUNO\" supercomputer or other)"
   echo "10. Nr of reruns (bigger than 1 only useful if you add timers to it)"
   echo ""
   echo "The following files in the \"full_config/\" dir  may be edited:"
   echo "\"inputfile.txt\" for the variables in the namelist that need to be changed across runs"
   echo "\"inputfile.txt\" for the MPI grid-decomposition that need to be changed across runs"
   SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
   echo ""
   full_configdir=$(realpath $SCRIPT_DIR"/..")
   echo "the \"full_config\" dir should be located in $full_configdir"
   exit 1
}


function error_print {
   echo "" >&2
   echo "" >&2
   echo "-------------------------" >&2
   echo "-------------------------" >&2
   echo "ERROR ENCOUNTERED IN $0" >&2
   echo "-------------------------" >&2
   echo "-------------------------" >&2
}

##########################################

nrargs_ncsry=10

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
######## SETUP and  DIRS ########
#################################
TESTNR="$1"
RERUN="${10}"

#Set general dirs
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
USR_DIR=$(realpath "$SCRIPT_DIR/../full_config")"/$5"
INPUTFILE="$USR_DIR/$4"
TESTS_DIR=$( realpath $SCRIPT_DIR"/../tests")"/$1"
RESULTS_DIR=$( realpath $SCRIPT_DIR"/../results")"/$2"
mkdir -p -v $RESULTS_DIR

#Set dirs for plotting
mkdir -p $SCRIPT_DIR"/../plots/$2"
PLOTS_DIR=$( realpath $SCRIPT_DIR"/../plots")"/$2/$3"
PLOTS_DIR_DIFF="$PLOTS_DIR/Diff_Xaxes"
PLOTS_DIR_CORR="$PLOTS_DIR/NC_Corr_Xaxes"
mkdir -p "$PLOTS_DIR_DIFF" "$PLOTS_DIR_CORR"
cp "$INPUTFILE" "$PLOTS_DIR/"

#Lists for later
RESULT_DIRS=()
PROCS_DIRS=()

#Read which files are OceanVar output for which differences need to be calculated
mapfile -d " " OUTPUTFILES <<< $7
OUTPUTFILES_2=$(IFS='|'; echo "${OUTPUTFILES[*]}" | sed -e "s/ //g")


#Check if python file exists and is complete
PYTHON_FILE=$(realpath "$SCRIPT_DIR/utils")"/gen_plots.py"
PYTHON_FILE_NEW="$PLOTS_DIR/gen_plots.py"
PFN_exists=true

if [[ ! -f $PYTHON_FILE_NEW ]]; then
   PFN_exists=false
elif grep -q "\@LIST_HERE\@" "$PYTHON_FILE_NEW"; then
   PFN_exists=false
fi

#If does not exist, setup new one
if ! $PFN_exists ; then

   #echo "copying Python file"
   cp -v $PYTHON_FILE $PYTHON_FILE_NEW

   #Set errortype="max" or "sum
   perl -i -pe "s/\@ERRTYPE\@/\"$6\"/g" $PYTHON_FILE_NEW
   for FILE in ${OUTPUTFILES[@]}; do
      if [[ $FILE == *".nc" ]] ; then
         cvar=${FILE/corr_}
         cvar=${cvar/.nc}
         perl -i -pe "s/\@file\@/"\'$cvar\'",\@file@/g" "$PYTHON_FILE_NEW"
      fi
   done
   perl -i -pe "s/,\@file\@//g" "$PYTHON_FILE_NEW"

fi



#####################################
######## READ INPUT FROM FILE ######
#####################################
echo ""
echo "------------------------"
echo "------------------------"
echo "----- READING INPUT -----"
echo "------------------------"
echo "------------------------"
source "$SCRIPT_DIR/utils/Read_input.sh" "after_test" $4 $PYTHON_FILE
#Reset SCRIPT_DIR, may have been changed in the sourced script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )




echo ""
echo "------------------------"
echo "------------------------"
echo "--- DONE READING INPUT ---"
echo "------------------------"
echo "------------------------"

#################################
##### DONE READING INPUT ########
#################################
#################################
##### START CALCULATIONS ########
#################################
echo ""
echo "------------------------"
echo "------------------------"
echo "--- DOING CALCULATIONS ---"
echo "------------------------"
echo "------------------------"


START_TIME=$SECONDS
if  ! $PFN_exists ; then #only run if the python file has not been finished yet

   curcalc=0
   replace0='10e-24'

   FILE_dgnstcs="$PLOTS_DIR/diagnostics$6.txt"
   echo "Calculating and writing metadata to $FILE_dgnstcs"

   rm -f $FILE_dgnstcs
   touch "$FILE_dgnstcs"


   for ((i=0; i<${#brcktExp_names[@]};i++)); do #loop over fixed subspaces

      read -r -a varsubspace <<<$(eval echo ${brcktExp_names[$i]})
      read -r -a vardirnm    <<<$(eval echo ${dir_nm[$i]})
      totalcalc=$(( ${#brcktExp_names[@]}*(${#varsubspace[@]} - 1)*( ${#varsubspace[@]} )/2 ))

      for ((j=0;j<${#varsubspace[@]};j++)); do #loop over variables subspaces

         for ((k=0;k<${#varsubspace[@]};k++)); do

            #If j=k then difference is zero, because they are the same points
            if [ $j == $k ]; then
               perl -i -pe "s/\@LIST_HERE\@/[$replace0,$replace0,$replace0,$replace0,$replace0], \#diag${vardirnm[$k]}\n \@LIST_HERE\@/g" $PYTHON_FILE_NEW
               continue
            fi

            if [ $j -gt $k ]; then
               curcalc=$((curcalc+1))
            fi


            #Check if directories exist, otherwise exit
            cdir1="$TESTS_DIR/${varsubspace[$k]}"
            vdirnm1=${vardirnm[$k]}
            cdir2="$TESTS_DIR/${varsubspace[$j]}"
            vdirnm2=${vardirnm[$j]}

            if ! test -d $cdir1 ; then
               error_print
               echo "$cdir1 not found" >&2
               echo "Please make sure that the test results are present" >&2
               exit 1
            fi
            if !  test -d $cdir2 ; then
               error_print
               echo "$cdir1 not found" >&2
               echo "Please make sure that the test results are present" >&2
               exit 1
            fi


            if [ $j -ge $k ]; then
               cdirnm="$vdirnm1"_V_"$vdirnm2"
            else
               cdirnm="$vdirnm2"_V_"$vdirnm1"
            fi
            RESULTS_CURDIR="$RESULTS_DIR/$cdirnm"
            mkdir -p $RESULTS_CURDIR

            FILE_res="$RESULTS_CURDIR/L1$6.txt"
            curprct=$(bc <<< "scale=0; 100*$curcalc/$totalcalc")



            ########################################
            ######## STARTING CALCULATIONS #########
            ########################################
            if [ -f "$FILE_res" ]; then  #If file containing differences already exists, don't recalculate
               if [ $j -ge $k ]; then
                  echo "Skipping ${varsubspace[$k]} vs ${varsubspace[$j]}" | tee -a "$FILE_dgnstcs"
                  echo "File L1$6 already calculated" | tee -a "$FILE_dgnstcs"
               fi

            else #else calculate the differences
               echo "Calculating ${varsubspace[$k]} vs ${varsubspace[$j]}" | tee -a "$FILE_dgnstcs"

               #CALCULATE THE DIFFERENCES
               SET_TIME=$SECONDS
               bash "$SCRIPT_DIR/utils/CompareOutput.sh" "$cdir1" "$cdir2" "$RESULTS_CURDIR" "$replace0" "$6"\
                  "$OUTPUTFILES_2" "$FILE_dgnstcs" 'fast' '64b' "$9"
               CALC_TIME=$(( $SECONDS - SET_TIME ))
               echo "Check took time = $CALC_TIME" | tee -a "$FILE_dgnstcs"

            fi
            ########################################
            ######## FINISHED CALCULATIONS #########
            ########################################





            #Get differences from file
            mapfile -t values_cvar < <( awk '{ print $2 }' "$FILE_res" )
            ELAPSED_TIME=$(($SECONDS - $START_TIME))
            perl -i -pe "s/\@LIST_HERE\@/[$(IFS=,; echo "${values_cvar[*]}" )], \#$cdirnm\n    \@LIST_HERE\@/g" "$PYTHON_FILE_NEW"

            if [ $j -ge $k ]; then
               echo "with values ${values_cvar[@]}" | tee -a "$FILE_dgnstcs"

               echo "at $curprct % of calculations" | tee -a "$FILE_dgnstcs"

               echo "Time running calculations = $ELAPSED_TIME" | tee -a "$FILE_dgnstcs"
               echo "---------------" | tee -a $FILE_dgnstcs
            fi


         done
      done
   done
   perl -i -pe 's/.*\@LIST_HERE\@.*/]/g' "$PYTHON_FILE_NEW"
else
   if grep -q "\@file\@" "$PYTHON_FILE_NEW"; then
      for FILE in ${OUTPUTFILES[@]}; do
         if [[ $FILE == *".nc" ]] ; then
            cvar=${FILE/corr_}
            cvar=${cvar/.nc}
            perl -i -pe "s/\@file\@/"\'$cvar\'",\@file\@/g" "$PYTHON_FILE_NEW"
         fi
      done
      perl -i -pe "s/,\@file\@//g" "$PYTHON_FILE_NEW"
   fi
   echo "---"
   echo "!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!"
   echo "ALL results seem to already have been calculated, skipping all calculations"
   echo "-------------- PLEASE READ THE FOLLOWING LINE ------------"
   echo "If calculations need to be redone, please (re)move the file $PYTHON_FILE_NEW and some of the files \"L1$6.txt\" in the \"results/$2/\" subdirectories which need to be recalculated"
   echo "!!!!!!!!!!!!!!!!END OF WARNING!!!!!!!!!!!!!!!!!!!"
   echo "---"
fi
if [ $RERUN -gt 1 ]; then
   #Timing.sh should not be in the released version
   source "$SCRIPT_DIR/utils/timing/set_timing.sh"
fi

ELAPSED_TIME=$(($SECONDS - $START_TIME))
echo ""
echo "------------------------"
echo "------------------------"
echo "--- DONE CALCULATIONS ---"
echo "------------------------"
echo "------------------------"
echo "Total time running calculations = $ELAPSED_TIME"
#################################
####### END CALCULATIONS ########
#################################


#################################
######## START PLOTTING #########
#################################


#if  ! $PFN_exists ; then
#   cp "$FILE_dgnstcs" "$PLOTS_DIR/diagnostics$6.txt"
#fi

#################################
#Run python file to get the plots
set +u
set +e
echo ""
echo "--------------------"
#echo "Running Python script to make plots"
if [ "$9" == "JUNO" ]; then
   source /juno/opt/anaconda/3-2022.10/etc/profile.d/conda.sh
fi
conda activate $8
echo "Initializing python script..."
python3 "$PYTHON_FILE_NEW" "$PLOTS_DIR_DIFF" "$PLOTS_DIR_CORR"
