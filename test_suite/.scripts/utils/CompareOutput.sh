#!/usr/bin/env bash
##!/bin/bash
#This file is part of the OceanVartesting suite:
#It checks if the output files are bitwise equal, if not then for a netcdf file differences are calculated
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a High-level overview of the testing suite

set -u
set -e



##############################################


function usage() {
   echo "-------"
   echo "Please pass the following arguments:"
   echo "1: Directory of first testresults (to be compared with second)"
   echo "2: Directory of second testresults (to be compared with first)"
   echo "3: Directory where calculations should be outputted to"
   #echo "4: Full path (including name) of Python file makin the plots"
   echo "4: Value of bitwise equal netcdf files (e.g. if netcdf in single precision and bitwise equal then difference can be set to
   smaller than single precision, but perhaps not to zero (for logarithmic plots))"
   echo "5: How differences should be calculate (July 2024: \"sum\" for sum over grid of L1 differences, otherwise maximum over
   grid of L1 differences"
   echo "6: The list of outputfiles which should be checked, seperator should be the pipe '|' symbol"
   echo "7: Diagnostics file to which write partial output"
   echo "8: Flag intermediate, cdo-output NetCDF files (differences, maximum) should be saved. \
      Will be saved if \$8=\"slow\", otherwise will be skipped (performance- and memory-wise it matters)"
   echo "9: Flag in whic precision cdo should run ('32b' = 32 bit, otherwise 64 bit)"
   echo "10. Architecture (\"JUNO\" supercomputer or other)"
   echo ""
   exit 1
}


##############################################

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
######## Reading Input ##########
#################################
#ToCheck=("corr_eta.nc" "corr_sal.nc" "corr_tem.nc" "corr_uvl.nc" "corr_vvl.nc" "iterate.dat" "sla_stat.txt" "arg_stat.txt" "obs.dat")
IFS='|' read -r -a ToCheck <<< $6


P1=${1%/}
P2=${2%/}
DIFFERENT=()

FL_dgnstcs="$3/diagnostics.txt"
touch "$FL_dgnstcs"
echo -n "" > "$FL_dgnstcs"
echo "COMPARING OUTPUT OF DIRECTORIES" | tee -a  $FL_dgnstcs >> $7
echo $P1 | tee -a  $FL_dgnstcs >> $7
echo $P2 | tee -a  $FL_dgnstcs >> $7

TIME_CDO=0
TIME_NCDUMP=0
TIME_CALC=0

#################################
####### End Reading Input #######
#################################








##################################
####### Start Check different ####
##################################

echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "-CHECKING IF FILES ARE DIFFERENT-" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "" >> $7
for FILE in "${ToCheck[@]}"
do
   lcontinue=false
   FILE1=$P1/$FILE
   FILE2=$P2/$FILE

   #Check if files exist
   if ! test -f $FILE1 ; then
      echo "$FILE1 does not exist" >&2
      lcontinue=true
   fi
   if ! test -f $FILE2; then
      echo "$FILE2 does not exist" >&2
      lcontinue=true
   fi
   if $lcontinue; then
      #continue
      exit 1
   fi


   #Check for bitwise equality of files
   if [ $FILE == "iterate.dat" ]; then #if iterate.dat than remove first 20 line
      if cmp -s <( tail -n +20 $FILE1 ) <( tail -n +20 $FILE2 ) ; then
         echo "files $FILE are bitwise iequal from 20th line onwards" | tee -a  $FL_dgnstcs >> $7
      else
         echo "files $FILE are bitwise different from 20th line onwards" | tee -a  $FL_dgnstcs >> $7
         DIFFERENT+=($FILE)
      fi
   else #if not FILE=iterate.dat
      if cmp -s "$FILE1" "$FILE2"; then
         echo "files $FILE are bitwise equal" | tee -a  $FL_dgnstcs >> $7
      else
         echo "files $FILE are bitwise different" | tee -a  $FL_dgnstcs >> $7
         DIFFERENT+=($FILE)
      fi
   fi
done
##################################
####### END Check different ######
##################################



#########################################################
if [ "${10}" == "JUNO" ]; then
   module load --auto intel-2021.6.0/cdo-threadsafe/2.1.1-lyjsw > "/dev/null"
fi
########################################################




##########################################
####### Start calculate differences ######
##########################################
echo "" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "INFO ON NETCDF WHICH ARE BITWISE DIFFERENT" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "---------------------------------" | tee -a  $FL_dgnstcs >> $7
echo "" | tee -a  $FL_dgnstcs >> $7
echo "Maximum errors over 2D/3D grids are:" | tee -a  $FL_dgnstcs >> $7






#Value to bitwise equal (0) error with
#Since netcdf precision seems to be around 10^-24 it usually set to 10e-24
replace0=$4

#CDO option for precision
if [ "$9" == '32b' ]; then
   prc_opt='F32'
else
   prc_opt='F64'
fi

#CDO option for sum or max
if [ "$5" == 'sum' ]; then
   calc_opt='sum'
else
   calc_opt='max'
fi

#Check all the differences
FILE_res="$3/L1$5.txt"

##########################################
####### Start calculate differences ######
##########################################


touch "$FILE_res"
for FILE in "${ToCheck[@]}"
do
   if [[ ${DIFFERENT[@]} ==  *"$FILE"* ]]; then #if the file is different
      START_TIME_CALC=$SECONDS
      echo "File $FILE is different" | tee -a  $FL_dgnstcs >> $7
      FILE1=$P1/$FILE
      FILE2=$P2/$FILE

      if [[ $FILE == *".nc" ]] ; then #if netcdf file then compute the max/sum of differences
         #cvar = "eta" or "sal" or ...
         cvar="${FILE/corr_}"
         cvar="${cvar/.nc}"

         FILE_dif="$3/L1diff$FILE"
         FILE_typ="$3/L1$5$FILE"

         if [ "$8" == 'slow' ]; then
            if [ ! -f "$FILE_typ" ]; then
               echo "Making $FILE_typ" > "$FL_dgnstcs"

               #check if differences are already calculated
               if [ ! -f "$FILE_dif" ]; then
                  echo "Making $FILE_dif" >> "$FL_dgnstcs"
                  SET_TIME=$SECONDS
                  #Turn of 'u' flag of debugger since problematic with cdo
                  set +u
                  cdo -b $prc_opt sub "$FILE1" "$FILE2" "$FILE_dif"
                  set -u
                  TIME_CDO=$(( TIME_CDO + $SECONDS - SET_TIME ))
               fi

               #calculate the max or sum over grid of L1 distances
               SET_TIME=$SECONDS
               set +u
               cdo -b "$prc_opt" -vert"$calc_opt" -fld"$calc_opt" -abs "$FILE_dif" "$FILE_typ"
               set -u
               TIME_CDO=$(( TIME_CDO + $SECONDS - SET_TIME ))

               #Next line extracts the max/sum of differences from the netcdf file
               #TODO: Possibly the next line of code could be made a lot faster, there is a cdo operator which directly outputs. Can it be used?
               SET_TIME=$SECONDS
               set +u
               value_cvar=$(ncdump "$FILE_typ" | sed -z -e 's/.* '"$cvar"'.*=\(.*\)/\1/' | grep "\w" | sed -e "s/ //; s/;//")
               set -u
               TIME_NCDUMP=$(( TIME_NCDUMP + $SECONDS - SET_TIME ))
               TIME_CALC=$(( TIME_CALC + $SECONDS - START_TIME_CALC  ))

            fi
         else  #if $8 != slow
            SET_TIME=$SECONDS
            set +u
            value_cvar=$(cdo -b "$prc_opt" -output -vert"$calc_opt" -fld"$calc_opt" -abs -sub "$FILE1" "$FILE2")
            set -u
            TIME_CDO=$(( TIME_CDO + $SECONDS - SET_TIME ))
         fi


         if [ $value_cvar == 0 ]; then
            value_cvar=$replace0
         fi

         echo "$cvar: $value_cvar" | tee -a "$FILE_res" | tee -a $7 | tee -a  $FL_dgnstcs


      fi
   else
      if [[ $FILE == *".nc" ]] ; then  #if netcdf file not different then put value $replace0
         cvar="${FILE/corr_}"
         cvar="${cvar/.nc}"
         value_cvar=$replace0
         echo "$cvar: $replace0" | tee -a "$FILE_res" | tee -a $7 | tee -a  $FL_dgnstcs
      fi
   fi

done
#Output time
echo "Current time calculating: $TIME_CALC" | tee -a $7 | tee -a  $FL_dgnstcs
if [ $TIME_CALC -ne 0 ]; then
   echo "Time using cdo (calculations): $TIME_CDO. Relative time: $(bc <<< "scale=0; 100*$TIME_CDO/$TIME_CALC")%" |\
      tee -a $7 | tee -a $FL_dgnstcs
   echo "Time reading cdo's NetCDF output: $TIME_NCDUMP. Relative time: $(bc <<< "scale=0; 100*$TIME_NCDUMP/$TIME_CALC")%" |\
      tee -a $7 | tee -a  $FL_dgnstcs
fi
