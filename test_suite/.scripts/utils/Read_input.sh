#!/usr/bin/env bash
#This file is part of the OceanVar testing suite
#This file parses the settings files
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a High-level overview of the testing suite
set -u
set -e



function usage() {
   echo "-------"
   echo "Please pass the following arguments:"
   echo "1: Mode of reading (July 2024: either \"Start\" which is the first read, or not \"START\""
   echo "2: Full path (including name) of file to be considered"
   echo "3: TWO CASES:"
   echo "In case of \$1 == \'Start\' 3.1: name/reference of test which will be run"
   echo "In case of \$1 != \'Start\' 3.2: full path (including name) of python file TEMPLATE to compute plots"
   #Note 3.2 actually unused at the moment
   echo ""
   exit 1
}

##########################################

nrargs_ncsry=3

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


#Remove trailing whitespace from inputfile and copying to test directory
#sed -i 's/[ \t]*$//' "$INPUTFILE"
echo ""
echo "copying inputfile"
if [ $1 == 'Start' ]; then
   cp -v $INPUTFILE $TEST_DIR/$3
else
   cp -v $INPUTFILE $RESULTS_DIR/$2
fi




################ START: READ INPUT FILE ###############


declare -a TOT_var=()
declare -a TOT_fxd=()


#Explanation of variables:
#  "*_var": if differences should be calculated (i.e. in $RESULTSSETTINGS_NAME, variable set to "var").
#           For example, if two executables are given and one is interested in the differences of the
#           netcdf files between runs with these different executables, then EXE_var='var'.
#           If one is not interested in these differences then it is not equal to 'var' (but usually
#           'fixed')
# "*_NM" :  The reference names used for building the directories and other naming
# "*_PT" :  Path where executables or input files can be found
# "NML_ID": For namelist used identify where in the namelist a value should be substituted
#           In particular, "@$NML_ID@" will be substituted by the corresponding value in "NML_SB"
# "NML_SB": Contains values that are to be substituted in the namelist at the corresponding "@NML_ID@"
#           position
# "PROC_{X,Y}: About the X or Y decomposition of the structured grid in OceanVar
# "PROC_XY: Total nr of processes


#######################################
################ EXE ##################
#######################################

nexe=$(grep "nexe.*=" "$INPUTFILE" |  cut -d '=' -f2 ) #read number of executables considered in testing
line_nexe=$(grep -n "nexe.*=" "$INPUTFILE" | cut -d : -f 1) #orientation in file, get linenumber of previously read line

if [ $1 != 'Start' ]; then #read executable and check if "fixed" or "var"
   EXE_var=$( sed -n -e $(( line_nexe + 1 ))p "$INPUTFILE" )
   line_nexe=$(( line_nexe+1 ))
   [[ $EXE_var == 'var' ]] && TOT_var+=($nexe) || TOT_fxd+=($nexe)
fi

mapfile -t EXE_NM < <( sed -n -e $(( line_nexe + 1 )),$((line_nexe+nexe ))p "$INPUTFILE" | awk '{print $1}')
if [ $1 == 'Start' ]; then
   mapfile -t EXE_PT < <( sed -n -e $(( line_nexe + 1 )),$((line_nexe+nexe ))p "$INPUTFILE" | awk '{print $2}')
   declare -a EXE_RPT=()
   for ((i=0;i<${#EXE_PT[@]};i++));do
      if [[ ${EXE_PT[$i]} != /* ]]; then
         EXE_RPT+=( "$(realpath $PWD/${EXE_PT[$i]})" )
      else
         EXE_RPT+=( ${EXE_PT[$i]} )
      fi
   done
fi




#######################################
############### FILES #################
#######################################

ninp=$(grep "ninput.*=" "$INPUTFILE" |  cut -d '=' -f2 )
line_ninp=$(grep -n "ninput.*=" "$INPUTFILE" | cut -d : -f 1)

if [ $1 != 'Start' ]; then
   INP_var=$( sed -n -e $(( line_ninp + 1 ))p "$INPUTFILE" )
   line_ninp=$(( line_ninp+1 ))
   [[ $INP_var == 'var' ]] && TOT_var+=($ninp) || TOT_fxd+=($ninp)

fi

mapfile -t INP_NM < <( sed -n -e $(( line_ninp + 1 )),$((line_ninp+ninp ))p "$INPUTFILE" | awk '{print $1}')
if [ $1 == 'Start' ]; then
   mapfile -t INP_PT < <( sed -n -e $(( line_ninp + 1 )),$((line_ninp+ninp ))p "$INPUTFILE" | awk '{print $2}')
   declare -a INP_RPT=()
   for ((i=0;i<${#INP_PT[@]};i++));do
      if [[ ${INP_PT[$i]} != /* ]]; then
         INP_RPT+=( "$(realpath $PWD/${INP_PT[$i]})" )
      else
         INP_RPT+=( ${INP_PT[$i]} )
      fi
   done

fi


#######################################
############# NAMELIST ################
#######################################
read -r -a nvar < <(grep "varpergroup.*=" "$INPUTFILE" | cut -d' ' -f2-)

line_nml=$(grep -n "varper.*=" "$INPUTFILE" | cut -d : -f 1)
line_nml=$((line_nml+1))

CAP=()            #CAP counts the number of namelist configurations for each variable e.g. 3 if testing 3 filters
STAT=()
cidx=1            #cidx keeps track of which linenumbers should be read next
varidx=0          #varidx helps count for CAP
declare -a NML_ID=()
declare -a NML_ID2=()
declare -a NML_DR=()
declare -a NML_NM=()
declare -a NML_SB=()
declare -a NML_var=()

#echo "nvar = ${nvar[@]}"

j=0
#echo "nvar = ${nvar[@]}"
for i in "${nvar[@]}"; do
   IFS=' ' read -r temparr <<< "$( sed -n -e $(( line_nml + cidx  ))p "$INPUTFILE" )"
   if [ $1 == 'Start' ]; then
      NML_ID+=(${temparr[@]})
      NML_DR+=( $( awk '{print $1}' <<< ${temparr[0]} ) )
      c_varname=$( sed -n -e $(( line_nml + cidx+1  ))p "$INPUTFILE" )
   else
      NML_ID+=( $( awk '{print $1}' <<< ${temparr[0]} ) )
      NML_var+=( $( sed -n -e $(( line_nml + cidx+1 ))p "$INPUTFILE" ) )
      c_varname=$( sed -n -e $(( line_nml + cidx+2  ))p "$INPUTFILE" )
   fi
   NML_NM+=("$c_varname")
   #echo "NML_NM = ${NML_NM[@]}"
   #echo "NML_DR= ${NML_DR[@]}"

   varidx=$( wc -w <<<"$c_varname")
   #echo "$varidx = varidx"
   CAP+=( $varidx )
   STAT+=( 1 )
   #echo "CAP = ${CAP[@]}"
   #ccapi=$(( $ccapi + $varidx ))
   #CAPP+=($ccapi)

   read -r -a check_names <<<  $c_varname
   if [ $1 == 'Start' ]; then
      for ((j=1;j<=i;j++));do
         c_varnamesgrp=$( sed -n -e $(( line_nml + cidx+1+j  ))p "$INPUTFILE" )

         #### Start small test to check if input ok ######
         IFS='|' read -r -a check_grp <<<  $c_varnamesgrp
         if [ $varidx != ${#check_grp[@]} ]; then
            echo "-----------------">&2
            echo "ERROR in $0">&2
            echo "------------">&2
            #echo "FOR VARIABLE $c_varname"
            echo "Substitute argument for ID = ${NML_ID[-1]}, and name = ${check_names[$j]} equal to:">&2
            echo "$c_varnamesgrp">&2
            echo "Thus, number of substitute arguments given is ${#check_grp[@]}">&2
            echo "According to the amount of names given: \"$c_varname\" it should be $varidx">&2
            echo "------------">&2
            echo "Please check your TESTSETTINGS file">&2
            echo "END OF ERROR">&2
            exit 1
         fi
         #### END small test to check if input ok" ######

         NML_SB+=("$c_varnamesgrp")
      done
   else
      [[ $NML_var[$j] == 'var' ]] && TOT_var+=(${CAP[$j]}) || TOT_fxd+=(${CAP[$j]})
      j=$((j+1))
   fi

   if [ $1 == 'Start' ]; then
      cidx=$(( cidx+3+i ))
   else
      cidx=$(( cidx+3+1 ))
   fi
done


#######################################
################ MPI ##################
#######################################
nproc=$(grep "nproc.*="  "$INPUTFILE" | awk '{print $2}')

line_nprc=$(grep -n "nproc.*=" "$INPUTFILE" | cut -d : -f 1)
line_nprc=$((line_nprc+1))

if [ $1 != 'Start' ]; then
   PRC_var=$( sed -n -e $(( line_nprc ))p "$INPUTFILE" )
   line_nprc=$((line_nprc+1))

   [[ $PRC_var == 'var' ]] && TOT_var+=($nproc) || TOT_fxd+=($nproc)
fi

declare -a PROCX=()
declare -a PROCY=()
declare -a PROCXY=()
for ((i=0;i<nproc;i++)); do
   #curprc=$(sed -n -e $(( line_nprc + i  ))p "$INPUTFILE" )
   #   echo "curprc="$curprc
   curpx=$(sed -n -e $(( line_nprc + i ))p "$INPUTFILE" | awk -F x '{ print $1 }' )
   curpy=$(sed -n -e $(( line_nprc + i ))p "$INPUTFILE" | awk -F x '{ print $2 }' )
   PROCX+=("$curpx")
   PROCY+=("$curpy")
   PROCXY+=("$curpx"x"$curpy")
done


################ END READ IN PUT FILE ###############


################ START PREP FOR PYTHON_FILE ###############
#Mainly for the conversion of all names to strings

#echo "New python file = $PYTHON_FILE_NEW"
if [ $1 != 'Start' ]; then
   declare -a EXE_NM_str=()
   declare -a INP_NM_str=()
   declare -a NML_NM_str=()
   declare -a PROCXY_str=()

   for line in "${EXE_NM[@]}"; do
      EXE_NM_str+=("\"$line"\")
   done
   for line in "${INP_NM[@]}"; do
      INP_NM_str+=("\"$line"\")
   done
   for ((j=0;j<${#CAP[@]};j++)); do
      #echo "NML_NM[$j] = ${NML_NM[$j]}"
      read -r -a NML_nm_str_cur <<< "${NML_NM[$j]}"
      #echo "NML_nm_str_cur = ${NML_nm_str_cur[@]}"
      for ((i=0;i<${#NML_nm_str_cur[@]};i++)); do
         NML_nm_str_cur[$i]="\"${NML_nm_str_cur[$i]}\""
      done
      NML_NM_str+=("${NML_nm_str_cur[*]}")
   done
   #for ((j=0;j<${#CAP[@]};j++)); do
   #   #echo "NML_NM_str[$j] = ${NML_NM_str[$j]}"
   #done
   for line in "${PROCXY[@]}"; do
      PROCXY_str+=("\"$line"\")
   done

   ################ END PREP FOR PYTHON_FILE ###############



   ################ START SETUP PYTHON FILE ###############

   #echo $(IFS=,; printf "${EXE_NM[*]}" >>> mapfileq   )
   #perl -i -pe "s/@IDS@/["$(IFS=,; printf "%q " "${EXE_NM[@]}" )"], @IDS@/g" "$PYTHON_FILE_NEW"
   #exit 0
   perl -i -pe "s/\@IDS\@/["$(IFS=,; echo "${EXE_NM_str[*]}" )"], \@IDS\@/g" "$PYTHON_FILE_NEW"
   perl -i -pe "s/\@IDS\@/["$(IFS=,; echo "${INP_NM_str[*]}" )"], \@IDS\@/g" "$PYTHON_FILE_NEW"
   perl -i -pe "s/\@TYPES\@/\"$EXE_var\", \@TYPES\@/g" "$PYTHON_FILE_NEW"
   perl -i -pe "s/\@TYPES\@/\"$INP_var\", \@TYPES\@/g" "$PYTHON_FILE_NEW"
   for ((j=0;j<${#CAP[@]};j++)); do
      perl -i -pe "s/\@IDS\@/["$(sed -e 's/ /,/g' <<<"${NML_NM_str[$j]}" )"], \@IDS\@/g" "$PYTHON_FILE_NEW"
      perl -i -pe "s/\@TYPES\@/\"${NML_var[$j]}\", \@TYPES\@/g" "$PYTHON_FILE_NEW"
   done
   perl -i -pe "s/\@IDS\@/["$(IFS=,; echo "${PROCXY_str[*]}" )"]/g" "$PYTHON_FILE_NEW"
   perl -i -pe "s/\@TYPES\@/\"$PRC_var\"/g" "$PYTHON_FILE_NEW"
   line_inpIds=$(grep -m 1 -n "inp_Ids=" "$PYTHON_FILE_NEW" | cut -d : -f 1)
   perl -i -pe "$line_inpIds""s/\w.*,\|]/\"&\"/g" $PYTHON_FILE_NEW
   #perl -i -pe "$line_inpIds""s/\w.*,\|]/\"&\"/g" $PYTHON_FILE_NEW

   ################ END SETUP PYTHON FILE ###############






   ################ START: BRACKET EXP VARIABLE  ##############
   # Bracket expansion is used to get the right loop iterations

   ################ EXE ##################
   if [[ ${#EXE_NM[@]} -gt 1 ]]; then
      [[ $EXE_var == var ]] &&\
         brcktExp_names="exe_{"$(IFS="-"; echo "${EXE_NM[*]}" )"}/" ||\
         brcktExp_names="exe_{"$(IFS=,; echo "${EXE_NM[*]}" )"}/"

   else
      brcktExp_names="exe_${EXE_NM[0]}/"
   fi

   ############### FILES #################
   if [[ ${#INP_NM[@]} -gt 1 ]]; then
      [[ $INP_var == var ]] &&\
         brcktExp_names=$brcktExp_names"inp_{"$(IFS="-"; echo "${INP_NM[*]}" )"}/" ||\
         brcktExp_names=$brcktExp_names"inp_{"$(IFS=,; echo "${INP_NM[*]}" )"}/"
   else
      brcktExp_names="$brcktExp_names""inp_${INP_NM[0]}/"
   fi

   ############# NAMELIST ################
   #echo "len CAP = ${#CAP[@]}"
   #echo "CAP = ${CAP[@]}"


   for ((j=0;j<${#CAP[@]};j++)); do
      #echo "j = $j"
      if [[ "${NML_NM[$j]}" == *" "* ]]; then
         [[ ${NML_var[$j]} == var ]] &&\
            brcktExp_names="$brcktExp_names""${NML_ID[$j]}{$(sed -e 's/ /-/g' <<<"${NML_NM[$j]}" )}_" ||
            brcktExp_names="$brcktExp_names""${NML_ID[$j]}{$(sed -e 's/ /,/g' <<<"${NML_NM[$j]}" )}_"
      else
         brcktExp_names="$brcktExp_names""${NML_ID[$j]}${NML_NM[$j]}_"
      fi
   done
   brcktExp_names="${brcktExp_names::-1}/"

   ################ MPI ##################
   if [[ ${#PROCXY[@]} -gt 1 ]]; then
      [[ $PRC_var == var ]] &&\
         brcktExp_names="$brcktExp_names{"$(IFS=-; echo "${PROCXY[*]}" )"}/" ||\
         brcktExp_names="$brcktExp_names{"$(IFS=,; echo "${PROCXY[*]}" )"}/"
   else
      brcktExp_names="$brcktExp_names""${PROCXY[0]}/"
   fi

   #Set names of result directories
   dir_nm=$( sed -e 's|exe_||g' <<< $brcktExp_names | sed -e 's|inp_||g' | sed -e 's|/||g' \
      | sed -e 's|_||g'  )
   for ((j=0;j<${#CAP[@]};j++)); do
      dir_nm=$(sed -e "s|${NML_ID[$j]}||g" <<<$dir_nm)
   done

   #mapfile -t brcktExp_names < <( eval echo $brcktExp_names | sed -e 's/-/,/g')
   #mapfile -t dir_nm < <( eval echo $dir_nm | sed -e 's/-/,/g')
   #echo "dir_nm= $dir_nm"

   #Expand a first time, over the 'fixed' variables
   read -r -a brcktExp_names <<<$( eval echo $brcktExp_names | sed -e 's/-/,/g')
   read -r -a dir_nm <<<$( eval echo $dir_nm | sed -e 's/-/,/g')

   #echo "brcktexp = ${brcktExp_names[@]}"
   #echo "len brcktex[= ${#brcktExp_names[@]}"
   #echo "len dir_nm = ${#dir_nm[@]}"
   #echo "dir_nm= ${dir_nm[@]}"
   ################ END: BRACKET EXP VARIABLE  ###############





   ################ START ECHO RESULTS ################
   # diagnostics written for good practice, to improve reproducibiliity,
   #  structure and understanding of (older) tests
else
   echo "--------------------------"
   echo "Start writing diagnostics"
   echo "--------------------------"
   source .scripts/utils/WriteDiagnostics.sh $1 $2 $3
   echo "--------------------------"
   echo "End of writing diagnostics"
   echo "--------------------------"
fi

#If timing set Calc_timing.py file
if [ $RERUN -gt 1 ]; then
   echo "Executing timing"
   source .scripts/utils/timing/Read_timing.sh
fi
