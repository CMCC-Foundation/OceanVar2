#!/usr/bin/env bash
#This file is part of the OceanVar testing suite: in particular writing a diagnostics file
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a High-level overview of the testing suite
#Not sure if this file is completely functional for the case !='Start"


function usage() {
   echo "-------"
   echo "Please pass the following arguments:"
   echo "1: Mode of writing diagnostics (July 2024: either \"Start\" which is the first read, or not which is the second read"
   echo "2: Full path (including name) of metadata/input file to be considered"
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

#######################################
################ EXE ##################
#######################################
echo ""
echo "-----------------"
echo "INPUT: executables"
echo "-----------------"

echo "Number of executables = $nexe"
echo ""



#Write to console,stdout
echo "Executable reference names are:"
for((i=0;i<${#EXE_NM[@]};i++)); do
   echo "   $i: ${EXE_NM[$i]}" 
done
echo ""
int_b=false
if [ $1 == 'Start' ]; then
   echo "Executable reference paths are: "
   for((i=0;i<${#EXE_PT[@]};i++)); do
      echo "   $i: ${EXE_PT[$i]}"
      if [[ ${EXE_PT[$i]} != /* ]]; then
         int_b=true
      fi
   done
   if $int_b; then
      echo "   Interpreted as"
      for((i=0;i<${#EXE_PT[@]};i++)); do
         if [[ ${EXE_PT[$i]} != /* ]]; then
            echo "   $i: $(realpath $PWD/${EXE_PT[$i]})"
         fi
      done
   fi
fi

#Write to diagnostics file
if [ $1 == 'Start' ]; then #write to metadata file
timestamp=$(date)
echo "This file contains metadata for test with testname $3 executed on" $timestamp>$2
int_b=false
   echo "###########################">>$2
   echo "########## EXE ############">>$2
   echo "###########################">>$2
   echo "Executable reference names are:">>$2
   for((i=0;i<${#EXE_NM[@]};i++)); do
      echo "   $i: ${EXE_NM[$i]}" >> $2
   done
   echo ""
   echo "Executable reference paths are:">>$2
   for((i=0;i<${#EXE_PT[@]};i++)); do
      echo "   $i: ${EXE_PT[$i]}">>$2
      if [[ ${EXE_PT[$i]} != /* ]]; then
         int_b=true
      fi
   done
   if $int_b; then
      echo "   Interpreted as">>$2
      for((i=0;i<${#EXE_PT[@]};i++)); do
         if [[ ${EXE_PT[$i]} != /* ]]; then
            echo "   $i: $(realpath $PWD/${EXE_PT[$i]})">>$2
         fi
      done
   fi
fi

#######################################
############### FILES #################
#######################################
echo ""
echo "-----------------"
echo "INPUT: input files"
echo "-----------------"
echo "Number of directories with input files = $nexe"
echo ""



int_b=false
echo "Input directory reference names are:"
for((i=0;i<${#INP_NM[@]};i++)); do
   echo "   $i: ${INP_NM[$i]}"
done
echo ""
if [ $1 == 'Start' ]; then
   echo "Input directory paths are:"
   for((i=0;i<${#INP_PT[@]};i++)); do
      echo "   $i: ${INP_PT[$i]}"
      if [[ ${INP_PT[@]} != /* ]]; then
         int_b=true
      fi
   done

   if  $int_b ; then
      echo "   Interpreted as"
      for((i=0;i<${#INP_PT[@]};i++)); do
         if [[ ${INP_PT[@]} != /* ]]; then
            echo "   $i: $(realpath $PWD/${INP_PT[$i]})"
         fi
      done
   fi
fi

int_b=false
if [ $1 == 'Start' ]; then #write to metadata file
   echo "###########################">>$2
   echo "########## INP ############">>$2
   echo "###########################">>$2
   echo "Input reference names are:">>$2
   for((i=0;i<${#INP_NM[@]};i++)); do
      echo "   $i: ${INP_NM[$i]}">>$2
   done
   echo ""
   echo "Input reference paths are:">>$2
   for((i=0;i<${#INP_PT[@]};i++)); do
      echo "   $i: ${INP_PT[$i]}">>$2
      if [[ ${INP_PT[@]} != /* ]]; then
         int_b=true
      fi
   done

   if  $int_b ; then
      echo "   Interpreted as">>$2
      for((i=0;i<${#INP_PT[@]};i++)); do
         if [[ ${INP_PT[@]} != /* ]]; then
            echo "   $i: $(realpath $PWD/${INP_PT[$i]})">>$2
         fi
      done
   fi
fi

#######################################
############# NAMELIST ################
#######################################
echo ""
echo "-----------------"
echo "INPUT: namelist configurations"
echo "-----------------"
echo "Number of configurations = ${nvar[@]}"


echo "Namelist configuration Identities (as to be replaced in the namelist) are:"
for((i=0;i<${#NML_ID[@]};i++)); do
   #echo "NML_ID[$i] = ${NML_ID[$i]}"
   echo "   $i: ${NML_ID[$i]}"
done
echo "Namelist configuration reference names are:"
for((i=0;i<${#NML_NM[@]};i++)); do
   #echo "NML_NM[$i] = "${NML_NM[$i]}
   echo "   $i: ${NML_NM[$i]}"
done
if [ $1 == 'Start' ]; then
   echo "Namelist configuration substitutes (the acutal values replacing the configuration identities) are:"
   for((i=0;i<${#NML_SB[@]};i++)); do
      #echo "NML_SB[$i] = "${NML_SB[$i]}
      echo "   $i :${NML_SB[$i]}"
   done
fi


if [ $1 == 'Start' ]; then #write to metadata file
   echo "###########################">>$2
   echo "########## NML ############">>$2
   echo "###########################">>$2
   echo "Namelist configuration identities:">>$2
   for((i=0;i<${#NML_ID[@]};i++)); do
      echo "   $i: ${NML_ID[$i]}">>$2
   done
   echo "Namelist configuration reference names:">>$2
   for((i=0;i<${#NML_NM[@]};i++)); do
      echo "   $i: ${NML_NM[$i]}">>$2
   done
   echo "Namelist configuration substitutes:">>$2
   for((i=0;i<${#NML_SB[@]};i++)); do
      echo "   $i: ${NML_SB[$i]}">>$2
   done
fi


#######################################
################ MPI ##################
#######################################
echo ""
echo "-----------------"
echo "INPUT: MPI decomposition"
echo "-----------------"

#Get linenumbers of previous information
echo "Number of MPI decompositions = $nproc"



echo "MPI decompositions are:"
for((i=0;i<${#PROCX[@]};i++)); do
   echo "${PROCX[$i]} x ${PROCY[$i]}"
done

if [ $1 == 'Start' ]; then #write to metadata file
   echo "###########################">>$2
   echo "########## MPI ############">>$2
   echo "###########################">>$2
   echo "MPI decompositions are:">>$2
   for((i=0;i<${#PROCX[@]};i++)); do
      echo "   $i: ${PROCX[$i]} x ${PROCY[$i]}">>$2
   done
fi
