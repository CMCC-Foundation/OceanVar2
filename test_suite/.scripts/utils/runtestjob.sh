#!/usr/bin/env bash
#This file is part of the OceanVar testing suite
#File for submitting jobs for automated testing of OceanVar
#Created 18 nov 2024 13:28 by F. Carere (francesco.carere@cmcc.it)

set -u
set -e


#TODO FCARERE20/05/2024: Write perl script instead??

#Loop over and prepare all testdirectories and submit tests
#Nested loop as in configuration file, loops over
#  executables
#  directories with input files
#  namelist options
#  mpi process decompositions
for((idexe=0;idexe<nexe;idexe++));do
   #copy executable to right directory
   curfulldir="$OUTT_DIR/exe_${EXE_NM[$idexe]}"
   mkdir -p "$curfulldir"
   cur_exe=$( realpath "$curfulldir")"/exe_${EXE_NM[$idexe]}"
   cp "${EXE_PT[$idexe]}" "$cur_exe"

   for((idinp=0;idinp<ninp;idinp++));do
      #copy input to right directory
      curfulldir="$OUTT_DIR/exe_${EXE_NM[$idexe]}/inp_${INP_NM[$idinp]}"
      mkdir -p "$curfulldir"
      #Input is heavy, so let's not copy it
      #cp -r "${INP_PT[$idinp]}" "$curfulldir"

      while true ; do
         #condition checked at end of loop
         #Namelist condition harder than simple for loop

         ##################################
         ### Setting the directory name ###
         ##################################
         cur_dirnm=""
         jobname="$1_${EXE_NM[$idexe]}${INP_NM[$idinp]}"
         for ((j=0;j<${#STAT[@]};j++));do
            curvargrp=${NML_NM[ $j ]}
            curgrp=$( awk -v idt="${STAT[$j]}"  ' { print $idt }' <<< "${NML_NM[$j]}" )
            jobname="$jobname"$( awk -v idt="${STAT[$j]}"  ' { print $idt }' <<< "${NML_NM[$j]}" )

            cur_dirnm=$cur_dirnm"_"${NML_DR[$j]}$curgrp
            #curvarnmgrp=${NML_SB[$j]}
         done

         cur_dirnm=${cur_dirnm:1}
         curfulldir="$OUTT_DIR/exe_${EXE_NM[$idexe]}/inp_${INP_NM[$idinp]}/$cur_dirnm"
         mkdir -v -p "$curfulldir"

         ###################################
         ## END Setting the directory name #
         ###################################

         #####################################
         ########## SET NML AND BSUB #########
         #####################################
         kk=0 #kk keeps track of how many values are substituted
         curnml="$curfulldir/OceanVar_nml"
         touch "$curnml"

         #Substitute values in the namelist
         sed -e 's,@INPDIR@,'"${INP_RPT[$idinp]},g" "$NMLFILE" > "$curnml"
         for ((j=0;j<${#STAT[@]};j++));do
            #cp $NMLFILE $curnml
            for ((l=0;l<${nvar[$j]};l++));do
               idtt=$((l + kk ))
               curnmlval=${NML_ID[$idtt]}
               rplnmlval=$(awk -F '|' -v rn="${STAT[$j]}" '{print $rn}' <<< "${NML_SB[$idtt]}")


               #Error check if variable in inputfile is also set in namelist template
               if ! grep -q "@${curnmlval^^}@" "$NMLFILE"; then
                  echo "------------------">&2
                  echo "ERROR in $0">&2
                  echo "The variable @${curnmlval^^}@ is non-existent in the namelist, but exists in the settings file $INPUTFILE">&2
                  echo "Please add it to the namelist template file $NMLFILE">&2
                  exit 1
               fi
               #End error check

               #If error check okay, replace given value in namelist template
               perl -i -pe "s/\@${curnmlval^^}\@/$rplnmlval/g" "$curnml"

            done
            kk=$(( kk + ${nvar[$j]}))

         done

         #Set the processors decomposition in the namelist and in the LSF job file
         for ((i=0;i<nproc;i++));do
            cprocx=${PROCX[$i]}
            cprocy=${PROCY[$i]}
            cproc=$(( cprocx * cprocy ))
            curpcdir="$curfulldir/$cprocx"x"$cprocy"
            jobnameproc="$jobname$cprocx$cprocy"
            nmlnm=$curpcdir"/OceanVar_nml"
            runnm=$curpcdir"/run_test.sh"

            #Update processes in namelist
            mkdir -p "$curpcdir"
            sed -e "s/@PX@/$cprocx/g" "$curnml" | sed -e "s/@PY@/$cprocy/g" > "$nmlnm"

            #Error check for namelist if there are any "@[...]@" variables unset
            if grep -q -e "@.*@" "$nmlnm"; then
               echo "-----------------">&2
               echo "ERROR IN $0">&2

               errvar=$( grep "@.*@" "$nmlnm" )
               echo "Variables $errvar is still in the namelist after substituting all values">&2
               echo "Perhaps you forgot to remove this in the namelist template $NMLFILE">&2

               exit 1
            fi
            bsubq="s"
            if [[ $cproc -gt 1 ]]; then
               bsubq="p"
            fi
            echo "---"



            #Update bsub script
            sed -e "s/@NPROC@/$cproc/g" "$RUNTEST" |  \
               sed -e "s/@NPROCX@/$(printf '%02d' $cprocx)/g" | \
               sed -e "s/@NPROCY@/$(printf '%02d' $cprocy)/g" | \
               sed -e "s/@Q@/$bsubq/g" | \
               sed -e "s,@EXE@,$cur_exe,g" | \
               sed -e "s/@JOBNAME@/$jobnameproc""0/g" |\
               sed -e "s,@SCRIPTDIR@,$SCRIPT_DIR,g" > "$runnm"
            #sed -e "s,@EXE@,${EXE_RPT[$idexe]},g" | \

            #BSUB THE SCRIPT
            olddir=$PWD
            cd "$curpcdir" || (echo "error $curpcdir does not exist"; exit 1)
            #check if all output files exist already, otherwise run
            for OUTFILE in ${OUTPUTFILES[@]}; do
               #echo "checking if $OUTFILE exists in $curpcdir"
               if [ ! -f "$OUTFILE" ]; then
                  echo "Submitting ${curpcdir/$OUTT_DIR/}"
                  echo "Submitted since test is not done: output $OUTFILE is missing"
                  if [ "$ARCH" == "JUNO" ]; then
                     perl -i -pe "s,tempidx=.*,tempidx=0,g" "$runnm"
                     bsub < $runnm
                     perl -i -pe "s,$jobnameproc""0,$jobnameproc,g" "$runnm"
                     #submit extra jobs for timing (in released version of OceanVar this never happens)
                     for((ii=1;ii<$RERUN;ii++)) do
                        perl -i -pe "s,$jobnameproc,$jobnameproc$ii,g" "$runnm"
                        perl -i -pe "s,tempidx=.*,tempidx=$ii,g" "$runnm"
                        bsub -w "done($jobnameproc$((ii-1)))" -ti < $runnm
                        perl -i -pe "s,$jobnameproc$ii,$jobnameproc,g" "$runnm"
                     done
                  else
                     tty_n=$(tty)
                     echo "Starting run $runnm" | tee "$tty_n"
                     #echo "Timing nr $ii" | tee "$tty_n"
                     #submit nr of jobs for timing (in released version of OceanVar RERUN is equal to one)
                     for((ii=0;ii<$RERUN;ii++)) do
                        bash $runnm
                        mv "profile_$cprocx_$cproct.txt" "profile_$cprocx_$cproct_$(($ii-1)).txt" 2>/dev/null
                     done
                     echo "Finished run $runnm" | tee "$tty_n"
                  fi
                  break
               else
                  echo "$OUTFILE exists in $curpcdir"
               fi
            done
            #if compgen -G "$curpcdir/*.nc" > /dev/null; then
            #   echo "Not submitting script in $curpcdir since output files seem to be present:"
            #   compgen -G *.sh
            #else
            #   bsub < $runnm
            ##echo "bsub < ${runnm/$OUTT_DIR}"
            #fi
            #echo "bsub with name $jobnameproc"
            cd "$olddir" || exit


         done
         rm $curnml
         #########################################
         ########## END SET NML AND BSUB #########
         #########################################


         #####################################
         ## Updating end of loop condition ##
         #####################################
         k=0
         while [ ${STAT[$k]} == ${CAP[$k]} ]; do
            STAT[$k]=1
            k=$(( k + 1 ))
            #####################################
            ## Checking end of loop condition ##
            #####################################
            if [ $k -ge ${#STAT[@]}  ]; then
               break 2
            fi
         done
         STAT[$k]=$(( ${STAT[$k]} + 1 ))
         ########################################
         ## END Updating end of loop condition ##
         ########################################
      done
   done
done
