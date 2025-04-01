#!/usr/bin/env bash
#THIS FILE IS NOT FIT FOR RELEASE
#This file is part of the OceanVar testing suite
#File for creating diagnostics of tested jobs
#Created 28 nov 2024 10:20 by F. Carere (francesco.carere@cmcc.it)

set -u
set -e


###############################
echo ""
echo "------------------------"
echo "------------------------"
echo "--- STARTING TIMING ---"
echo "------------------------"
echo "------------------------"
for ((i=0; i<${#brcktExp_names[@]};i++)); do #loop over fixed subspaces
   read -r -a varsubspace <<<$(eval echo ${brcktExp_names[$i]})
   read -r -a vardirnm    <<<$(eval echo ${dir_nm[$i]})
   #for((j=0;j<${#varsubspace[@]};j++)); do
   #   echo "varsubspace[$j] = ${varsubspace[$j]}"
   #done
   #for((j=0;j<${#vardirnm[@]};j++)); do
   #   echo "vardirnm[$j] = ${vardirnm[$j]}"
   #done
   for ((j=0;j<${#varsubspace[@]};j++)); do #loop over variables subspaces
      for ((k=0;k<${#varsubspace[@]};k++)); do
         #Check if directories exist, otherwise exit
         cdir1="$TESTS_DIR/${varsubspace[$k]}"
         #cdir2="$TESTS_DIR/${varsubspace[$j]}"
         echo "cdir1 = $cdir1"
         #echo "cdir2 = $cdir2"
         #If timing on, then do timing routine
         echo "doing timing results"
         TIMING_DIR="$SCRIPT_DIR/../timing/$TESTNR"
         proc1=$(echo $cdir1 | perl -w -pe 's|.*\/(\w*x\w*)\/|$1|g')
         tdir1=${cdir1/$proc1\///}
         tdir1="${tdir1/tests/timing}"
         tdir1="${tdir1/\/\//\/}"
         mkdir -p "$tdir1"
         echo "tdir1 = $tdir1"
         #proc2=$(echo $cdir2 | perl -w -pe 's|.*\/(\w*x\w*)\/|$1|g')
         #tdir2=${cdir2/$proc2\///}
         #tdir2="${tdir2/tests/timing}"
         #tdir2="${tdir2/\/\//\/}"
         #mkdir -p "$tdir2"
         bash "$SCRIPT_DIR/utils/calcavg_std.sh" "$cdir1"/profile_* > "$tdir1""profileavg$proc1.txt"
         #bash "$SCRIPT_DIR/utils/calcavg_std.sh" "$cdir2"/profile_* > "$tdir2""profileavg$proc2.txt"
      done
   done
done
echo ""
echo "------------------------"
echo "------------------------"
echo "--- END TIMING ---"
echo "------------------------"
echo "------------------------"
