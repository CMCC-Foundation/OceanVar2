#!/usr/bin/env bash
#This file is part of the OceanVar testing suite
#THIS FILE SHOULD NOT BE RELEASED
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a High-level overview of the testing suite
set -u
set -e

echo "Executing timing"

thisSCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
TIMING_DIR="$thisSCRIPT_DIR/../../timing/$TESTNR"
CALCTIMEFILE="$TIMING_DIR/Calc_timing.py"
mkdir -p "$TIMING_DIR/globaltime"

cp "$SCRIPT_DIR/utils/timing/Calc_timing.py"  
perl -i -pe "s/\@RERUN\@/$RERUN/g"  > "$CALCTIMEFILE"
#if test -d "$SCRIPT_DIR/../timing/$TESTNR/"; then
perl -i -pe "s/\@EXE\@/$(IFS=,; echo "${EXE_NM[*]@Q}")/g" "$CALCTIMEFILE"
perl -i -pe "s/\@INP\@/$(IFS=,; echo "${INP_NM[*]@Q}" )/g" "$CALCTIMEFILE"
jj=0
for i in "${nvar[@]}"; do
   perl -i -pe "s/\@NMLID\@/[\@NMLID\@/g" "$CALCTIMEFILE"
   for((kk=$jj;kk<$i+$jj;kk++)); do
      perl -i -pe "s/\@NMLID\@/$(IFS=,; echo "${NML_ID[kk]@Q}" ),\@NMLID\@/g" "$CALCTIMEFILE"
   done
   jj=$jj+$i
   perl -i -pe "s/,\@NMLID\@/],\@NMLID\@/g" "$CALCTIMEFILE"
done
perl -i -pe "s/,\@NMLID\@//g" "$CALCTIMEFILE"
#for((kk=0;kk<${#NML_NM[@]};kk++)); do
#   perl -i -pe "s/\@NMLNM\@/[$(IFS=,; echo "${NML_NM[kk]@Q}" )],\@NMLNM\@/g" "$CALCTIMEFILE"
#done
#perl -i -pe "s/,\@NMLNM\@//g" "$CALCTIMEFILE"
perl -i -pe "s/\@NMLNM\@/$(IFS=,; echo "${NML_NM[*]@Q}" )/g" "$CALCTIMEFILE"
perl -i -pe "s/\@PROC\@/$(IFS=,; echo "${PROCXY[*]@Q}" )/g" "$CALCTIMEFILE"
