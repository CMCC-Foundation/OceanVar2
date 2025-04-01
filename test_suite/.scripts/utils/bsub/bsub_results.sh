#!/usr/bin/env bash

#BSUB -q s_short
#BSUB -J ResultsOceanVar
#BSUB -n 1
#BSUB -R "rusage[mem=5G]"
#BSUB -o out/logresultsTest@POSTFIX@.%J.txt
#BSUB -e err/ResultsTest@POSTFIX@.%J.txt
#BSUB -P R000
@ADD_DEPENDENCY@

#This file is part of the OceanVar testing suite: this file bsubs the calculations
#and the plotting
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a high-level overview of the testing suite
mkdir -p err out

module load --auto intel-2021.6.0/2021.6.0
module load --auto impi-2021.6.0/2021.6.0

export I_MPI_HYDRA_BOOTSTRAP=lsf
export I_MPI_HYDRA_BRANCH_COUNT=2
export I_MPI_HYDRA_COLLECTIVE_LAUNCH=1


mpiexec.hydra -l @PATH@/.scripts/ResultsAfterStart.sh \
   $1 `#TEST_NAME` \
   $2 `#RESULTS_NAME`\
   $3 `#PLOTS_NAME`\
   $4 `#RESULTS_CONFFILE_NAME`\
   $5 `#FULL_CONFIG_SUBDIR`\
   $6 `#TYPE_RESULTS`\
   $7 `#FILES_OUTPUT`\
   $8 `#CONDA_ENVIRONMENT`\
   $9 `#ARCHITECTURE JUNO or LOCAL`\
   $(10) `#RERUN amount, used only for timing`
