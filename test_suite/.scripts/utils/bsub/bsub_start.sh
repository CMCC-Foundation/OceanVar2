#!/usr/bin/env bash

#BSUB -q s_short
#BSUB -J TestOceanVar
#BSUB -n 1
#BSUB -R "rusage[mem=1G]"
#BSUB -o out/logStartTest@POSTFIX@.%J.txt
#BSUB -e err/StartTest@POSTFIX@.%J.txt
#BSUB -P R000
#BSUB -K

#This file is part of the OceanVar testing suite: this file submits the startjob of the test
#Written by F. Carere (francesco.carere@cmcc.it) in June-July 2024
#Please see the documentation for a high-level overview of the testing suite

mkdir -p err out

module load --auto intel-2021.6.0/2021.6.0
module load --auto impi-2021.6.0/2021.6.0

export I_MPI_HYDRA_BOOTSTRAP=lsf
export I_MPI_HYDRA_BRANCH_COUNT=2
export I_MPI_HYDRA_COLLECTIVE_LAUNCH=1

mpiexec.hydra -l @PATH@/.scripts/StartTestOceanVar.sh\
   $1 `#TEST_NAME`\
   $2 `#TEST_CONFFILE_NAME`\
   $3 `#NAMELIST_TEMPLATE_NAME`\
   $4 `#RUNTEST_TEMPLATE_NAME`\
   $5 `#FULL_CONFIG_SUBDIR`\
   $6 `#FILES_OUTPUT`\
   $7 `#ARCHITECTURE`\
   $8 `#NR OF RERUN, only used for timing`
