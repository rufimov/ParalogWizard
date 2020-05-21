#!/bin/bash
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=2gb
#PBS -N Convert_HP_HPM
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

echo 'Transfering data from HybPiper to HybPhyloMaker'
echo
echo

#Copy file with settings from home and set variables from settings.cfg
echo 'Setting variables'
cp ${PBS_O_WORKDIR}/settings.cfg ${PBS_O_WORKDIR}/settings_transfer.cfg .
. settings.cfg
. settings_transfer.cfg
path_HP=/storage/${server_HP}/home/${LOGNAME}/${data_HybPiper}
path_HPM=/storage/${server}/home/${LOGNAME}/${data}
source=/storage/${server}/home/${LOGNAME}/HybSeqSource
path_to_data_HP="${data_HybPiper}"
path_to_data_HPM="${data}"
probe_HP="${probe_HP_one_repr}"

echo

echo 'Going to scratch'
#Move to scratch
cd ${SCRATCHDIR}
#Add necessary modules
module add blast+-2.8.0a

echo

echo 'Copying data to scratch'
#Copy data to scratch
mkdir ${path_to_data_HP}
cp -r ${path_HP}/*.dedup ${path_to_data_HP}

#Copy script and reference to scratch
cp ${source}/${probe_HP} .
cp ${source}/convert_HP_HPM_with_paralogs.py .

echo

echo 'Running script'

python3 convert_HP_HPM_with_paralogs.py ${path_to_data_HP} ${path_to_data_HPM} ${probe_HP_one_repr} | tee convert_HP_HPM.log
echo

echo 'Copying results back to working directory'

#Copy results back
mkdir ${path_HPM}
cp -r ${path_to_data_HPM}/* ${path_HPM}
cp convert_HP_HPM.log ${PBS_O_WORKDIR}

echo
echo

echo 'Transfering finished!'
