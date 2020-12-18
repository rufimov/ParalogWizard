#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=6:mem=1gb:scratch_local=2gb
#PBS -N ParalogWizard-Analyze
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Analyzing data and searching for paralogs...'
env echo
env echo

env echo 'Going to scratch'
cd "${SCRATCHDIR}" || exit 1

env echo


#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp -f "${PBS_O_WORKDIR}"/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
block_list=(${blocklist})
path_HPM=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
path_to_data_HPM="${data}"

#Add necessary modules
module add python-3.6.2-gcc
module add python36-modules-gcc
module add mafft-7.453
module add fasttree-2.1.8

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${path_to_data_HPM}"/exons
cp "${path_HPM}"/exons/all_hits.txt "${SCRATCHDIR}"/"${path_to_data_HPM}"/exons

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy scripts and reference to scratch
cp "${source}"/ParalogWizard.py .

env echo

env echo 'Running script...'
env echo


python3 ParalogWizard.py cast_analyze -d "${path_to_data_HPM}" -b "${block_list[@]}" -nc 6 || exit 1
env echo

env echo 'Copying results back to working directory'

#Copy results back
cp -r "${path_to_data_HPM}"/exons/aln_orth_par "${path_HPM}"/exons
cp *.log "${PBS_O_WORKDIR}"/

env echo
env echo

env echo 'Analyzing finished!'

exit 0
