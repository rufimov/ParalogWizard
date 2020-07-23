#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=6:mem=1gb:scratch_local=2gb
#PBS -N HybWizard
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Analyzing data and searching for paralogs...'
env echo
env echo

#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp "${PBS_O_WORKDIR}"/settings.cfg "${PBS_O_WORKDIR}"/HybWizard_Settings.cfg .
. settings.cfg
. HybWizard_Settings.cfg
path_HPM=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
path_to_data_HPM="${data}"

env echo

env echo 'Going to scratch'

#Add necessary modules
module add python-3.6.2-gcc
module add python36-modules-gcc

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${path_to_data_HPM}"/exons
cp "${path_HPM}"/exons/all_hits.txt "${SCRATCHDIR}"/"${path_to_data_HPM}"/exons

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy scripts and reference to scratch
cp "${source}"/HybWizard_CastAnalyze.py .
cp "${source}"/HybWizard_Functions.py .

env echo

env echo 'Running script...'
env echo


python3 HybWizard_CastAnalyze.py "${path_to_data_HPM}" "${blacklist}"  || exit 1
env echo

env echo 'Copying results back to working directory'

#Copy results back
cp -r "${path_to_data_HPM}"/exons/new_reference_for_HybPhyloMaker.fas "${path_HPM}"/exons

env echo
env echo

env echo 'Analyzing finished!'

exit 0
