#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=6:mem=1gb:scratch_local=2gb
#PBS -N ParalogWizard-Convert
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Transferring data from HybPiper to HybPhyloMaker'
env echo
env echo

env echo 'Going to scratch'
cd "${SCRATCHDIR}" || exit 1

env echo

#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp "${PBS_O_WORKDIR}"/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
path_HP=/storage/"${server_HP}/home/${LOGNAME}/${data_HybPiper}"
path_HPM=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
path_to_data_HPM="${data}"

#Add necessary modules
module add blast+-2.8.0a
module add python-3.6.2-gcc
module add python36-modules-gcc

env echo

#Copy scripts and reference to scratch
grep "^[^>].\{${exon_length}\}" -B1 --no-group-separator "${source}/${probe_HP_exons_split}" > "${probe_HP_exons_split}"
cp "${source}"/ParalogWizard.py .


env echo

#Copy data to scratch
env echo 'Copying data to scratch'
if [[ "$collect_contigs" =~ "yes" ]]; then
  python3 ParalogWizard.py cast_collect -c "${path_HP}" -d "${path_HPM}"
fi

mkdir -p "${SCRATCHDIR}/${data}"
cp -r /storage/"${server}"/home/"${LOGNAME}/${data}"/HybPiper_contigs "${SCRATCHDIR}/${data}"


env echo

env echo 'Running script...'
env echo


python3 ParalogWizard.py cast_retrieve -d "${data}" -pe "${probe_HP_exons_split}" -l "${length_cut}" -s "${spades_cover_cut}" -nc 6 || exit 1
env echo

env echo 'Copying results back to working directory'

#Copy results back
mkdir -p "${path_HPM}"
cp -r "${path_to_data_HPM}"/* "${path_HPM}"
cp *.log "${PBS_O_WORKDIR}"/

env echo
env echo

env echo 'Transferring finished!'

exit 0
