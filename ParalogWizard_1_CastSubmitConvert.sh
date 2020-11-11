#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=4:mem=1gb:scratch_local=2gb
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
cp "${source}"/ParalogWizard_CastConvert.py .
cp "${source}"/ParalogWizard_Functions.py .
cp "${source}"/ParalogWizard_CastCollect.py .

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
if [[ "$collect_contigs" =~ "yes" ]]; then
  python3 ParalogWizard_CastCollect.py "${path_HP}" "${path_HPM}"
fi

mkdir -p "${SCRATCHDIR}/${data}"
cp -r /storage/"${server}"/home/"${LOGNAME}/${data}"/HybPiper_contigs "${SCRATCHDIR}/${data}"


env echo

env echo 'Running script...'
env echo


python3 ParalogWizard_CastConvert.py "${data}" "${probe_HP_exons_split}" "${length_cut}" "${spades_cover_cut}"  || exit 1
env echo

env echo 'Copying results back to working directory'

#Copy results back
mkdir -p "${path_HPM}"
cp -r "${path_to_data_HPM}"/* "${path_HPM}"

env echo
env echo

env echo 'Transferring finished!'

exit 0
