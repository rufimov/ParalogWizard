#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=6:mem=1gb:scratch_local=2gb
#PBS -N ParalogWizard-Retrieve
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
cpu=$TORQUE_RESC_PROC

path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource


#Add necessary modules
module add blast+-2.8.0a
module add python-3.6.2-gcc
module add python36-modules-gcc

env echo

#Copy scripts and reference to scratch
grep "^[^>].\{${exon_length}\}" -B1 --no-group-separator "${source}/${probe_exons_split}" > "${probe_exons_split}"
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .

env echo

#Copy data to scratch
mkdir -p "${SCRATCHDIR}/${data}"

env echo 'Copying data to scratch'
if [[ "$collect_contigs" =~ "yes" ]]; then
  cp -r "${path_to_data}/20assemblies" "${SCRATCHDIR}/${data}"
else
  cp -r "${path_to_data}/30raw_contigs" "${SCRATCHDIR}/${data}"
fi

env echo

env echo 'Running script...'
env echo

if [[ "$collect_contigs" =~ "yes" ]]; then
  python3 ParalogWizard.py cast_retrieve -d "${data}" -c -pe "${probe_exons_split}" -l "${length_cut}" -s "${spades_cover_cut}" -nc 6 || exit 1

else
  python3 ParalogWizard.py cast_retrieve -d "${data}" -pe "${probe_exons_split}" -l "${length_cut}" -s "${spades_cover_cut}" -nc 6 || exit 1

fi

env echo

env echo 'Copying results back to working directory'

#Copy results back
if [[ "$collect_contigs" =~ "yes" ]]; then
  cp -r "${data}/30raw_contigs" "${path_to_data}"
  cp -r "${data}/31exonic_contigs" "${path_to_data}"
else
  cp -r "${data}/31exonic_contigs" "${path_to_data}"
fi
cp *.log "${PBS_O_WORKDIR}"/

env echo
env echo

env echo 'Transferring finished!'

exit 0
