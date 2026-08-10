#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=32gb:scratch_local=50gb
#PBS -N ParalogWizard-Analyze
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
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
cpu=$TORQUE_RESC_PROC


#Add necessary modules
module unload python
module add python/3.7.7-intel-19.0.4-mgiwa7z
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH
module add mafft
module add fasttree

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${data}"/31exonic_contigs
cp "${path_to_data}"/31exonic_contigs/all_hits.tsv "${SCRATCHDIR}"/"${data}"/31exonic_contigs/

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy scripts to scratch
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .

env echo

env echo 'Running script...'
env echo

if [[ -z "${block_list}"  ]]; then
  python3 ParalogWizard.py cast_analyze -d "${data}" -nc "${cpu}" || exit 1
  env echo
else
  python3 ParalogWizard.py cast_analyze -d "${data}" -b "${block_list[@]}" -nc "${cpu}"
  env echo
fi

env echo 'Copying results back to working directory'

#Copy results back
cp -r "${data}"/40aln_orth_par "${path_to_data}"/
cp *.log "${PBS_O_WORKDIR}"/

env echo
env echo

env echo 'Analyzing finished!'

exit 0
