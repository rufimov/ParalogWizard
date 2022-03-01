#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=12:mem=8gb:scratch_local=2gb
#PBS -N ParalogWizard-Detect
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Creating new reference'
env echo
env echo

env echo 'Going to scratch'
cd "${SCRATCHDIR}" || exit 1

env echo

#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp "${PBS_O_WORKDIR}"/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
block_list=(${blocklist})
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
cpu=$TORQUE_RESC_PROC
# shellcheck disable=SC2003
#cpu=$( expr "$TORQUE_RESC_PROC" '*' 6 )

#Add necessary modules
module add python-3.6.2-gcc
module add python36-modules-gcc

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${data}"/40aln_orth_par/
mkdir -p "${SCRATCHDIR}/${data}"/31exonic_contigs/
cp "${path_to_data}"/31exonic_contigs/all_hits.tsv "${SCRATCHDIR}"/"${data}"/31exonic_contigs

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy scripts and reference to scratch
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .
#grep "^[^>].\{${exon_length}\}" -B1 --no-group-separator "${source}/${probe_exons_split}" > "${probe_exons_split}"
cp "${source}"/"${probe_exons_split}" .

env echo

env echo 'Running script...'
env echo

env echo 'Copying data to scratch'
if [[ "${paralogs}" =~ "yes" ]]; then
  cp "${path_to_data}"/40aln_orth_par/pairwise_distances.tsv "${SCRATCHDIR}"/"${data}"/40aln_orth_par/
  if [[ -z "${block_list}"  ]]; then
    python3 ParalogWizard.py cast_detect -d "${data}" -p -mi "${paralog_min_divergence}" -ma "${paralog_max_divergence}" -pe "${probe_exons_split}" -nc "${cpu}"|| exit 1
    env echo
  else
    python3 ParalogWizard.py cast_detect -d "${data}" -b "${block_list[@]}" -p -mi "${paralog_min_divergence}" -ma "${paralog_max_divergence}" -pe "${probe_exons_split}" -nc "${cpu}" || exit 1
    env echo
  fi
  env echo 'Copying results back to working directory'

  #Copy results back
  cp -r "${data}"/41detected_par "${path_to_data}"
  cp *.log "${PBS_O_WORKDIR}"/

  env echo
  env echo
else
  if [[ -z "${block_list}"  ]]; then
    python3 ParalogWizard.py cast_detect -d "${data}" -pe "${probe_exons_split}" -nc "${cpu}" || exit 1
    env echo
  else
    python3 ParalogWizard.py cast_detect -d "${data}" -b "${block_list[@]}" -pe "${probe_exons_split}" -nc "${cpu}" || exit 1
    env echo
  fi

  env echo 'Copying results back to working directory'

  #Copy results back
  cp -r "${data}"/41without_par "${path_to_data}"
  cp *.log "${PBS_O_WORKDIR}"/
  env echo
  env echo
fi



env echo 'New reference created!'

exit 0
