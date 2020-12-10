#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=2gb
#PBS -N ParalogWizard-Create
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
path_HPM=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
path_to_data_HPM="${data}"

#Add necessary modules
module add python-3.6.2-gcc
module add python36-modules-gcc

env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${path_to_data_HPM}"/exons/aln_orth_par/
cp "${path_HPM}"/exons/all_hits.txt "${SCRATCHDIR}"/"${path_to_data_HPM}"/exons
cp "${path_HPM}"/exons/aln_orth_par/pairwise_distances.txt "${SCRATCHDIR}"/"${path_to_data_HPM}"/exons/aln_orth_par/

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy scripts and reference to scratch
cp "${source}"/ParalogWizard.py .
grep "^[^>].\{${exon_length}\}" -B1 --no-group-separator "${source}/${probe_HP_exons_split}" > "${probe_HP_exons_split}"

env echo

env echo 'Running script...'
env echo

env echo 'Copying data to scratch'
if [[ "${paralogs}" =~ "yes" ]]; then
  python3 ParalogWizard.py cast_create -d "${path_to_data_HPM}" -b "${blocklist}" -p -mi "${paralog_min_divergence}" -ma "${paralog_max_divergence}" -pe "${probe_HP_exons_split}" || exit 1
else
  python3 ParalogWizard.py cast_create -d "${path_to_data_HPM}" -b "${blocklist}" -pe "${probe_HP_exons_split}" || exit 1
fi

env echo

env echo 'Copying results back to working directory'

#Copy results back
mkdir -p "${path_HPM}"
cp -r "${path_to_data_HPM}"/exons/new_reference_for_HybPhyloMaker*.fas "${path_HPM}"/exons/
if [[ "$paralogs" =~ "yes" ]]; then
  cp -r "${path_to_data_HPM}"/exons/paralog_statistics*.tsv "${path_HPM}"/exons/
  cp -r "${path_to_data_HPM}"/exons/locus_statistics*.tsv "${path_HPM}"/exons/
  cp -r "${path_to_data_HPM}"/exons/refined* "${path_HPM}"/exons/
  cp -r "${path_to_data_HPM}"/exons/warnings.txt "${path_HPM}"/exons/
else
  cp -r "${path_to_data_HPM}"/exons/new_reference_for_HybPiper*.fas "${path_HPM}"/exons/
fi
cp *.log "${PBS_O_WORKDIR}"/


env echo
env echo

env echo 'New reference created!'

exit 0
