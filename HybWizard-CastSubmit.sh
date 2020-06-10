#!/bin/bash
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=2gb
#PBS -N HybWizard
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Transferring data from HybPiper to HybPhyloMaker'
env echo
env echo

#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp "${PBS_O_WORKDIR}"/settings.cfg "${PBS_O_WORKDIR}"/HybWizard-Settings.cfg .
. settings.cfg
. HybWizard-Settings.cfg
path_HP=/storage/"${server_HP}/home/${LOGNAME}/${data_HybPiper}"
path_HPM=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
path_to_data_HP="${data_HybPiper}"
path_to_data_HPM="${data}"
probe_HP="${probe_HP_one_repr}"

env echo

env echo 'Going to scratch'

#Add necessary modules
module add blast+-2.8.0a

env echo

env echo 'Copying data to scratch'
#Copy data to scratch
mkdir -p "${SCRATCHDIR}/${path_to_data_HP}"
cd "${path_HP}" || exit 1
for folder in $(find . -maxdepth 1 -type d | sed 's/.\///' | tail -n +2); do
  cd "${folder}" || exit 1
  env echo "Processing ${folder}"
  for gene in $(find . -maxdepth 1 -type d | sed 's/.\///' | tail -n +2); do
     locus=$gene
     if test -f "$locus/${locus}_contigs.fasta"; then
     sed "s/>/>${locus}_/g" < "${locus}"/"${locus}"_contigs.fasta >> "${SCRATCHDIR}/${path_to_data_HP}/${folder}"_contigs.fasta
     fi          
   done
   cd ..
done

#Move to scratch
cd "${SCRATCHDIR}" || exit 1

#Copy script and reference to scratch
cp "${source}/${probe_HP}" .
cp "${source}"/HybWizard-CastConvert.py .

env echo

env echo 'Running script...'
env echo

# shellcheck disable=SC2086
python3 HybWizard-CastConvert.py "${path_to_data_HP}" "${path_to_data_HPM}" "${probe_HP_one_repr}" "${length_cut}" "${spades_cover_cut}" "${new_reference_bool}" || exit 1
env echo

env echo 'Copying results back to working directory'

#Copy results back
mkdir "${path_HPM}"
cp -r "${path_to_data_HPM}"/* "${path_HPM}"
cp HybWizard-CastConvert.log "${PBS_O_WORKDIR}"

env echo
env echo

env echo 'Transferring finished!'

exit 0
