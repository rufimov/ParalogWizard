#!/bin/bash
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -N ParalogWizard-Separate
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

#Move to scratch
cd "${SCRATCHDIR}"

#Copy file with settings from home and set variables from settings.cfg
cp -f $PBS_O_WORKDIR/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
red_list=(${redlist})
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource

#Add necessary modules
module add blat-suite-34
module add python-3.6.2-gcc
module add python36-modules-gcc
module add mafft-7.453


#Copy fasta from home folder to scratch, reference, script for generating and correcting pslx files
mkdir -p "${data}"
cp -r "${path_to_data}"/31exonic_contigs "${SCRATCHDIR}"/"${data}"
cp "${source}/${customized_probes}" .
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .



#Run script
python3 ParalogWizard.py cast_separate -d "${data}" -pc "${customized_probes}" -i "${minident}" -r "${red_list[@]}"

#Copy results back
cp -r "${data}"/50pslx "${path_to_data}"
cp -r "${data}"/60mafft "${path_to_data}"
cp -r "${data}"/70concatenated_exon_alignments "${path_to_data}"
cp *.log "${PBS_O_WORKDIR}"/

exit 0
