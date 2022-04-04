#!/bin/bash
#PBS -l walltime=8:0:0
#PBS -l select=1:ncpus=6:mem=1gb:scratch_local=2gb
#PBS -N ParalogWizard-Assemble
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

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



module add python36-modules-gcc || exit 1 # biopython
module add spades-3.14.0 || exit 1 # spades.py
module add parallel-20200322 || exit 1 # parallel
module add bwa-0.7.17 || exit 1 # bwa
module add samtools-1.10 || exit 1 # samtools
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .
cp "${source}/${probe_exons_concat}" .
#Copy scripts and reference to scratch



#Copy data to scratch
mkdir -p "${SCRATCHDIR}/${data}"/10deduplicated_reads

cp -r "${path_to_data}"/10deduplicated_reads/"$sample"*.fastq "${SCRATCHDIR}/${data}"/10deduplicated_reads
echo "${sample}" > "${data}"/10deduplicated_reads/samples_list.txt
python3 ParalogWizard.py cast_assemble -d "${data}" -pr "${probe_exons_concat}" -nc "${cpu}"

cp -r "${data}"/20assemblies/* "${path_to_data}"/20assemblies/
cp *.log "${PBS_O_WORKDIR}"/



