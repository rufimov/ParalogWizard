#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=6:mem=16gb:scratch_local=10gb
#PBS -N ParalogWizard-Extend
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
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
cpu=$TORQUE_RESC_PROC

#Add necessary modules
module add parallel-20200322
module add exonerate-2.2.0
module add python-3.6.2-gcc
module add python36-modules-gcc
module add mafft-7.453


#Copy fasta from home folder to scratch, reference, script for generating and correcting pslx files
mkdir -p "${data}"
cp -r "${path_to_data}"/30raw_contigs "${SCRATCHDIR}"/"${data}"
cp -r "${path_to_data}"/41detected_par "${SCRATCHDIR}"/"${data}"
cp "${source}/${probe_exons_concat}" .
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .



#Run script

python3 ParalogWizard.py cast_extend -d "${data}" -pr "${probe_exons_concat}" -nc "${cpu}"


#Copy results back
cp -r "${data}"/21supercontigs "${path_to_data}"
cp *.log "${PBS_O_WORKDIR}"/

exit 0
