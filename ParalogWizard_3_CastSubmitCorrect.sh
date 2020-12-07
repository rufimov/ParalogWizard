#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -N ParalogWizard-Correct
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
path="/storage/${server}/home/${LOGNAME}/${data}"
source=/"storage/${server}/home/${LOGNAME}/HybSeqSource"

#Add necessary modules
module add blat-suite-34
module add python-3.6.2-gcc
module add python36-modules-gcc


#Copy fasta from home folder to scratch, reference, script for generating and correcting pslx files
mkdir -p "${data}"/exons/40contigs
cp -r "${path}"/exons/40contigs "${data}"/exons/
cp -r "${source}/${probes}" .
cp -r "${source}"/ParalogWizard.py .


#Make a new folder for results
mkdir -p "${path}"/exons/50pslx

#Run script
python3 ParalogWizard.py cast_correct -d "${data}" -pp "${probes}" -i "${minident}" -r "${redlist}"
cp -r "${data}"/exons/50pslx/* "${path}"/exons/50pslx/

