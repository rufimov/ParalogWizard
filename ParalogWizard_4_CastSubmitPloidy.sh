#!/bin/bash
#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=12:mem=128gb:scratch_local=1Tb
#PBS -N ParalogWizard-Ploidy
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
module unload python
module add python/3.7.7-intel-19.0.4-mgiwa7z
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH
module add samtools


#Copy data to scratch
mkdir -p "${SCRATCHDIR}/${data}"/10deduplicated_reads
mkdir -p "${SCRATCHDIR}/${data}"/100remapped

cp "${path_to_data}"/10deduplicated_reads/samples_list.txt "${SCRATCHDIR}/${data}"/10deduplicated_reads/
cp -r "${path_to_data}"/100remapped/exons* "${SCRATCHDIR}/${data}"/100remapped/
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .
cp -r "${source}"/nQuire .
chmod +x nQuire

PATH="$(pwd):$PATH"

python3 ParalogWizard.py cast_ploidy -d "${data}" -nc "${cpu}" || exit 1

cp -r "${data}"/102ploidy "${path_to_data}"/
cp *.log "${PBS_O_WORKDIR}"/
