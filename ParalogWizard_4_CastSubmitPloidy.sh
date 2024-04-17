#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=16gb:scratch_local=10gb
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
module add python/python-3.7.7-gcc-8.3.0-t4loj4a
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37/cvmfs/software.metacentrum.cz/spack1/software/python/linux-debian10-x86_64/3.7.7-gcc-t4loj4
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH
module add bwa-0.7.17 || exit 1 # bwa
module add mafft-7.453
module add samtools-1.10 || exit 1 # samtools


#Copy data to scratch
mkdir -p "${SCRATCHDIR}/${data}"/10deduplicated_reads
mkdir -p "${SCRATCHDIR}/${data}"/41detected_par

cp -r "${path_to_data}"/10deduplicated_reads/* "${SCRATCHDIR}/${data}"/10deduplicated_reads
cp "${path_to_data}"/41detected_par/all_paralogs_for_reference.tsv "${SCRATCHDIR}/${data}"/41detected_par
cp "${source}/${customized_probes}" .
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .
cp -r "${source}"/nQuire .
chmod +x nQuire

python3 ParalogWizard.py cast_ploidy -d "${data}" -pc "${customized_probes}"  -nc "${cpu}" -e 300

cp -r "${data}"/ploidy/ "${path_to_data}"/
cp *.log "${PBS_O_WORKDIR}"/