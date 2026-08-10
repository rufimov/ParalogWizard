#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=6:mem=64gb:ompthreads=6:scratch_local=1Tb
#PBS -N ParalogWizard-Assemble
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



# module unload python
# module add python/3.7.7-intel-19.0.4-mgiwa7z
# export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
# export PATH=$PYTHONUSERBASE/bin:$PATH
# export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH
# # module add python36-modules/
# # module unload perl
# # module add perl/5.30.3-intel-19.0.4-ww4zelp
# module add bwa
# module add samtools/1.11
# # module add spades/3.14.0bin 

# module add spades/3.1

# module add parallel

# module add spades-3.14.0 || exit 1 # spades.py
# module add parallel-20200322 || exit 1 # parallel

# module unload python
# module add  python/3.9.12-gcc-10.2.1-rg2lpmk
# export PYTHONUSERBASE=/storage/brno12-cerit/home/rufimov/python39/cvmfs/software.metacentrum.cz/spack18/software/linux-debian11-x86_64_v2/gcc-10.2.1/python-3.9.12-rg2lpmkxpcq423gx5gmedbyam7eibwtc
# export PATH=$PYTHONUSERBASE/bin:$PATH
# export PYTHONPATH=$PYTHONUSERBASE/lib/python3.9/site-packages:$PYTHONPATH

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

cp -r "${path_to_data}"/10deduplicated_reads/"$sample"*.fastq.gz .
gunzip *gz
mv *fastq "${SCRATCHDIR}/${data}"/10deduplicated_reads
echo "${sample}" > "${data}"/10deduplicated_reads/samples_list.txt
python3 ParalogWizard.py cast_assemble -d "${data}" -pr "${probe_exons_concat}" -nc "${cpu}"

cp -r "${data}"/20assemblies/* "${path_to_data}"/20assemblies/
cp *.log "${PBS_O_WORKDIR}"/



