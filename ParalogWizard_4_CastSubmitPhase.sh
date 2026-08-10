#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=2:mem=32gb:scratch_local=50Gb
#PBS -N ParalogWizard-Phase
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

#Move to scratch
cd "${SCRATCHDIR}"

#Copy file with settings from home and set variables from settings.cfg
cp -f $PBS_O_WORKDIR/ParalogWizard_Settings.cfg .
cp -f /auto/brno12-cerit/nfs4/home/rufimov/Crataegus_phylogeny_all/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
cpu=$TORQUE_RESC_PROC

#Add necessary modules
module unload python
module add python/python-3.7.7-gcc-8.3.0-t4loj4a
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH


module add bcftools/1.11
module add samtools/1.11
module add raxml/8.2.12


#mkdir -p "${SCRATCHDIR}/${data}"/100remapped

#for dir in "${path_to_data}"/100remapped/exons*; do
#  [ -d "$dir" ] || continue
  
#  subfolder=$(basename "$dir")
  
#  mkdir -p "${SCRATCHDIR}/${data}"/100remapped/"${subfolder}"
  
#  cp "$dir"/reference* "${SCRATCHDIR}/${data}"/100remapped/"${subfolder}" 

#done

#cp "${path_to_data}"/100remapped/all_variants* "${SCRATCHDIR}/${data}"/100remapped/
#cp "${path_to_data}"/100remapped/samples_list.txt "${SCRATCHDIR}/${data}"/100remapped/

mkdir -p "${SCRATCHDIR}/${data}"/101phased
#cp "${path_to_data}"/101phased/* "${SCRATCHDIR}/${data}"/101phased
cp "${path_to_data}"/101phased/combined* "${SCRATCHDIR}/${data}"/101phased

cp "${source}/${customized_probes}" .
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .



#Run script

python3 ParalogWizard.py cast_phase -d ${data} -nc ${cpu}

#Copy results back
cp -r "${data}"/101phased "${path_to_data}"

cp *.log "${PBS_O_WORKDIR}"/

exit 0
