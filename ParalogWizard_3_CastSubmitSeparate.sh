#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=64gb:scratch_local=50gb
#PBS -N ParalogWizard-Separate
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
cpu=$TORQUE_RESC_PROC

#Add necessary modules
module add blat-suite
module add mafft

module unload python
module add python/3.7.7-intel-19.0.4-mgiwa7z
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH


#Copy fasta from home folder to scratch, reference, script for generating and correcting pslx files
mkdir -p "${data}"
cp -r "${path_to_data}"/31exonic_contigs "${SCRATCHDIR}"/"${data}"
cp -r "${path_to_data}"/41detected_par "${SCRATCHDIR}"/"${data}"
cp -r "${path_to_data}"/60mafft "${SCRATCHDIR}"/"${data}"
cp "${source}/${customized_probes}" .
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .



#Run script
if [[ -z "${red_list}"  ]]; then
  python3 ParalogWizard.py cast_separate -d "${data}" -pc "${customized_probes}" -i "${minident}" -nc "${cpu}"
else
  python3 ParalogWizard.py cast_separate -d "${data}" -pc "${customized_probes}" -i "${minident}" -r "${red_list[@]}" -nc "${cpu}"
fi

#Copy results back
cp -r "${data}"/50pslx "${path_to_data}"
cp -r "${data}"/60mafft "${path_to_data}"
cp -r "${data}"/70concatenated_exon_alignments "${path_to_data}"
cp *.log "${PBS_O_WORKDIR}"/

exit 0
