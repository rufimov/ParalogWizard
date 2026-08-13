#!/bin/bash
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=12:mem=16gb:scratch_local=100gb
#PBS -N ParalogWizard-Polyploid
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

env echo 'Phasing polyploids (cast_polyploid)...'
env echo
env echo

env echo 'Going to scratch'
cd "${SCRATCHDIR}" || exit 1

env echo

#Copy file with settings from home and set variables from settings.cfg
env echo 'Setting variables'
cp -f "${PBS_O_WORKDIR}"/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource
cpu=$TORQUE_RESC_PROC
ns="${polyploid_n_nearest:-5}"

env echo "data=${data}"
env echo "path_to_data=${path_to_data}"
env echo "source=${source}"
env echo "cpu=${cpu}"
env echo "n_nearest=${ns}"
env echo

# Python logs live in the PBS stdout (unbuffered) and in *.log copied home
export PYTHONUNBUFFERED=1

#Add necessary modules
env echo 'Loading modules'
module unload python
module add python/3.7.7-intel-19.0.4-mgiwa7z
export PYTHONUSERBASE=/storage/${server}/home/${LOGNAME}/python37
export PATH=$PYTHONUSERBASE/bin:$PATH
export PYTHONPATH=$PYTHONUSERBASE/lib/python3.7/site-packages:$PYTHONPATH
module add mafft
module add iqtree
env echo "python3=$(command -v python3)"
env echo "mafft=$(command -v mafft || true)"
env echo "iqtree=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)"
env echo

#Copy data to scratch
env echo 'Copying data to scratch'
mkdir -p "${SCRATCHDIR}/${data}"/41detected_par
mkdir -p "${SCRATCHDIR}/${data}"/10deduplicated_reads

cp "${path_to_data}"/41detected_par/all_paralogs_for_reference.tsv \
  "${SCRATCHDIR}/${data}"/41detected_par/ || exit 1
cp "${path_to_data}"/10deduplicated_reads/samples_list.txt \
  "${SCRATCHDIR}/${data}"/10deduplicated_reads/ || exit 1

# Backbone / polyploid lists: HybSeqSource, else submit cwd, else data folder
if [[ -f "${source}/backbone_list.txt" ]]; then
  cp -f "${source}/backbone_list.txt" .
elif [[ -f "${PBS_O_WORKDIR}/backbone_list.txt" ]]; then
  cp -f "${PBS_O_WORKDIR}/backbone_list.txt" .
elif [[ -f "${path_to_data}/backbone_list.txt" ]]; then
  cp -f "${path_to_data}/backbone_list.txt" .
else
  env echo "ERROR: backbone_list.txt not found"
  exit 1
fi

if [[ -f "${source}/polyploid_list.txt" ]]; then
  cp -f "${source}/polyploid_list.txt" .
elif [[ -f "${PBS_O_WORKDIR}/polyploid_list.txt" ]]; then
  cp -f "${PBS_O_WORKDIR}/polyploid_list.txt" .
elif [[ -f "${path_to_data}/polyploid_list.txt" ]]; then
  cp -f "${path_to_data}/polyploid_list.txt" .
else
  env echo "ERROR: polyploid_list.txt not found"
  exit 1
fi

#Copy scripts to scratch
cp "${source}"/ParalogWizard.py .
cp -r "${source}"/ParalogWizard .

# Prefer highest IQ-TREE binary shipped with HybSeqSource
for _iq in iqtree3 iqtree2 iqtree; do
  if [[ -x "${source}/${_iq}" ]]; then
    cp "${source}/${_iq}" .
    chmod +x "${_iq}"
    PATH="$(pwd):$PATH"
    env echo "Staged ${_iq} from HybSeqSource onto PATH"
    break
  fi
done
unset _iq

env echo
env echo 'backbone_list.txt:'
cat backbone_list.txt
env echo
env echo 'polyploid_list.txt:'
cat polyploid_list.txt
env echo

copy_logs() {
  env echo 'Copying log files back'
  cp -f ./*.log "${PBS_O_WORKDIR}"/ 2>/dev/null || true
}

env echo 'Running script...'
env echo "python3 ParalogWizard.py cast_polyploid -d ${data} -bb backbone_list.txt -pl polyploid_list.txt -ns ${ns} -nc ${cpu}"
env echo

python3 ParalogWizard.py cast_polyploid \
  -d "${data}" \
  -bb backbone_list.txt \
  -pl polyploid_list.txt \
  -ns "${ns}" \
  -nc "${cpu}"
rc=$?

env echo
if [[ "${rc}" -ne 0 ]]; then
  env echo "cast_polyploid failed with exit ${rc}"
  copy_logs
  exit "${rc}"
fi

env echo 'Copying results back to working directory'
cp -r "${data}"/101polyploid "${path_to_data}"/ || exit 1
copy_logs

env echo
env echo
env echo 'cast_polyploid finished!'

exit 0
