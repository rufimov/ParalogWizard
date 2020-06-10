#!/bin/bash
#PBS -l walltime=4:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=2gb
#PBS -N HybWizzard
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

echo 'Transfering data from HybPiper to HybPhyloMaker'
echo
echo

#Copy file with settings from home and set variables from settings.cfg
echo 'Setting variables'
cp ${PBS_O_WORKDIR}/settings.cfg ${PBS_O_WORKDIR}/settings_transfer.cfg .
. settings.cfg
. HybWizzard-Settings.cfg
path_HP=/storage/${server_HP}/home/${LOGNAME}/${data_HybPiper}
path_HPM=/storage/${server}/home/${LOGNAME}/${data}
source=/storage/${server}/home/${LOGNAME}/HybSeqSource
path_to_data_HP="${data_HybPiper}"
path_to_data_HPM="${data}"
probe_HP="${probe_HP_one_repr}"

echo

echo 'Going to scratch'

#Add necessary modules
module add blast+-2.8.0a

echo

echo 'Copying data to scratch'
#Copy data to scratch
mkdir -p ${SCRATCHDIR}/${path_to_data_HP}
cd $path_HP
for folder in $(find . -maxdepth 1 -type d | sed 's/.\///' | tail -n +2); do
  cd $folder
  echo "Processing ${folder}"
  for gene in $(find . -maxdepth 1 -type d | sed 's/.\///' | tail -n +2); do
     locus=$gene
     if test -f "$locus/${locus}_contigs.fasta"; then
     cat $locus/${locus}_contigs.fasta | sed "s/>/>${locus}_/g" >> ${SCRATCHDIR}/${path_to_data_HP}/${folder}_contigs.fasta
     fi          
   done
   cd ..
done

#Move to scratch
cd ${SCRATCHDIR}

#Copy script and reference to scratch
cp ${source}/${probe_HP} .
cp ${source}/HybWizzard-CastConvert.py .

echo

echo 'Running script...'
echo

python3 HybWizzard-CastConvert.py ${path_to_data_HP} ${path_to_data_HPM} ${probe_HP_one_repr} ${length_cut} ${spades_cover_cut} | tee HybWizzard-CastConvert.log || exit 1
echo

echo 'Copying results back to working directory'

#Copy results back
mkdir ${path_HPM}
cp -r ${path_to_data_HPM}/* ${path_HPM}
cp HybWizzard-CastConvert.log ${PBS_O_WORKDIR}

echo
echo

echo 'Transfering finished!'

exit
