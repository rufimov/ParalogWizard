#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker3_generate_pslx
#PBS -m abe

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker3 is running on MetaCentrum..."
	#Move to scratch
	cd $SCRATCHDIR

	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	cp -f $PBS_O_WORKDIR/ParalogWizard_Settings.cfg .
	. settings.cfg
	. ParalogWizard_Settings.cfg
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx

	#Add necessary modules
	module add blat-suite-34
	module add python-3.6.2-gcc
  module add python36-modules-gcc

else
	echo -e "\nHybPhyloMaker3 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. ParalogWizard_Settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir03
	cd workdir03
fi

#Setting for the case when working with cpDNA
echo -e "Working with exons\n"
type="exons"


#Copy fasta from home folder to scratch/workdir
cp -r $path/$type/40contigs/* .

#-----------------------BLAT ASSEMBLIES TO REFERENCE-----------------------
echo -e "Generating pslx files using BLAT...\n"

#Copy reference
cp -r $source/$probes .


#Copy script for correcting pslx files
cp -r $source/ParalogWizard_CastCorrect.py .

#Make a new folder for results
mkdir $path/$type/50pslx
mkdir -p $data/exons/50pslx

#Make a list of all files with contigs
ls *.fas > contig_names.txt

#A loop to process all contig files specified in contig_names.txt
for contigfile in $(cat contig_names.txt)
do
	echo -e "\nProcessing $contigfile..."
	blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $probes $contigfile ${contigfile}.pslx
	cp $contigfile.pslx $path/$type/50pslx
	cp $contigfile.pslx $data/exons/50pslx
done

#Correct pslx files and copy results
python3 ParalogWizard_CastCorrect.py $data $probes $whitelist
cp -r $data/exons/50pslx/corrected $path/$type/50pslx/

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r workdir03
fi

echo -e "\nScript HybPhyloMaker3 finished...\n"
