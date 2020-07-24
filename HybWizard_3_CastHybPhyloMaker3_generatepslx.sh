#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=8gb
#PBS -j oe
#PBS -N HybPhyloMaker3_generate_pslx
#PBS -m abe
#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=4G,h_data=4G,h_vmem=4G
#$ -cwd
#$ -j y
#$ -N HybPhyloMaker3_generate_pslx
#$ -o HybPhyloMaker3_generate_pslx.log


# ********************************************************************************
# *    HybPhyloMaker - Pipeline for Hyb-Seq data processing and tree building    *
# *                  https://github.com/tomas-fer/HybPhyloMaker                  *
# *         Script 03 - Process consensus after mapping, make pslx files         *
# *                                   v.1.6.4                                    *
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2018 *
# * tomas.fer@natur.cuni.cz                                                      *
# * based on Weitemier et al. (2014), Applications in Plant Science 2(9): 1400042*
# ********************************************************************************

# Input: Consensus sequences from HybPhyloMaker2 or Geneious: must be named consensus.fasta or consensus_cpDNA.fasta 
# it is multiple fasta file with names (in case of Geneious mapping):
# e.g., Camptandra-latifolia_S4-all-no-dups_assembled_to_Curcuma_exons_reference_400Ns__consensus_sequence)
# prepared in /storage/$server/home/$LOGNAME/data/30consensus/

#Complete path and set configuration for selected location
if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\nHybPhyloMaker3 is running on MetaCentrum..."
	#settings for MetaCentrum
	#Copy file with settings from home and set variables from settings.cfg
	cp -f $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	othersourcepath=/storage/$server/home/$LOGNAME/$othersource
	otherpslxpath=/storage/$server/home/$LOGNAME/$otherpslx
	#Move to scratch
	cd $SCRATCHDIR
	#Add necessary modules
	module add blat-suite-34
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\nHybPhyloMaker3 is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	othersourcepath=../$othersource
	otherpslxpath=../$otherpslx
	#Make and enter work directory
	mkdir -p workdir03
	cd workdir03
	#Add necessary modules
	module load bioinformatics/blat/36x1
else
	echo -e "\nHybPhyloMaker3 is running locally..."
	#settings for local run
	#set variables from settings.cfg
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
if [[ $cp =~ "yes" ]]; then
	echo -e "Working with cpDNA\n"
	type="cp"
else
	echo -e "Working with exons\n"
	type="exons"
fi

#Copy fasta from home folder to scratch/workdir
cp -r $path/$type/40contigs/* .

#-----------------------BLAT ASSEMBLIES TO REFERENCE-----------------------
echo -e "Generating pslx files using BLAT...\n"
#Copy other transcriptome/genome data from home to scratch/workdir (must be named with suffix *.fas)
if [ "$othersource" != "" ] && [ "$othersource" != "NO" ]; then 
	cp -r $othersourcepath/* .
fi

if [[ $cp =~ "yes" ]]; then
	#Copy cpDNA reference
	cp -r $source/$cpDNACDS .
else
	#Copy reference
	cp -r $source/$probes .
fi

#Copy script for correcting pslx files
cp -r $source/HybWizard-CastCorrect.py .

#Make a new folder for results
mkdir $path/$type/50pslx
mkdir /storage/$server/home/$LOGNAME/$otherpslx

#Make a list of all files with contigs
ls *.fas > contig_names.txt

#A loop to process all contig files specified in contig_names.txt
for contigfile in $(cat contig_names.txt)
do
	echo -e "\nProcessing $contigfile..."
	if [[ $cp =~ "yes" ]]; then
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $cpDNACDS $contigfile ${contigfile}.pslx
	else
		blat -t=DNA -q=DNA -out=pslx -minIdentity=$minident $probes $contigfile ${contigfile}.pslx
		#Modify PSLX if the consensus was called using 'ConsensusFixer'
		#BLAT change all ambiguous bases to 'n', following lines put ambiguities back
		if [[ $conscall =~ "consensusfixer" ]]; then
			echo -e "Modifying PSLX..."
			head -5 ${contigfile}.pslx > pslx_header.txt #get first 5 lines as a header
			sed -i '1,5d' ${contigfile}.pslx #delete header (first 5 linees)
			cut -f22 < ${contigfile}.pslx > pslx_sequences.txt #extract column 22, i.e. column with samples-specific sequences
			sed -i "s/,$/\n/" pslx_sequences.txt #replace ',' at thee end of line by EOL (separate lines by empty line)
			sed -i "s/,/\n/g" pslx_sequences.txt #replace each ',' by EOL (put each fragment to a separate line)
			sed -i 's/n/\[ACGTRYSWKMBDHVN\]/g' pslx_sequences.txt #replace any 'n' by every possible DNA letter incl. ambiguities (as regex)
			cat pslx_sequences.txt | while read line; do
				if [ ! -z "$line" ]; then
					grep -io -m 1 $line $contigfile | head -1 >> pslx_sequences_ambig.txt
				else
					echo >> pslx_sequences_ambig.txt
				fi
			done
			#convert output to lowercase
			tr A-Z a-z < pslx_sequences_ambig.txt > pslx_sequences_ambigLower.txt
			#format back to conform PSLX standards
			sed -i '/^$/!s/$/,/' pslx_sequences_ambigLower.txt #add ',' on non-empty lines
			tr '\n' 'Q' < pslx_sequences_ambigLower.txt | sed "s/,QQ/,\n/g" | sed 's/Q//g' > pslx_sequences_ambigLowerModif.txt
			#combine with original PSLX
			cat ${contigfile}.pslx | cut -f 1-21 > first.txt #extract first 21 columns
			cat ${contigfile}.pslx | cut -f 23 > second.txt #extract 23rd column
			paste first.txt pslx_sequences_ambigLowerModif.txt second.txt > final.txt #add the modified 22nd column
			cat pslx_header.txt final.txt > ${contigfile}.pslx #add original header
			rm pslx_header.txt pslx_sequences.txt pslx_sequences_ambig.txt pslx_sequences_ambigLower.txt pslx_sequences_ambigLowerModif.txt first.txt second.txt final.txt
		fi
	fi
	cp $contigfile.pslx $path/$type/50pslx
done

#Correct pslx files and copy results
python3 HybWizard-CastCorrect.py $probes
cp -r corrected/* /storage/$server/home/$LOGNAME/$otherpslx

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
