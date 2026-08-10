#!/bin/bash
#PBS -l walltime=162:0:0
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=500gb
#PBS -N ParalogWizard-ABBA_BABA
#PBS -m abe
#PBS -j oe

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'clean_scratch' TERM

cd "${SCRATCHDIR}"
cp "${PBS_O_WORKDIR}"/ParalogWizard_Settings.cfg .
. ParalogWizard_Settings.cfg
path_to_data=/storage/"${server}/home/${LOGNAME}/${data}"
source=/storage/"${server}/home/${LOGNAME}"/HybSeqSource


module unload python-2.6.6-gcc #avoid conflict with the next module
module add python36-modules-gcc
module add newick-utils-13042016



mkdir -p "${data}"/110abba_baba

cd "${data}"/110abba_baba

# Prefer compressed cast_call output; fall back to legacy uncompressed VCF.
if [[ -f "${path_to_data}"/100remapped/all_variants.vcf.gz ]]; then
  gzip -dc "${path_to_data}"/100remapped/all_variants.vcf.gz \
    | sed "s/${OUTGROUP}/Outgroup/" \
    | gzip > all_variants.vcf.gz
elif [[ -f "${path_to_data}"/100remapped/all_variants.vcf ]]; then
  cp "${path_to_data}"/100remapped/all_variants.vcf .
  sed -i "s/${OUTGROUP}/Outgroup/" all_variants.vcf
  gzip all_variants.vcf
else
  echo "ERROR: neither all_variants.vcf.gz nor all_variants.vcf found in ${path_to_data}/100remapped" >&2
  exit 1
fi

# Modify species tree
cp "${path_to_data}"/72trees70_75/RAxML/species_trees/Astral/Astral_70_75.tre .
# replace ' ' back to '-' and '_'
sed 's/ \([^ ]*\) / \1_/g' Astral_70_75.tre > sptree.tre #replace every second occurrence of ' ' by '_'
sed -i 's/ /-/g' sptree.tre #replace all spaces by '-'
# remove bootstrap values
nw_topology -bI sptree.tre > tmp && mv tmp sptree.tre
# rename ${OUTGROUP} to 'Outgroup' (requirement of Dsuite)
sed -i "s/$OUTGROUP/Outgroup/g" sptree.tre

# Create species set (each sample is a separate species)
nw_labels -I sptree.tre > SpeciesSet.txt
paste SpeciesSet.txt SpeciesSet.txt > tmp && mv tmp SpeciesSet.txt

cp "${source}"/DtriosParallel .
cp "${source}"/Dsuite .
cp "${source}"/dtools.py .
cp "${source}"/plot_d.rb .
cp "${source}"/plot_f4ratio.rb .

./Dsuite Dtrios -c -n gene_flow -t sptree.tre --ABBAclustering all_variants.vcf.gz SpeciesSet.txt
./Dsuite Fbranch sptree.tre SpeciesSet_gene_flow_tree.txt > SpeciesSet_gene_flow_Fbranch.txt
./dtools.py -n gene_flow SpeciesSet_gene_flow_Fbranch.txt sptree.tre
cairosvg gene_flow.svg -o gene_flow.pdf

cut -f 2 SpeciesSet.txt | uniq > plot_order.txt
ruby ./plot_d.rb SpeciesSet_gene_flow_BBAA.txt plot_order.txt 0.7 SpeciesSet_gene_flow_BBAA_D.svg
cairosvg SpeciesSet_gene_flow_BBAA_D.svg -o SpeciesSet_gene_flow_BBAA_D.pdf
ruby ./plot_f4ratio.rb SpeciesSet_gene_flow_BBAA.txt plot_order.txt 0.2 SpeciesSet_gene_flow_BBAA_f4ratio.svg
cairosvg SpeciesSet_gene_flow_BBAA_f4ratio.svg -o SpeciesSet_gene_flow_BBAA_f4ratio.pdf
rm plot_order.txt

cp -r ../110abba_baba "${path_to_data}"/

