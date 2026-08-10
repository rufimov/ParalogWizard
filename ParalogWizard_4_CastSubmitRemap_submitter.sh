#!/bin/bash

. ParalogWizard_Settings.cfg

while read -r sample;
do
  echo "Submitting ${sample}"
  qsub -v sample="${sample}" -N ParalogWizard-Remap_"${sample}" ParalogWizard_4_CastSubmitRemap.sh
done < /storage/"${server}/home/${LOGNAME}/${data}"/10deduplicated_reads/samples_list.txt
