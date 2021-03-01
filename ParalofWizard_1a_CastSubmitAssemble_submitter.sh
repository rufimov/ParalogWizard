#!/bin/bash

. ParalogWizard_Settings.cfg

while read -r sample;
do
  echo "Submitting ${sample}"
  qsub ParalofWizard_1a_CastSubmitAssemble.sh -v sample="${sample}"
done < /storage/"${server}/home/${LOGNAME}/${data}"/10deduplicated_reads/samples_list.txt
