#!/bin/bash

while read d; do
echo $d
DONOR=$d

bsub -o log -e log -M 8000 -R "rusage[mem=8000]" sh run_complete_pipeline.sh $DONOR

done </nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/list_of_singlecell_endodiff_donors.tsv

