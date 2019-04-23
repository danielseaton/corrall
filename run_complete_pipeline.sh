#!/bin/bash

DONOR=$1

#DONORMAPPINGFILE='/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_qc_filt_mnn_correct_slalom_metadata_20180310.subset.tsv'
#DONORMAPPINGFILE='/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_afterqc_beforefilt_metadata_20180618.linemapping.reformatted.tsv'
DONORMAPPINGFILE='/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.mapping_donors_to_ase_files.tsv'


OUTDIR='/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/ase_aggregated_by_donor/'

TOTAL=$OUTDIR$DONOR'.ase.lowthresh.totalcount.tsv'
ALT=$OUTDIR$DONOR'.ase.lowthresh.altcount.tsv'
SNPINFO=$OUTDIR$DONOR'.ase.lowthresh.snp_info.tsv'

TOTAL2=$OUTDIR$DONOR'.ase.lowthresh.totalcount.phased.tsv'
ALT2=$OUTDIR$DONOR'.ase.lowthresh.chrBcount.phased.tsv'

GENEBED='/nfs/leia/research/stegle/dseaton/genomes/hg19/annotation/Homo_sapiens.GRCh37.75.exons.genenames.bed'

python aggregate_data_by_donor.py $DONOR $DONORMAPPINGFILE $OUTDIR
source activate cyvcf2_env && python phase_ase_quantifications.py $ALT $TOTAL $SNPINFO 
source activate py27 && python aggregate_ase_by_gene.py $ALT2 $TOTAL2 $SNPINFO $GENEBED genelevel

