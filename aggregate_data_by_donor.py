from __future__ import print_function
import os
import pandas as pd
import numpy as np
import glob
import re
import collections
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(description='Collect data for all samples from a particular individual.')
parser.add_argument('donor_id')
parser.add_argument('donor_id_mapping_file')
args = parser.parse_args()

filepattern = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_raw/scrnaseq/{run_id}/ase/low_thresh/{sample_id}.ase.lowthresh.tsv'
#filepattern = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_raw/scrnaseq/{run_id}/ase/high_thresh/{sample_id}.ase.highthresh.tsv'

out_dir = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/ase_aggregated_by_donor'

#outfile_template = os.path.join(out_dir,'{donor}.ase.lowthresh.tsv')
outfile_template = os.path.join(out_dir,'{donor_id}.ase.highthresh.tsv')

donor_id = args.donor_id
donor_id_mapping_file = args.donor_id_mapping_file

donor_df = pd.read_csv(donor_id_mapping_file, sep='\t')
donors = donor_df['donor_long_id'].unique().tolist()
sub_donor_df = donor_df.query('donor_long_id==@donor_id')

filelist = [filepattern.format(run_id=row['run_id'],sample_id=row['sample_id']) for idx,row in sub_donor_df.iterrows()]
filelist = [x for x in filelist if os.path.exists(x)]

print('{} files for donor {}'.format(len(filelist), donor_id))

ase_outfile = outfile_template.format(donor_id=donor_id)
snp_outfile = ase_outfile.replace('.tsv','.snp_info.tsv')

#threshold on number of cells in which a specific snp was quantified
count_threshold = len(filelist) // 100
#    count_threshold = 15
#    count_threshold = 2 ##for testing only
# threshold to look for heterozygosity, while allowing for ase bias
lower_mean_threshold = 0.02
upper_mean_threshold = 1.0 - lower_mean_threshold


#select snps
counter_dict = collections.Counter()
sum_dict = collections.defaultdict(lambda :0)

selected_snps = []

for filename in filelist[:]:
    df = pd.read_csv(filename, sep='\t')
    variant_ids = df['variantID'].tolist()
    values = (df['refCount']/df['totalCount']).tolist()
    
    for vid,value in zip(variant_ids,values):
        sum_dict[vid] += value
        counter_dict[vid] += 1
            
selected_snps = [x for x in counter_dict if counter_dict[x]>count_threshold]
selected_snps = [x for x in selected_snps if (sum_dict[x]/counter_dict[x]>=lower_mean_threshold 
                                            and sum_dict[x]/counter_dict[x]<=upper_mean_threshold)]

print('{} SNPs selected of {}'.format(len(selected_snps),len(counter_dict.keys())))

#make tables
sample_ids = [os.path.basename(x).split('.')[0] for x in filelist]
ase_df = pd.DataFrame(index=selected_snps,columns=sample_ids)

snp_df_columns= ['contig','position','variantID','refAllele','altAllele']
snp_df = pd.DataFrame(np.nan,index=selected_snps,columns=snp_df_columns)
    
for idx,filename in enumerate(filelist):
    sample_id = os.path.basename(filename).split('.')[0]
    
    df = pd.read_csv(filename, sep='\t', index_col=2)
    
    df['ratio'] = df['refCount']/df['totalCount']
    
    snp_subset = list(set(selected_snps) & set(df.index))
    
    ase_df.loc[snp_subset,sample_id] = df.loc[snp_subset,'ratio']
    snp_df.loc[snp_subset,snp_df_columns] = df.loc[snp_subset,snp_df_columns]

ase_df.to_csv(ase_outfile,sep='\t')
snp_df.to_csv(snp_outfile,sep='\t')
