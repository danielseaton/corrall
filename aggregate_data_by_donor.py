from __future__ import print_function
import os
import pandas as pd
import numpy as np
import glob
import re
import collections
from datetime import datetime
import argparse
#import plink_utils
import vcf_utils

parser = argparse.ArgumentParser(description='Collect data for all samples from a particular individual.')
parser.add_argument('donor_id', help='donor ID, to match format from the donor_id column in donor_id_mapping_file')
parser.add_argument('donor_id_mapping_file', help='tab-separated file with \"donor_id\" and \"filepath\" columns, mapping donors to ase quantification files')
parser.add_argument('out_dir', help='output directory')
args = parser.parse_args()

donor_id = args.donor_id
donor_id_mapping_file = args.donor_id_mapping_file
out_dir = args.out_dir

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

outfile_template = os.path.join(out_dir,'{donor_id}.ase.lowthresh.{filetype}.tsv')

donor_df = pd.read_csv(donor_id_mapping_file, sep='\t')

sub_donor_df = donor_df.query('donor_id==@donor_id')
filelist = sub_donor_df['filepath'].tolist()
n_files_expected = len(filelist)
missing_filelist = [x for x in filelist if not os.path.exists(x)]
filelist = [x for x in filelist if os.path.exists(x)]
n_files = len(filelist)

print('{} files for donor {} ({} missing)'.format(n_files, donor_id,n_files_expected - n_files))

if len(missing_filelist)>0:
    print('Missing:')
    print([os.path.basename(x) for x in missing_filelist])

alt_outfile = outfile_template.format(donor_id=donor_id, filetype='altcount')
total_outfile = outfile_template.format(donor_id=donor_id, filetype='totalcount')
snp_outfile = outfile_template.format(donor_id=donor_id, filetype='snp_info')

#threshold on number of cells in which a specific snp was quantified
count_threshold = len(filelist) // 100

#select snps
counter_dict = collections.Counter()
sum_dict = collections.defaultdict(lambda :0)

selected_snps = []

for filename in filelist[:]:
    df = pd.read_csv(filename, sep='\t')
    variant_ids = df['variantID'].tolist()
    values = (df['altCount']/df['totalCount']).tolist()
    
    for vid,value in zip(variant_ids,values):
        sum_dict[vid] += value
        counter_dict[vid] += 1
            
selected_snps = [x for x in counter_dict if counter_dict[x]>count_threshold]
selected_snps = vcf_utils.filter_to_het_snps(selected_snps, donor_id)

print('{} SNPs selected of {}'.format(len(selected_snps),len(counter_dict.keys())))

#make tables
sample_ids = [os.path.basename(x).split('.')[0] for x in filelist]
alt_df = pd.DataFrame(index=selected_snps,columns=sample_ids)
total_df = pd.DataFrame(index=selected_snps,columns=sample_ids)

snp_df_columns= ['contig','position','variantID','refAllele','altAllele']
snp_df = pd.DataFrame(np.nan,index=selected_snps,columns=snp_df_columns)
    
for idx,filename in enumerate(filelist):
    sample_id = os.path.basename(filename).split('.')[0]
    
    df = pd.read_csv(filename, sep='\t', index_col=2)
    
    #df['ratio'] = df['refCount']/df['totalCount']
    
    snp_subset = list(set(selected_snps) & set(df.index))
    
    alt_df.loc[snp_subset,sample_id] = df.loc[snp_subset,'altCount']
    total_df.loc[snp_subset,sample_id] = df.loc[snp_subset,'totalCount']
    snp_df.loc[snp_subset,snp_df_columns] = df.loc[snp_subset,snp_df_columns]

alt_df.to_csv(alt_outfile,sep='\t')
total_df.to_csv(total_outfile,sep='\t')
snp_df.to_csv(snp_outfile,sep='\t')
