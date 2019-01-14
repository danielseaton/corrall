from __future__ import print_function
import pandas as pd
import numpy as np
import os
import glob
import cyvcf2
import scipy.stats
import random
import argparse
import sys
import corrall

metadata_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'
metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)

# test_df = pd.read_csv('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/ase_env_interactions/combined_donor_tests/PC1.pearsonr.logcounts_all_107don_expt_defendo_leads.testrun.tsv', sep='\t', index_col=0)
# test_df = test_df.query('n_cells>2').sort_values(by='pval')
# gene_list = test_df.head(1000)['gene_id'].tolist()



# qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/qtl_results/all_results_combined.for_snp_annotation.tsv'
# qtl_file_short_name = 'all_leads_one_per_gene'
# qtl_df = pd.read_csv(qtl_filename,sep='\t')
# qtl_df['ensembl_gene_id'] = qtl_df['feature'].apply(lambda x: x.split('_')[0])
# #select lowest p val for each gene
# qtl_df = qtl_df.groupby('feature').apply(lambda x: x.sort_values(by='empirical_feature_p_value').iloc[0,:])
# #output to similar file to record gene-snp relationships
# qtl_df.to_csv(qtl_filename.replace('.tsv','.subset_for_ase_phasing.tsv'), sep='\t',index=False)
# #process for het phasing
# qtl_df = qtl_df.set_index('ensembl_gene_id', drop=False)
# qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# #qtl_df = qtl_df.query('empirical_feature_p_value < 0.00001')
# # sort low to high by p-value
# qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
# qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id'])


allelic_datafile = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/all_donors.ase.lowthresh.chrBcount.phased.genelevel.tsv'
total_datafile = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/all_donors.ase.lowthresh.totalcount.phased.genelevel.tsv'

qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/qtl_results/all_results_combined.for_snp_annotation.tsv'
qtl_file_short_name = 'all_leads'
outfile_prefix = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/ase/complete_ase_phased'
qtl_df = pd.read_csv(qtl_filename,sep='\t')
qtl_df['ensembl_gene_id'] = qtl_df['feature'].apply(lambda x: x.split('_')[0])
#process for het phasing
qtl_df = qtl_df.set_index(['ensembl_gene_id','snp_id'], drop=False)
qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# sort low to high by p-value
qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
#qtl_df = qtl_df.head(5)


# allelic_datafile = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/test_subset_of_donors.ase.lowthresh.chrBcount.phased.genelevel.tsv'
# total_datafile = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/test_subset_of_donors.ase.lowthresh.totalcount.phased.genelevel.tsv'

# qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/qtl_results/all_results_combined.for_snp_annotation.tsv'
# qtl_file_short_name = 'testing_example'
# outfile_prefix = './tests/data/complete_ase_phased'
# qtl_df = pd.read_csv(qtl_filename,sep='\t')
# qtl_df['ensembl_gene_id'] = qtl_df['feature'].apply(lambda x: x.split('_')[0])
# #process for het phasing
# qtl_df = qtl_df.set_index(['ensembl_gene_id','snp_id'], drop=False)
# qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# # sort low to high by p-value
# qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
# qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
# qtl_df = qtl_df.head(50)


# qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/qtl_analysis/leads_to_test/top_qtl_results_all_primary_global_fdr0.05.txt'
# qtl_file_short_name = 'all_hipsci_ipsc_bulk_leads'
# qtl_df = pd.read_csv(qtl_filename, sep='\t')
# qtl_df['ensembl_gene_id'] = qtl_df['feature_id'].apply(lambda x: x.split('_')[0])
# # process for het phasing
# qtl_df = qtl_df.set_index(['ensembl_gene_id','snp_id'], drop=False)
# qtl_df = qtl_df.rename(columns={'snp_chromosome': 'chrom', 'snp_position': 'pos'})
# # sort low to high by p-value
# qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
# qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
# #qtl_df = qtl_df.head(10)


# output file name
outfile_template = '{prefix}.{datatype}.{qtl_file_short_name}.tsv'
outfile_allelic_fractions = outfile_template.format(prefix=outfile_prefix,datatype='allelic_fractions',qtl_file_short_name=qtl_file_short_name)
outfile_allelic_counts = outfile_template.format(prefix=outfile_prefix,datatype='allelic_counts',qtl_file_short_name=qtl_file_short_name)
outfile_total_counts = outfile_template.format(prefix=outfile_prefix,datatype='total_counts',qtl_file_short_name=qtl_file_short_name)

outfile_type = 'all'

#get list of genes to be evaluated

gene_list = qtl_df['ensembl_gene_id'].drop_duplicates().tolist()


allelic_df = pd.read_csv(allelic_datafile, sep='\t', index_col=0)
allelic_df = allelic_df.loc[gene_list, :]
total_df = pd.read_csv(total_datafile, sep='\t', index_col=0)
total_df = total_df.loc[gene_list, :]


ase_df = allelic_df/total_df

all_cells = list(set(ase_df.columns) & set(metadata_df.index))
metadata_df = metadata_df.loc[all_cells, :]
donor_list = metadata_df['donor_long_id'].drop_duplicates().tolist()
donor2cell_dict = dict([(x,metadata_df.query('donor_long_id==@x').index) for x in donor_list])

# gene list is reduced to whatever we actually have data for
gene_list = list(ase_df.index)

snp_df = qtl_df.query('ensembl_gene_id in @gene_list')



#snp_df = snp_df.head(20)
snp_phase_df = corrall.vcf_utils.get_het_snp_phase_dataframe(snp_df, donor_list)



def phase_allelic_data(gene_name, snp_id, datatype):
    if datatype not in ['alleliccount', 'totalcount']:
        raise ValueError
    phasing_data = snp_phase_df.loc[(gene_name,snp_id), :].dropna()

    #subset to heterozygous donors
    het_donors = list(phasing_data.index)
    #map from heterozygous donors to cells
    selected_cells = [donor2cell_dict[x] for x in het_donors]
    selected_cells = [x for y in selected_cells for x in y]
    
    allelic_data = allelic_df.loc[gene_name, all_cells]
    total_data = total_df.loc[gene_name, all_cells]

    # set ASE for cells from homozygous donors to NaN
    cells_to_drop = [x for x in all_cells if x not in selected_cells]
    allelic_data.loc[cells_to_drop] = np.nan
    total_data.loc[cells_to_drop] = np.nan

    # flip data for the relevant cells
    for donor in het_donors:
        cells = donor2cell_dict[donor]
        phase = phasing_data.loc[donor]
        if phase == 0:
            # phase ASE to the selected SNP by swapping phase
            allelic_data.loc[cells] = total_data.loc[cells] - allelic_data.loc[cells]
        elif phase == 1:
            # ASE is already phased to the selected SNP
            pass

    if datatype=='alleliccount':
        return allelic_data
    elif datatype=='totalcount':
        return total_data


new_phased_allelic_df = snp_df[['ensembl_gene_id','snp_id']].apply(lambda x: phase_allelic_data(x['ensembl_gene_id'], x['snp_id'], 'alleliccount'), axis=1)

new_phased_total_df = snp_df[['ensembl_gene_id','snp_id']].apply(lambda x: phase_allelic_data(x['ensembl_gene_id'], x['snp_id'], 'totalcount'), axis=1)

new_phased_ase_df = new_phased_allelic_df/new_phased_total_df

if outfile_type in ['counts', 'all']:
    new_phased_total_df.to_csv(outfile_total_counts, sep='\t')
    new_phased_allelic_df.to_csv(outfile_allelic_counts, sep='\t')
if outfile_type in ['fractions', 'all']:
    new_phased_ase_df.to_csv(outfile_allelic_fractions, sep='\t')
