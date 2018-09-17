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




qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/qtl_results/all_results_combined.for_snp_annotation.tsv'
qtl_file_short_name = 'all_leads'
qtl_df = pd.read_csv(qtl_filename,sep='\t')
qtl_df['ensembl_gene_id'] = qtl_df['feature'].apply(lambda x: x.split('_')[0])
#process for het phasing
qtl_df = qtl_df.set_index(['ensembl_gene_id','snp_id'], drop=False)
qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# sort low to high by p-value
qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
# #qtl_df = qtl_df.head(50)



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


# qtl_filename = '/nfs/leia/research/stegle/acuomo/mean_day3_chr1/logcounts_all_107don_expt_defendo_leads/top_qtl_results_all.txt'
# qtl_file_short_name = os.path.basename(os.path.dirname(qtl_filename))
# qtl_df = pd.read_csv(qtl_filename,sep='\t')
# qtl_df = qtl_df.set_index('ensembl_gene_id', drop=False)
# qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# #qtl_df = qtl_df.query('empirical_feature_p_value < 0.00001')
# # sort low to high by p-value
# qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
# qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id'])

#qtl_df = qtl_df.loc[gene_list,:].dropna(how='all')



#get list of genes to be evaluated

gene_list = qtl_df['ensembl_gene_id'].tolist()

#subset donors

n_extra_donors = 3000

donor_list = ['HPSI0514i-puie_5','HPSI0214i-poih_4','HPSI0514i-letw_1','HPSI0813i-guss_1','HPSI0413i-nudd_1','HPSI1014i-sehl_6','HPSI0114i-joxm_1']
with open('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/list_of_singlecell_endodiff_donors.tsv', 'r') as f:
    extra_donors = [x.strip() for x in f.readlines()][:n_extra_donors]
    donor_list.extend(extra_donors)

donor_list = list(set(donor_list))[:]


snp_df = qtl_df

# ase_df gives data for each gene, using phased SNP info to give proportion of expression from chrB
#ase_df

#donor2cell_dict

donor2cell_dict = dict()


for count_type in ['chrBcount','totalcount']:
    df_list = []
    for donor in donor_list[:]:
        try:
            filename = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/ase_aggregated_by_donor/{donor}.ase.lowthresh.{count_type}.phased.genelevel.tsv'.format(donor=donor,count_type=count_type)
            df = pd.read_csv(filename,sep='\t',index_col=0)
            df = df.reindex(gene_list)
            donor2cell_dict[donor] = list(df.columns)
            df_list.append(df)
        except:
            print('No file for cell line {}'.format(donor))
            pass
    print('Combining data from {} cell lines'.format(len(df_list)))
    df = pd.concat(df_list, axis=1)
    if count_type=='totalcount':
        total_df = df
    else:
        allelic_df = df

# donor list is reduced to whatever we actually have data for
donor_list = donor2cell_dict.keys()

ase_df = allelic_df/total_df

# re-index ase_df to match gene-snp pairs in snp_df
gene_to_snp_mapping_df = snp_df[['ensembl_gene_id','snp_id']].set_index('ensembl_gene_id')

ase_df = ase_df.join(gene_to_snp_mapping_df, how='inner')
ase_df['ensembl_gene_id'] = ase_df.index
ase_df = ase_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
ase_df = ase_df.set_index(['ensembl_gene_id','snp_id'])


all_cells = list(ase_df.columns)

#snp_df = snp_df.head(20)
snp_phase_df = corrall.vcf_utils.get_het_snp_phase_dataframe(snp_df, donor_list)



def phase_data_for_gene(gene_name):
    phasing_data = snp_phase_df.loc[gene_name, :].dropna()

    #subset to heterozygous donors
    het_donors = list(phasing_data.index)
    #map from heterozygous donors to cells
    selected_cells = [donor2cell_dict[x] for x in het_donors]
    selected_cells = [x for y in selected_cells for x in y]
    

    ase_data = ase_df.loc[gene_name, all_cells]

    # set ASE for cells from homozygous donors to NaN
    cells_to_drop = [x for x in all_cells if x not in selected_cells]
    ase_data.loc[cells_to_drop] = np.nan

    # flip data for the relevant cells
    for donor in het_donors:
        cells = donor2cell_dict[donor]
        phase = phasing_data.loc[donor]
        if phase == 0:
            # phase ASE to the selected SNP
#            pass
            ase_data.loc[cells] = ase_data.loc[cells].apply(lambda x: 1.0-x)
        elif phase == 1:
            # ASE is already phased to the selected SNP
            pass

    return ase_data

# gene_list = snp_df['ensembl_gene_id'].tolist()
# gene_df = pd.DataFrame(gene_list, columns=['gene_id'], index=gene_list)
# new_phased_ase_df = gene_df['gene_id'].apply(phase_data_for_gene)

index_list = list(snp_df.index)
new_phased_ase_data = [phase_data_for_gene(x) for x in index_list]
new_phased_ase_df = pd.DataFrame(new_phased_ase_data)

new_phased_ase_df.to_csv('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/ase/complete_ase_phased.{}.tsv'.format(qtl_file_short_name), sep='\t')
