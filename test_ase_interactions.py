import pandas as pd
import numpy as np
import argparse
import statistical_models
import patsy
import limix
import loader_utils

# parser = argparse.ArgumentParser()
# parser.add_argument('--factor', '-f',
#                     help='Name of factor to test.')
# parser.add_argument('--covariates', '-cv', nargs='+', default=[],
#                     help='List of covariates to be accounted for.')
# parser.add_argument('--random_effect', '-cv', default=None,
#                     help='Name of random effect to be modelled.')
# parser.add_argument('--gene_id', '-g', type=str)
# parser.add_argument('--snp_id', '-s', type=str)
# parser.add_argument('--outfile_prefix', '-o', type=str)
# args = parser.parse_args()



gene_id, snp_id = 'ENSG00000170291','17_7151111_A_C'

#datafiles
# total counts
# allelic counts
# ase fractions
# factor dataframe
# cell metadata

#parameters
# factor name
# covariate factor names
# min_n_cells
# SNP/gene name(s)/feature variant filter

min_n_cells = 50

n_genes = 20

ase_short_name = 'all_leads'
n_index_cols = 2


working_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff'

alleliccounts_file = working_dir + '/data/ase/complete_ase_phased.allelic_counts.all_leads.tsv'
allelicfractions_file = working_dir + '/data/ase/complete_ase_phased.allelic_fractions.all_leads.tsv'
totalcounts_file = working_dir + '/data/ase/complete_ase_phased.total_counts.all_leads.tsv'

factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_pseudotimeandmodules_20180618.tsv'
outfile_pattern = working_dir + '/data/ase_env_interactions/pseudotimeandmodules.{test_type}.'+ase_short_name+'.tsv'
metadata_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'


metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
#allelic_df = loader_utils.read_csv_by_idx(alleliccounts_file, sep='\t', idx_list=[(gene_id, snp_id)], index_col=list(range(n_index_cols)), nrows=n_genes)
#total_df = loader_utils.read_csv_by_idx(totalcounts_file, sep='\t', idx_list=[(gene_id, snp_id)], index_col=list(range(n_index_cols)), nrows=n_genes)
ase_df = loader_utils.read_csv_by_idx(allelicfractions_file, idx_list=[(gene_id, snp_id)], sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
factor_df = pd.read_csv(factor_file, sep='\t', index_col=0)


# qtl_list = [('ENSG00000233927', '19_8387207_G_A')]
# allelic_df = allelic_df.loc[qtl_list, :]
# total_df = total_df.loc[qtl_list, :]
# ase_df = ase_df.loc[qtl_list, :]

qtl_list = [(gene_id, snp_id)]
#qtl_list = list(ase_df.index)[10:11]

if qtl_list is None:
    qtl_list = ase_df.index

factor = 'pseudotime'
covariate_factors = []

# factor = 'respiration'
# covariate_factors = []

# factor = 'respiration'
# covariate_factors = ['pseudotime']

list_of_outputs = []


random_effect = None
#random_effect = 'donor_long_id'


for ase_idx in qtl_list[:]:
    ase_ds = ase_df.loc[ase_idx, :].dropna()
    
    if ase_ds.shape[0]< min_n_cells:
        continue

    y_ds, candidate_ds, covariate_df = statistical_models.get_model_dataframes(ase_ds, factor_df, factor, covariate_factors)
    
    if random_effect is not None:
        dummy_df = pd.get_dummies(metadata_df.loc[y_ds.index, random_effect])
        if dummy_df.shape[1]<2:
            K = None
        K = np.dot(dummy_df.values, dummy_df.values.T)
    else:
        K = None
    
    output_ds = statistical_models.test_limix_lmm(y_ds, candidate_ds, covariate_df, K)
    output_ds.name = ase_idx
    
    list_of_outputs.append(output_ds)

output_df = pd.concat(list_of_outputs, axis=1).transpose()
output_df['factor'] = factor
output_df['covariates'] = ', '.join(covariate_factors)
output_df['random_effect'] = random_effect
