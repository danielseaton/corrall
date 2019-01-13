import pandas as pd
import numpy as np
import argparse
import statistical_models
import patsy
import limix

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

n_genes = 5

ase_short_name = 'all_leads'
n_index_cols = 2

min_n_cells = 50

working_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff'

alleliccounts_file = working_dir + '/data/ase/subset_complete_ase_phased.allelic_counts.all_leads.tsv'
allelicfractions_file = working_dir + '/data/ase/subset_complete_ase_phased.allelic_fractions.all_leads.tsv'
totalcounts_file = working_dir + '/data/ase/subset_complete_ase_phased.total_counts.all_leads.tsv'

factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_pseudotimeandmodules_20180618.tsv'
outfile_pattern = working_dir + '/data/ase_env_interactions/pseudotimeandmodules.{test_type}.'+ase_short_name+'.tsv'
metadata_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'


allelic_df = pd.read_csv(alleliccounts_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
total_df = pd.read_csv(totalcounts_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
ase_df = pd.read_csv(allelicfractions_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
factor_df = pd.read_csv(factor_file, sep='\t', index_col=0)


# qtl_list = [('ENSG00000233927', '19_8387207_G_A')]
# allelic_df = allelic_df.loc[qtl_list, :]
# total_df = total_df.loc[qtl_list, :]
# ase_df = ase_df.loc[qtl_list, :]

if qtl_list is None:
    qtl_list = ase_df.index

factor = 'pseudotime'
covariate_factors = []

# factor = 'respiration'
# covariate_factors = []

# factor = 'respiration'
# covariate_factors = ['pseudotime']

list_of_outputs = []

for ase_idx in qtl_list:
    ase_ds = ase_df.loc[ase_idx, :].dropna()

    if ase_ds.shape[0]< min_n_cells:
        continue
    
    y_ds, candidate_ds, covariate_df = statistical_models.get_model_dataframes(ase_ds, factor_df, factor, covariate_factors)
    
    output_ds = statistical_models.test_limix_lmm(y_ds, candidate_ds, covariate_df)
    output_ds.name = ase_idx

    list_of_outputs.append(output_ds)

output_df = pd.concat(list_of_outputs, axis=1).transpose()
output_df['factor'] = factor
output_df['covariates'] = ', '.join(covariate_factors)
