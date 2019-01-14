import pandas as pd
import numpy as np
import argparse
import statistical_models
import patsy
import limix
import loader_utils

parser = argparse.ArgumentParser()
parser.add_argument('--factor', '-f',
                    help='Name of factor to test.')
parser.add_argument('--covariates', '-cv', nargs='+', default=[],
                    help='List of covariates to be accounted for.')
parser.add_argument('--random_effect', '-re', default=None,
                    help='Name of random effect to be modelled.')
parser.add_argument('--gene_id', '-g', type=str)
parser.add_argument('--snp_id', '-s', type=str)
parser.add_argument('--outfile_prefix', '-o', type=str)
args = parser.parse_args()

factor = args.factor
covariates = args.covariates
random_effect = args.random_effect
gene_id = args.gene_id
snp_id = args.snp_id
outfile_prefix = args.outfile_prefix

min_n_cells = 50

ase_short_name = 'all_leads'
n_index_cols = 2


qtl_list = [(gene_id, snp_id)]


working_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff'

alleliccounts_file = working_dir + '/data/ase/complete_ase_phased.allelic_counts.all_leads.tsv'
allelicfractions_file = working_dir + '/data/ase/complete_ase_phased.allelic_fractions.all_leads.tsv'
totalcounts_file = working_dir + '/data/ase/complete_ase_phased.total_counts.all_leads.tsv'

factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_pseudotimeandmodules_20180618.tsv'
outfile_name = '{outfile_prefix}.ase_interaction_results.{gene_id}.{snp_id}.tsv'.format(outfile_prefix=outfile_prefix,
                                                                                        gene_id=gene_id,
                                                                                        snp_id=snp_id)
metadata_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'

metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
#allelic_df = loader_utils.read_csv_by_idx(alleliccounts_file, sep='\t', idx_list=qtl_list, index_col=list(range(n_index_cols)))
#total_df = loader_utils.read_csv_by_idx(totalcounts_file, sep='\t', idx_list=qtl_list, index_col=list(range(n_index_cols)))
ase_df = loader_utils.read_csv_by_idx(allelicfractions_file, idx_list=qtl_list, sep='\t', index_col=list(range(n_index_cols)))
factor_df = pd.read_csv(factor_file, sep='\t', index_col=0)

list_of_outputs = []

for ase_idx in qtl_list[:]:
    ase_ds = ase_df.loc[ase_idx, :].dropna()
    
    if ase_ds.shape[0]< min_n_cells:
        continue

    y_ds, candidate_ds, covariate_df = statistical_models.get_model_dataframes(ase_ds, factor_df, factor, covariates)
    
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

if len(list_of_outputs)>0:
    output_df = pd.concat(list_of_outputs, axis=1).transpose()
    output_df['factor'] = factor
    output_df['covariates'] = ', '.join(covariates)
    output_df['random_effect'] = random_effect
    output_df.to_csv(outfile_name, sep='\t')
else:
    with open(outfile_name, 'w') as f:
        f.write('')
