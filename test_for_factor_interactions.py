import re
import pandas as pd
import numpy as np
import os
import scipy.stats
import random
import argparse
import statsmodels.api as sm
from sklearn import decomposition
from itertools import combinations


permute = False
#permute = True

working_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff'

#ase_file = working_dir + '/data/ase/complete_ase_phased_to_defendo_leads_old.tsv'
#ase_file = working_dir + '/data/ase/complete_ase_phased.all_leads_one_per_gene.tsv'
#ase_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_exprs_20180618.tsv'
#n_index_cols = 1
#ase_short_name = os.path.basename(ase_file).replace('.tsv','')


ase_file = working_dir + '/data/ase/complete_ase_phased.all_leads.tsv'
#ase_file = working_dir + '/data/ase/complete_ase_phased.all_hipsci_ipsc_bulk_leads.tsv'
ase_short_name = os.path.basename(ase_file).split('.')[1]
n_index_cols = 2
#
# factor_file = working_dir + '/data/clustering/sce_merged_afterqc_filt_allexpts_exprs_20180618/cluster_means_eDF_all_data.transposed.tsv'
# outfile_pattern = working_dir + '/data/ase_env_interactions/cluster_means.{test_type}.'+ase_short_name+'.tsv'
# #outfile_pattern = working_dir + '/data/ase_env_interactions/cluster_means.{test_type}.'+ase_short_name+'.codetest.tsv'
# factors_to_test = ['0','6','10','28','30','59']
# combinations_to_test = list(combinations(factors_to_test, 2))


factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts.PCA.cluster_means.tsv'
outfile_pattern = working_dir + '/data/ase_env_interactions/PCA_and_cluster_means.{test_type}.'+ase_short_name+'.tsv'
factors_to_test = ['PC1','0','10','28','30','51','54','59']
combinations_to_test = list(combinations(factors_to_test, 2))

#factors_to_test = []
#combinations_to_test = []


# factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_PCA_20180618.tsv'
# factors_to_test = ['PC{}'.format(x) for x in range(1,11)]
# combinations_to_test = list(combinations(factors_to_test, 2))


factor_df = pd.read_csv(factor_file, sep='\t', index_col=0)

###phased dataframe
ase_df = pd.read_csv(ase_file, sep='\t', index_col=list(range(n_index_cols)))

if permute:
    outfile_pattern = outfile_pattern.replace('.tsv','.permuted.tsv')



# #in_dir = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25419_5245STDY7352551_hg19-1_2_0'

# in_dir = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25528_5245STDY7386983_hg19-1_2_0'
# donor = 'HPSI0115i-hecn_6'
# #donor = 'HPSI0514i-uenn_3'

# in_dir = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25528_5245STDY7386984_hg19-1_2_0'
# #donor = 'HPSI0214i-pelm_1'
# donor = 'HPSI0715i-meue_5'

# e_file = os.path.join(in_dir, 'expression_data.vargenefiltered.tsv')
# ase_filepattern = os.path.join(in_dir,'ase/{donor}.ase.{filetype}.tsv')

# e_df = pd.read_csv(e_file, sep='\t', index_col=0)
# pca_mat = decomposition.PCA(n_components=10).fit_transform(e_df.transpose().values)
# factor_df = pd.DataFrame(data=pca_mat, index=e_df.columns, columns = ['PC{}'.format(x) for x in range(1,pca_mat.shape[1]+1)])

# #ase_filepattern = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25419_5245STDY7352549_hg19-1_2_0/ase/HPSI0115i-paim_1.ase.{filetype}.tsv'
# #ase_filepattern = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25419_5245STDY7352549_hg19-1_2_0/ase/HPSI0714i-iudw_4.ase.{filetype}.tsv'
# #ase_filepattern = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25419_5245STDY7352551_hg19-1_2_0/ase/HPSI1113i-podx_1.ase.{filetype}.tsv'
# #ase_filepattern = working_dir + '/../singlecell_neuroseq/data/data_processed/cellranger211_count_25419_5245STDY7352551_hg19-1_2_0/ase/HPSI0114i-eipl_1.ase.{filetype}.tsv'

# alt_df = pd.read_csv(ase_filepattern.format(donor=donor, filetype='altcount'), sep='\t', index_col=0)
# total_df = pd.read_csv(ase_filepattern.format(donor=donor, filetype='totalcount'), sep='\t', index_col=0)
# ase_df = alt_df/total_df
# ase_df.columns = [x.replace('-1','') for x in ase_df.columns]

# outfile = ase_filepattern.format(donor=donor, filetype='ase_pc_interaction_tests')
# if permute:
#     outfile = outfile.replace('.tsv','.permuted.tsv')


def test_multiple_interactions_lm(ase_data, env_factor, permute=False):
    '''Rotate the env_factor to test the last one while including all others.'''
    output = pd.Series(index=['coef','pval','n_cells'])
    cells = list(set(env_factor.index) & set(ase_data.dropna().index))
    data1 = ase_data.loc[cells]
    data2 = env_factor.loc[cells]
    if permute:
        permuted_cells = [x for x in cells]
        random.shuffle(permuted_cells)
        data2 = env_factor.loc[permuted_cells]
    try:
        data1 = data1.values
        data2 = sm.add_constant(data2.values, prepend=True)
        res = sm.OLS(data1, data2).fit()
        # last column is the column to test
        pval = res.pvalues[-1]
        coef = res.params[-1]
        output['n_cells'] = len(cells)
        output['coef'] = coef
        output['pval'] = pval
    except:
        pass
    return output



def test_single_interaction_lm(ase_data, env_factor, permute=False):
    assert(env_factor.shape[1]==1)
    output = pd.Series(index=['coef','pval','n_cells'])
    cells = list(set(env_factor.index) & set(ase_data.dropna().index))
    data1 = ase_data.loc[cells]
    data2 = env_factor.loc[cells]
    if permute:
        permuted_cells = [x for x in cells]
        random.shuffle(permuted_cells)
        data2 = env_factor.loc[permuted_cells]
    try:
        data1 = data1.values
        data2 = sm.add_constant(data2.values, prepend=False)
        res = sm.OLS(data1, data2).fit()
        pval = res.pvalues[0]
        coef = res.params[0]
        output['n_cells'] = len(cells)
        output['coef'] = coef
        output['pval'] = pval
    except:
        pass
    return output



def test_quadratic_interaction_lm(ase_data, env_factor, permute=False):
    assert(env_factor.shape[1]==1)
    output = pd.Series(index=['coef','pval','n_cells'])
    cells = list(set(env_factor.index) & set(ase_data.dropna().index))
    data1 = ase_data.loc[cells]
    data2 = env_factor.loc[cells]
    data2['quadratic_term'] = data2.iloc[:,0]**2
    if permute:
        permuted_cells = [x for x in cells]
        random.shuffle(permuted_cells)
        data2 = env_factor.loc[permuted_cells]
    try:
        data1 = data1.values
        data2 = sm.add_constant(data2.values, prepend=False)
        res = sm.OLS(data1, data2).fit()
        pval = res.pvalues[1]
        coef = res.params[1]
        output['n_cells'] = len(cells)
        output['coef'] = coef
        output['pval'] = pval
    except:
        pass
    return output


def test_combined_interaction_lm(ase_data, env_factor, permute=False):
    assert(env_factor.shape[1]==3)
    output = pd.Series(index=['coef','pval','n_cells'])
    cells = list(set(env_factor.index) & set(ase_data.dropna().index))
    data1 = ase_data.loc[cells]
    data2 = env_factor.loc[cells]
    if permute:
        permuted_cells = [x for x in cells]
        random.shuffle(permuted_cells)
        data2 = env_factor.loc[permuted_cells]

    try:
        data1 = data1.values
        data2 = sm.add_constant(data2.values, prepend=False)
        res = sm.OLS(data1, data2).fit()
        pval = res.pvalues[:3]
        coef = res.params[:3]
        output['n_cells'] = len(cells)
        output['coef'] = coef[2]
#        output['pval_factor1'] = pval[0]
#        output['pval_factor2'] = pval[1]
        output['pval'] = pval[2]
    except:
        pass
    return output


def run_tests(env_factor, test_fcn, permute=False):
    test_columns = ['coef','pval','n_cells']
    test_df = pd.DataFrame(columns=test_columns, index=ase_df.index)
    test_df['index'] = test_df.index
    test_df[test_columns] = test_df['index'].apply(lambda x : test_fcn(ase_df.loc[x,:], env_factor, permute=permute))
    test_df['mean_ase'] = test_df['index'].apply(lambda x : ase_df.loc[x,:].mean())
    test_df = test_df.sort_values(by='pval')
    return test_df




#### With pseudotime as only a squared term (centred)

### NOTE - NEEDS TO BE CENTERED. WORKS FOR PC1, BUT NOT FOR CLUSTER MEANS
factor = 'PC1'
env_factor = factor_df[[factor]].applymap(lambda x : x**2)
df = run_tests(env_factor, test_single_interaction_lm, permute=permute)
df['factor'] = factor

quadratic_only_test_df = df.sort_values(by='pval')

quadratic_only_test_df.to_csv(outfile_pattern.format(test_type='quadratic_only_test'), sep='\t')
coef_df = quadratic_only_test_df.pivot(columns='factor')['coef']
coef_df.to_csv(outfile_pattern.format(test_type='quadratic_only_test.coef_matrix'), sep='\t')


# for each factor individually

list_of_dfs = []

for factor in factors_to_test:
    env_factor = factor_df[[factor]]
    df = run_tests(env_factor, test_single_interaction_lm, permute=permute)
    df['factor'] = factor
    list_of_dfs.append(df)

single_test_df = pd.concat(list_of_dfs).sort_values(by='pval')

single_test_df.to_csv(outfile_pattern.format(test_type='single_factor_test'), sep='\t')
coef_df = single_test_df.pivot(columns='factor')['coef']
coef_df.to_csv(outfile_pattern.format(test_type='single_factor_test.coef_matrix'), sep='\t')



list_of_dfs = []

for factor in factors_to_test:
    env_factor = factor_df[[factor]]
    df = run_tests(env_factor, test_quadratic_interaction_lm, permute=permute)
    df['factor'] = factor
    list_of_dfs.append(df)

quadratic_test_df = pd.concat(list_of_dfs).sort_values(by='pval')

quadratic_test_df.to_csv(outfile_pattern.format(test_type='quadratic_factor_test'), sep='\t')
coef_df = quadratic_test_df.pivot(columns='factor')['coef']
coef_df.to_csv(outfile_pattern.format(test_type='quadratic_factor_test.coef_matrix'), sep='\t')




list_of_dfs = []

for factor in factors_to_test:
    other_factors = [x for x in factors_to_test if x!=factor]
    env_factor = factor_df[other_factors + [factor]]
    df = run_tests(env_factor, test_multiple_interactions_lm, permute=permute)
    df['factor'] = factor
    list_of_dfs.append(df)

multi_test_df = pd.concat(list_of_dfs).sort_values(by='pval')

multi_test_df.to_csv(outfile_pattern.format(test_type='multi_factor_test'), sep='\t')
coef_df = multi_test_df.pivot(columns='factor')['coef']
coef_df.to_csv(outfile_pattern.format(test_type='multi_factor_test.coef_matrix'), sep='\t')


#### With pseudotime as the only covariate
list_of_dfs = []

covariate_factor = 'PC1'

for factor in [x for x in factors_to_test if x!=covariate_factor]:
    other_factors = [covariate_factor]
    env_factor = factor_df[other_factors + [factor]]
    df = run_tests(env_factor, test_multiple_interactions_lm, permute=permute)
    df['factor'] = factor
    list_of_dfs.append(df)

single_covariate_test_df = pd.concat(list_of_dfs).sort_values(by='pval')

single_covariate_test_df.to_csv(outfile_pattern.format(test_type='pseudotime_covariate_test'), sep='\t')
coef_df = multi_test_df.pivot(columns='factor')['coef']
coef_df.to_csv(outfile_pattern.format(test_type='pseudotime_covariate_test.coef_matrix'), sep='\t')






# for each pair of factors as a combination

list_of_dfs = []

if len(combinations_to_test)>0:
    for factor1,factor2 in combinations_to_test:
        env_factor = factor_df[[factor1,factor2]]
        env_factor['interaction'] = env_factor[factor1] * env_factor[factor2]
        df = run_tests(env_factor, test_combined_interaction_lm, permute=permute)
        df['factor1'] = factor1
        df['factor2'] = factor2
        list_of_dfs.append(df)
    comb_test_df = pd.concat(list_of_dfs).sort_values(by='pval')


comb_test_df.to_csv(outfile_pattern.format(test_type='combined_factor_test'), sep='\t')
