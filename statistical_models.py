import re
import pandas as pd
import numpy as np
import os
import scipy.stats
import random
import argparse
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn import decomposition
from itertools import combinations
import limix

def test_interaction_limix_glmm(alt_data, total_data, env_factor, permute=False):
    output = pd.Series(index=['coeff','pval','n_cells'])
    cells = list(set(env_factor.index) & set(alt_data.dropna().index))
    successes = alt_data.loc[cells].values
    trials = total_data.loc[cells].values
    y = (successes, trials)
    exog = env_factor.loc[cells].values.reshape((len(cells),1))
    if permute:
        permuted_cells = [x for x in cells]
        random.shuffle(permuted_cells)
        exog = env_factor.loc[permuted_cells].values.reshape((len(permuted_cells),1))
#    try:
    glm = limix.qtl.glmm.qtl_test_glmm(exog, y, lik='binomial', K=np.identity(len(cells)), verbose=False)
    pval = glm.getPv()[0]
    coeff = glm.getBetaSNP()[0]
    output['n_cells'] = len(cells)
    output['coeff'] = coeff
    output['pval'] = pval
#    except:
#        pass
    return output


def test_single_interaction_lmm(ase_data, env_factor, permute=False):
    ### TODO introduce permutations
    assert(not permute)
    output = pd.Series(index=['coef','pval','n_cells'])
    ase_data.name = 'ASE'
    env_factor.columns = ['factor']
    try:
        data = pd.concat([env_factor,donor_data,ase_data],axis=1,join='inner').dropna()
        n_cells = data.shape[0]
        md = smf.mixedlm("ASE ~ factor", data, groups=data["donor"])
        mdf = md.fit()
        pval = mdf.pvalues['factor']
        coef = mdf.params['factor']
        if mdf.converged:
            output['n_cells'] = n_cells
            output['coef'] = coef
            output['pval'] = pval
    except:
        pass
    return output


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



def run_tests(ase_df, env_factor, test_fcn, permute=False):
    test_columns = ['coef','pval','n_cells']
    test_df = pd.DataFrame(columns=test_columns, index=ase_df.index)
    test_df['index'] = test_df.index
    test_df[test_columns] = test_df['index'].apply(lambda x : test_fcn(ase_df.loc[x,:], env_factor, permute=permute))
    test_df['mean_ase'] = test_df['index'].apply(lambda x : ase_df.loc[x,:].mean())
    test_df = test_df.sort_values(by='pval')
    return test_df
