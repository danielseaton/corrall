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
import patsy


def get_model_dataframes(ase_ds, factor_df, factor, covariate_factors):
    '''This function takes allelic data and a dataframe of factors/covariates, 
    along with strings specifying the factor to be tested and the covariates 
    to be accounted for as fixed effects.'''
    
    for element in covariate_factors + [factor]:
        if '*' in element:
            raise(ValueError(''' Error in specifying input factors and covariates:\n'''
                             ''' Use "a:b" notation rather than "a*b" notation to denote interaction terms, \n'''
                             ''' and include any individual terms as desired (e.g. "a:b + a + b" for the \n'''
                             ''' interactions and both individual terms). '''))
    
    ase_ds.name = 'ASE'
    
    # convert arguments into a formula, and convert that into a dataframe for general use
    formula = ' ASE ~ {}'.format(' + '.join([factor]+covariate_factors))
    
    # join dataframes, aligning indices
    df = factor_df.join(ase_ds, how='inner')
    
    y_ds, model_df = patsy.dmatrices(formula, df, return_type='dataframe')
    candidate_ds = model_df[factor]
    covariate_df = model_df[[x for x in model_df.columns if x!=factor]]
    
    return y_ds, candidate_ds, covariate_df

def test_limix_lmm(y_ds, candidate_ds, covariate_df, random_effect_df=None):
    model = limix.qtl.st_scan(candidate_ds, y_ds, lik='normal', K=random_effect_df, M=covariate_df, verbose=False)    
    pval = model.variant_pvalues[0]
    coeff = model.variant_effsizes[0]
    
    output = pd.Series(index=['coeff','pval','n_cells'])
    output['n_cells'] = candidate_ds.shape[0]
    output['coeff'] = coeff
    output['pval'] = pval
    return output


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
    model = limix.qtl.st_scan(exog, successes, ('binomial', trials), K=None, verbose=False)    
    pval = model.variant_pvalues[0]
    coeff = model.variant_effsizes[0]
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
