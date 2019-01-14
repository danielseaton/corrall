import pandas as pd

def read_csv_by_idx(*args, chunksize=500, idx_list=None, **kwargs):
    ''' Read in a dataframe in chunks and keep only a subset 
    of rows. Takes arguments as usually passed to pandas.read_csv,
    plus idx_list (required) and chunksize (optional - default is
    500). '''
    if idx_list is None:
        raise(ValueError(''' idx_list must be specified. '''))
    list_of_dfs = []
    for df in pd.read_csv(*args, **kwargs, chunksize=chunksize):
        idx_subset = list(set(idx_list) & set(df.index))
        if len(idx_subset) > 0:
            df = df.loc[idx_subset,:]
            list_of_dfs.append(df)
    out_df = pd.concat(list_of_dfs, axis=0)
    return out_df



working_dir = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff'

alleliccounts_file = working_dir + '/data/ase/subset_complete_ase_phased.allelic_counts.all_leads.tsv'
allelicfractions_file = working_dir + '/data/ase/subset_complete_ase_phased.allelic_fractions.all_leads.tsv'
totalcounts_file = working_dir + '/data/ase/subset_complete_ase_phased.total_counts.all_leads.tsv'

factor_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_pseudotimeandmodules_20180618.tsv'
metadata_file = working_dir + '/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'


# metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
# allelic_df = pd.read_csv(alleliccounts_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
# total_df = pd.read_csv(totalcounts_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
# ase_df = pd.read_csv(allelicfractions_file, sep='\t', index_col=list(range(n_index_cols)), nrows=n_genes)
factor_df = pd.read_csv(factor_file, sep='\t', index_col=0)

selected_rows = list(factor_df.index)[::200]

df1 = factor_df.loc[selected_rows, :]

df2 = read_csv_by_idx(factor_file, idx_list=selected_rows, chunksize=500, sep='\t', index_col=0)
df2 = df2.loc[df1.index, :]
