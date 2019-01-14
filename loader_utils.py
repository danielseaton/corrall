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
