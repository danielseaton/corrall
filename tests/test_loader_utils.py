import pandas as pd
from context import loader_utils


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

df2 = loader_utils.read_csv_by_idx(factor_file, idx_list=selected_rows, chunksize=500, sep='\t', index_col=0)
df2 = df2.loc[df1.index, :]

assert((df1==df2).all().all())
