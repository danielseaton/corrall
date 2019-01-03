import numpy as np
import pandas as pd
import corrall

exprs_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_afterqc_filt_allexpts_exprs_20180618.tsv'

qtl_filename = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/qtl_results/all_results_combined.for_snp_annotation.tsv'
qtl_file_short_name = 'all_leads'

metadata_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/sce_merged_afterqc_filt_allexpts_metadata_20180618.tsv'

assert(exprs_file.endswith('.tsv'))
outfile = exprs_file.replace('.tsv','.subset_to_het_donors.{}.tsv'.format(qtl_file_short_name))


# expression data
e_df = pd.read_csv(exprs_file, sep='\t')
e_df.index = [x.split('_')[0] for x in e_df.index]

# metadata
metadata_df = pd.read_csv(metadata_file, sep='\t')

# snp_df
qtl_df = pd.read_csv(qtl_filename,sep='\t')
qtl_df['ensembl_gene_id'] = qtl_df['feature'].apply(lambda x: x.split('_')[0])
#process for het phasing
qtl_df = qtl_df.set_index(['ensembl_gene_id','snp_id'], drop=False)
qtl_df = qtl_df.rename(columns={'snp_chromosome':'chrom','snp_position':'pos'})
# sort low to high by p-value
qtl_df = qtl_df.sort_values(by='empirical_feature_p_value')
qtl_df = qtl_df.drop_duplicates(subset=['ensembl_gene_id','snp_id'])
# qtl_df = qtl_df.head(40)
gene_list = list(e_df.index)
qtl_df = qtl_df.query('ensembl_gene_id in @gene_list')

snp_df = qtl_df


# expand expression dataframe to include snp in the index
e_df = e_df.loc[snp_df['ensembl_gene_id']]
e_df.index = snp_df.index

# donor list
with open('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/list_of_singlecell_endodiff_donors.tsv', 'r') as f:
    donor_list = [x.strip() for x in f.readlines()]

#identify heterozygous donors (this is a roundabout way of doing this)
snp_phase_df = corrall.vcf_utils.get_het_snp_phase_dataframe(snp_df, donor_list)

#for each gene,snp, mask irrelevant cells.
for index in list(set(snp_df.index) & set(e_df.index)):
    homozygous_donors = list(set(donor_list) - set(snp_phase_df.loc[index].dropna().index))
    cells_to_mask = metadata_df.query('donor_long_id in @homozygous_donors')['cell_name']
    e_df.loc[index, cells_to_mask] = np.nan


e_df.to_csv(outfile, sep='\t')
