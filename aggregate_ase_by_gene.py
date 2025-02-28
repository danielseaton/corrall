import pandas as pd
import pybedtools
import argparse

parser = argparse.ArgumentParser(description='Aggregate phased allele-specific expression quantifications at the gene level.')
parser.add_argument('altcounts_file')
parser.add_argument('totalcounts_file')
parser.add_argument('snp_info_file')
parser.add_argument('gene_bedfile')
parser.add_argument('outfile_suffix', default='genelevel')
args = parser.parse_args()

altcounts_file = args.altcounts_file
totalcounts_file = args.totalcounts_file
snp_info_file = args.snp_info_file
gene_bedfile = args.gene_bedfile
outfile_suffix = args.outfile_suffix

# infer other filenames based on input file
snp_info_bedfile = snp_info_file.replace('.snp_info.tsv', '.snp_info.bed')
snp_gene_intersection_file = snp_info_file.replace('.snp_info.tsv', '.snp_gene_intersection.tsv')
alt_file_out = altcounts_file.replace('.tsv', '.{}.tsv'.format(outfile_suffix))
total_file_out = totalcounts_file.replace('.tsv', '.{}.tsv'.format(outfile_suffix))

# input ASE quantifications
alt_df = pd.read_csv(altcounts_file, sep='\t', index_col=0)
total_df = pd.read_csv(totalcounts_file, sep='\t', index_col=0)

# gene location file
gene_bed = pybedtools.BedTool(gene_bedfile)

# convert SNP info table into a bed file
snp_df = pd.read_csv(snp_info_file, sep='\t', index_col=0)
snp_df = snp_df.loc[alt_df.index, :]
# convert chr names
snp_df['chromosome'] = snp_df['contig'].astype(str).apply(lambda x: x.replace('chr', '').split('.')[0])
# snp_df['chromosome'] = snp_df['contig']
snp_df['start'] = snp_df['position'].apply(int)
snp_df['end'] = snp_df['position'].apply(int)
snp_df['score'] = '.'
snp_df['strand'] = '.'
snp_df['name'] = snp_df.index

snp_df = snp_df[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
snp_df = snp_df.sort_values(by=['chromosome', 'start'])
snp_df.to_csv(snp_info_bedfile, sep='\t', index=False, header=False)

snp_bed = pybedtools.BedTool(snp_info_bedfile)


# find intersection of snps and genes, write to file

joint_bed = snp_bed.intersect(gene_bed, wa=True, wb=True)
joint_bed.saveas(snp_gene_intersection_file)


# process SNP mapping to genes, filtering out multimapping SNPs

map_df = pd.read_csv(snp_gene_intersection_file, sep='\t', header=None)
map_df = map_df[[3, 9]]
map_df.columns = ['snp_id', 'gene_id']
map_df = map_df.drop_duplicates()
map_df = map_df.drop_duplicates(subset=['snp_id'], keep=False)
map_df = map_df.set_index('gene_id', drop=False)

# map from genes to sets of SNPs

gene_list = map_df['gene_id'].drop_duplicates().tolist()


def aggregate_counts_to_gene(df):
    gene_df = pd.DataFrame(index=gene_list, columns=df.columns)
    for gene in gene_list[:]:
        snps = map_df.loc[[gene], 'snp_id']
        
        data = df.loc[snps, :]
        gene_df.loc[gene, :] = data.sum()
    return gene_df

gene_alt_df = aggregate_counts_to_gene(alt_df)
gene_total_df = aggregate_counts_to_gene(total_df)

gene_alt_df.to_csv(alt_file_out, sep='\t')
gene_total_df.to_csv(total_file_out, sep='\t')
