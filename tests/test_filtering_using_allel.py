import pandas as pd
import allel


metadata_df = pd.read_csv('./data/mapping_donors_to_ase_files.test.small.tsv', sep='\t')

test_idx = 0
filepath = metadata_df['filepath'].iloc[test_idx]
sample = metadata_df['donor_long_id'].iloc[test_idx]

aseout_df = pd.read_csv(filepath, sep='\t')
snp_list = aseout_df['variantID'].tolist()

#chromosome_list = [str(x) for x in range(1,23)] + ['X']
chromosome_list = [str(x) for x in range(1,23)]

output_snp_list = []

for chromosome in chromosome_list:
    vcf_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Renamed/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.{chromosome}.norm.renamed.vcf.gz'.format(chromosome=chromosome)
    
    vcf_data = allel.read_vcf(vcf_file, samples=[sample])

    gen_mat = allel.GenotypeArray(vcf_data['calldata/GT']).to_n_alt()
    gen_df = pd.DataFrame(data=gen_mat, columns=['genotype'], index=vcf_data['variants/ID'])
    gen_df['snp_id'] = gen_df.index
    
    snp_subset = list(set(snp_list) & set(gen_df.index))
    gen_df = gen_df.loc[snp_subset, :]
    
    gen_df = gen_df.query('genotype==1')
    output_snp_list += gen_df['snp_id'].tolist()

print('Total number of SNPs retained: {}'.format(len(output_snp_list)))
print('First 5 SNPs (sorted):')
print(sorted(output_snp_list)[:5])

# bim['snp'].duplicated(keep=False)
# snp_idx = bim.loc[, 'i'].tolist()

# genotypes = bed[snp_idx,sample_idx].compute()

# gen_df = pd.DataFrame(data={'genotype':genotypes,'i':snp_idx})
