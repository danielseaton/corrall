import pandas as pd
import limix


metadata_df = pd.read_csv('./data/mapping_donors_to_ase_files.test.small.tsv', sep='\t')

test_idx = 0
filepath = metadata_df['filepath'].iloc[test_idx]
donor = metadata_df['donor_long_id'].iloc[test_idx]

aseout_df = pd.read_csv(filepath, sep='\t')
snp_list = aseout_df['variantID'].tolist()

#snp_list = [x for x in snp_list if x.startswith('20')]

geno_prefix = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Plink/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.norm.renamed.recode.vcf.gz'
bim,fam,bed = limix.io.plink.read(geno_prefix,verbose=False)
bim = bim.set_index('snp', drop=False)
fam = fam.set_index('iid')


snp_idx = bim.loc[snp_list,'i'].tolist()
sample_idx = fam.loc[donor, 'i']

genotypes = bed[snp_idx,sample_idx].compute()

gen_df = pd.DataFrame(data={'genotype':genotypes,'i':snp_idx,'snp_id':snp_list})
gen_df = gen_df.query('genotype==1.0')

output_snp_list = gen_df['snp_id'].tolist()

print('Total number of SNPs retained: {}'.format(len(output_snp_list)))
print('First 5 SNPs (sorted):')
print(sorted(output_snp_list)[:5])

# bim['snp'].duplicated(keep=False)
# snp_idx = bim.loc[, 'i'].tolist()

# genotypes = bed[snp_idx,sample_idx].compute()

# gen_df = pd.DataFrame(data={'genotype':genotypes,'i':snp_idx})
