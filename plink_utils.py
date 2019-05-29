import pandas as pd
import limix

def filter_to_het_snps(snp_list, sample):
    '''Takes a list of SNPs and a sample ID, and returns the subset of SNPs heterozygous in that sample'''

    # load genotype data
    geno_prefix = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Plink/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.norm.renamed.recode.vcf.gz'
    bim,fam,bed = limix.io.plink.read(geno_prefix,verbose=False)
    bim = bim.drop_duplicates(subset='snp').set_index('snp')
    fam = fam.set_index('iid')

    
    # restrict to SNPs for which genotypes are available
    snp_subset = list(set(snp_list) & set(bim.index))
    n_dropped = len(snp_list) - len(snp_subset)
    print('{} SNPs dropped - not found in PLINK file'.format(n_dropped))

    # collect indices, and then use indices to get genotypes
    snp_idx = bim.loc[snp_subset,'i'].tolist()
    sample_idx = fam.loc[sample, 'i']    
    genotypes = bed[snp_idx,sample_idx].compute()
    
    gen_df = pd.DataFrame(data={'genotype':genotypes,'i':snp_idx,'snp_id':snp_subset})
    gen_df = gen_df.query('genotype==1.0')
    
    output_snp_list = gen_df['snp_id'].tolist()

    return output_snp_list
