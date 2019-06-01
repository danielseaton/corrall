from __future__ import print_function
import cyvcf2
import allel
import pandas as pd
import numpy as np

def get_snp_genotypes(chromosome, position, samples=None):
    '''Returns a pandas DataFrame of genotypes, along with phasing status, for a specific SNP.

    >>> samples = ['HPSI0516i-pebf_2', 'HPSI0516i-zujs_5', 'HPSI1116pf-peru']
    >>> df = get_snp_genotypes(1,714439,samples)
    >>> df
                      chrA  chrB  phased
    HPSI0516i-pebf_2     0     0    True
    HPSI0516i-zujs_5     0     0    True
    HPSI1116pf-peru      0     0    True
    '''

    vcf_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.{chromosome}.norm.renamed.recode.vcf.gz'.format(chromosome=chromosome)

    if samples is None:
        vcf = cyvcf2.VCF(vcf_file)
        samples = vcf.samples
    else:
        if len(samples) != len(set(samples)):
            raise(ValueError('Duplicated samples in input list'))
        vcf = cyvcf2.VCF(vcf_file, samples=samples)
        for sample in samples:
            if sample not in vcf.samples:
                raise(KeyError('{} not in vcf'.format(sample)))
        # reorder samples to match order in which they'll be give by the vcf object
        samples = vcf.samples
    
    query_string = '{chromosome}:{position}-{position}'.format(chromosome=chromosome, position=position)

    variants = [x for x in vcf(query_string)]

    # only keep SNPs
    variants = [x for x in variants if x.is_snp]

    if len(variants) > 1:
        error_message = '''Input vcf file contains more than one SNP at position {chromosome}:{position}'''.format(chromosome=chromosome,
                                                                                                               position=position)
        raise(ValueError(error_message))
    if len(variants) == 0:
        error_message = '''Input vcf file has no SNP at position {chromosome}:{position}'''.format(chromosome=chromosome,
                                                                                                   position=position)
        raise(ValueError(error_message))
    
    var = variants[0]
    
    genotype_df = pd.DataFrame(index=samples, columns=['chrA','chrB','phased'], data=var.genotypes)

    return genotype_df


def get_het_snp_phase_dataframe(snp_df, samples):
    '''Get phase information for SNPs across samples.
    
    Each heterozygous SNP is encoded as 0 (chrA) or 1 (chrB). Homozygous samples
    are encoded as NaN.'''
    
    for colname in ['chrom','pos']:
        if colname not in snp_df.columns:
            raise(ValueError('{} not present in the input SNP dataframe.'.format(colname)))
    
    if any(snp_df.index.duplicated()):
        print('Warning: duplicated indices in snp_df to get_het_snp_phase_dataframe.')
        snp_df = snp_df[~snp_df.index.duplicated(keep='first')]

    out_df = pd.DataFrame(index=snp_df.index, columns=samples)
    
    failed = []
    
    for idx in snp_df.index:
        try:
            chrom, pos = snp_df.loc[idx, ['chrom', 'pos']]
            # get genotypes across donors
            genotype_df = get_snp_genotypes(chrom, pos, samples=samples)
            # encode phasing as NaN if not heterozygous, otherwise 1 if chrB, 0 if chrA
            results = genotype_df.apply(lambda x: np.nan if x[['chrA','chrB']].sum()!=1 else x['chrB'],axis=1)
            out_df.loc[idx,samples] = results.loc[samples]
        except Exception as e:
            print('{} failed. {}'.format(idx, e))
            failed.append(idx)
            pass

    return out_df

def filter_to_het_snps(snp_list, sample):
    '''Takes a list of SNPs and a sample ID, and returns the subset of SNPs heterozygous in that sample'''
    
    chromosome_list = [str(x) for x in range(1,23)] + ['X']
    output_snp_list = []
    
    for chromosome in chromosome_list:
        vcf_file = _get_vcf_filepath(chromosome)
        
        vcf_data = allel.read_vcf(vcf_file, samples=[sample])

        gen_mat = allel.GenotypeArray(vcf_data['calldata/GT']).to_n_alt()
        gen_df = pd.DataFrame(data=gen_mat, columns=['genotype'], index=vcf_data['variants/ID'])
        gen_df['snp_id'] = gen_df.index
        
        snp_subset = list(set(snp_list) & set(gen_df.index))
        gen_df = gen_df.loc[snp_subset, :]
        
        gen_df = gen_df.query('genotype==1')
        output_snp_list += gen_df['snp_id'].tolist()
    return output_snp_list


def _get_vcf_filepath(chromosome):
    if chromosome=='X':
        vcf_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/dseaton_genotype_processing/REL-2018-01/Renamed/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.{chromosome}.norm.renamed.vcf.gz'.format(chromosome=chromosome)
    else:
        vcf_file = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Renamed/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.{chromosome}.norm.renamed.vcf.gz'.format(chromosome=chromosome)
    return vcf_file

if __name__ == "__main__":
    import doctest
    doctest.testmod()
