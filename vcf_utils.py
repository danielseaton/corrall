from __future__ import print_function
import cyvcf2
import pandas as pd

def get_snp_genotypes(chromosome, position, samples=None):
    '''Returns a pandas DataFrame of genotypes, along with phasing status, for a specific SNP.'''

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
    
    query_string = '{chromosome}:{position}-{position}'.format(chromosome=chromosome, position=position)

    variants = [x for x in vcf(query_string)]

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

#samples = ['HPSI0516i-pebf_2', 'HPSI0516i-zujs_5', 'HPSI1116pf-peru']
#df = get_snp_genotypes(1,714439,samples)
