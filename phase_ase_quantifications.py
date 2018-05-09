from __future__ import print_function
import cyvcf2
import pandas as pd
import glob
import os
import re
import argparse

parser = argparse.ArgumentParser(description='Convert allele-specific expression from REF vs ALT to chromosome A vs chromosome B, based on phased genotypes.')
parser.add_argument('altcounts_file')
parser.add_argument('totalcounts_file')
parser.add_argument('snp_info_file')
args = parser.parse_args()

altcounts_file = args.altcounts_file
totalcounts_file = args.totalcounts_file
snp_info_file = args.snp_info_file

alt_outfile = altcounts_file.replace('altcount.tsv','chrBcounts.phased.tsv')
total_outfile = altcounts_file.replace('altcount.tsv','totalcounts.phased.tsv')

print('Phasing ASE quantifications from {}'.format(altcounts_file))

vcffile_pattern = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.*.norm.renamed.recode.vcf.gz'
vcf_files = glob.glob(vcffile_pattern)


alt_df = pd.read_csv(altcounts_file, sep='\t', index_col=0)
total_df = pd.read_csv(totalcounts_file, sep='\t', index_col=0)
snp_df = pd.read_csv(snp_info_file, sep='\t', index_col=0)
sorted_snp_list = snp_df.sort_values(by=['contig','position']).index

alt_df = alt_df.loc[sorted_snp_list, :]
total_df = total_df.loc[sorted_snp_list, :]

excluded_snps = []
flipped_snps = []

#get sample list
vcf = cyvcf2.VCF(vcf_files[0])
available_samples = vcf.samples

#match donor with sample
donor_id = os.path.basename(altcounts_file).split('.')[0]
sample = donor_id
if sample not in available_samples:
    raise(ValueError('Sample {} not present in input vcf file.'.format(sample)))

for vcffile in vcf_files:
    
    vcf =  cyvcf2.VCF(vcffile, samples=[sample])
    
    for var in vcf:
                    
        ID = var.ID
        if ID not in alt_df.index:
            #not a relevant SNP
            continue
            
        genotypes = var.genotypes[0]
        if not genotypes[2] or not (sum(genotypes[:2])==1):
            #no phase information, or homozygous SNP - add SNP to excluded list
            excluded_snps.append(var.ID)
            continue
            
        if genotypes[0]==1:
            # 1|0 - flip quantification so it refers to chrB
            flipped_snps.append(ID)

print('{} SNPs filtered out'.format(len(excluded_snps)))
print('{} SNPs with switched ASE ratios to match phases across SNPs'.format(len(flipped_snps)))

total_df = total_df.drop(excluded_snps)
alt_df = alt_df.drop(excluded_snps)
alt_df.loc[flipped_snps,:] = total_df.loc[flipped_snps,:] - alt_df.loc[flipped_snps,:]

alt_df.to_csv(alt_outfile, sep='\t')
total_df.to_csv(total_outfile, sep='\t')
