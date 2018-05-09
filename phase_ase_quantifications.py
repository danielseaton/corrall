from __future__ import print_function
import cyvcf2
import pandas as pd
import glob
import os
import re
import argparse

parser = argparse.ArgumentParser(description='Convert allele-specific expression from REF vs ALT to chromosome A vs chromosome B, based on phased genotypes.')
parser.add_argument('ase_file')
args = parser.parse_args()


#ase_file = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/ase_aggregated_by_donor/ffdm.ase.lowthresh.tsv'
ase_file = args.ase_file
snp_info_file = ase_file.replace('.tsv','.snp_info.tsv')

outfile = ase_file.replace('thresh.tsv','thresh.phased.tsv')

print('Phasing ASE quantifications from {}'.format(ase_file))

vcffile_pattern = '/hps/nobackup/hipsci/scratch/genotypes/imputed/REL-2018-01/Full_Filtered/hipsci.wec.gtarray.HumanCoreExome.imputed_phased.20180102.genotypes.chr.*.norm.renamed.recode.vcf.gz'
vcf_files = glob.glob(vcffile_pattern)


ase_df = pd.read_csv(ase_file, sep='\t', index_col=0)
snp_df = pd.read_csv(snp_info_file, sep='\t', index_col=0)
sorted_snp_list = snp_df.sort_values(by=['contig','position']).index
ase_df = ase_df.loc[sorted_snp_list, :]

excluded_snps = []
flipped_snps = []

#get sample list
vcf = cyvcf2.VCF(vcf_files[0])
samples = vcf.samples

#match donor with sample
donor_id = os.path.basename(ase_file).split('.')[0]
regex_string = 'HPSI[0-9]{{4}}i-{donor}_[0-9]+'.format(donor=donor_id)
matching_samples = [x for x in samples if re.search(regex_string, x)]
sample = matching_samples[0]


for vcffile in vcf_files:
    
    vcf =  cyvcf2.VCF(vcffile, samples=[sample])
    
    for var in vcf:
                    
        ID = var.ID
        if ID not in ase_df.index:
            #not a relevant SNP
            continue
            
        genotypes = var.genotypes[0]
        if not genotypes[2] or not (sum(genotypes[:2])==1):
            #no phase information, or homozygous SNP - add SNP to excluded list
            excluded_snps.append(var.ID)
            continue
            
        if genotypes[0]==1:
            # 1|0 - swap phase so that all alternate alleles are on one chromosome
            flipped_snps.append(ID)

print('{} SNPs filtered out'.format(len(excluded_snps)))
print('{} SNPs with switched ASE ratios to match phases across SNPs'.format(len(flipped_snps)))

ase_df = ase_df.drop(excluded_snps)
ase_df.loc[flipped_snps,:] = 1.0 - ase_df.loc[flipped_snps,:]

ase_df.to_csv(outfile, sep='\t')
