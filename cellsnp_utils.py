import pandas as pd
import allel

vcf_fields = ['variants/REF', 'variants/FILTER_PASS', 'variants/POS', 'variants/ID', 'variants/ALT', 'variants/QUAL', 'samples', 'variants/CHROM', 'calldata/GT', 'calldata/DP', 'calldata/AD']

#DP is depth in each sample (REF+ALT), AD is number of ALT counts

vcf = allel.read_vcf('./tests/data/cellsnp_output_example.vcf', fields=vcf_fields)

total_count_mat = vcf['calldata/DP']
total_count_mat[total_count_mat == -1] = 0


#this needs to be massively reshaped
alt_count_mat = vcf['calldata/AD']
alt_count_mat = alt_count_mat[:,:,0] # I think this is what I want
#alt_count_mat = alt_count_mat[alt_count_mat != -1]

alt_count_mat[alt_count_mat == -1] = 0

#samples = vcf['samples']
