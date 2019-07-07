import pandas as pd
import allel

def read_cellsnp_vcf(cellsnp_vcf_file):
    #DP is depth in each sample (REF+ALT), AD is number of ALT counts
    vcf_fields = ['variants/REF', 'variants/FILTER_PASS', 'variants/POS', 'variants/ID', 'variants/ALT', 'variants/QUAL', 'samples', 'variants/CHROM', 'calldata/GT', 'calldata/DP', 'calldata/AD']    
    
    vcf = allel.read_vcf(cellsnp_vcf_file, fields=vcf_fields)
    
    total_count_mat = vcf['calldata/DP']
    total_count_mat[total_count_mat == -1] = 0
    
    alt_count_mat = vcf['calldata/AD']
    alt_count_mat = alt_count_mat[:,:,0]
    alt_count_mat[alt_count_mat == -1] = 0

    variants = vcf['variants/ID']
    samples = vcf['samples']

    total_count_df = pd.DataFrame(total_count_mat, index=variants, columns=samples)
    alt_count_df = pd.DataFrame(alt_count_mat, index=variants, columns=samples)

    return alt_count_df, total_count_df

