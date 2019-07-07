from context import cellsnp_utils

cellsnp_vcf_file = './data/cellsnp_output_example.vcf'

alt_df, total_df = cellsnp_utils.read_cellsnp_vcf(cellsnp_vcf_file)

