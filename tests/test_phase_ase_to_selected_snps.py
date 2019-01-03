import numpy as np
import pandas as pd

new_output = './data/complete_ase_phased.allelic_fractions.testing_example.tsv'
saved_output = new_output+'.saved'

allow_partial_index_match = False

new_df = pd.read_csv(new_output, sep='\t', index_col=[0,1])
saved_df = pd.read_csv(saved_output, sep='\t', index_col=[0,1])


if allow_partial_index_match:
    shared_cells = list(set(new_df.columns) & set(saved_df.columns))
    saved_df = saved_df[shared_cells]
    new_df = new_df[shared_cells]

max_diff = (new_df - saved_df).abs().max().max()
print('Max difference = {}'.format(max_diff))
assert(max_diff < 1e-12)
