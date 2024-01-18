import random
import numpy as np
import pandas as pd

def read_hg38_chrom_sizes(hg38_chrom_sizes_file='hg38.chrom.sizes'):
    chrom_sizes = pd.read_csv(hg38_chrom_sizes_file, header=None, sep='\t')
    chrom_sizes.columns = ['chrom', 'size']
    chrom_list = [f'chr{i}' for i in list(range(23)) + ['X', 'Y']]
    chrom_sizes = chrom_sizes[[x in chrom_list for x in chrom_sizes['chrom']]].copy()
    return dict(zip(chrom_sizes['chrom'], chrom_sizes['size']))

def get_sample_freqs(chrom_sizes: dict):
    total_size = sum(list(chrom_sizes.values()))
    sample_freqs = {}
    for k, v in chrom_sizes.items():
        sample_freqs[k] = v / total_size
    return sample_freqs

def sample_random_ranges(chrom_sizes: dict, sample_freqs: dict, n=1, seq_len=1000):
    chrom_choices = np.random.choice(list(sample_freqs.keys()), n, list(sample_freqs.values()))
    upper_bounds = [chrom_sizes[chrom] - seq_len for chrom in chrom_choices]
    starts = np.random.randint(0, upper_bounds)
    ends = starts + seq_len
    return list(zip(chrom_choices, starts, ends))

# Example usage:
# chrom_sizes = read_hg38_chrom_sizes()
# sample_freqs = get_sample_freqs(chrom_sizes)
# sample_random_ranges(chrom_sizes, sample_freqs, 100, 18000)
