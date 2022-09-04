# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Downloading data

# https://malariagen.github.io/vector-data/ag3/download.html
# !mkdir -p data/AG1000G-AO/
# !gsutil -m rsync -r \
#         -x '.*/calldata/(AD|GQ|MQ)/.*' \
#         gs://vo_agam_release/v3/snp_genotypes/all/AG1000G-AO/ \
#         data/AG1000G-AO/ > /dev/null

# !mkdir -p data/metadata/
# !gsutil -m rsync -r gs://vo_agam_release/v3/metadata/ data/metadata/

# # BLA

# +
import numpy as np
import zarr

mosquito = zarr.open('data/AG1000G-AO')
print(mosquito.tree())

gt_2l = mosquito['/2L/calldata/GT']
gt_2l.info

dir(gt_2l)
gt_2l.shape[0]

# +
from math import ceil
from multiprocessing import Pool


def calc_stats(my_chunk):
    num_miss = np.sum(np.equal(my_chunk[:,:,0], -1), axis=1)
    num_anc_hom = np.sum(
        np.all([
            np.equal(my_chunk[:,:,0], 0),
            np.equal(my_chunk[:,:,0], my_chunk[:,:,1])], axis=0), axis=1)
    num_het = np.sum(
        np.not_equal(
            my_chunk[:,:,0],
            my_chunk[:,:,1]), axis=1)
    return num_miss, num_anc_hom, num_het


chunk_pos_size = gt_2l.chunks[0]
max_pos = gt_2l.shape[0]


intervals = []
for chunk_pos in range(ceil(max_pos / chunk_pos_size)):
    start_pos = chunk_pos * chunk_pos_size
    end_pos = min(max_pos + 1, (chunk_pos + 1) * chunk_pos_size)
    intervals.append((start_pos, end_pos))


def compute_interval(interval):
    start_pos, end_pos = interval
    my_chunk = gt_2l[start_pos:end_pos, :, :]
    num_samples = my_chunk.shape[1]
    num_miss, num_anc_hom, num_het = calc_stats(my_chunk)
    chunk_complete_data = np.sum(np.equal(num_miss, 0))
    chunk_more_anc_hom = np.sum(num_anc_hom > num_het)
    return chunk_complete_data, chunk_more_anc_hom


with Pool() as p:
    print(p)
    chunk_returns = p.map(compute_interval, intervals)
    complete_data = sum(map(lambda x: x[0], chunk_returns))
    more_anc_hom = sum(map(lambda x: x[1], chunk_returns))
    
    print(complete_data, more_anc_hom)
# -


