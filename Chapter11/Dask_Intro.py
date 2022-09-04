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

import zarr

mosquito = zarr.open('data/AG1000G-AO/2L/calldata/GT')
mosquito
zarr.array(mosquito, chunks=(1 + 48525747 // 4, 81, 2), store='data/rechunk')

mosquito = zarr.open('data/rechunk')
mosquito.chunks

# +
import numpy as np
import dask.array as da

mosquito = da.from_zarr('data/rechunk')
#mosquito = da.from_zarr('data/AG1000G-AO/2L/calldata/GT')
# ^^^ load array
# -

mosquito

print(mosquito[0])

mosquito[0].compute()

mosquito.visualize(rankdir='TB')


def calc_stats(variant):
    variant = variant.reshape(variant.shape[0] // 2, 2)
    num_misses = np.sum(np.equal(variant, -1)) // 2
    return num_misses


mosquito_2d = mosquito.reshape(mosquito.shape[0], mosquito.shape[1] * mosquito.shape[2])
mosquito_2d.visualize(rankdir='TB')

mosquito_2d

max_pos = 10000000
stats = da.apply_along_axis(
    calc_stats, 1, mosquito_2d[:max_pos,:],
    shape=(max_pos,), dtype=np.int64)

stats.visualize('x.png',rankdir='TB')

a = stats.compute()

a


