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

# +
#import dask
#from dask.base import get_scheduler
#import dask.array as da
#
#mosquito = da.from_zarr('data/AG1000G-AO/2L/calldata/GT')
#print(get_scheduler(collections=[mosquito]).__module__) 

# +
import zarr
import dask.dataframe as dd
from dask.distributed import Client

#client = Client('127.0.0.1:8786')
client = Client()
client

# +
import numpy as np
import dask.array as da

mosquito = da.from_zarr('data/AG1000G-AO/2L/calldata/GT')
# -

mosquito

mosquito.shape[0]

mosquito = mosquito.rechunk((mosquito.shape[0]//8, 81, 2))

mosquito = mosquito.persist()

mosquito.visualize()

mosquito

mosquito.chunks


def calc_stats(my_chunk):
    num_miss = np.sum(np.equal(my_chunk[0][0][:,:,0], -1), axis=1)
    return num_miss


stats = da.blockwise(calc_stats, 'i', mosquito, 'ijk', dtype=np.uint8)

stats.visualize()

stat_results = stats.compute()

stat_results


