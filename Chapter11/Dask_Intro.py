# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import zarr

mosquito = zarr.open('data/AG1000G-AO/2L/calldata/GT')
mosquito
zarr.array(mosquito, chunks=(1 + 48525747 // 4, 81, 2), store='data/rechunk')
# -

mosquito = zarr.open('data/rechunk')
mosquito.chunks

# +
import numpy as np
import dask.array as da

mosquito = da.from_zarr('data/rechunk')
# ^^^ load array
# -

mosquito

print(mosquito[0])

mosquito[0].compute()

mosquito.visualize()


def calc_stats(variant):
    # vvv NumPy
    variant = variant.reshape(variant.shape[0] // 2, 2)
    misses = np.equal(variant, -1)
    hets = np.not_equal(
            variant[:,0],
            variant[:,1])
    return misses
    #return num_miss, num_het



mosquito_2d = mosquito.reshape(mosquito.shape[0], mosquito.shape[1] * mosquito.shape[2])
mosquito_2d.shape
x = da.apply_along_axis(
    calc_stats, 1, mosquito_2d,
    shape=(5, 81), dtype=bool)
x


cs = da.gufunc(calc_stats, signature="(n,m)->(n)", output_dtypes=float, vectorize=True)


cs(mosquito_2d)
# mosquito_2d


help(np.sum)

