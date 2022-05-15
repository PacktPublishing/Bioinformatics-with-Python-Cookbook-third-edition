# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# ## Loading HapMap data

# +
import numpy as np
import xarray as xr
import sgkit as sg
from sgkit.io import plink

data = plink.read_plink(path='hapmap10_auto_noofs_ld', fam_sep='\t')
# -

data

print(data.dims)

variant_stats = sg.variant_stats(data)
variant_stats

variant_stats.variant_call_rate.to_series().describe()

print(type(variant_stats.variant_call_rate.to_series()))

sample_stats = sg.sample_stats(data)
sample_stats

sample_stats.sample_call_rate.to_series().hist()

data['sample_cohort'] = xr.DataArray(
    np.zeros(data.dims['samples'], dtype=np.int64),
    dims='samples')
# data["sample_cohort"] = xr.DataArray(np.repeat([0, 1], data.dims["samples"] // 2), dims="samples")

sg.cohort_allele_frequencies(data)['cohort_allele_frequency'][:,:,0].values

sg.cohort_allele_frequencies(data)['cohort_allele_frequency'][:,:,0].to_series().hist()



# # maf

cohort_allele_frequency = sg.cohort_allele_frequencies(data)['cohort_allele_frequency'].values

min_freqs = map(
    lambda x: x if x < 0.5 else 1 - x,
    filter(
        lambda x: x not in [0, 1],
        cohort_allele_frequency[:, 0, 0]))


