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

# ## Loading HapMap meta-data

# +
from collections import defaultdict
from pprint import pprint
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import xarray as xr
import sgkit as sg
from sgkit.io import plink

data = plink.read_plink(path='hapmap10_auto_noofs_ld', fam_sep='\t')
# -

data

f = open('relationships_w_pops_041510.txt')
pop_ind = defaultdict(list)
f.readline()  # header
for line in f:
    toks = line.rstrip().split('\t')
    fam_id = toks[0]
    ind_id = toks[1]
    pop = toks[-1]
    pop_ind[pop].append((fam_id, ind_id))

pops = list(pop_ind.keys())


def assign_cohort(pops, pop_ind, sample_family_id, sample_id):
    cohort = []
    for fid, sid in zip(sample_family_id, sample_id):
        processed = False
        for i, pop in enumerate(pops):
            if (fid, sid) in pop_ind[pop]:
                processed = True
                cohort.append(i)
                break
        if not processed:
            raise Exception(f'Not processed {fid}, {sid}')
    return cohort


cohort = assign_cohort(pops, pop_ind, data.sample_family_id.values, data.sample_id.values)

data['sample_cohort'] = xr.DataArray(
    cohort, dims='samples')

# # monomorphic positions per pop

cohort_allele_frequency = sg.cohort_allele_frequencies(data)['cohort_allele_frequency'].values

monom = {}
for i, pop in enumerate(pops):
    monom[pop] = len(list(filter(lambda x: x, np.isin(cohort_allele_frequency[:, i, 0], [0, 1]))))
pprint(monom)

# # MAF

mafs = {}
for i, pop in enumerate(pops):
    min_freqs = map(
        lambda x: x if x < 0.5 else 1 - x,
        filter(
            lambda x: x not in [0, 1],
            cohort_allele_frequency[:, i, 0]))
    mafs[pop] = pd.Series(min_freqs)

maf_plot, maf_ax = plt.subplots(nrows=2, sharey=True)
mafs['YRI'].hist(ax=maf_ax[0], bins=50)
maf_ax[0].set_title('*YRI*')
mafs['JPT'].hist(ax=maf_ax[1], bins=50)
maf_ax[1].set_title('*JPT*')
maf_ax[1].set_xlabel('MAF')

# # Fst

fst = sg.Fst(data)

fst = fst.assign_coords({"cohorts_0": pops, "cohorts_1": pops})

remove_nan = lambda data: filter(lambda x: not np.isnan(x), data)
ceu_chb = pd.Series(remove_nan(fst.stat_Fst.sel(cohorts_0='CEU', cohorts_1='CHB').values))
chb_chd = pd.Series(remove_nan(fst.stat_Fst.sel(cohorts_0='CHB', cohorts_1='CHD').values))

ceu_chb.describe()

chb_chd.describe()

mean_fst = {}
for i, pop_i in enumerate(pops):
    for j, pop_j in enumerate(pops):
        if j <= i:
            continue
        pair_fst = pd.Series(remove_nan(fst.stat_Fst.sel(cohorts_0=pop_i, cohorts_1=pop_j).values))
        mean = pair_fst.mean()
        mean_fst[(pop_i, pop_j)] = mean

min_pair = min(mean_fst.values())
max_pair = max(mean_fst.values())

sns.set_style("white")
num_pops = len(pops)
arr = np.ones((num_pops - 1, num_pops - 1, 3), dtype=float)
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111)
for row in range(num_pops - 1):
    pop_i = pops[row]
    for col in range(row + 1, num_pops):
        pop_j = pops[col]
        val = mean_fst[(pop_i, pop_j)]
        norm_val = (val - min_pair) / (max_pair - min_pair)
        ax.text(col - 1, row, '%.3f' % val, ha='center')
        if norm_val == 0.0:
            arr[row, col - 1, 0] = 1
            arr[row, col - 1, 1] = 1
            arr[row, col - 1, 2] = 0
        elif norm_val == 1.0:
            arr[row, col - 1, 0] = 1
            arr[row, col - 1, 1] = 0
            arr[row, col - 1, 2] = 1
        else:
            arr[row, col - 1, 0] = 1 - norm_val
            arr[row, col - 1, 1] = 1
            arr[row, col - 1, 2] = 1
ax.imshow(arr, interpolation='none')
ax.set_title('Multilocus Pairwise FST')
ax.set_xticks(range(num_pops - 1))
ax.set_xticklabels(pops[1:])
ax.set_yticks(range(num_pops - 1))
ax.set_yticklabels(pops[:-1])
