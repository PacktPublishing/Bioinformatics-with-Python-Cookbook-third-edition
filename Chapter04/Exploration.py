# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import gzip
import pickle
import random

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import scatter_matrix

# %matplotlib inline
# -

fit = np.load(gzip.open('balanced_fit.npy.gz', 'rb'))
ordered_features = np.load(open('ordered_features', 'rb'), allow_pickle=True)
num_features = len(ordered_features)
fit_df = pd.DataFrame(fit, columns=ordered_features + ['pos', 'error'])
num_samples = 80
del fit

fig,ax = plt.subplots(figsize=(16,9))
_ = fit_df.hist(column=ordered_features, ax=ax)

fit_df['MeanDP'] = fit_df['DP'] / 80
fig, ax = plt.subplots()
_ = ax.hist(fit_df[fit_df['MeanDP']<50]['MeanDP'], bins=100)

errors_df = fit_df[fit_df['error'] == 1]
ok_df = fit_df[fit_df['error'] == 0]

ok_qual_above_df = ok_df[ok_df['QUAL']>0.005]
errors_qual_above_df = errors_df[errors_df['QUAL']>0.005]
print(ok_df.size, errors_df.size, ok_qual_above_df.size, errors_qual_above_df.size)
print(ok_qual_above_df.size / ok_df.size, errors_qual_above_df.size / errors_df.size)

ok_qd_above_df = ok_df[ok_df['QD']>0.05]
errors_qd_above_df = errors_df[errors_df['QD']>0.05]
print(ok_df.size, errors_df.size, ok_qd_above_df.size, errors_qd_above_df.size)
print(ok_qd_above_df.size / ok_df.size, errors_qd_above_df.size / errors_df.size)

not_bad_area_errors_df = errors_df[(errors_df['QUAL']<0.005)&(errors_df['QD']<0.05)]
_ = scatter_matrix(not_bad_area_errors_df[['FS', 'ReadPosRankSum', 'MQ', 'HRun']], diagonal='kde', figsize=(16, 9), alpha=0.02)

not_bad_area_ok_df = ok_df[(ok_df['QUAL']<0.005)&(ok_df['QD']<0.05)]
_ = scatter_matrix(not_bad_area_ok_df[['FS', 'ReadPosRankSum', 'MQ', 'HRun']], diagonal='kde', figsize=(16, 9), alpha=0.02)

all_fit_df = pd.DataFrame(np.load(gzip.open('feature_fit.npy.gz', 'rb')), columns=ordered_features + ['pos', 'error'])
potentially_good_corner_df = all_fit_df[(all_fit_df['QUAL']<0.005)&(all_fit_df['QD']<0.05)]
all_errors_df=all_fit_df[all_fit_df['error'] == 1]
print(len(all_fit_df), len(all_errors_df), len(all_errors_df) / len(all_fit_df))

potentially_good_corner_errors_df = potentially_good_corner_df[potentially_good_corner_df['error'] == 1]
print(len(potentially_good_corner_df), len(potentially_good_corner_errors_df), len(potentially_good_corner_errors_df) / len(potentially_good_corner_df))
print(len(potentially_good_corner_df)/len(all_fit_df))


