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
# %matplotlib inline

from collections import defaultdict
import gzip

import numpy as np
import matplotlib.pylab as plt
# -

num_parents = 8
dp_2L = np.load(gzip.open('DP_2L.npy.gz', 'rb'))
dp_2L.shape

for i in range(num_parents):
    print(np.median(dp_2L[:,i]), np.median(dp_2L[50000:150000,i]))

window_size = 200000
parent_DP_windows = [defaultdict(list) for i in range(num_parents)]


# +
def insert_in_window(row):
    for parent in range(num_parents):
        parent_DP_windows[parent][row[-1] // window_size].append(row[parent])

insert_in_window_v = np.vectorize(insert_in_window, signature='(n)->()')
_ = insert_in_window_v(dp_2L)
# -

fig, axs = plt.subplots(2, num_parents // 2, figsize=(16, 9), sharex=True, sharey=True, squeeze=True)
for parent in range(num_parents):
    ax = axs[parent // 4][parent % 4]
    parent_data = parent_DP_windows[parent]
    ax.set_ylim(10, 40)
    ax.plot(*zip(*[(win*window_size, np.mean(lst)) for win, lst in parent_data.items()]), '.')


