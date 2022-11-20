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

# + jupyter={"outputs_hidden": false}
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import tree

# + [markdown] jupyter={"outputs_hidden": false}
# http://archive.ics.uci.edu/ml/datasets/breast+cancer+wisconsin+%28diagnostic%29

# + jupyter={"outputs_hidden": false}
# !wget http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data
# !wget http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.names
# -

# ## With scikit-learn

# + jupyter={"outputs_hidden": false}
f = open('breast-cancer-wisconsin.data')
w = open('clean.data', 'w')
for line in f:
    if line.find('?') > -1:
        continue
    w.write(line)
f.close()
w.close()

# + jupyter={"outputs_hidden": false}
column_names = [
    'sample_id', 'clump_thickness', 'uniformity_cell_size',
    'uniformity_cell shape', 'marginal_adhesion',
    'single_epithelial_cell_size', 'bare_nuclei',
    'bland_chromatin', 'normal_nucleoli', 'mitoses',
    'class'
]
samples = pd.read_csv('clean.data', header=None, names=column_names, index_col=0)
samples

# + jupyter={"outputs_hidden": false}
training_input = samples.iloc[:,:-1]
target = samples.iloc[:,-1].apply(lambda x: 0 if x == 2 else 1)

# + jupyter={"outputs_hidden": false}
clf = tree.DecisionTreeClassifier(max_depth=3)

# + jupyter={"outputs_hidden": false}
clf.fit(training_input, target)

# + jupyter={"outputs_hidden": false}
importances = pd.Series(
    clf.feature_importances_ * 100,
    index=training_input.columns).sort_values(ascending=False)
importances

# + jupyter={"outputs_hidden": false}
100 * clf.score(training_input, target)

# + jupyter={"outputs_hidden": false}
fig, ax = plt.subplots(1, dpi=300)
tree.plot_tree(clf,ax=ax, feature_names=training_input.columns, class_names=['Benign', 'Malignant'])
# -




