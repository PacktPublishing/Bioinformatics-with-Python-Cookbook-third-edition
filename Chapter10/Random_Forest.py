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
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.tree import export_graphviz

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
trainning_input = samples.iloc[:,:-1]
target = samples.iloc[:,-1]

# + jupyter={"outputs_hidden": false}
clf = RandomForestClassifier(max_depth=3, n_estimators=200)

# + jupyter={"outputs_hidden": false}
clf.fit(trainning_input, target)

# + jupyter={"outputs_hidden": false}
importances = pd.Series(
    clf.feature_importances_ * 100,
    index=trainning_input.columns).sort_values(ascending=False)
importances
# -

100 * clf.score(trainning_input, target)



for test_size in [0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 0.99]:
    X_train, X_test, y_train, y_test = train_test_split(
        trainning_input, target, test_size=test_size)
    tclf = RandomForestClassifier(max_depth=3)
    tclf.fit(X_train, y_train)
    score = tclf.score(X_test, y_test)
    print(f'{1 - test_size:.1%} {score:.2%}')
# Random number generator


