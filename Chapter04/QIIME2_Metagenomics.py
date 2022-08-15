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

# # Important: Read this!
#
# This recipe does not work with the standard conda environment.
#
# If you are in the standard environment, do this:
#
# 1. Stop Jupyter
# 2. Activate QIIME2 environment on conda
# 3. Do `jupyter serverextension enable --py qiime2 --sys-prefix`
# 4. Start Jupyter inside QIIME2 environment
#
# Note that other recipes will not work inside this environment. 

# # Check this out!
#
# This is based on on [QIIME2 Fecal Microbiota Transpant example](https://docs.qiime2.org/2018.8/tutorials/fmt/) (for the command line). You are strongly advised to read it before proceeding.
#
# There is an [amazing example](http://nbviewer.jupyter.org/gist/tkosciol/29de5198a4be81559a075756c2490fde) of using the Artifact API using the "Moving Pictures" tutorial of QIIME 2 produced by Tomasz Kościółek. I use a more convoluted approach than Tomasz's in order to go a little deeper in terms of understanding of the Python internals. That is more of a learning experience on the internals than a practical recommendatin. **My recommendation is to use Tomasz's dialect, not mine**.
#

# # Getting the data

# !wget https://data.qiime2.org/2018.8/tutorials/fmt/sample_metadata.tsv
# !wget https://data.qiime2.org/2018.8/tutorials/fmt/fmt-tutorial-demux-1-10p.qza
# !wget https://data.qiime2.org/2018.8/tutorials/fmt/fmt-tutorial-demux-2-10p.qza

# # The recipe

# +
import pandas as pd

from qiime2.metadata.metadata import Metadata
from qiime2.metadata.metadata import CategoricalMetadataColumn
from qiime2.sdk import Artifact
from qiime2.sdk import PluginManager
from qiime2.sdk import Result
# -

pm = PluginManager()
demux_plugin = pm.plugins['demux']
#demux_emp_single = demux_plugin.actions['emp_single']
demux_summarize = demux_plugin.actions['summarize']
pm.plugins

print(demux_summarize.description)
demux_summarize_signature = demux_summarize.signature
print(demux_summarize_signature.inputs)
print(demux_summarize_signature.parameters)
print(demux_summarize_signature.outputs)

# +
seqs1 = Result.load('fmt-tutorial-demux-1-10p.qza')
sum_data1 = demux_summarize(seqs1)

sum_data1.visualization

# +
seqs2 = Result.load('fmt-tutorial-demux-2-10p.qza')
sum_data2 = demux_summarize(seqs2)

print(dir(sum_data2))
print(type(sum_data2.visualization))
print(dir(sum_data2.visualization))
sum_data2.visualization
# -

#Quality control
dada2_plugin = pm.plugins['dada2']
dada2_denoise_single = dada2_plugin.actions['denoise_single']
qual_control1 = dada2_denoise_single(demultiplexed_seqs=seqs1,
                                    trunc_len=150, trim_left=13)

qual_control2 = dada2_denoise_single(demultiplexed_seqs=seqs2,
                                    trunc_len=150, trim_left=13)

metadata_plugin = pm.plugins['metadata']
metadata_tabulate = metadata_plugin.actions['tabulate']
stats_meta1 = metadata_tabulate(input=qual_control1.denoising_stats.view(Metadata))
stats_meta1.visualization

stats_meta2 = metadata_tabulate(input=qual_control2.denoising_stats.view(Metadata))
stats_meta2.visualization

# +
ft_plugin = pm.plugins['feature-table']
ft_merge = ft_plugin.actions['merge']
ft_merge_seqs = ft_plugin.actions['merge_seqs']
ft_summarize = ft_plugin.actions['summarize']
ft_tab_seqs = ft_plugin.actions['tabulate_seqs']

table_merge = ft_merge(tables=[qual_control1.table, qual_control2.table])
seqs_merge = ft_merge_seqs(data=[qual_control1.representative_sequences, qual_control2.representative_sequences])
# -

ft_sum = ft_summarize(table=table_merge.merged_table)
ft_sum.visualization

tab_seqs = ft_tab_seqs(data=seqs_merge.merged_data)
tab_seqs.visualization


