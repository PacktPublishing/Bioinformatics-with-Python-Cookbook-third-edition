# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Getting the necessary data

# You just need to do this only once

# !rm -f NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam 2>/dev/null
# !rm -f NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam.bai 2>/dev/null
# !wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam
# !wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam.bai

# # The recipe

# +
#pip install pysam
from collections import defaultdict

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

import pysam
# -

bam = pysam.AlignmentFile('NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam', 'rb')

headers = bam.header
for record_type, records in headers.items():
    print (record_type)
    for i, record in enumerate(records):
        if type(record) == dict:
            print('\t%d' % (i + 1))
            for field, value in record.items():
                print('\t\t%s\t%s' % (field, value))
        else:
            print('\t\t%s' % record)

#0-based
for rec in bam:
    if rec.cigarstring.find('M') > -1 and rec.cigarstring.find('S') > -1 and not rec.is_unmapped and not rec.mate_is_unmapped:
        break
print(rec.query_name, rec.reference_id, bam.getrname(rec.reference_id), rec.reference_start, rec.reference_end)
print(rec.cigarstring)
print(rec.query_alignment_start, rec.query_alignment_end, rec.query_alignment_length)
print(rec.next_reference_id, rec.next_reference_start, rec.template_length)
print(rec.is_paired, rec.is_proper_pair, rec.is_unmapped, rec.mapping_quality)
print(rec.query_qualities)
print(rec.query_alignment_qualities)
print(rec.query_sequence)

counts = [0] * 76
for n, rec in enumerate(bam.fetch('20', 0, 10000000)):
    for i in range(rec.query_alignment_start, rec.query_alignment_end):
        counts[i] += 1
freqs = [100 * x / (n + 1) for x in counts]
fig, ax = plt.subplots(figsize=(16,9), dpi=300, tight_layout=True)
ax.plot(range(1, 77), freqs)
ax.set_xlabel('Read distance', fontsize='xx-large')
ax.set_ylabel('PHRED score', fontsize='xx-large')
fig.suptitle('Percentage of mapped calls as a function of the position from the start of the sequencer read', fontsize='xx-large')
fig.savefig('map_perc.png')

phreds = defaultdict(list)
for rec in bam.fetch('20', 0, None):
    for i in range(rec.query_alignment_start, rec.query_alignment_end):
        phreds[i].append(rec.query_qualities[i])

maxs = [max(phreds[i]) for i in range(76)]
tops = [np.percentile(phreds[i], 95) for i in range(76)]
medians = [np.percentile(phreds[i], 50) for i in range(76)]
bottoms = [np.percentile(phreds[i], 5) for i in range(76)]
medians_fig = [x - y for x, y in zip(medians, bottoms)]
tops_fig = [x - y for x, y in zip(tops, medians)]
maxs_fig = [x - y for x, y in zip(maxs, tops)]

fig, ax = plt.subplots(figsize=(16,9),dpi=300, tight_layout=True)
ax.stackplot(range(1, 77), (bottoms, medians_fig, tops_fig, maxs_fig))
ax.plot(range(1, 77), maxs, 'k-')
ax.set_xlabel('Read distance', fontsize='xx-large')
ax.set_ylabel('PHRED score', fontsize='xx-large')
fig.suptitle('Distribution of PHRED scores as a function of the position in the read', fontsize='xx-large')
fig.savefig('phred2.png')
