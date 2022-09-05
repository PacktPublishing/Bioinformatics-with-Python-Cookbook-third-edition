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

# +
import gzip

import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO, SeqUtils
# -

# !rm -f atroparvus.fa.gz gambiae.fa.gz 2>/dev/null
# !wget https://vectorbase.org/common/downloads/Current_Release/AgambiaePEST/fasta/data/VectorBase-57_AgambiaePEST_Genome.fasta -O gambiae.fa
# !gzip -9 gambiae.fa
# !wget https://vectorbase.org/common/downloads/Current_Release/AatroparvusEBRO/fasta/data/VectorBase-57_AatroparvusEBRO_Genome.fasta -O atroparvus.fa
# !gzip -9 atroparvus.fa

gambiae_name = 'gambiae.fa.gz'
atroparvus_name = 'atroparvus.fa.gz'

recs = SeqIO.parse(gzip.open(gambiae_name, 'rt', encoding='utf-8'), 'fasta')
for rec in recs:
    print(rec.description)
#Do not do this with atroparvus

recs = SeqIO.parse(gzip.open(gambiae_name, 'rt', encoding='utf-8'), 'fasta')
chrom_Ns = {}
chrom_sizes = {}
for rec in recs:
    if rec.description.find('supercontig') > -1:
        continue
    print(rec.description, rec.id, rec)
    chrom = rec.id.split('_')[1]
    if chrom in ['UNKN']:#, 'Y_unplaced']:
        continue
    chrom_Ns[chrom] = []
    on_N = False
    curr_size = 0
    for pos, nuc in enumerate(rec.seq):
        if nuc in ['N', 'n']:
            curr_size += 1
            on_N = True
        else:
            if on_N:
                chrom_Ns[chrom].append(curr_size)
                curr_size = 0
            on_N = False
    if on_N:
        chrom_Ns[chrom].append(curr_size)
    chrom_sizes[chrom] = len(rec.seq)

for chrom, Ns in chrom_Ns.items():
    size = chrom_sizes[chrom]
    if len(Ns) > 0:
        max_Ns = max(Ns)
    else:
        max_Ns = 'NA'
    print(f'{chrom} ({size}): %Ns ({round(100 * sum(Ns) / size, 1)}), num Ns: {len(Ns)}, max N: {max_Ns}')

# ## Atroparvus super-contigs

recs = SeqIO.parse(gzip.open(atroparvus_name, 'rt', encoding='utf-8'), 'fasta')
sizes = []
size_N = []
for rec in recs:
    size = len(rec.seq)
    sizes.append(size)
    count_N = 0
    for nuc in rec.seq:
        if nuc in ['n', 'N']:
            count_N += 1
    size_N.append((size, count_N / size))

print(len(sizes), np.median(sizes), np.mean(sizes), max(sizes), min(sizes),
      np.percentile(sizes, 10), np.percentile(sizes, 90))

small_split = 4800
large_split = 540000
fig, axs = plt.subplots(1, 3, figsize=(16, 9), dpi=300, squeeze=False, sharey=True)
xs, ys = zip(*[(x, 100 * y) for x, y in size_N if x <= small_split])
axs[0, 0].plot(xs, ys, '.')
xs, ys = zip(*[(x, 100 * y) for x, y in size_N if x > small_split and x <= large_split])
axs[0, 1].plot(xs, ys, '.')
axs[0, 1].set_xlim(small_split, large_split)
xs, ys = zip(*[(x, 100 * y) for x, y in size_N if x > large_split])
axs[0, 2].plot(xs, ys, '.')
axs[0, 0].set_ylabel('Fraction of Ns', fontsize=12)
axs[0, 1].set_xlabel('Contig size', fontsize=12)
fig.suptitle('Fraction of Ns per contig size', fontsize=26)
fig.savefig('frac.png')


