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

# You just need to download this ~28 MB file only once

# !rm -f SRR003265.filt.fastq.gz 2>/dev/null
# !wget -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265.filt.fastq.gz

# # The recipe

# +
from collections import defaultdict
import gzip

import seaborn as sns
import matplotlib.pyplot as plt

from Bio import SeqIO
# -

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
rec = next(recs)
print(rec.id, rec.description, rec.seq)
print(rec.letter_annotations)

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
cnt = defaultdict(int)
for rec in recs:
    for letter in rec.seq:
        cnt[letter] += 1
tot = sum(cnt.values())
for letter, cnt in cnt.items():
    print('%s: %.2f %d' % (letter, 100 * cnt / tot, cnt))

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='UTF-8'), 'fastq')
n_cnt = defaultdict(int)
for rec in recs:
    for i, letter in enumerate(rec.seq):
        pos = i + 1
        if letter == 'N':
            n_cnt[pos] += 1
seq_len = max(n_cnt.keys())
positions = range(1, seq_len + 1)
fig, ax = plt.subplots(figsize=(16, 9), tight_layout=True, dpi=300)
fig.suptitle('Number of N calls as a function of the distance from the start of the sequencer read', fontsize='xx-large')
ax.plot(positions, [n_cnt[x] for x in positions])
ax.set_xlim(1, seq_len)
ax.set_xlabel('Read distance', fontsize='xx-large')
ax.set_ylabel('Number of N Calls', fontsize='xx-large')
fig.savefig('n_calls.png')

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
cnt_qual = defaultdict(int)
for rec in recs:
    for i, qual in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25:
            continue
        cnt_qual[qual] += 1
tot = sum(cnt_qual.values())
for qual, cnt in cnt_qual.items():
    print('%d: %.2f %d' % (qual, 100. * cnt / tot, cnt))

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
qual_pos = defaultdict(list)
for rec in recs:
    for i, qual in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25 or qual == 40:
            continue
        pos = i + 1
        qual_pos[pos].append(qual)
vps = []
poses = list(qual_pos.keys())
poses.sort()
for pos in poses:
    vps.append(qual_pos[pos])
fig, ax = plt.subplots(figsize=(16,9), dpi=300, tight_layout=True)
sns.boxplot(data=vps, ax=ax)
ax.set_xticklabels([str(x) for x in range(26, max(qual_pos.keys()) + 1)])
ax.set_xlabel('Read distance', fontsize='xx-large')
ax.set_ylabel('PHRED score', fontsize='xx-large')
fig.suptitle('Distribution of PHRED scores as a function of read distance', fontsize='xx-large')
fig.savefig('phred.png')

# # There is more...

# ## Do this to download the paired end data

# Be careful as this will be 1GB of data (and fully optional)

# !rm -f SRR003265_1.filt.fastq.gz 2>/dev/null
# !rm -f SRR003265_2.filt.fastq.gz 2>/dev/null
# !wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265_1.filt.fastq.gz
# !wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/sequence_read/SRR003265_2.filt.fastq.gz

# +
f1 = gzip.open('SRR003265_1.filt.fastq.gz', 'rt', encoding='utf8')
f2 = gzip.open('SRR003265_2.filt.fastq.gz', 'rt', encoding='utf8')
recs1 = SeqIO.parse(f1, 'fastq')
recs2 = SeqIO.parse(f2, 'fastq')
cnt = 0
for rec1, rec2 in zip(recs1, recs2):
    cnt +=1

print('Number of pairs: %d' % cnt)
# -


