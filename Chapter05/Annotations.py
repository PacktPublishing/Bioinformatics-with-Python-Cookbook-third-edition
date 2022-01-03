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
#pip install gffutils
from collections import defaultdict

import gffutils
import sqlite3
# -

# !rm -f ag.db
# !wget https://vectorbase.org/common/downloads/release-55/AgambiaePEST/gff/data/VectorBase-55_AgambiaePEST.gff -O gambiae.gff
# !gzip -9 gambiae.gff

try:
    db = gffutils.create_db('gambiae.gff.gz', 'ag.db')
except sqlite3.OperationalError:
    db = gffutils.FeatureDB('ag.db')

print(list(db.featuretypes()))
for feat_type in db.featuretypes():
    print(feat_type, db.count_features_of_type(feat_type))

seqids = set()
for e in db.all_features():
    seqids.add(e.seqid)
for seqid in seqids:
    print(seqid)

num_mRNAs = defaultdict(int)
num_exons = defaultdict(int)
max_exons = 0
max_span = 0
for seqid in seqids:
    cnt = 0
    for gene in db.region(seqid=seqid, featuretype='protein_coding_gene'):
        cnt += 1
        span = abs(gene.start - gene.end) # strand
        if span > max_span:
            max_span = span
            max_span_gene = gene
        my_mRNAs = list(db.children(gene, featuretype='mRNA'))
        num_mRNAs[len(my_mRNAs)] += 1
        if len(my_mRNAs) == 0:
            exon_check = [gene]
        else:
            exon_check = my_mRNAs
        for check in exon_check:
            my_exons = list(db.children(check, featuretype='exon'))
            num_exons[len(my_exons)] += 1
            if len(my_exons) > max_exons:
                max_exons = len(my_exons)
                max_exons_gene = gene
    print(f'seqid {seqid}, number of genes {cnt}')
print('Max number of exons: %s (%d)' % (max_exons_gene.id, max_exons))
print('Max span: %s (%d)' % (max_span_gene.id, max_span))
print(num_mRNAs)
print(num_exons)


