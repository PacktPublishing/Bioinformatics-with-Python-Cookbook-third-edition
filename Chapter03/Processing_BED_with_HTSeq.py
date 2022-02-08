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

from collections import defaultdict
import re
import HTSeq

lct_bed = HTSeq.BED_Reader('LCT.bed')

# +
feature_types = defaultdict(int)

for rec in lct_bed:
    last_rec = rec
    feature_types[re.search('([A-Z]+)', rec.name).group(0)] += 1

print(feature_types)

#Code specific to this dataset, document
# -

print(last_rec)
print(last_rec.name)
print(type(last_rec))
interval = last_rec.iv
print(interval)
print(type(interval))

# +
print(interval.chrom, interval.start, interval.end)
print(interval.strand)
print(interval.length)
print(interval.start_d)
print(interval.start_as_pos)
print(type(interval.start_as_pos))

#talk about overlaps

# -

exon_start = None
exon_end = None
sizes = []
for rec in lct_bed:
    if not rec.name.startswith('CCDS'):
        continue
    interval = rec.iv
    exon_start = min(interval.start, exon_start or interval.start)
    exon_end = max(interval.length, exon_end or interval.end)
    sizes.append(interval.length)
sizes.sort()
print("Num exons: %d / Begin: %d / End %d" % (len(sizes), exon_start, exon_end))
print("Smaller exon: %d / Larger exon: %d / Mean size: %.1f" % (sizes[0], sizes[-1], sum(sizes)/len(sizes)))


