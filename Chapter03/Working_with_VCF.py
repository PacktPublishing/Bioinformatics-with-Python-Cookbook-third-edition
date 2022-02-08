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

# !rm -f genotypes.vcf.gz 2>/dev/null
# !tabix -fh ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/vcf_with_sample_level_annotation/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5_extra_anno.20130502.genotypes.vcf.gz 22:1-17000000|bgzip -c > genotypes.vcf.gz
# !tabix -p vcf genotypes.vcf.gz

# +
from collections import defaultdict

import seaborn as sns
import matplotlib.pyplot as plt

from cyvcf2 import VCF

# +
v = VCF('genotypes.vcf.gz')
rec = next(v)
print('Variant Level information')
info = rec.INFO
for info in rec.INFO:
    print(info)

print('Sample Level information')
for fmt in rec.FORMAT:
    print(fmt)

# +
v = VCF('genotypes.vcf.gz')
samples = v.samples
print(len(samples))  # Order change

variant = next(v)
print(variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT, variant.QUAL, variant.FILTER)
print(variant.INFO)
print(variant.FORMAT)
print(variant.is_snp)

#rec.format('DP')
#rec.format('GT')

str_alleles = variant.gt_bases[0]
alleles = variant.genotypes[0][0:2]
is_phased = variant.genotypes[0][2]
print(str_alleles, alleles, is_phased)
print(variant.format('DP')[0])

# +
f = VCF('genotypes.vcf.gz')

my_type = defaultdict(int)
num_alts = defaultdict(int)

for variant in f:
    my_type[variant.var_type, variant.var_subtype] += 1
    if variant.var_type == 'snp':
        num_alts[len(variant.ALT)] += 1
print(my_type)
print(num_alts)

# +
f = VCF('genotypes.vcf.gz')

sample_dp = defaultdict(int)
for variant in f:
    if not variant.is_snp or len(variant.ALT) != 1:
        continue
    for dp in variant.format('DP'):
        #dp = int(dp)
        sample_dp[dp] += 1
# -

dps = list(sample_dp.keys())
dps.sort()
dp_dist = [sample_dp[x] for x in dps]
fig, ax = plt.subplots(figsize=(16, 9))
ax.plot(dp_dist[:50], 'r')
ax.axvline(dp_dist.index(max(dp_dist)))


