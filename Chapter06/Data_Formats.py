# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# ## Data download

# +
# !wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz
# !wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz

# !wget https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3/relationships_w_pops_041510.txt
# -

# !gzip -d hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz
# !gzip -d hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz

# # Preparation

import os
from collections import defaultdict

# ## Loading HapMap meta-data

f = open('relationships_w_pops_041510.txt')
pop_ind = defaultdict(list)
f.readline()  # header
offspring = []
for l in f:
    toks = l.rstrip().split('\t')
    fam_id = toks[0]
    ind_id = toks[1]
    mom = toks[2]
    dad = toks[3]
    if mom != '0' or dad != '0':
        offspring.append((fam_id, ind_id))
    pop = toks[-1]
    pop_ind[pop].append((fam_id, ind_id))
f.close()

# ## Sub-sampling

os.system('plink2 --pedmap hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap10 --thin 0.1 --geno 0.1 --export ped')
os.system('plink2 --pedmap hapmap3_r3_b36_fwd.consensus.qc.poly --out hapmap1 --thin 0.01 --geno 0.1 --export ped')


# ## Getting only autosomal data

def get_non_auto_SNPs(map_file, exclude_file):
    f = open(map_file)
    w = open(exclude_file, 'w')
    for l in f:
        toks = l.rstrip().split('\t')
        try:
            chrom = int(toks[0])
        except ValueError:
            rs = toks[1]
            w.write('%s\n' % rs)
    w.close()


get_non_auto_SNPs('hapmap10.map', 'exclude10.txt')
get_non_auto_SNPs('hapmap1.map', 'exclude1.txt')

# !plink2 --pedmap hapmap10 --out hapmap10_auto --exclude exclude10.txt --export ped
# !plink2 --pedmap hapmap1 --out hapmap1_auto --exclude exclude1.txt --export ped


# ## Removing offspring

# !plink2 --pedmap hapmap10_auto --filter-founders --out hapmap10_auto_noofs --export ped

# ## LD-prunning

# !plink2 --pedmap hapmap10_auto_noofs --indep-pairwise 50 10 0.1 --out keep --export ped
# !plink2 --pedmap hapmap10_auto_noofs --extract keep.prune.in --out hapmap10_auto_noofs_ld --export ped

# ## Different encoding

# !plink2 --pedmap hapmap10_auto_noofs_ld --out hapmap10_auto_noofs_ld_12 --export ped 12
# !plink2 --make-bed --pedmap hapmap10_auto_noofs_ld --out hapmap10_auto_noofs_ld

# ## Single chromosome

# !plink2 --pedmap hapmap10_auto_noofs --chr 2 --out hapmap10_auto_noofs_2 --export ped
