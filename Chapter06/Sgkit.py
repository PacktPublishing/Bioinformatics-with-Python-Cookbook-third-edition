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

import os
from collections import defaultdict

# ## Loading HapMap data

# +
import numpy as np
from sgkit.io import plink

data = plink.read_plink(path='hapmap10_auto_noofs_ld', fam_sep='\t')
# -

data

print(data.dims)

print(len(data.sample_id.values))
print(data.sample_id.values)
print(data.sample_family_id.values)
print(data.sample_sex.values)

print(data.contigs)

print(len(data.variant_contig.values))
print(data.variant_contig.values)
print(data.variant_position.values)
print(data.variant_allele.values)
print(data.variant_id.values)

data.call_genotype

call_genotype = data.call_genotype.values
print(call_genotype.shape)
first_individual = call_genotype[:,0,:]
first_variant = call_genotype[0,:,:]
first_variant_of_first_individual = call_genotype[0,0,:]
print(first_variant_of_first_individual)
print(data.sample_family_id.values[0], data.sample_id.values[0])
print(data.variant_allele.values[0])


