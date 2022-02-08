# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import os
from collections import OrderedDict

import numpy as np
import pandas as pd

import dendropy
from dendropy.calculate import popgenstat
# -

# ## Genes

# +
genes_species = OrderedDict()
my_species = ['RESTV', 'SUDV']
my_genes = ['NP', 'L', 'VP35', 'VP40']


for name in my_genes:
    gene_name = name.split('.')[0]
    char_mat = dendropy.DnaCharacterMatrix.get_from_path('%s_align.fasta' % name, 'fasta')
    genes_species[gene_name] = {}
    
    for species in my_species:
        genes_species[gene_name][species] = dendropy.DnaCharacterMatrix()
    for taxon, char_map in char_mat.items():
        species = taxon.label.split('_')[0]
        if species in my_species:
            genes_species[gene_name][species].taxon_namespace.add_taxon(taxon)
            genes_species[gene_name][species][taxon] = char_map
# -

summary = np.ndarray(shape=(len(genes_species), 4 * len(my_species)))
stats = ['seg_sites', 'nuc_div', 'taj_d', 'wat_theta']
for row, (gene, species_data) in enumerate(genes_species.items()):
    for col_base, species in enumerate(my_species):
        summary[row, col_base * 4] = popgenstat.num_segregating_sites(species_data[species])
        summary[row, col_base * 4 + 1] = popgenstat.nucleotide_diversity(species_data[species])
        summary[row, col_base * 4 + 2] = popgenstat.tajimas_d(species_data[species])
        summary[row, col_base * 4 + 3] = popgenstat.wattersons_theta(species_data[species])
columns = []
for species in my_species:
    columns.extend(['%s (%s)' % (stat, species) for stat in stats])
df = pd.DataFrame(summary, index=genes_species.keys(), columns=columns)
df # vs print(df)


# ## Genomes

def do_basic_popgen(seqs):
    num_seg_sites = popgenstat.num_segregating_sites(seqs)
    avg_pair = popgenstat.average_number_of_pairwise_differences(seqs)
    nuc_div = popgenstat.nucleotide_diversity(seqs)
    print('Segregating sites: %d, Avg pairwise diffs: %.2f, Nucleotide diversity %.6f' % (num_seg_sites, avg_pair, nuc_div))
    print("Watterson's theta: %s" % popgenstat.wattersons_theta(seqs))
    print("Tajima's D: %s" % popgenstat.tajimas_d(seqs))


#XXX change
ebov_seqs = dendropy.DnaCharacterMatrix.get_from_path(
    'trim.fasta', schema='fasta', data_type='dna')
sl_2014 = []
drc_2007 = []
ebov2007_set = dendropy.DnaCharacterMatrix()
ebov2014_set = dendropy.DnaCharacterMatrix()
for taxon, char_map in ebov_seqs.items():
    print(taxon.label)
    if taxon.label.startswith('EBOV_2014') and len(sl_2014) < 8:
        sl_2014.append(char_map)
        ebov2014_set.taxon_namespace.add_taxon(taxon)
        ebov2014_set[taxon] = char_map
    elif taxon.label.startswith('EBOV_2007'):
        drc_2007.append(char_map)
        ebov2007_set.taxon_namespace.add_taxon(taxon)
        ebov2007_set[taxon] = char_map
        #ebov2007_set.extend_map({taxon: char_map})
del ebov_seqs

# +
print('2007 outbreak:')
print('Number of individuals: %s' % len(ebov2007_set.taxon_namespace))
do_basic_popgen(ebov2007_set)

print('\n2014 outbreak:')
print('Number of individuals: %s' % len(ebov2014_set.taxon_namespace))
do_basic_popgen(ebov2014_set)
# -

print(len(sl_2014))
print(len(drc_2007))

pair_stats = popgenstat.PopulationPairSummaryStatistics(sl_2014, drc_2007)

print('Average number of pairwise differences irrespective of population: %.2f' %
      pair_stats.average_number_of_pairwise_differences)
print('Average number of pairwise differences between populations: %.2f' %
      pair_stats.average_number_of_pairwise_differences_between)
print('Average number of pairwise differences within populations: %.2f' %
      pair_stats.average_number_of_pairwise_differences_within)
print('Average number of net pairwise differences : %.2f' %
      pair_stats.average_number_of_pairwise_differences_net)
print('Number of segregating sites: %d' %
      pair_stats.num_segregating_sites)
print("Watterson's theta: %.2f" %
      pair_stats.wattersons_theta)
print("Wakeley's Psi: %.3f" % pair_stats.wakeleys_psi)
print("Tajima's D: %.2f" % pair_stats.tajimas_d)


