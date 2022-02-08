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

import dendropy
from dendropy.interop import genbank


# ## Getting the data

# +
def get_ebov_2014_sources():
    #EBOV_2014
    #yield 'EBOV_2014', genbank.GenBankDna(id_range=(233036, 233118), prefix='KM')
    yield 'EBOV_2014', genbank.GenBankDna(id_range=(34549, 34563), prefix='KM0')
    
def get_other_ebov_sources():
    #EBOV other
    yield 'EBOV_1976', genbank.GenBankDna(ids=['AF272001', 'KC242801'])
    yield 'EBOV_1995', genbank.GenBankDna(ids=['KC242796', 'KC242799'])
    yield 'EBOV_2007', genbank.GenBankDna(id_range=(84, 90), prefix='KC2427')
    
def get_other_ebolavirus_sources():
    #BDBV
    yield 'BDBV', genbank.GenBankDna(id_range=(3, 6), prefix='KC54539')
    yield 'BDBV', genbank.GenBankDna(ids=['FJ217161'])

    #RESTV
    yield 'RESTV', genbank.GenBankDna(ids=['AB050936', 'JX477165', 'JX477166', 'FJ621583', 'FJ621584', 'FJ621585']) 

    #SUDV
    yield 'SUDV', genbank.GenBankDna(ids=['KC242783', 'AY729654', 'EU338380',
                                          'JN638998', 'FJ968794', 'KC589025', 'JN638998'])
    #yield 'SUDV', genbank.GenBankDna(id_range=(89, 92), prefix='KC5453')    

    #TAFV
    yield 'TAFV', genbank.GenBankDna(ids=['FJ217162'])


# +
other = open('other.fasta', 'w')
sampled = open('sample.fasta', 'w')

for species, recs in get_other_ebolavirus_sources():
    tn = dendropy.TaxonNamespace()
    char_mat = recs.generate_char_matrix(taxon_namespace=tn,
        gb_to_taxon_fn=lambda gb: tn.require_taxon(label='%s_%s' % (species, gb.accession)))
    char_mat.write_to_stream(other, 'fasta')
    char_mat.write_to_stream(sampled, 'fasta')
other.close()
ebov_2014 = open('ebov_2014.fasta', 'w')
ebov = open('ebov.fasta', 'w')
for species, recs in get_ebov_2014_sources():
    tn = dendropy.TaxonNamespace()
    char_mat = recs.generate_char_matrix(taxon_namespace=tn,
        gb_to_taxon_fn=lambda gb: tn.require_taxon(label='EBOV_2014_%s' % gb.accession))
    char_mat.write_to_stream(ebov_2014, 'fasta')
    char_mat.write_to_stream(sampled, 'fasta')
    char_mat.write_to_stream(ebov, 'fasta')
ebov_2014.close()

ebov_2007 = open('ebov_2007.fasta', 'w')
for species, recs in get_other_ebov_sources():
    tn = dendropy.TaxonNamespace()
    char_mat = recs.generate_char_matrix(taxon_namespace=tn,
        gb_to_taxon_fn=lambda gb: tn.require_taxon(label='%s_%s' % (species, gb.accession)))
    char_mat.write_to_stream(ebov, 'fasta')
    char_mat.write_to_stream(sampled, 'fasta')
    if species == 'EBOV_2007':
        char_mat.write_to_stream(ebov_2007, 'fasta')

ebov.close()
ebov_2007.close()
sampled.close()
# -

# ## Genes

# +
my_genes = ['NP', 'L', 'VP35', 'VP40']

def dump_genes(species, recs, g_dls, p_hdls):
    for rec in recs:

        for feature in rec.feature_table:
                    if feature.key == 'CDS':
                        gene_name = None
                        for qual in feature.qualifiers:
                            if qual.name == 'gene':
                                if qual.value in my_genes:
                                    gene_name = qual.value
                            elif qual.name == 'translation':
                                protein_translation = qual.value
                        if gene_name is not None:
                            locs = feature.location.split('.')
                            start, end = int(locs[0]), int(locs[-1])
                            g_hdls[gene_name].write('>%s_%s\n' % (species, rec.accession))
                            p_hdls[gene_name].write('>%s_%s\n' % (species, rec.accession))
                            g_hdls[gene_name].write('%s\n' % rec.sequence_text[start - 1 : end])
                            p_hdls[gene_name].write('%s\n' % protein_translation)

g_hdls = {}
p_hdls = {}
for gene in my_genes:
    g_hdls[gene] = open('%s.fasta' % gene, 'w')
    p_hdls[gene] = open('%s_P.fasta' % gene, 'w')
for species, recs in get_other_ebolavirus_sources():
    if species in ['RESTV', 'SUDV']:
        dump_genes(species, recs, g_hdls, p_hdls)
for gene in my_genes:
    g_hdls[gene].close()
    p_hdls[gene].close()


# -

# ## Genome exploration

def describe_seqs(seqs):
    print('Number of sequences: %d' % len(seqs.taxon_namespace))
    print('First 10 taxon sets: %s' % ' '.join([taxon.label for taxon in seqs.taxon_namespace[:10]]))
    lens = []
    for tax, seq in seqs.items():
        lens.append(len([x for x in seq.symbols_as_list() if x != '-']))
    print('Genome length: min %d, mean %.1f, max %d' % (min(lens), sum(lens) / len(lens), max(lens)))


ebov_seqs = dendropy.DnaCharacterMatrix.get_from_path('ebov.fasta', schema='fasta', data_type='dna')
print('EBOV')
describe_seqs(ebov_seqs)
del ebov_seqs

print('ebolavirus sequences')
ebolav_seqs = dendropy.DnaCharacterMatrix.get_from_path('other.fasta', schema='fasta', data_type='dna')
describe_seqs(ebolav_seqs)
from collections import defaultdict
species = defaultdict(int)
for taxon in ebolav_seqs.taxon_namespace:
    toks = taxon.label.split('_')
    my_species = toks[0]
    if my_species == 'EBOV':
        ident = '%s (%s)' % (my_species, toks[1])
    else:
        ident = my_species
    species[ident] += 1
for my_species, cnt in species.items():
    print("%20s: %d" % (my_species, cnt))
del ebolav_seqs

# ## Genes

# +
import os
gene_length = {}
my_genes = ['NP', 'L', 'VP35', 'VP40']

for name in my_genes:
    gene_name = name.split('.')[0]
    seqs = dendropy.DnaCharacterMatrix.get_from_path('%s.fasta' % name, schema='fasta', data_type='dna')
    gene_length[gene_name] = []
    for tax, seq in seqs.items():
        gene_length[gene_name].append(len([x for x in seq.symbols_as_list() if x != '-']))
for gene, lens in gene_length.items():
    print ('%6s: %d' % (gene, sum(lens) / len(lens)))
# -


