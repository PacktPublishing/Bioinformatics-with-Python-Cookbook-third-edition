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

import gffutils
import gzip
from Bio import Seq, SeqIO

# ## Retrieving data

# !rm -f ag.db
# !wget https://vectorbase.org/common/downloads/release-55/AgambiaePEST/gff/data/VectorBase-55_AgambiaePEST.gff -O gambiae.gff
# !gzip -9 gambiae.gff

db = gffutils.FeatureDB('ag.db')

# # Getting a gene

gene_id = 'AGAP004707'

gene = db[gene_id]

print(gene)
print(gene.seqid, gene.strand)

recs = SeqIO.parse(gzip.open('gambiae.fa.gz', 'rt', encoding='utf-8'), 'fasta')
for rec in recs:
    print(rec.description)
    if rec.id == gene.seqid:
        my_seq = rec.seq
        break


# +
def get_sequence(chrom_seq, CDSs, strand):
    seq = Seq.Seq('')
    for CDS in CDSs:
        # #FRAME???
        my_cds = Seq.Seq(str(chrom_seq[CDS.start - 1: CDS.end]))
        seq += my_cds
    return seq if strand == '+' else seq.reverse_complement()


# +
mRNAs = db.children(gene, featuretype='mRNA')
for mRNA in mRNAs:
    print(mRNA.id)
    if mRNA.id.endswith('RA'):
        break

CDSs = db.children(mRNA, featuretype='CDS', order_by='start')
gene_seq = get_sequence(my_seq, CDSs, gene.strand)

print(len(gene_seq), gene_seq)
prot = gene_seq.translate()
print(len(prot), prot)
# -

# # Reverse strand

reverse_transcript_id = 'AGAP004708-RA'

# +
reverse_CDSs = db.children(reverse_transcript_id, featuretype='CDS', order_by='start')
reverse_seq = get_sequence(my_seq, reverse_CDSs, '-')

print(len(reverse_seq), reverse_seq)
reverse_prot = reverse_seq.translate()
print(len(reverse_prot), reverse_prot)
# -


