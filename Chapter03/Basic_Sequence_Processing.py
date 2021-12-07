# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from Bio import Entrez, Seq, SeqIO
from Bio.Alphabet import IUPAC

Entrez.email = "put@your_email.here" 
hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='fasta')  # Lactase gene
#for l in hdl:
#    print l
seq = SeqIO.read(hdl, 'fasta')

w_seq = seq[11:5795]
w_seq

w_hdl = open('example.fasta', 'w')
SeqIO.write([w_seq], w_hdl, 'fasta')
w_hdl.close()

recs = SeqIO.parse('example.fasta', 'fasta')
for rec in recs:
    seq = rec.seq
    print(rec.description)
    print(seq[:10])
    print(seq.alphabet)

seq = Seq.Seq(str(seq), IUPAC.unambiguous_dna)
seq

print((seq[:12], seq[-12:]))
rna = seq.transcribe()
rna

prot = seq.translate()
prot


