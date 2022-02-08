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

from Bio import Entrez, Seq, SeqIO, SeqRecord

Entrez.email = "put@your_email.here" 
hdl = Entrez.efetch(db='nucleotide', id=['NM_002299'], rettype='gb')  # Lactase gene
#for l in hdl:
#    print l
gb_rec = SeqIO.read(hdl, 'gb')

for feature in gb_rec.features:
    if feature.type == 'CDS':
        location = feature.location  # Note translation existing
cds = SeqRecord.SeqRecord(gb_rec.seq[location.start:location.end], 'NM_002299', description='LCT CDS only')

w_hdl = open('example.fasta', 'w')
SeqIO.write([cds], w_hdl, 'fasta')
w_hdl.close()

recs = SeqIO.parse('example.fasta', 'fasta')
for rec in recs:
    seq = rec.seq
    print(rec.description)
    print(seq[:10])

print((seq[:12], seq[-12:]))
rna = seq.transcribe()
rna

prot = seq.translate()
prot


