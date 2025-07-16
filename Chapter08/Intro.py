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

# +
from collections import defaultdict

import requests

from Bio import ExPASy, SwissProt
# -

#explain why not biopython
server = 'https://rest.uniprot.org/uniprotkb/search'
def do_request(server, **kwargs):
    params = ''
    req = requests.get(server, params=kwargs)
    if not req.ok:
        req.raise_for_status()
    return req


req = do_request(server,
    # 1. Filtering human p53, reviewed entries
    query='gene:p53 AND reviewed:true AND organism_id:9606',
    format='tsv',
    # 2. Specifying output columns with REST API field names
    fields='accession,id,protein_name,gene_names,organism_name,length',
    size=50
)
print(req.text)

#We might revisit this for KEGG

# +
#XXX - stringio
import pandas as pd
import io

uniprot_list = pd.read_table(io.StringIO(req.text))
uniprot_list.rename(columns={'Organism ID': 'ID'}, 
inplace=True)
print(uniprot_list)
# -

p53_human = uniprot_list[
    (uniprot_list.Entry == 'P04637') &
    (uniprot_list['Entry Name'].str.contains('P53_HUMAN'))]['Entry'].iloc[0]

handle = ExPASy.get_sprot_raw(p53_human)

sp_rec = SwissProt.read(handle)

print(sp_rec.entry_name, sp_rec.sequence_length, sp_rec.gene_name)
print(sp_rec.description)
print(sp_rec.organism, sp_rec.seqinfo)
print(sp_rec.sequence)

print(sp_rec.comments)
print(sp_rec.keywords)

help(sp_rec)

done_features = set()
print('Total features:', len(sp_rec.features))
for feature in sp_rec.features:
    if feature in done_features:
        continue
    else:
        done_features.add(feature)
        print(feature)
print('Cross references: ',len(sp_rec.cross_references))
per_source = defaultdict(list)
for xref in sp_rec.cross_references:
    source = xref[0]
    per_source[source].append(xref[1:])
print(per_source.keys())
done_GOs = set()
print('Annotation SOURCES:', len(per_source['GO']))
for annot in per_source['GO']:
    if annot[1][0] in done_GOs:
        continue
    else:
        done_GOs.add(annot[1][0])
        print(annot)


