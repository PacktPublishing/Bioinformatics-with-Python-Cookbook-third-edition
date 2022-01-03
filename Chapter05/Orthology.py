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

# +
import requests
 
ensembl_server = 'http://rest.ensembl.org'

def do_request(server, service, *args, **kwargs):
    url_params = ''
    for a in args:
        if a is not None:
            url_params += '/' + a
    req = requests.get('%s/%s%s' % (server, service, url_params),
                       params=kwargs,
                       headers={'Content-Type': 'application/json'})
 
    if not req.ok:
        req.raise_for_status()
    return req.json()


# -

answer = do_request(ensembl_server, 'info/species')
for i, sp in enumerate(answer['species']):
    print(i, sp['name'])

ext_dbs = do_request(ensembl_server, 'info/external_dbs', 'homo_sapiens', filter='HGNC%')
print(ext_dbs)

answer = do_request(ensembl_server, 'lookup/symbol', 'homo_sapiens', 'LCT')
print(answer)
lct_id = answer['id']

lct_seq = do_request(ensembl_server, 'sequence/id', lct_id)
print(lct_seq)

lct_xrefs = do_request(ensembl_server, 'xrefs/id', lct_id)
for xref in lct_xrefs:
    print(xref['db_display_name'])
    print(xref)

refs = do_request(ensembl_server, 'xrefs/id', lct_id, external_db='GO', all_levels='1')
print(lct_id, refs)

hom_response = do_request(ensembl_server, 'homology/id', lct_id, type='orthologues', sequence='none')
#print(hom_response['data'][0]['homologies'])
homologies = hom_response['data'][0]['homologies']
for homology in homologies:
    print(homology['target']['species'])
    if homology['target']['species'] != 'equus_caballus':
        continue
    print(homology)
    print(homology['taxonomy_level'])
    horse_id = homology['target']['id']

horse_req = do_request(ensembl_server, 'lookup/id', horse_id)
print(horse_req)

# +
#maybe synteny of MCM6 and LCT with caballus and gorilla
