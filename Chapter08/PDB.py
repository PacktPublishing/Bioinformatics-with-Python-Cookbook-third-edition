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

from Bio import PDB

repository = PDB.PDBList()
repository.retrieve_pdb_file('1TUP', pdir='.', file_format='pdb')
repository.retrieve_pdb_file('1OLG', pdir='.', file_format='pdb')
repository.retrieve_pdb_file('1YCQ', pdir='.', file_format='pdb')

parser = PDB.PDBParser()
p53_1tup = parser.get_structure('P 53 - DNA Binding', 'pdb1tup.ent')
p53_1olg = parser.get_structure('P 53 - Tetramerization', 'pdb1olg.ent')
p53_1ycq = parser.get_structure('P 53 - Transactivation', 'pdb1ycq.ent')


# +
def print_pdb_headers(headers, indent=0):
    ind_text = ' ' * indent
    for header, content in headers.items():
        if type(content) == dict:
            print('\n%s%20s:' % (ind_text, header))
            print_pdb_headers(content, indent + 4)
            print()
        elif type(content) == list:
            print('%s%20s:' % (ind_text, header))
            for elem in content:
                print('%s%21s %s' % (ind_text, '->', elem))
        else:
            print('%s%20s: %s' % (ind_text, header, content))

print_pdb_headers(p53_1tup.header)
# -

print(p53_1tup.header['compound'])
print(p53_1olg.header['compound'])
print(p53_1ycq.header['compound'])


def describe_model(name, pdb):
    print()
    for model in pdb:
        for chain in model:
            print('%s - Chain: %s. Number of residues: %d. Number of atoms: %d.' %
                  (name, chain.id, len(chain), len(list(chain.get_atoms()))))
describe_model('1TUP', p53_1tup)
describe_model('1OLG', p53_1olg)
describe_model('1YCQ', p53_1ycq)
#will go deep in a next recipe (bottom up)

for residue in p53_1tup.get_residues():
    if residue.id[0] in [' ', 'W']:
        continue
    print(residue.id)

res = next(p53_1tup[0]['A'].get_residues())
print(res)
for atom in res:
    print(atom, atom.serial_number, atom.element)
print(p53_1tup[0]['A'][94]['CA'])

# +
from Bio.SeqIO import PdbIO, FastaIO
from Bio import SeqIO

def get_fasta(pdb_file, fasta_file, transfer_ids=None):
    records = list(PdbIO.PdbSeqresIterator(pdb_file))
    if transfer_ids is not None:
        records = [rec for rec in records if rec.id in transfer_ids and len(rec.seq) > 0]
    else:
        records = [rec for rec in records if len(rec.seq) > 0]
    
    with open(fasta_file, 'w') as out_handle:
        SeqIO.write(records, out_handle, 'fasta')
    for rec in records:
       print(rec.id, rec.seq, len(rec.seq))
        
        
get_fasta('pdb1tup.ent', '1tup.fasta', transfer_ids=['1TUP:B'])
get_fasta('pdb1olg.ent', '1olg.fasta', transfer_ids=['1OLG:B'])
get_fasta('pdb1ycq.ent', '1ycq.fasta', transfer_ids=['1YCQ:B'])
# -


