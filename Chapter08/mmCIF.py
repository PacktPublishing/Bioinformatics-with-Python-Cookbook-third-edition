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

from Bio import PDB

# !rm -f 1tup.cif 2>/dev/null
# !wget "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=cif&compression=NO&structureId=1TUP" -O 1tup.cif

parser = PDB.MMCIFParser()
p53_1tup = parser.get_structure('P53_HUMAN', '1tup.cif')


def describe_model(name, pdb):
    print()
    for model in p53_1tup:
        for chain in model:
            print('%s - Chain: %s. Number of residues: %d. Number of atoms: %d.' %
                  (name, chain.id, len(chain), len(list(chain.get_atoms()))))
describe_model('1TUP', p53_1tup)

done_chain = set()
for residue in p53_1tup.get_residues():
    chain = residue.parent
    if chain.id in done_chain:
        continue
    done_chain.add(chain.id)
    print(chain.id, residue.id)

mmcif_dict = PDB.MMCIF2Dict.MMCIF2Dict('1tup.cif')

for k, v in mmcif_dict.items():
    print(k, v)
    print()


