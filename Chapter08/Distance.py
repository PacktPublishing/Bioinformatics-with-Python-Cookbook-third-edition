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
import math
import timeit

from Bio import PDB
# -

repository = PDB.PDBList()
parser = PDB.PDBParser()
repository.retrieve_pdb_file('1TUP', file_format='pdb', pdir='.')  # XXX
p53_1tup = parser.get_structure('P 53', 'pdb1tup.ent')

zns = []
for atom in p53_1tup.get_atoms():
    if atom.element == 'ZN':
        #print(atom, dir(atom), atom.mass, atom.element, atom.coord[0])
        zns.append(atom)
for zn in zns:
        print(zn, zn.coord)


# +
#Suggest a pymol viewing
# -

#Try this in numba?
def get_closest_atoms(pdb_struct, ref_atom, distance):
    atoms = {}
    rx, ry, rz = ref_atom.coord
    for atom in pdb_struct.get_atoms():
        if atom == ref_atom:
            continue
        x, y, z = atom.coord
        my_dist = math.sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2) 
        if my_dist < distance:
            atoms[atom] = my_dist
    return atoms


for zn in zns:
    print()
    print(zn.coord)
    atoms = get_closest_atoms(p53_1tup, zn, 4)
    for atom, distance in atoms.items():
        print(atom.element, distance, atom.coord)

for distance in [1, 2, 4, 8, 16, 32, 64, 128]:
    my_atoms = []
    for zn in zns:
        atoms = get_closest_atoms(p53_1tup, zn, distance)
        my_atoms.append(len(atoms))
    print(distance, my_atoms)

nexecs = 10
print(timeit.timeit('get_closest_atoms(p53_1tup, zns[0], 4.0)',
                    'from __main__ import get_closest_atoms, p53_1tup, zns',
                    number=nexecs) / nexecs * 1000)


def get_closest_alternative(pdb_struct, ref_atom, distance):
    atoms = {}
    rx, ry, rz = ref_atom.coord
    for atom in pdb_struct.get_atoms():
        if atom == ref_atom:
            continue
        x, y, z = atom.coord
        if abs(x - rx) > distance or abs(y - ry) > distance or abs(z - rz) > distance:
            continue
        my_dist = math.sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2) 
        if my_dist < distance:
            atoms[atom] = my_dist
    return atoms


print(timeit.timeit('get_closest_alternative(p53_1tup, zns[0], 4.0)',
                    'from __main__ import get_closest_alternative, p53_1tup, zns',
                    number=nexecs) / nexecs * 1000)

print('Standard')
for distance in [1, 4, 16, 64, 128]:
    print(timeit.timeit('get_closest_atoms(p53_1tup, zns[0], distance)',
                        'from __main__ import get_closest_atoms, p53_1tup, zns, distance',
                        number=nexecs) / nexecs * 1000)
print('Optimized')
for distance in [1, 4, 16, 64, 128]:
    print(timeit.timeit('get_closest_alternative(p53_1tup, zns[0], distance)',
                        'from __main__ import get_closest_alternative, p53_1tup, zns, distance',
                        number=nexecs) / nexecs * 1000)



# +
#for interesting distances
