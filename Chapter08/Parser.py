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

#XXX
repository = PDB.PDBList()
repository.retrieve_pdb_file('1TUP', pdir='.', file_format='pdb')

# +
rec_types = {
    #single line
    'HEADER': [(str, 11, 49), (str, 50, 58), (str, 62, 65)],
    #multi_line
    'SOURCE': [(int, 7, 9), (str, 10, 78)],
    #multi_rec
    'LINK' : [(str, 12, 15), (str, 16, 16), (str, 17, 19), (str, 21, 21), (int, 22, 25),
              (str, 26, 26), (str, 42, 45), (str, 46, 46), (str, 47, 49), (str, 51, 51),
              (int, 52, 55), (str, 56, 56), (str, 59, 64), (str, 66, 71), (float, 73, 77)],
    'HELIX': [(int, 7, 9), (str, 11, 13), (str, 15, 17), (str, 19, 19), (int, 21, 24),
              (str, 25, 25), (str, 27, 29), (str, 31, 31),
              (int, 33, 36), (str, 37 ,37), (int, 38, 39), (str, 40, 69), (int, 71, 75)],
    'SHEET': [(int, 7, 9), (str, 11, 13), (int, 14, 15), (str, 17, 19), (str, 21, 21),
              (int, 22, 24), (str, 26, 26), (str, 28, 30),
              (str, 32, 32), (int, 33, 36), (str, 37, 37), (int, 38, 39), (str, 41, 44),
              (str, 45, 47), (str, 49, 49), (int, 50, 53), (str, 54, 54), (str, 56, 59),
              (str, 60, 62), (str, 64, 64), (int, 65, 68), (str, 69, 69)],
}

def parse_pdb(hdl):
    for line in hdl:
        line = line[:-1]  # remove \n but not other whitespace
        toks = []
        for section, elements in rec_types.items():
            if line.startswith(section):
                for fun, start, end in elements:
                    try:
                        toks.append(fun(line[start: end + 1]))
                    except ValueError:
                        toks.append(None)  # eg continuation
                yield (section, toks)
        if len(toks) == 0:
            yield ('UNKNOWN', line)
                


# -

hdl = open('pdb1tup.ent')
done_rec = set()
for rec in parse_pdb(hdl):
    if rec[0] == 'UNKNOWN' or rec[0] in done_rec:
        continue
    print(rec)
    done_rec.add(rec[0])

# +
multi_lines = ['SOURCE']

#assume multi is just a string
def process_multi_lines(hdl):
    current_multi = ''
    current_multi_name = None
    for rec_type, toks in parse_pdb(hdl):
        if current_multi_name is not None and current_multi_name != rec_type:
            yield current_multi_name, [current_multi]
            current_multi = ''
            current_multi_name = None
        if rec_type in multi_lines:
            current_multi += toks[1].strip().rstrip() + ' '
            current_multi_name = rec_type
        else:
            if len(current_multi) != 0:
                yield current_multi_name, [current_multi]
                current_multi = ''
                current_multi_name = None                
            yield rec_type, toks
    if len(current_multi) != 0:
        yield current_multi_name, [current_multi]


# -

hdl = open('pdb1tup.ent')
done_rec = set()
for rec in process_multi_lines(hdl):
    if rec[0] == 'UNKNOWN' or rec[0] in done_rec:
        continue
    print(rec)
    done_rec.add(rec[0])


# +
def get_spec_list(my_str):
    #ignoring escape characters
    spec_list = {}
    elems = my_str.strip().strip().split(';')
    for elem in elems:
        toks = elem.split(':')
        spec_list[toks[0].strip()] = toks[1].strip()
    return spec_list

struct_types = {
    'SOURCE': [get_spec_list] 
}

def process_struct_types(hdl):
    for rec_type, toks in process_multi_lines(hdl):
        if rec_type in struct_types.keys():
            funs = struct_types[rec_type]
            struct_toks = []
            for tok, fun in zip(toks, funs):
                struct_toks.append(fun(tok))
            yield rec_type, struct_toks
        else:
            yield rec_type, toks


# -

hdl = open('pdb1tup.ent')
for rec in process_struct_types(hdl):
    if rec[0] != 'SOURCE':
        continue
    print(rec)


