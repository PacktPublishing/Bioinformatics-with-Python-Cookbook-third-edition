import shutil
import pandas as pd


def restore_db(file_name):
    shutil.copyfile(f'{file_name}.base', file_name)


def load(file_name):
    df = pd.read_csv(file_name).set_index('gene')
    return dict(df['count'])


def save(dict_db, file_name):
    pd.Series(dict_db).to_csv(
        file_name, index_label='gene', header=['count'])


def add_sample_dict(dict_db, gene_list):
    for gene in gene_list:
        dict_db[gene] = dict_db.get(0) + 1


def add_sample_new_dict(dict_db, gene_list):
    my_dict_db = dict(dict_db)  # next recipe
    for gene in gene_list:
        my_dict_db[gene] = my_dict_db.get(0) + 1
    return my_dict_db


gene_count = load('my_genes.csv')

add_sample_dict(gene_count, ['DEPP'])

new_gene_count = add_sample_new_dict(gene_count, ['DEPP'])
