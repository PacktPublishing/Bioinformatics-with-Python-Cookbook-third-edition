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


def add_sample_new_dict(dict_db, gene_list):
    my_dict_db = dict(dict_db)  # next recipe
    for gene in gene_list:
        my_dict_db[gene] = my_dict_db.get(0) + 1
    return my_dict_db


restore_db('my_genes.csv')

gene_count = load('my_genes.csv')
gene_count = add_sample_new_dict(gene_count, ['MC4R', 'TYR'])
gene_count = add_sample_new_dict(gene_count, ['LCT', 'HLA-A'])
gene_count = add_sample_new_dict(gene_count, ['HLA-B', 'HLA-C'])
save(gene_count, 'my_genes.csv')
