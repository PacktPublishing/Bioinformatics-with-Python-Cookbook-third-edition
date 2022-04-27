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


def add_sample_csv(gene_list):
    gene_count = load('my_genes.csv')
    for gene in gene_list:
        gene_count[gene] = gene_count.get(0) + 1
    save(gene_count, 'my_genes.csv')


restore_db('my_genes.csv')

add_sample_csv(['MC4R', 'TYR'])
add_sample_csv(['LCT', 'HLA-A'])
add_sample_csv(['HLA-B', 'HLA-C'])
