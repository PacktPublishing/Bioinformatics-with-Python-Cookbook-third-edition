import pandas as pd


def load(file_name):
    df = pd.read_csv(file_name).set_index('gene')
    return dict(df['count'])


def get_min_reads(all_data, min_reads):
    return {
        gene: count
        for gene, count in all_data.items()
        if count >= min_reads
    }


def has_min_observations(subset_data, min_observations):
    return len(subset_data) >= min_observations


print(has_min_observations(
    get_min_reads(
        load('my_genes.csv'), 4
    ), 3))


def get_rec(file_name):
    with open(file_name) as f:
        f.readline()  # header
        for line in f:
            toks = line.strip().split(',')
            yield toks[0], int(toks[1])


def gene_min_reads(source, min_reads):
    for gene, count in source:
        if count >= min_reads:
            yield gene


def gene_min_observations(subset_source, min_observations):
    my_observations = 0
    for gene in subset_source:
        my_observations += 1
        if my_observations == min_observations:
            return True
    return False


print(gene_min_observations(
    gene_min_reads(
        get_rec('my_genes.csv'), 4
    ), 2))
