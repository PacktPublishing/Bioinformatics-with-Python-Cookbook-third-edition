import pandas as pd

eigen_fname = snakemake.input[0] if snakemake.input[0].endswith('eigenvec') else snakemake.input[1]
pca_df = pd.read_csv(eigen_fname, sep='\t') 
ax = pca_df.plot.scatter(x=2, y=3, figsize=(16, 9))
ax.figure.savefig(snakemake.output[0]) 
