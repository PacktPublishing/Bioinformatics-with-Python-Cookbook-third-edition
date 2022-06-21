nextflow.enable.dsl=2

download_root = "https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3_r3"


process plink_download {
  output:
  path 'hapmap.map.gz'//, emit: mapgz
  path 'hapmap.ped.gz'//, emit: pedgz
 
  script:
  """
  wget $download_root/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.map.gz -O hapmap.map.gz
  wget $download_root/plink_format/hapmap3_r3_b36_fwd.consensus.qc.poly.ped.gz -O hapmap.ped.gz
   """
}


process uncompress_plink {
  publishDir 'data', glob: '*', mode: 'copy'
  
  input:
  path mapgz
  path pedgz

  output:
  path 'hapmap.map'
  path 'hapmap.ped'

  script:
  """
  gzip -dc $mapgz > hapmap.map
  gzip -dc $pedgz > hapmap.ped
  """
}

//DSL 2 and docs
//conda

process subsample_1p {
  input:
  path 'hapmap.map'
  path 'hapmap.ped'

  output:
  path 'hapmap1.map'
  path 'hapmap1.ped'

  script:
  """
  plink2 --pedmap hapmap --out hapmap1 --thin 0.01 --geno 0.1 --export ped
  """
}

process plink_pca {
  input:
  path 'hapmap.map'
  path 'hapmap.ped'

  output:
  path 'hapmap.eigenvec'
  path 'hapmap.eigenval'

  script:
  """
  plink2 --pca --pedmap hapmap -out hapmap
  """
}


process plot_pca {
  publishDir '.', glob: '*', mode: 'copy'

  input:
  path 'hapmap.eigenvec'
  path 'hapmap.eigenval'

  output:
  path 'pca.png'

  script:
  """
  #!/usr/bin/env python
  import pandas as pd

  pca_df = pd.read_csv('hapmap.eigenvec', sep='\t') 
  ax = pca_df.plot.scatter(x=2, y=3, figsize=(16, 9))
  ax.figure.savefig('pca.png')
  """
}


/*
workflow {
    plink_download | uncompress_plink
}
*/


/*
workflow {
    ped_file = file('data/hapmap.ped')
    map_file = file('data/hapmap.map')
    if (!ped_file.exists() | !map_file.exists()) {
        plink_download | uncompress_plink
    }
}
*/


workflow {
    ped_file = file('data/hapmap.ped')
    map_file = file('data/hapmap.map')
    if (!ped_file.exists() | !map_file.exists()) {
        plink_download | uncompress_plink | subsample_1p | plink_pca | plot_pca
    }
    else {
        subsample_1p(
            Channel.fromPath('data/hapmap.map'),
            Channel.fromPath('data/hapmap.ped')) | plink_pca | plot_pca
    }
}
