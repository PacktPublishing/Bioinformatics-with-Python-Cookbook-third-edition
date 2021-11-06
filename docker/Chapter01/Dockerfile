FROM tiagoantao/bio3
MAINTAINER Tiago Antao <tiago@tiago.org>
RUN conda create -n bioinformatics-r --clone bionformatics_base

RUN conda init bash
RUN conda activate bioinformatics-r; conda install r-ggplot2=3.3.5 r-lazyeval=0.2.2 r-gridextra=2.3 rpy2
CMD conda init bash; conda activate bioinformatics-r; jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --port=9875 --NotebookApp.token='' --NotebookApp.password=''
