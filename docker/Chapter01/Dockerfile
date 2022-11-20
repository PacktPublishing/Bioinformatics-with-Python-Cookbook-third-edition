FROM tiagoantao/bio3
MAINTAINER Tiago Antao <tiago@tiago.org>
# RUN conda create -n bioinformatics_r --clone bioinformatics_base

#RUN conda install -n bioinformatics_r r-base=4.1.3 r-ggplot2=3.3.5 r-lazyeval=0.2.2 r-gridextra=2.3 rpy2
RUN conda create -n bioinformatics_r jupyterlab jupytext pandas
RUN conda install -n bioinformatics_r r-base r-ggplot2 r-lazyeval r-gridextra rpy2
CMD conda run --no-capture-output -n bioinformatics_r jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --port=9875 --NotebookApp.token='' --NotebookApp.password=''
