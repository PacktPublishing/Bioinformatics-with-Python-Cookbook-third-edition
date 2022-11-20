FROM continuumio/anaconda3:2021.05
MAINTAINER Tiago Antao <tiago@tiago.org>
#ENV DEBIAN_FRONTEND noninteractive

#RUN apt-get update && apt-get upgrade -y && apt-get install -y git wget build-essential unzip graphviz libgraphviz-dev pkg-config swig libx11-dev libgsl0-dev libopenblas-dev liblapacke-dev
#RUN apt-get install -y samtools mafft muscle raxml tabix

RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-third-Edition.git

#RUN conda upgrade -n base conda
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda create -n bioinformatics_base --file /Bioinformatics-with-Python-Cookbook-third-Edition/Chapter01/bioinformatics_base.txt
RUN pip install pyarrow==8.0.0
RUN conda init bash

EXPOSE 9875

WORKDIR /Bioinformatics-with-Python-Cookbook-third-Edition

RUN echo setterm -foreground magenta >> /etc/bash.bashrc
CMD conda run --no-capture-output -n bioinformatics_base jupyter-lab --ip=0.0.0.0 --no-browser --allow-root --port=9875 --NotebookApp.token='' --NotebookApp.password=''
