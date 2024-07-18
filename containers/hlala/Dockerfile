FROM continuumio/anaconda3:2022.05 

LABEL maintainer="Anne Marie Noronha (noronhaa@mskcc.org)" \
    version.base="4.10.3" \
    version.image="0.0.1" 

ENV envName "hlala"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
RUN conda update -n base -c defaults conda
RUN conda env create --name ${envName} -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/${envName}/bin:$PATH
RUN echo "export PATH=$PATH" > /etc/environment

