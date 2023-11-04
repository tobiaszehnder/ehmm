FROM bioconductor/bioconductor_docker:RELEASE_3_11

COPY . /install

WORKDIR /install

RUN Rscript install.R