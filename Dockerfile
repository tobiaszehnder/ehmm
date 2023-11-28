FROM bioconductor/bioconductor_docker:RELEASE_3_11

COPY dependencies.R /install/dependencies.R
WORKDIR /install
RUN Rscript dependencies.R

COPY . /install
RUN Rscript install.R