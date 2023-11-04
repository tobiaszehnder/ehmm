FROM rocker/r-ver:4.3.1

# Install devtools
RUN Rscript -e "install.packages('devtools')"

COPY . /tmp

# Install package
RUN Rscript -e "devtools::install('/tmp')"