# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## copy files
COPY /install_packages.R /install_packages.R
COPY /CIRC_NEIGH_PARALLEL_SNP_Khaya.R /CIRC_NEIGH_PARALLEL_SNP_Khaya.R

## install R-packages
RUN Rscript /install_packages.R

# run the Rscript automatically
CMD Rscript /CIRC_NEIGH_PARALLEL_SNP_Khaya.R

# build image:
# docker build -t myname/myimage .

# start container to test the myimage:
# docker run -it --rm -v ~/"R-Script in Docker" myname/myimage