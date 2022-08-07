# Base image https://hub.docker.com/u/rocker/
FROM rocker/verse:4.1.3

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>"

RUN /rocker_scripts/install_geospatial.sh

## create directories
RUN mkdir -p /code
RUN mkdir -p /output

## copy files
COPY /code/CIRC_NEIGH_PARALLEL_SNP_Khaya_2.R /code/CIRC_NEIGH_PARALLEL_SNP_Khaya_2.R
COPY /code/bnut_filtered.csv /code/bnut_filtered.csv

RUN R -e "install.packages( c( 'install2.r', 'parallel', 'adegenet', 'plyr', 'stringr', 'pacman'), dependencies = TRUE )"

# WORKDIR /payload/

## run the script
CMD Rscript /code/CIRC_NEIGH_PARALLEL_SNP_Khaya_2.R

# build image:
# docker build -t gabrielforest/ic .

# start container to test the myimage:
# docker run -it --rm -v ~/"R-Script in Docker"/output:/output gabrielforest/ic
# docker run -it -v ~/output:/output gabrielforest/ic

#The password is set to apohkohhaiqu6IaJ
#If you want to set your own password, set the PASSWORD environment variable. e.g. run with:
#docker run -e PASSWORD=<YOUR_PASS> -p 8787:8787 rocker/rstudio