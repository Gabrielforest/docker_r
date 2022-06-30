# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

## create directories
RUN mkdir -p /code
RUN mkdir -p /output

## copy files
COPY /code/nice.R /code/nice.R

## install R-packages
RUN R -e "install.packages(c('raster', 'parallel', 'adegenet', 'plyr', 'stringr', 'maptools'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

## run the script
CMD Rscript /code/nice.R

# build image:
# docker build -t gabriel/bnut .

# start container to test the myimage:
# docker run -it --rm -v ~/"R-Script in Docker"/output:/output gabriel/bnut