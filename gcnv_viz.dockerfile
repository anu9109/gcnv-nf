FROM rocker/tidyverse:4.3.2

RUN R -e "install.packages(c('data.table', 'patchwork'))"
