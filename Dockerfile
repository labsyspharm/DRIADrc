FROM rocker/tidyverse:3.6.2

RUN apt-get update && \
    apt-get install -y texlive-latex-recommended

RUN R -e "install.packages(c('here', 'glue', 'openxlsx', 'toOrdinal', 'numform', \
          'gridExtra', 'gtable', 'cowplot', 'RColorBrewer', \
          'ggridges', 'ggbeeswarm', 'ggthemes', 'ggrepel', 'egg', \
          'pheatmap', 'grConvert', 'grImport2'))"

COPY *.tex Makefile /
COPY figures/*.R /figures/
COPY results /results/
COPY schematics /schematics/
