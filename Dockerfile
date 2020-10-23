FROM rocker/tidyverse:3.6.2

RUN apt-get update
RUN apt-get install -y texlive-latex-recommended
RUN apt-get install -y libcairo2-dev libspectre-dev librsvg2-dev \
    libpoppler-glib-dev r-base-dev

RUN R -e "install.packages(c('devtools','here', 'glue', 'openxlsx', 'toOrdinal', 'numform'))"
RUN R -e "install.packages(c('gridExtra', 'gtable', 'cowplot', 'RColorBrewer'))"
RUN R -e "install.packages(c('ggridges', 'ggbeeswarm', 'ggthemes', 'ggrepel', 'egg'))"
RUN R -e "install.packages(c('pheatmap', 'grImport2'))"
RUN R -e "devtools::install_github('sjp/grConvert')"

COPY *.tex Makefile /
COPY figures/*.R /figures/
COPY results /results/
COPY schematics /schematics/
