## Currently installed packages
ipkgs <- rownames(installed.packages())

## CRAN
cran <- c("tidyverse", "devtools", "here", "glue", "openxlsx",
          "toOrdinal", "numform",
          "gridExtra", "gtable", "cowplot", "RColorBrewer",
          "ggridges", "ggbeeswarm", "ggthemes", "ggrepel", "egg",
          "pheatmap", "grConvert", "grImport2" )

install.packages( setdiff(cran, ipkgs) )

## Synapser
if( !("synapser" %in% ipkgs) )
    install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))    

## SynExtra
if( !("synExtra" %in% ipkgs) )
    devtools::install_github( "ArtemSokolov/synExtra" )
