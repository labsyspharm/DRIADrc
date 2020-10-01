library( tidyverse )
library( DRIAD )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/rebuttal", ifcollision="overwrite.local" )

## Set up the prediction task
XY <- prepareTask( "~/data/ROSMAP/rosmap.tsv.gz", "AC" )

## Load relevant data
## Prefilter to the common space of genes
read_csv2 <- function( ... ) { read_csv( ..., col_types=cols() ) %>% filter(Gene %in% colnames(XY)) }
dfx1 <- syn( "syn18145776" ) %>% read_csv2
dfx2 <- syn( "syn17167348" ) %>% read_csv2

## Count the number of differentially expressed genes for each drug
nn1 <- dfx1 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )
nn2 <- dfx2 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))

## Panel A
ggplot( bind_rows(nn1, nn2), aes(x=nn) ) + theme_bw() +
    geom_histogram(binwidth=100, fill="steelblue") +
    geom_vline(xintercept=300) +
    scale_x_continuous( breaks=c(0, 300, 1000, 2000, 3000) ) +
    xlab( "Number of diff. expressed genes with FDR < 0.05" ) +
    ylab( "Number of drugs" ) +
    ggsave( "R4-4.png", width=8, height=4 )
