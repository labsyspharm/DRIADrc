library( tidyverse )
library( DRIAD )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/rebuttal", ifcollision="overwrite.local" )

## Set up the prediction task
XY <- prepareTask( "~/data/ROSMAP/rosmap.tsv.gz", "AC" )
lP <- preparePairs(XY)

## Load relevant data
## Prefilter to the common space of genes
read_csv2 <- function( ... ) { read_csv( ..., col_types=cols() ) %>% filter(Gene %in% colnames(XY)) }
dfx1 <- syn( "syn18145776" ) %>% read_csv2
dfx2 <- syn( "syn17167348" ) %>% read_csv2

## Count the number of differentially expressed genes for each drug
##nn1 <- dfx1 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )
##nn2 <- dfx2 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )
##bind_rows(nn1, nn2) %>% ggplot( aes(x=nn) ) + geom_histogram(binwidth=100) + ggsave( "test.pdf" )

## Identify drugs with a large number of differentially expressed genes
d1 <- dfx1 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% filter( n() > 300 ) %>% pull(Drug) %>% unique
d2 <- dfx2 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% filter( n() > 300 ) %>% pull(Drug) %>% unique

## Selects the top k differentially expressed genes for drugs d
topDFX <- function( .dfx, d, k )
    .dfx %>% filter(Drug %in% d) %>% group_by(Drug) %>% arrange(FDR) %>% slice(1:k) %>%
        summarize_at("Gene", list) %>% with( set_names(Gene, Drug) )

## Select increasingly large sets of differentially expressed genes for these drugs
v <- seq(100, 900, by=100)
DFX1 <- set_names( v, str_c("L",v) ) %>% map( ~topDFX(dfx1, d1, .x) )
DFX2 <- set_names( v, str_c("L",v) ) %>% map( ~topDFX(dfx2, d2, .x) )
DFX <- map2( DFX1, DFX2, c )

## Evalute all gene sets
##res <- imap( DFX, ~{cat("Evaluating ",.y,"\n"); evalGeneSets(.x, XY, lP)} )
##saveRDS( res,"driad-geq300.rds" )
res <- readRDS("driad-geq1000.rds" )

RR <- bind_rows( res, .id="nn" ) %>% mutate( nFeats = map_int(Feats, length) ) %>%
    select( nn, nFeats, Drug=Set, AUC ) %>% filter( nFeats %% 100 == 0 ) %>%
    arrange(Drug) %>% group_by( nn ) %>% summarize_at("AUC",list)
CM <- expand_grid(
    rename_all( RR, str_c, "x" ),
    rename_all( RR, str_c, "y" )) %>%
    mutate(cr = map2_dbl(AUCx, AUCy, cor, method="sp"),
           lbl = as.character(round(cr,2)))

pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
ggplot( CM, aes(x=nnx, y=nny, fill=cr) ) +
    geom_tile() + theme_bw() +
    geom_text(aes(label=lbl)) +
    scale_fill_gradientn( colors=pal, limits=c(-1,1) ) +
    ggsave("cr.pdf")

