library( tidyverse )
library( DRIAD )
library( here )

##future::plan( future::multiprocess )

## Set up the prediction task
XY <- prepareTask( "~/data/ROSMAP/rosmap.tsv.gz", "AC" )
lP <- preparePairs(XY)

read_csv2 <- function( ... ) { read_csv(..., col_types=cols()) %>% filter(Gene %in% colnames(XY)) }
DFX <- bind_rows(read_csv2( "~/data/DRIAD/rebuttal/DGE1-dfx.csv.gz" ),
                 read_csv2( "~/data/DRIAD/rebuttal/DGE2-dfx.csv.gz" ))

gsi <- DFX %>% filter( FDR < 0.05 ) %>% group_by( Drug ) %>% summarize_at( "Gene", list ) %>%
    with( set_names(Gene, Drug) ) %>% keep( ~length(.x) >= 10 )

## Run individual methods and store intermediates
##RR1 <- evalGeneSets( gsi, XY, lP, 0, "lgr" )
##saveRDS( RR1, "pred/lgr-all.rds" )
RR1 <- readRDS( "pred/lgr-all.rds" )

##RR2 <- evalGeneSets( gsi, XY, lP, 0, "svm" )
##saveRDS( RR2, "pred/svm-all.rds" )
RR2 <- readRDS( "pred/svm-all.rds" )

##RR3 <- evalGeneSets( gsi, XY, lP, 0, "xgb" )
##saveRDS( RR3, "pred/xgb-all.rds" )
RR3 <- readRDS( "pred/xgb-all.rds" )

##RR4 <- evalGeneSets( gsi, XY, lP, 0, "nn" )
##saveRDS( RR4, "pred/nn-all.rds" )
RR4 <- readRDS( "pred/nn-all.rds" )

## Combine everything into a common data frame
RR <- list( lgr = RR1, svm = RR2, xgb = RR3, nn = RR4 ) %>%
    bind_rows( .id="Method" ) %>% select( -BK, -pval )
saveRDS( RR, here("results/allpred.rds") )
