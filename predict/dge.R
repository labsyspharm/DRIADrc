source("~/projects/DRIAD/R/lpocv.R")
suppressMessages(library( purrr ))
future::plan( future::multiprocess )

## Evaluates a gene set from a GMT file
## fnGMT  - filename of gene set in a .gmt format
## fnData - filename of a wrangled AMP-AD dataset
## task   - one of {AB, AC, BC}
evalGMTset <- function( fnGMT, fnData, task )
{
    stopifnot( task %in% c("AB","AC","BC") )
    cat( "Evaluating", fnGMT, "on", fnData, "for", task, "\n" )
    
    ## Download the gene sets
    gsis <- read_gmt( fnGMT )

    ## Set up the prediction task
    XY <- prepareTask( fnData, task )
    lP <- preparePairs(XY)

    ## Evaluate all gene sets
    RR <- evalGeneSets( gsis, XY, lP, 1000 )
    RR
}

## Runs both DGE experiments on all datasets
vGS <- c( DGE1 = "genesets/DGE1.gmt",
         DGE2 = "genesets/DGE2.gmt" )
vDS <- c( ROSMAP = "~/data/amp-ad/rosmap/rosmap.tsv.gz",
         MSBB10 = "~/data/amp-ad/msbb/msbb10.tsv.gz",
         MSBB22 = "~/data/amp-ad/msbb/msbb22.tsv.gz",
         MSBB36 = "~/data/amp-ad/msbb/msbb36.tsv.gz",
         MSBB44 = "~/data/amp-ad/msbb/msbb44.tsv.gz" )

X <- tidyr::crossing( DGE=names(vGS), Dataset=names(vDS), Task=c("AB","BC","AC") ) %>%
    dplyr::mutate( fnGMT=vGS[DGE], fnData=vDS[Dataset] )

allRes <- X %>% dplyr::mutate( Results=purrr::pmap(list(fnGMT, fnData, Task), evalGMTset) ) %>%
    dplyr::select( -fnGMT, -fnData ) %>% tidyr::unnest() %>% dplyr::rename( Plate=DGE, Drug=Set )

save( allRes, file=stringr::str_c("results/results-",Sys.Date(),".RData") )
