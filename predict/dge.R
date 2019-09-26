source("~/projects/DRIAD/R/lpocv.R")

## Parses a .gmt file and puts it into the list format
## iName - index of the column containing pathway names
##    (This is typically 1 for Broad MSigDB sets, and 2 for PathwayCommons sets)
read_gmt <- function( fn, iName=1 )
{
    readr::read_lines(fn) %>% stringr::str_split( "\\t" ) %>%
        set_names( purrr::map_chr(., dplyr::nth, iName) ) %>%
        purrr::map( ~.x[-2:-1] )
}

## Evaluates a gene set that lives on Synapse
## Saves results into a local .RData
## fnGMT   - filename of gene set in a .gmt format
## fnData  - filename of a wrangled AMP-AD dataset
## fnOut   - filename where output will be written
evalGMTset <- function( fnGMT, fnData, fnOut )
{
    ## Download the gene sets
    gsis <- read_gmt( fnGMT )

    ## Set up the prediction task
    XY <- prepareTask( fnData, "AC" )
    lP <- preparePairs(XY)

    ## Evaluate all gene sets
    RR <- evalGeneSets( gsis, XY, lP, 100 )
    save( RR, file=fnOut )
}

## Runs both DGE experiments on all datasets
vGS <- c( DGE1 = "genesets/DGE1.gmt",
         DGE2 = "genesets/DGE2.gmt" )
vDS <- c( ROSMAP = "~/data/amp-ad/rosmap/rosmap.tsv.gz",
         MSBB10 = "~/data/amp-ad/msbb/msbb10.tsv.gz",
         MSBB22 = "~/data/amp-ad/msbb/msbb22.tsv.gz",
         MSBB36 = "~/data/amp-ad/msbb/msbb36.tsv.gz",
         MSBB44 = "~/data/amp-ad/msbb/msbb44.tsv.gz" )

X <- tidyr::crossing( DGE=names(vGS), Dataset=names(vDS) ) %>%
    dplyr::mutate( fnGMT=vGS[DGE], fnData=vDS[Dataset],
                  fnOut=stringr::str_c("results/", DGE, "-", Dataset, ".RData") )

library( purrr )
future::plan( future::multiprocess )
R <- furrr::future_pmap( list(X$fnGMT, X$fnData, X$fnOut), evalGMTset )
