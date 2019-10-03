## Evaluating gene sets from past AMP-AD analyses in the literature
##
## by Artem Sokolov

source("~/projects/DRIAD/R/lpocv.R")
suppressMessages(library( purrr ))
future::plan( future::multiprocess )

myeval <- function( fnData )
{
    cat( "Evaluating on", fnData, "\n" )
    gsis <- read_gmt("genesets/prev-ampad.gmt" )

    ## Set up the prediction task
    XY <- prepareTask( fnData, "AC" )
    lP <- preparePairs(XY)

    ## Evaluate all gene sets
    evalGeneSets( gsis, XY, lP, 100 )
}

## Run evaluation on all datasets
vDS <- c( ROSMAP = "~/data/amp-ad/rosmap/rosmap.tsv.gz",
         MSBB10 = "~/data/amp-ad/msbb/msbb10.tsv.gz",
         MSBB22 = "~/data/amp-ad/msbb/msbb22.tsv.gz",
         MSBB36 = "~/data/amp-ad/msbb/msbb36.tsv.gz",
         MSBB44 = "~/data/amp-ad/msbb/msbb44.tsv.gz" )

litRes <- tibble::enframe( vDS, "Dataset", "fnData" ) %>%
    dplyr::mutate( Results=purrr::map(fnData, myeval) ) %>%
    dplyr::select( -fnData ) %>% tidyr::unnest()

save( litRes, file=stringr::str_c("results/lit-",Sys.Date(),".RData") )
