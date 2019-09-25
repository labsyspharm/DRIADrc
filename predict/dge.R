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
## synid - synapse ID of the gene set in a .gmt format
## fnData - filename of a wrangled AMP-AD dataset
## fnOut - filename of the output .RData
evalSynSet <- function( synid, fnData, fnOut )
{
    ## Download the gene sets
    synapser::synLogin()
    fn <- synapser::synGet( synid, downloadLocation="." )$path
    gsiROSMAP <- read_gmt(fn)

    ## Set up the prediction task
    XY <- prepareTask( fnData, "AC" )
    lP <- preparePairs(XY)

    ## Evaluate all gene sets
    RR <- evalGeneSets( gsiROSMAP, XY, lP, 100 )
    save( RR, file=fnOut )
}

