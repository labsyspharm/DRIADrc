## Gene set composition for DGE experiments
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/AMP-AD", ifcollision="overwrite.local" )

## Given a differential expression matrix, composes gene sets for each drug that
##   consist of [significantly] differentially-expressed genes
## Caps the number at 300 (approximate size of the largest mined set)
genesets_dfexp <- function( DX )
{
    ## Given a differential expression matrix (for a single drug), identifies and returns
    ## genes with FDR < thresh. If there are fewer than atLeast such genes, returns
    ## genes with PValue < thresh instead
    top_genes <- function( .df, thresh = 0.05, atLeast=10 )
    {
        vFDR  <- with( .df, set_names(FDR, Gene) ) %>% keep( ~ .x < thresh ) %>% names
        vpval <- with( .df, set_names(PValue, Gene) ) %>% keep( ~ .x < thresh ) %>% names
        if( length(vFDR) < atLeast ) return(vpval)
        vFDR
    }

    GS <- DX %>% mutate_at( "Drug", str_to_lower ) %>% nest( -Drug, .key="DFX" ) %>%
        mutate_at( "DFX", map, top_genes ) %>%
        mutate_at( "DFX", map_if, ~(length(.x) > 300), ~.x[1:300] )

    ## Annotate with metadata
    M <- syn("syn11801537") %>% read_csv(col_types=cols()) %>%
        mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, URL = link, Drug = name )
    stopifnot( all(GS$Drug %in% M$Drug) )

    inner_join( M, GS )
}

## Takes the output of genesets_dfexp() and concatenates it into tab-delimited
##   strings in preparation for writing to .gmt files
prep4gmt <- function( GS )
{
    GS %>% mutate( Content = map_chr(DFX, str_flatten, "\t"),
                  GMT = str_c(LINCSID, "\t", URL, "\t", Content) ) %>%
        with( set_names(GMT, LINCSID) )
}

## Process each DGE experiment separately
main <- function()
{
    ## Load the raw differential expression matrix
    DX <- c( DGE1="syn18145776", DGE2="syn17167348" ) %>% map(syn) %>%
        map( read_csv, col_types=cols() )

    ## Compute dfx-based gene sets and concatenate into tab-delimited strings
    GSdfx <- map( DX, genesets_dfexp ) %>% map( prep4gmt )
    imap( GSdfx, ~cat(.x, sep="\n", file=str_c(.y,".gmt")) )
}

