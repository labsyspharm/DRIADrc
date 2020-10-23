## Pulls relevant results slices needed for figures into here/results
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
library( here )

synapser::synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/DRIAD/figs", ifcollision="overwrite.local" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Raw TAS values
TASvalues <- function()
{
    synapser::synGet( "syn20830941", version=2, ifcollision="overwrite.local",
                     downloadLocation="~/data/DRIAD/figs" ) %>%
        pluck( "path" ) %>% read_rds() %>% chuck("data", 2 ) %>%
        write_csv( here("results","TAS-values.csv.gz") )
}

## The composite score is defined as the harmonic mean p value
##   p-values below 0.01 are thresholded at 0.005 to avoid "zero" issues
## Annotates a results data frame with drug, target and toxicity information
DGEcomposite <- function( task="AC" )
{
    ## Harmonic mean
    hmean <- function(v) {length(v)/sum(1/v)}
    
    load( here("results","results-2019-10-06.RData") )
    
    ## Load the associated LINCS metadata (drug names)
    M <- syn_csv( "syn11801537" ) %>% mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, Drug = name, Target = target_name )

    ## Retrieve additional drug approval status information
    DBA <- syn_csv( "syn19042441" ) %>%
        select( Drug = `Drug Name`, Approval = `FDA Approval Stage` ) %>%
        mutate_at( "Drug", str_to_lower )

    ## Load toxicity information (threshold at 2,200 nuclei count for toxicity)
    ## MG-132 is known to be toxic from previous work (add it by hand)
    TOX <- syn_csv("syn18496657") %>% mutate( Drug = str_to_lower(`Fluid name`) ) %>%
        group_by( Drug ) %>% summarize_at( "Nuclei counts", mean ) %>%
        mutate( IsToxic = as.integer( `Nuclei counts` < 2200 ) ) %>%
        bind_rows( list(Drug="mg-132", IsToxic=1) ) %>%
        select( Drug, IsToxic )

    ## Combine everything into a single dataframe
    R <- allRes %>% filter( Task == task ) %>%
        select( Plate, Dataset, LINCSID=Drug, pval ) %>%
        inner_join( M, by="LINCSID" ) %>% 
        left_join( DBA, by="Drug" ) %>%
        mutate_at( "Approval", replace_na, "experimental" ) %>%
        left_join( TOX, by="Drug" )

    ## Compute the composite score
    R %>% mutate_at( "pval", pmax, 0.0005 ) %>%
        spread( Dataset, pval ) %>%
        mutate( HMP = pmap_dbl(list(MSBB10, MSBB22, MSBB36, MSBB44, ROSMAP), lift_vd(hmean)) ) %>%
        arrange( HMP ) %>%
        write_csv( here("results", "DGE-composite.csv") )
}
