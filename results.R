## Scripts that simplify loading of results for all figures
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/AMP-AD/figs" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Returns all available results files
synResults <- function()
{
    ## Given a synapse ID, retrieves settings_md5 annotation associated with it
    getSettingsMd5 <- function( synid )
    {
        s <- synapser::synGet( synid, downloadFile=FALSE )
        s$annotations$get("settings_md5")
    }

    ## List of valid types
    vTypes = c("background_predictions", "hypothesis_predictions", "score", "stats")
    
    cat( "Cataloguing available runs...\n" )
    R <- c( "MAYO" = "syn12180241", "ROSMAP" = "syn15589860", "MSBB" = "syn15588043" ) %>%
        map( synExtra::synChildren ) %>% map( enframe, "name", "synid" ) %>%
        bind_rows( .id = "Dataset" ) %>%
        mutate( md5 = map_chr( synid, getSettingsMd5 ),
               chunks = str_split(str_sub(name, 1, -6), "\\_"),
               Region = map_chr(chunks, nth, 3),
               Task = map_chr(chunks,6) ) %>% select( -name, -synid, -chunks ) %>%
        mutate_at( "Region", recode, TempCortex = "TCX", cerebellum = "CBE" )

    cat( "Identifying files in each run...\n" )
    R %>% mutate( synRun = map2_chr(Dataset, md5, ~synExtra::synPluck("syn15590460",
                                                                      str_c(.x, "pc"), .y)),
                 Type = rep( list(vTypes), n() ), md5=NULL ) %>% unnest() %>%
        mutate( synType = map2_chr(synRun, Type, synExtra::synPluck),
               Files = map(synType, synExtra::synChildren) ) %>%
        mutate_at( "Files", map, enframe, "fileName", "fileId" ) %>%
        select( -synRun, -synType )
}


## Pre-defined index of background performances, built by the following code:
##   synResults( "score", Region != "cerebellum" ) %>% filter( grepl( "background", name ) ) %>%
##     select( -name, -parentName ) %>% dput()
## Dataset and Region names are then adjusted by hand:
## 1. Removal of pc suffix from dataset names
## 2. Recoding TempCortex as TCX
indexBackground <- function()
{
    structure(list(Dataset = c("ROSMAP", "ROSMAP", "ROSMAP", "MAYO", "MAYO", "MAYO",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB",
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB"),
                   Region = c("DLPFC", "DLPFC", "DLPFC", "TCX", "TCX",
                              "TCX", "BM10", "BM36", "BM44", "BM22", "BM36", 
                              "BM10", "BM22", "BM10", "BM44", "BM36", "BM44", "BM22"),
                   Task = c("AB", "AC", "BC", "BC", "AB", "AC", "AC", "AB", "BC",
                            "BC", "AC", "AB", "AB", "BC", "AC", "BC", "AB", "AC"),
                   id = c("syn15589822", "syn15589816", "syn15589810", "syn15572112",
                          "syn15570670", "syn15572002", "syn15584696", "syn15585576",
                          "syn15582756", "syn15585289", "syn15581358", "syn15582188", 
                          "syn15586016", "syn15585432", "syn15583710", "syn15584001",
                          "syn15583564", "syn15584563")),
              row.names = c(NA, -18L), class = c("tbl_df", "tbl", "data.frame"))
}

## Fetches all background matrices provided by the given index
bkFromIndex <- function( IDX = indexBackground() )
{
    IDX %>% mutate( BK = map(id, syn_csv), id=NULL ) %>%
        mutate_at( "BK", map, ~gather(.x, Size, AUC) ) %>%
        mutate_at( "BK", map, ~mutate_at( .x, "Size", as.numeric ) )
}

## Pre-defined index of DGE results, constructed using the following code:
##    synResults() %>% filter( Type == "stats", Region != "CBE" ) %>%
##        unnest() %>% mutate( Plate = str_sub(fileName, 1, 4) ) %>%
##        select( -Type, -fileName ) %>% rename( id = fileId ) %>% dput()
indexDGE <- function()
{
    structure(list(Dataset = c("MAYO", "MAYO", "MAYO", "MAYO", "MAYO", 
                               "MAYO", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", "ROSMAP", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", 
                               "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB", "MSBB"),
                   Region = c("TCX", "TCX", "TCX", "TCX", "TCX", "TCX", "DLPFC", 
                              "DLPFC", "DLPFC", "DLPFC", "DLPFC", "DLPFC", "BM10", "BM10", 
                              "BM10", "BM10", "BM10", "BM10", "BM22", "BM22", "BM22", "BM22", 
                              "BM22", "BM22", "BM36", "BM36", "BM36", "BM36", "BM36", "BM36", 
                              "BM44", "BM44", "BM44", "BM44", "BM44", "BM44"),
                   Task = c("AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", "BC", "BC", 
                            "AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", "BC", 
                            "BC", "AB", "AB", "AC", "AC", "BC", "BC", "AB", "AB", "AC", "AC", 
                            "BC", "BC"),
                   id = c("syn18201042", "syn18143076", "syn18201045", 
                          "syn18143082", "syn18201026", "syn18143046", "syn18201029", "syn18143052", 
                          "syn18201066", "syn18143125", "syn18201051", "syn18143094", "syn18201060", 
                          "syn18143112", "syn18201063", "syn18143119", "syn18201017", "syn18143028", 
                          "syn18201032", "syn18143058", "syn18201014", "syn18143022", "syn18201038", 
                          "syn18143070", "syn18201075", "syn18143143", "syn18201072", "syn18143137", 
                          "syn18201023", "syn18143040", "syn18201020", "syn18143034", "syn18201048", 
                          "syn18143088", "syn18201035", "syn18143064"),
                   Plate = c("DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", 
                             "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", 
                             "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", 
                             "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2", "DGE1", 
                             "DGE2", "DGE1", "DGE2")),
              class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -36L))
}

## Annotates a results data frame with drug, target and toxicity information
annotateResults <- function( R )
{
    ## Load the associated LINCS metadata (drug names)
    M <- syn_csv( "syn11801537" ) %>% mutate_at( "name", str_to_lower ) %>%
        select( LINCSID = lincs_id, Drug = name, Target = target_name, URL=link )

    ## Retrieve additional drug approval status information
    DBA <- syn( "syn19042441" ) %>% read_csv( col_types=cols() ) %>%
        select( Drug = `Drug Name`, Approval = `FDA Approval Stage` ) %>%
        mutate_at( "Drug", str_to_lower )

    ## Load the list of approved drugs, taken from DrugBank
    ## vAppr <- syn_csv("syn16932412") %>% with( str_to_lower(Name) )

    ## Load toxicity information (threshold at 2,200 nuclei count for toxicity)
    ## MG-132 is known to be toxic from previous work (add it by hand)
    TOX <- syn_csv("syn18496657") %>% mutate( Drug = str_to_lower(`Fluid name`) ) %>%
        group_by( Drug ) %>% summarize_at( "Nuclei counts", mean ) %>%
        mutate( IsToxic = as.integer( `Nuclei counts` < 2200 ) ) %>%
        bind_rows( list(Drug="mg-132", IsToxic=1) )
    
    R %>% rename( URL = id ) %>% inner_join( M, by="URL" ) %>% 
        left_join( DBA, by="Drug" ) %>% mutate_at( "Approval", replace_na, "experimental" ) %>%
        left_join( TOX, by="Drug" ) %>%
        select( LINCSID, Approval, IsToxic, Drug, Target, Size = intersect, AUC, p_value )
}

## Fetches all results matrices associated with a given results index
##   Annotates results with relevant tags
resFromIndex <- function( IDX = indexDGE() )
{
    ## Fetch all the relevant results values
    ## Generate a tag for each Dataset / Region / Task triplet
    IDX %>% mutate( Results = map(id, syn_csv) ) %>% select( -id ) %>%
        mutate_at( "Results", map, annotateResults )
}

## Given a vector of matching Dataset and Region names, converts them into a single tag
tagDataset <- function( vDataset, vRegion )
{
    recode( vRegion, DLPFC="", TCX="", CBE="" ) %>%
        str_sub( 3, 5 ) %>% str_c( vDataset, . )
}

## Gives more explicit names to classification tasks
tagTask <- function( vTask )
    c( AB="A-vs-B", AC="A-vs-C", BC="B-vs-C" )[vTask]

## Retags a dataset according to tagDataset() and tagTask()
retag <- function( .df )
{.df %>% mutate( Dataset = tagDataset(Dataset, Region) ) %>%
     mutate_at( "Task", tagTask ) %>% select( -Region )}

## Retrieves a slice of DGE results relevant to the requested task
## Removes MAYO results, due to the observed batch effect
DGEslice <- function( task="AC" )
{
    indexDGE() %>% filter( Task == task, Dataset != "MAYO" ) %>%
        resFromIndex() %>% retag() %>% unnest() %>%
        select( Dataset, Plate, LINCSID, Approval, IsToxic, Drug, Target, p_value )
}

## The composite score is defined by the geometric average of p-values
##   each of the MSBB regions is given a 0.25 weight to avoid MSBB dominating the score
##   p-values below 0.01 are thresholded at 0.005 to avoid "zero" issues
DGEcomposite <- function( task="AC" )
{
    ## Harmonic mean
    hmean <- function(v) {length(v)/sum(1/v)}
    
    DGEslice(task) %>% mutate_at( "p_value", pmax, 0.005 ) %>%
        spread( Dataset, p_value ) %>%
        mutate( HMP = pmap_dbl(list(MSBB10, MSBB22, MSBB36, MSBB44, ROSMAP), lift_vd(hmean)) ) %>%
        arrange( HMP )
}

## Pulls pre-computed DGEcomposite() output from Synapse
DGEcompositePre <- function( synID = "syn20617283" ) { syn_csv(synID) }

