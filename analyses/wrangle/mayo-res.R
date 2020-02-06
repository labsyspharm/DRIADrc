## Wrangles MAYO results from old BTR runs
##
## by Artem Sokolov

library( tidyverse )

synapser::synLogin()

synParent <- function( synid )
{
    s <- synapser::synGet(synid, downloadFile=FALSE)
    s$properties$parentId
}

syn_csv <- function( synid )
{
    gp <- synParent( synid ) %>% synParent() %>% synExtra::synName()
    p <- file.path( "~/data/DRIAD/BTR", gp )
    syn <- synExtra::synDownloader(p)
    syn(synid) %>% read_csv(col_types = cols())
}

indexDGE <- function()
{
    structure(list(Dataset = c("MAYO", "MAYO", "MAYO", "MAYO", "MAYO", "MAYO"),
                   Task = c("AB", "AB", "AC", "AC", "BC", "BC"),
                   id = c("syn18201042", "syn18143076", "syn18201045",
                          "syn18143082", "syn18201026", "syn18143046"),
                   Plate = c("DGE1", "DGE2", "DGE1", "DGE2", "DGE1", "DGE2")),
              class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -6L))
}

indexBK <- function()
{
    structure(list(Dataset = c("MAYO", "MAYO", "MAYO"),
                   Task = c("BC", "AB", "AC"),
                   id = c("syn15572112", "syn15570670", "syn15572002")),
              row.names = c(NA, -3L), class = c("tbl_df", "tbl", "data.frame"))
}

main <- function()
{
    ## Pull the background data off Synapse
    BK0 <- indexBK() %>% mutate( BK = map(id, syn_csv), id=NULL )

    ## Wrange background values
    BK <- BK0 %>% mutate_at( "BK", map, summarize_all, list ) %>% unnest() %>%
        gather( "Size", "BK", `10`:`1000` ) %>% mutate_at( "Size", as.integer )

    ## Load DGE results
    X0 <- indexDGE() %>% mutate( Results = map(id, syn_csv), id=NULL )

    ## Combine everything into a single data frame
    s <- rlang::exprs( Drug=description, nFeats=intersect, AUC )
    mayoRes <- X0 %>% mutate_at( "Results", map, select, !!!s ) %>% unnest() %>%
        mutate( Size=round(nFeats/10)*10 ) %>%
        inner_join(BK, c("Dataset", "Task", "Size")) %>%
        mutate( pval=map2_dbl(AUC, BK, ~mean(.y>=.x)), Size=NULL )

    save( mayoRes, file=stringr::str_c("MAYO-",Sys.Date(),".RData") )
}
