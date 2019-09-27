## Wrangling data from DGE1 plate that was re-sequenced
##  with the new protocol
##
## by Artem Sokolov

library( tidyverse )

## Set up the Synapse Downloader
synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/DGE1", ifcollision="overwrite.local" )

## Consolidates "Raw" concentration values in the origin spreadsheet into a set of common categories
map.conc <- function( v )
{
    setNames( c("0.3", "10", "1", "1", "10", "3", NA),
             c("0.299991", "9.99000001", "1_M", "0.99990001", "10_M", "2.932453597", "2_g/ml") )[v]
}

main <- function()
{
    ## Load the mapping of well barcodes to well IDs
    wmap <- syn( "syn17103737" ) %>% read_tsv( col_types=cols() ) %>%
        select( Well=well, Barcode=barcode )

    ## Load the mapping of ENSEMBL gene IDs to HUGO gene names
    gmap <- syn( "syn14236139" ) %>% read_csv( col_types=cols() ) %>%
        select( gene_id, gene_name, gene_biotype )

    ## Load the counts table (in long format)
    X <- syn( "syn18143723" ) %>% read_delim( " ", comment="%", col_types=cols() ) %>%
        rename( rowIndex=1, colIndex=2, Value=3 )

    ## Load the matching row and column names
    rn <- syn("syn18143724") %>% read_csv( col_names="gene_id", col_types=cols() ) %>%
        mutate( rowIndex=1:34947 )
    cn <- syn("syn18143725") %>% read_csv( col_names="Barcode", col_types=cols() ) %>%
        mutate( colIndex=1:386 ) %>% left_join( wmap, by="Barcode" ) %>% select( -Barcode )

    ## Data contains several duplicated gene_names
    ## For each gene name, we inspected the corresponding gene IDs and kept the one
    ##   with the highest number of total reads that map to it
    vDrop <- c( "ENSG00000213380", "ENSG00000268439", "ENSG00000234289", "ENSG00000284024",
               "ENSG00000280987", "ENSG00000284741", "ENSG00000258724", "ENSG00000206549",
               "ENSG00000130489", "ENSG00000158427" )
    
    ## Pull everything together into a single genes-by-wells matrix
    ## Map gene IDs to HUGO names and reduce to protein-coding space
    XX <- inner_join( X, rn, by="rowIndex" ) %>% inner_join( cn, by="colIndex" ) %>%
        select( Well, gene_id, Value ) %>% inner_join( gmap, by="gene_id" ) %>%
        filter( gene_biotype == "protein_coding", !(gene_id %in% vDrop) ) %>%
        select( -gene_biotype, -gene_id ) %>% spread( Well, Value, fill=0L ) %>%
        select( -P11_new2, -P11_old ) %>% rename( P11 = P11_new1, HUGO = gene_name )

    ## Load the metadata and adjust a couple of things by hand:
    ## - ruxolitinib+dsRNA is just ruxolitinib
    ## - metformin concentration in K16 is 0.3
    ## - paropanib is a typo; drug name is pazopanib
    Y <- syn( "syn11947020" ) %>% read_tsv( col_types=cols() ) %>%
        select( Well = `Dispensed well`, Drug = `Fluid name`, Concentration )%>%
        mutate_at( "Concentration", map.conc ) %>%
        mutate( Concentration = ifelse( Well == "K16", "0.3", Concentration ) ) %>%
        mutate_at( "Drug", recode, `ruxolitinib+dsRNA` = "ruxolitinib",
                  paropanib = "pazopanib", jq1 = "(+)-jq1" )
    
    ## Remove well P11 due to primer sequence issue
    ## Remove DMSO wells L09 and O14 because they cluster with Lipo controls
    vRemove <- c( "P11", "L09", "O14" )
    XX %>% select( -one_of(vRemove) ) %>% write_csv( "DGE1-counts.csv" )
    Y %>% filter( !(Well %in% vRemove) ) %>% write_csv( "DGE1-meta.csv" )
}
