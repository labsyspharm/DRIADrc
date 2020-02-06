## Wrangling of the second DGE experiment
##
## by Artem Sokolov

library( tidyverse )

## Set up the Synapse Downloader
synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/DGE2", ifcollision="overwrite.local" )

main <- function()
{
    ## Load the mapping of well barcodes to well IDs
    wmap <- syn( "syn17103737" ) %>% read_tsv( col_types=cols() ) %>% select( barcode, well )

    ## Load the mapping of ENSEMBL gene IDs to HUGO gene names
    gmap <- syn( "syn14236139" ) %>% read_csv( col_types=cols() ) %>%
        select( gene_id, gene_name, gene_biotype )
    
    ## Load the counts table (in long format)
    X <- syn("syn17103650") %>% read_delim( " ", comment="%", col_types=cols() ) %>%
        rename( rowIndex=1, colIndex=2, Value=3 )

    ## Load the matching row and column names
    rn <- syn("syn17103651") %>% read_csv( col_names="gene_id", col_types=cols() ) %>%
        mutate( rowIndex=1:34947 )
    cn <- syn("syn17103652") %>% read_csv( col_names="barcode", col_types=cols() ) %>%
        mutate( colIndex=1:386 ) %>% left_join( wmap, by="barcode" ) %>% select( -barcode )

    ## The following gene names map to multiple gene ids:
    ##   PINX1, TMSB15B, COG8, MATR3, HSPA14, SCO2, EMG1, H2BFS, PDE11A
    ## For each gene name, we inspected the corresponding gene IDs and kept the one
    ##   with the highest number of total reads that map to it
    vDrop <- c( "ENSG00000258724", "ENSG00000158427", "ENSG00000213380",
               "ENSG00000280987", "ENSG00000284024", "ENSG00000130489",
               "ENSG00000268439", "ENSG00000234289", "ENSG00000284741" )
    
    ## Combine everything into a common matrix
    ## Map gene IDs to gene names
    ## P11 is tagged as P11_new1 (P11_new2 and P11_old are not used and should have low counts)
    XX <- inner_join( X, rn, by="rowIndex" ) %>% inner_join( cn, by="colIndex" ) %>%
        select( Well = well, gene_id, Value ) %>% inner_join( gmap, by="gene_id" ) %>%
        filter( gene_biotype == "protein_coding", !(gene_id %in% vDrop) ) %>%
        select( -gene_biotype, -gene_id ) %>% spread( Well, Value, fill=0L ) %>%
        select( -P11_new2, -P11_old ) %>% rename( P11 = P11_new1, HUGO = gene_name )

    ## Wrangle the metadata matrix
    ## Fix kw2449 and AZD-1480 to match the LINCS database
    ## M03, E01 and M21 are dsRNA wells; label them so
    Y <- syn("syn17114439") %>% read_csv( col_types=cols() ) %>%
        rename( Well=1, Drug=2 ) %>% mutate_at( "Concentration", round, 1 ) %>%
        mutate_at( "Drug", replace_na, "DMSO" ) %>%
        mutate_at( "Drug", recode, `KW-2449` = "KW2449", AZD1480 = "AZD-1480" ) %>%
        mutate( Drug = ifelse(Well %in% c("E01","M21","M03"), "dsRNAmi", Drug) )

    ## Write to files
    Y %>% write_csv( "DGE2-meta.csv" )
    XX %>% write_csv( "DGE2-counts.csv" )
}
