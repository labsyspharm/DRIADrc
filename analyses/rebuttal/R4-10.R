library( tidyverse )
library( here )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/rebuttal", ifcollision="overwrite.local" )

## IDs of the 80 drugs profiled by DGE
load( here("results","results-2019-10-06.RData") )
idDrugs <- unique( allRes$Drug )

## Mapping of LINCS IDs to lspci_id and drug names
M1 <- syn( "syn11801537" ) %>% read_csv(col_types=cols()) %>%
    mutate_at( "name", str_to_lower ) %>%
    select( LINCSID = lincs_id, Drug = name ) %>%
    filter( LINCSID %in% idDrugs )

M <- syn("syn21094266") %>% read_csv() %>%
    select( LSPID = lspci_id, LINCSID = hmsl_id ) %>%
    inner_join(M1)

## Morgan fingerprints
FP <- syn("syn21614996") %>% read_csv() %>%
    filter( fingerprint_type == "morgan_normal" ) %>%
    select( LSPID = lspci_id, Morgan = fingerprint ) %>%
    inner_join( M, ., by="LSPID" ) %>%
    select( Drug, Morgan )

## Pairwise tanimoto similarity
TMT <- crossing(rename_all(FP, str_c, "1"),
                rename_all(FP, str_c, "2")) %>%
    mutate( Tanimoto = map2_dbl(Morgan1, Morgan2, morgancpp::tanimoto) )
write_csv( TMT, here("results", "tanimoto.csv") )
