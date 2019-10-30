source( "../figures/results.R" )
source( "ecdf.R" )

gmean <- function(x) {exp(mean(log(x)))}

## Load composite scores for each drug / plate combination
## Combine scores for each drug from the two plates
S <- DGEcomposite() %>% select( Drug, LINCSID, HMP ) %>%
    group_by( Drug ) %>% summarize( LINCSID=unique(LINCSID), HMP=gmean(HMP) ) %>%
    arrange( HMP ) %>% mutate( Rank = 1:n() )

## Combine with the TAS information
P <- syn("syn20830941") %>% read_rds() %>% chuck("data", 2) %>%
    select(LINCSID = compound_id, entrez_symbol, tas) %>%
    left_join( S, ., by="LINCSID" ) %>% mutate(binding = (tas!=10))

## Compute ecdfs for each target / tas combination
## Always include 1 and max rank for visualization
ECDF <- P %>% group_by( entrez_symbol, tas ) %>%
    summarize( Ranks=list(Rank) ) %>%
    mutate( Rank = map(Ranks, ~unique(c(1L, .x, nrow(S)))),
           CumProb = map2( Ranks, Rank, ~ecdf(.x)(.y) )) %>%
    select( -Ranks ) %>% unnest( c(Rank, CumProb) ) %>% ungroup()
