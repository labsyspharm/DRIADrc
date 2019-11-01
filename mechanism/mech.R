library( tidyverse )

## Synapse interface
synapser::synLogin()
syn <- synExtra::synDownloader("~/data/DRIAD/mech")

gmean <- function(x) {exp(mean(log(x)))}

calculate_ecdf <- function( drug_ranking, ...) {
    max_rank <- max( drug_ranking$Rank )
    as_tibble(drug_ranking) %>%
        group_by(...) %>%
        summarize( RR=list(Rank) ) %>%
        ## Always include 1 and max rank for visualization
        mutate(Rank = map(RR, ~unique(c(1L, .x, max_rank))),
               CumProb = map2( RR, Rank, ~ecdf(.x)(.y) )) %>%
        select( ..., Rank, CumProb ) %>%
        unnest( c(Rank, CumProb) ) %>%
        ungroup()
}

calculate_auc <- function(ecdf_df) {
    ecdf_df %>%
        arrange(Rank) %>%
        mutate(
            pos = scales::rescale(Rank, to = c(0, 1), from = c(min(Rank), max(Rank))),
            diff_pos = pos - lag(pos, 1),
            area = diff_pos*CumProb
        ) %>%
        with( sum(area, na.rm=TRUE) )
}

## Load composite scores for each drug / plate combination
## Combine scores for each drug from the two plates
S <- syn("syn20928503") %>% read_csv( col_types=cols() ) %>%
    select( Drug, LINCSID, HMP ) %>%
    group_by( Drug ) %>%
    summarize( LINCSID=unique(LINCSID), HMP=gmean(HMP) ) %>%
    arrange( HMP ) %>%
    mutate( Rank = 1:n() )

## Combine with the TAS information
P <- syn("syn20830941") %>% read_rds() %>%
    chuck("data", 2) %>%
    select(LINCSID = compound_id, Target=entrez_symbol, TAS=tas) %>%
    left_join( S, ., by="LINCSID" ) %>%
    mutate(binding = (TAS!=10))

## Compute correlations between performance ranking and binding affinities
## Some independent filtering may be in order to remove hopeless cases and improve
## FDR. Require three data points in all three binding TAS scores
cor.test2 <- function(...) {suppressWarnings( cor.test(...) )}
PAC <- P %>% filter( binding ) %>%
    group_by(Target) %>%
    filter(all(table(TAS)[c("1", "2", "3")] >= 3)) %>%
    summarize_at( c("HMP", "TAS"), list ) %>%
    transmute( Target=Target,
              COR = map2( HMP, TAS, cor.test2,
                         method="kendall", alternative="greater"),
              n_binding=map_int(HMP, length) ) %>%
    mutate_at( "COR", map, broom::tidy ) %>%
    mutate( p.value = map_dbl(COR, pull, "p.value"), COR=NULL )

## Compute ecdfs and areas under ecdfs for each target / tas combination
ECDF <- P %>% calculate_ecdf( Target, TAS ) %>%
    filter(TAS != "10") %>%
    mutate_at( "TAS", factor, levels = c("1", "2", "3", "10"),
              labels = c("1 (<100 nM)", "2 (100-999 nM)", "3 (1-10 uM)", "10 (>10 uM)")) %>%
    nest( ECDF = c(Rank, CumProb) ) %>%
    mutate( AUC = map_dbl(ECDF, calculate_auc) )

## Combine everything into a common results data frame
TAS_ECDF <- inner_join( ECDF, PAC, by="Target" )
save( TAS_ECDF, file="TAS-ecdfs.RData" )
