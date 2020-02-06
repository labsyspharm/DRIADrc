library(tidyverse)
library(here)
library(broom)
library(furrr)
library(IHW)

wd <- here("mechanism", "polypharmacology")
dir.create(wd, showWarnings = FALSE)

source( here("figures","results.R") )

gmean <- function(x) {exp(mean(log(x)))}

## Load composite scores for each drug / plate combination
## Combine scores for each drug from the two plates
S <- DGEcomposite() %>% 
  select( Drug, LINCSID, HMP ) %>%
  group_by( Drug ) %>%
  summarize( LINCSID=unique(LINCSID), HMP=gmean(HMP) ) %>%
  arrange( HMP ) %>%
  mutate( Rank = 1:n() )

## Combine with the TAS information
P <- read_rds( here("external", "tas_vector_annotated_long.rds") ) %>%
  chuck("data", 2) %>%
  select(LINCSID = compound_id, Target=entrez_symbol, TAS=tas) %>%
  left_join( S, ., by="LINCSID" ) %>%
  mutate(binding = (TAS!=10))

tas_binding_threshold <- 10L

generate_combos <- function(df, key) {
  ix_combs <- t(combn(nrow(df), 2))
  tibble(
    Target_1 = df$Target[ix_combs[, 1]],
    Target_2 = df$Target[ix_combs[, 2]],
    TAS_1 = df$TAS[ix_combs[, 1]],
    TAS_2 = df$TAS[ix_combs[, 2]]
  ) %>%
    mutate(
      binding_class = case_when(
        TAS_1 < tas_binding_threshold & TAS_2 < tas_binding_threshold ~ "T1_AND_T2",
        TAS_1 < tas_binding_threshold & TAS_2 >= tas_binding_threshold ~ "T1_AND_NOT_T2",
        TAS_1 >= tas_binding_threshold & TAS_2 < tas_binding_threshold ~ "T2_AND_NOT_T1",
        TAS_1 >= tas_binding_threshold & TAS_2 >= tas_binding_threshold ~ "NOT_T1_AND_NOT_T2",
        TRUE ~ "UNKNOWN"
      )
    )
}

## Find all pairwise combinations of targets for each drug
## and then aggregate this info per target combination
target_combinations <- P %>%
  arrange(LINCSID, Target) %>%
  group_by(LINCSID, HMP, Rank) %>%
  # Filter drugs with less than 2 annotated targets
  filter(n() >= 2) %>%
  group_modify(generate_combos) %>%
  ungroup() %>%
  group_nest(Target_1, Target_2) %>%
  mutate(
    data = map(data, ~split(.x, .x$binding_class))
  )

write_rds(
  target_combinations,
  file.path(wd, "target_combos.rds")
)

## Compare the composite score for two drug sets
## using the Wilcoxon rank-sum test
compare_drug_sets <- function(set_1, set_2) {
  n_1 <- nrow(set_1)
  n_2 <- nrow(set_2)
  # Require at least 3 drugs in both sets
  if (is.null(set_1) || is.null(set_2) || n_1 < 3 || n_2 < 3)
    return(NULL)
  wilcox.test(set_1$HMP, set_2$HMP, alternative = "two.sided", exact = TRUE, conf.int = TRUE, conf.level = 0.8) %>%
    tidy() %>%
    mutate(
      n_1 = n_1, n_2 = n_2, n = n_1 + n_2,
      HMP_median_1 = median(set_1$HMP), HMP_median_2 = median(set_2$HMP),
      set_1 = list(set_1), set_2 = list(set_2)
    )
}

## For each target pair, doing three comparisons
## Drugs that target one of the targets but not the other against
## drugs that target both.
## And drugs that target either vs drugs that target both.
process_target_pair <- function(Target_1, Target_2, data) {
  library(tidyverse)
  library(broom)
  res_and_not_1 <- list(
    "Target_1" = Target_1,
    "Target_2" = Target_2,
    "Comparison" = "T1_AND_NOT_T2",
    "Result" = list(compare_drug_sets(data[["T1_AND_NOT_T2"]], data[["T1_AND_T2"]]))
  )
  res_and_not_2 <- list(
    "Target_1" = Target_1,
    "Target_2" = Target_2,
    "Comparison" = "T2_AND_NOT_T1",
    "Result" = list(compare_drug_sets(data[["T2_AND_NOT_T1"]], data[["T1_AND_T2"]]))
  )
  res_xor <- list(
    "Target_1" = Target_1,
    "Target_2" = Target_2,
    "Comparison" = "T1_XOR_T2",
    "Result" = list(
        compare_drug_sets(
        bind_rows(data[["T1_AND_NOT_T2"]], data[["T2_AND_NOT_T1"]]),
        data[["T1_AND_T2"]]
      )
    )
  )
  res <- list(res_and_not_1, res_and_not_2, res_xor) %>%
    discard(~is.null(.x[["Result"]][[1]])) %>%
    bind_rows()
  if (nrow(res) > 0)
    return(unnest(res, Result))
  NULL
}

# process_target_pair(target_combinations$Target_1[[2]], target_combinations$Target_2[[2]], target_combinations$data[[2]])

plan(multisession(workers = 8))
target_combo_significance_raw <- future_pmap(
  target_combinations, process_target_pair
)

target_combo_significance <- bind_rows(target_combo_significance_raw) %>%
  mutate(
    padj = ihw(p.value, n, alpha = 0.1, covariate_type = "ordinal") %>%
      adj_pvalues()
  )

write_rds(
  target_combo_significance,
  file.path(wd, "target_combo_significance.rds")
)
