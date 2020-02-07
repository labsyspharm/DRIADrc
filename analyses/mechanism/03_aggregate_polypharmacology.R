library(tidyverse)
library(data.table)
library(philentropy)
library(magrittr)
library(furrr)
library(here)

wd <- here("results")

hmean <- function(v) {length(v)/sum(1/v)}

jak_members <- c(paste0("JAK", 1:3), "TYK2")
target_combo_significance <- read_rds(file.path(wd, "target_combo_signif-2020-02-07.rds"))
target_combos <- read_rds(file.path(wd, "target_combos-2020-02-07.rds"))

combo_direction_classes <- target_combo_significance %>%
  mutate(
    Class = case_when(
      p.value > 0.05 ~ "NS",
      estimate > 0 ~ "synergistic",
      estimate < 0 ~ "antagonistic",
      TRUE ~ "NA"
    )
  )

## Check if there are cases where a single target combination has both
## antagonistic and synergistic interactions
combo_direction_classes_sums <- combo_direction_classes %>%
  select(Target_1, Target_2, Comparison, Class) %>%
  spread(Comparison, Class) %>%
  group_by_at(vars(3:5)) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  arrange(desc(n))
## This only happens in very few cases, filtering those below

combo_combined <- combo_direction_classes %>%
  as.data.table() %>%
  ## Require that at least two drug sets out of three were compared against the
  ## AND drug set
  ## Also require that all significant effects go in one direction
  ## (only synergistic or only antagonistic, not both)
  .[
    ,
    if(!("synergistic" %in% Class & "antagonistic" %in% Class) & .N >= 2) .SD,
    keyby = list(Target_1, Target_2)
  ] %>%
  .[
    ,
    Class := sort(factor(Class, levels = c("synergistic", "antagonistic", "NS")))[[1]],
    keyby = list(Target_1, Target_2)
  ] %>%
  as_tibble() %>%
  select(Target_1, Target_2, Class, Comparison, p.value) %>%
  spread(Comparison, p.value) %>%
  ## Combine results of drug set tests against the AND drug set using
  ## the harmonic mean
  mutate(
    Combined = pmap_dbl(
      list(T1_AND_NOT_T2, T1_XOR_T2, T2_AND_NOT_T1),
      function(...) {
        hmean(na.omit(as.double(list(...))))
      }
    )
  )

write_rds(
  combo_combined,
  file.path(wd, paste0("tc_signif_agg-",Sys.Date(),".rds")),
  compress="gz"
)

calculate_jaccard <- function(df) {
  # browser()
  mat <- df %>%
    mutate(
      Dummy = 1L,
      Drugs = map2(
        set_1, set_2,
        ~union(.x$LINCSID, .y$LINCSID)
      )
    ) %>%
    unnest_longer(Drugs, values_to = "LINCSID") %>%
    select(LINCSID, Target_Combo, Dummy) %>%
    spread(LINCSID, Dummy, fill = 0L) %>%
    column_to_rownames("Target_Combo") %>%
    as.matrix() %>%
    is_greater_than(0L)
  jaccard = suppressMessages(distance(mat, method = "jaccard")) %>%
    set_rownames(rownames(mat)) %>%
    set_colnames(rownames(mat))
  # Return jaccard distance, instead of similarity
  1 - jaccard
}

# Aggregating the p-values of individual target-cotarget associations
# with a modification of Fisher's method called Brown's method, implemented in the
# R package [poolR](https://github.com/ozancinar/poolR/). As input, a covariance
# matrix between the pairwise comparisons is needed. Using the jaccard distance
# matrix measuring the overlap between drug sets for this purpose. Can use function
# mvnconv to convert it to covariance matrix.
calculate_pooled_p <- function(df, jaccard) {
  padj <- poolr::fisher(
     df$p.value,
     adjust = "generalized",
     seed = 1,
     type = 2,
     R = jaccard
  )
  p <- poolr::fisher(df$p.value, adjust = "none")
  tibble(
    n = nrow(df),
    padj = padj$p,
    p = p$p
  )
}

## Finding most significant co-targets
plan(multisession(workers = 8))
cotarget_significance <- target_combo_significance %>%
  mutate(Class = ifelse(estimate > 0, "synergistic", "antagonistic")) %>%
  ## Only include those that passed filters and only consider AND NOT tests here
  filter(Comparison != "T1_XOR_T2") %>%
  semi_join(combo_combined, by = c("Target_1", "Target_2")) %>%
  ## Make Target_1 always the included one and Target_2 the omitted one
  mutate(
    target = if_else(Comparison == "T2_AND_NOT_T1", Target_2, Target_1),
    cotarget = if_else(Comparison == "T2_AND_NOT_T1", Target_1, Target_2),
    Target_Combo = paste(target, cotarget, sep = "|")
    # Drugs = map2(
    #   set_1, set_2,
    #   ~union(.x$LINCSID, .y$LINCSID)
    # )
  ) %>%
  gather("Target_Class", "Symbol", target, cotarget) %>%
  group_nest(Class, Target_Class, Symbol, .key = "target_combinations", keep = TRUE) %>%
  filter(map_lgl(target_combinations, ~nrow(.x) >= 3)) %>%
  mutate(
    jaccard = future_map(target_combinations, calculate_jaccard, .progress = TRUE),
    pooled = map2(target_combinations, jaccard, calculate_pooled_p)
  ) %>%
  unnest(pooled)

cotarget_significance %>%
    select( Target_Class, Symbol, Class, n, p, padj ) %>%
    write_rds(
        file.path(wd, paste0("cotarget_signif-",Sys.Date(),".rds")),
        compress="gz"
    )
