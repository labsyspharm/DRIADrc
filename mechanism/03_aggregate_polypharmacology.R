library(tidyverse)
library(data.table)
library(here)

wd <- here("mechanism", "polypharmacology")
dir.create(wd, showWarnings = FALSE)


hmean <- function(v) {length(v)/sum(1/v)}

jak_members <- c(paste0("JAK", 1:3), "TYK2")
target_combo_significance <- read_rds(file.path(wd, "target_combo_significance.rds"))
target_combos <- read_rds(file.path(wd, "target_combos.rds"))

jak_combo_significance <- target_combo_significance %>%
  filter(Target_1 %in% jak_members | Target_2 %in% jak_members) %>%
  mutate(
    jak_member = map2(Target_1, Target_2, ~intersect(c(.x, .y), jak_members)),
    non_jak_member = map2(Target_1, Target_2, ~setdiff(c(.x, .y), jak_members))
  ) %>%
  filter(map_lgl(jak_member, ~length(.x) == 1)) %>%
  mutate_at(vars(jak_member, non_jak_member), as.character)

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
  file.path(wd, "target_combo_significance_aggregated.rds")
)
write_rds(
  jak_combo_significance,
  file.path(wd, "jak_combo_significance.rds")
)
