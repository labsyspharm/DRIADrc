library(tidyverse)
library(here)
library(callr)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader("~/data/AMP-AD/DGE")

# Download data ----------------------------------------------------------------
###############################################################################T

isg <- syn("syn11629935") %>%
  read_lines()

# Wrangle DGE metadata and counts ----------------------------------------------
###############################################################################T

r(
  function() {
    library(here)
    source(here("analyses", "wrangle", "DGE1.R"))
    main()
  }
)
r(
  function() {
    library(here)
    source(here("analyses", "wrangle", "DGE2.R"))
    main()
  }
)

meta_raw <- tribble(
  ~experiment, ~fn,
  "dge1", here("DGE1-meta.csv"),
  "dge2", here("DGE2-meta.csv")
) %>%
  mutate(data = map(fn, read_csv))

counts_raw <- tribble(
  ~experiment, ~fn,
  "dge1", here("DGE1-counts.csv"),
  "dge2", here("DGE2-counts.csv")
) %>%
  mutate(data = map(fn, read_csv))

meta <- meta_raw %>%
  select(-fn) %>%
  unnest(data) %>%
  mutate_at(vars(Drug), str_to_lower) %>%
  mutate(
    Sample = paste(experiment, Well, sep = "_"),
    Drug = recode(Drug, `drug control` = "dmso") %>%
      str_to_lower(),
    Concentration = if_else(Drug == "dmso", 0, Concentration),
    Condition = paste(Drug, Concentration, sep = "_")
  )

counts <- counts_raw %>%
  select(-fn) %>%
  mutate(
    data = map2(
      data, experiment,
      function(d, n) {
        d %>%
          rename_at(vars(-HUGO), ~paste0(n, "_", .x))
      }
    )
  ) %>%
  pull(data) %>%
  reduce(inner_join, by = "HUGO")

count_sums <- counts %>%
  summarize_at(vars(-HUGO), sum) %>%
  gather("sample", "count") %>%
  arrange(count)

write_csv(
  meta,
  here("results", "deseq_meta.csv")
)

# Normalizing all compounds together -------------------------------------------
###############################################################################T

library(DESeq2)

# Removing samples with <30k counts
passing_samples <- count_sums %>%
  filter(count >= 30000) %>%
  pull(sample) %>%
  sort()

deseq_all <- DESeqDataSetFromMatrix(
  counts %>%
    select(HUGO, one_of(passing_samples)) %>%
    column_to_rownames("HUGO") %>%
    as.matrix(),
  meta %>%
    filter(Sample %in% passing_samples) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample"),
  design = ~ Condition + experiment
)

deseq_all <- estimateSizeFactors(deseq_all)

counts_all <- counts(deseq_all, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  as_tibble() %>%
  gather("sample", "count", -gene_name)

write_csv(
  counts_all,
  here("results", "normalized_counts.csv.gz")
)

# Taking only interferon stimulated genes --------------------------------------
###############################################################################T

counts_isg <- counts_all %>%
  filter(gene_name %in% isg)

write_csv(
  counts_isg,
  here("results", "normalized_counts_isg.csv")
)
