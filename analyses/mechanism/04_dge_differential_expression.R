library(tidyverse)
library(here)
library(callr)
library(furrr)

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

# Differential expression per drug using rank-transformed dose -----------------
###############################################################################T

controls <- meta %>%
  filter(Drug == "dmso") %>%
  mutate(Dose = 0)

deseq_input_per_drug <- meta %>%
  filter(Drug != "dmso", !is.na(Concentration)) %>%
  group_nest(Drug, .key = "meta", keep = TRUE) %>%
  mutate(
    meta = map(
      meta,
      ~.x %>%
        mutate(
          Dose = factor(
            Concentration,
            levels = sort(unique(Concentration))
          ) %>%
            as.integer() %>%
            as.double()
        ) %>%
        bind_rows(controls) %>%
        arrange(Sample) %>%
        column_to_rownames("Sample")
    ),
    counts = map(
      meta,
      ~counts %>%
        select(HUGO, one_of(rownames(.x))) %>%
        column_to_rownames("HUGO") %>%
        as.matrix()
    )
  )

library(DESeq2)



deseq_objects_per_drug <- deseq_input_per_drug %>%
  mutate(
    data = map2(
      counts, meta,
      ~DESeqDataSetFromMatrix(
        .x, .y,
        design = ~ Dose + experiment
      )
    )
  )

plan(multisession(workers = 6))
deseq_de_per_drug <- deseq_objects_per_drug %>%
  transmute(
    Drug,
    data = future_map(
      data,
      DESeq,
      parallel = FALSE,
      .progress = TRUE
    )
  )

write_rds(
  deseq_de_per_drug,
  here("results", "deseq_de_per_drug.rds"),
  compress = "gz"
)

extract_result <- function(de, name) {
  res <- results(de, name = name, parallel = FALSE)
  shrunken <- lfcShrink(de, coef = name, res = res, type = "apeglm", parallel = FALSE)
  shrunken %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    as_tibble() %>%
    left_join(
      res %>%
        as.data.frame() %>%
        rownames_to_column("gene_symbol") %>%
        select(gene_symbol, log2FoldChange_MLE = log2FoldChange, lfcSE_MLE = lfcSE),
      by = "gene_symbol"
    ) %>%
    select(gene_symbol, everything()) %>%
    arrange(padj)
}

deseq_res_per_drug <- deseq_de_per_drug %>%
  mutate(
    data = future_map(
      data,
      extract_result,
      name = "Dose",
      .progress = TRUE
    )
  )

write_rds(
  deseq_res_per_drug,
  here("results", "deseq_res_per_drug.rds"),
  compress = "gz"
)

# Normalizing all compounds together -------------------------------------------
###############################################################################T

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

write_rds(
  counts_all,
  here("results", "normalized_counts.csv.gz")
)

