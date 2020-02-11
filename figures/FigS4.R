library(tidyverse)
library(here)
library(ggbeeswarm)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader("~/data/AMP-AD/DGE")

source(here("figures", "plot.R"))

# Download data ----------------------------------------------------------------
###############################################################################T

cmpd_name_map <- syn("syn21586544") %>%
  read_csv() %>%
  mutate_at(vars(name), str_to_lower)

tas <- syn("syn20830942") %>%
  read_csv() %>%
  filter(fp_name == "morgan_normal") %>%
  distinct(lspci_id, entrez_symbol, tas)

# Read and wrangle counts ------------------------------------------------------
###############################################################################T

meta <- here("results", "deseq_meta.csv") %>%
  read_csv()

counts <- here("results", "normalized_counts_isg.csv") %>%
  read_csv()

cmpds_mapped <- meta %>%
  distinct(Drug) %>%
  inner_join(cmpd_name_map, by = c("Drug" = "name")) %>%
  distinct()

cmpds_tas <- cmpds_mapped %>%
  inner_join(
    tas,
    by = "lspci_id"
  ) %>%
  mutate_at(vars(tas), factor, levels = c("1", "2", "3", "10"))

meta_mapped <- meta %>%
  inner_join(
    cmpds_mapped,
    by = "Drug"
  )

# Plot ISG gene expression -----------------------------------------------------
###############################################################################T

hand_picked_isg <- c("IFITM3", "IFI16", "IFI35", "PLSCR1")

tyk2_hand_picked_violin <- counts %>%
  filter(gene_name %in% hand_picked_isg) %>%
  inner_join(
    meta_mapped,
    by = c("sample" = "Sample")
  ) %>%
  inner_join(
    cmpds_tas %>%
      filter(entrez_symbol == "TYK2", tas != 10),
    by = "lspci_id"
  ) %>%
  ggplot(aes(tas, count, color = tas)) +
  geom_quasirandom() +
  facet_wrap(vars(gene_name), scale = "free_y") +
  scale_color_manual(
    values = c("1" = "#b2182b", "2" = "#ef8a62", "3" = "#fddbc7", "10" = "grey85"),
    na.value = "#d9d9d9",
  ) +
  labs(x = "TYK2 TAS", y = "Normalized count") +
  guides(color = FALSE) +
  theme_bold()

ggsave(
  here("figures", str_c("FigS4-", Sys.Date(), ".pdf")),
  tyk2_hand_picked_violin,
  width=6, height=4
)

