suppressMessages(library(tidyverse))
library(here)
library(ggbeeswarm)
library(synapser)
library(synExtra)

source(here("figures", "results.R"))
source(here("figures", "plot.R"))

# Compound name <-> LSP id map
cmpd_name_map <- syn("syn21586544") %>%
    read_csv( col_types=cols() ) %>%
    transmute( lspci_id, Drug = str_to_lower(name) )

# Counts and metadata
counts <- here("results", "normalized_counts_isg.csv") %>%
    read_csv(col_types=cols()) %>% rename(Sample=sample)
meta <- here("results", "deseq_meta.csv") %>%
    read_csv(col_types=cols()) %>%
    inner_join( cmpd_name_map, by="Drug" )

# TAS values
tas <- TASvalues() %>% distinct(lspci_id, entrez_symbol, tas) %>%
    filter( entrez_symbol == "TYK2", tas != 10 ) %>%
    mutate_at( "tas", factor, levels = c("1", "2", "3"))

# Plot ISG gene expression ----------------------------------------------------
###############################################################################

hand_picked_isg <- c("IFITM3", "IFI16", "IFI35", "PLSCR1")
X <- filter(counts, gene_name %in% hand_picked_isg) %>%
    inner_join(meta, by="Sample") %>%
    filter( Concentration == 10 ) %>%
    select( -Concentration ) %>%
    inner_join( tas, by="lspci_id" ) %>%
    mutate( Affinity = ifelse(tas == 1, "strong", "weak") )

pal <- c("1" = "#b2182b", "2" = "#ef8a62", "3" = "#fddbc7")
gg <- ggplot(X, aes(Affinity, count, color = tas)) + theme_bw() + theme_bold() +
    geom_quasirandom() + facet_wrap(vars(gene_name), scale = "free_y") +
    scale_color_manual( name = "TAS", values=pal ) +
    labs(x = "Binding affinity to TYK2", y = "Normalized count")

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Suppl4-", Sys.Date(), ".pdf") )
ggsave(fnOut, gg, width=7.25, height=4.75 )
