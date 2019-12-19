library(tidyverse)
library(here)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data/DRIAD/mech")

wd <- here("mechanism", "polypharmacology")

targets <- read_rds(file.path(wd, "cotarget_significance.rds")) %>%
    arrange(padj) %>% filter(Target_Class == "cotarget") %>%
    select( symbol = Symbol, direction = Class, n, p, padj )

write_csv( targets, "TableS1.csv" )

