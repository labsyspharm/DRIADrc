library(tidyverse)
library(here)

wd <- here("mechanism", "polypharmacology")

targets <- read_rds(file.path(wd, "cotarget_significance.rds")) %>%
    arrange(padj) %>% filter(Target_Class == "cotarget") %>%
    select( symbol = Symbol, direction = Class, n, p, padj )

write_csv( head(targets,10), "Table2.csv" )

openxlsx::write.xlsx( targets, "TableS1.xlsx" )


