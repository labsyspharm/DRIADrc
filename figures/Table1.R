## Table1: List of profiled drugs and their associated scores
##
## by Artem Sokolov

library( tidyverse )

read_csv(here("results","DGE-composite.csv"), col_types=cols()) %>%
    write_csv("Table1.csv")
