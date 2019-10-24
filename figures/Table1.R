## Table1: List of profiled drugs and their associated scores
##
## by Artem Sokolov

source( "results.R" )

DGEcomposite() %>% write_csv("Table1.csv")
