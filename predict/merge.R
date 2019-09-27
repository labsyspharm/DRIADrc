## Combines results from individual runs into a single data frame
##
## by Artem Sokolov

library( tidyverse )

readRData <- function(fn) {load(fn); RR}

## Load all the new files
allRes <- crossing( Plate = c("DGE1","DGE2"),
                   Dataset = c("ROSMAP", "MSBB10", "MSBB22", "MSBB36", "MSBB44") ) %>%
    mutate( fn = str_c("results/", Plate, "-", Dataset, ".RData"),
           Results = map(fn, readRData), fn=NULL ) %>%
    unnest() %>% rename( Drug = Set )

## Save the file into a new stand-alone .RData
save( allRes, file=str_c("results/results-",Sys.Date(),".RData") )

