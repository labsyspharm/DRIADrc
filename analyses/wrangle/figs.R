## Pulls relevant results slices needed for figures into here/results
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))
library( here )

synapser::synLogin()

## Loaders
syn <- synExtra::synDownloader( "~/data/DRIAD/figs", ifcollision="overwrite.local" )
syn_csv <- function( synid ) { syn(synid) %>% read_csv(col_types = cols()) }

## Raw TAS values
TASvalues <- function()
{
    synapser::synGet( "syn20830941", version=2, ifcollision="overwrite.local",
                     downloadLocation="~/data/DRIAD/figs" ) %>%
        pluck( "path" ) %>% read_rds() %>% chuck("data", 2 ) %>%
        write_csv( here("results","TAS-values.csv.gz") )
}
