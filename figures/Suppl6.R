## Relative performance of machine learning methods
##
## by Artem Sokolov

library( tidyverse )
library( here )

R <- readRDS( here("results/allpred.rds") ) %>%
    mutate( Mtd = recode(Method,
                         lgr = "Logistic\nRegression",
                         svm = "Support\nVector\nMachine",
                         xgb = "Boosted\nRandom\nForest",
                         nn  = "Neural\nNetwork") ) %>%
    mutate( Mtd = factor(Mtd, unique(Mtd)) )

gg <- ggplot( R, aes(x=Mtd, y=AUC) ) + theme_bw() +
    geom_boxplot() + xlab("") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Suppl5-", Sys.Date(), ".pdf") )
ggsave( fnOut, gg, width=6, height=4 )
