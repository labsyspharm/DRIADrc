library( tidyverse )
library( here )

R <- readRDS( here("results/allpred.rds") ) %>%
    mutate( Mtd = recode(Method,
                         lgr = "Logistic\nRegression",
                         svm = "Support\nVector\nMachine",
                         xgb = "Boosted\nRandom\nForest",
                         nn  = "Neural\nNetwork") ) %>%
    mutate( Mtd = factor(Mtd, unique(Mtd)) )

ggplot( R, aes(x=Mtd, y=AUC) ) + theme_bw() +
    geom_boxplot() + xlab("") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggsave( "test.pdf" )
