library( tidyverse )
library( here )

library( seriation )   # For optimal leaf reordering

## Plotting elements
source( here("figures","plot.R") )

## Load taimoto similarity data
TMT <- read_csv(here("results","tanimoto.csv"),
                col_types=cols()) %>%
    mutate_at( c("Drug1", "Drug2"), recode,
              "staurosporine aglycone" = "st. aglycone" )

## Perform hierarchicanl clustering with
##   optimal leaf reordering
DM <- TMT %>% select( Drug1, Drug2, Tanimoto ) %>%
    spread( Drug2, Tanimoto ) %>% as.data.frame() %>%
    column_to_rownames( "Drug1" ) %>% dist
lvl <- hclust(DM) %>% reorder(DM) %>% dendextend::order.hclust() %>% labels(DM)[.]

## Fix the order via factor levels
R <- TMT %>% mutate(Drug1 = factor(Drug1, lvl),
                    Drug2 = factor(Drug2, rev(lvl)))

colh <- "#B2182B"

## Plot the heatmap
hmp <- ggplot( R, aes(x=Drug1, y=Drug2, fill=Tanimoto) ) +
    geom_tile() + xlab("") + ylab("") +
    scale_fill_gradient(low="white", high=colh,
                        limits=c(0,1), guide=FALSE) +
    theme(axis.text.x = etxt(9, angle=90, vjust=0.5, hjust=1),
          axis.text.y = etxt(9))

## Create a legend
TMT0 <- TMT %>% filter( Drug1 != Drug2 )
LGD <- tibble(y = seq(0, 1, by=0.001))
ebl <- element_blank()
lgnd <- ggplot() + theme_bw() + theme_bold() +
    geom_hline( aes(yintercept=y, color=y), data=LGD ) +
    geom_boxplot( aes(y=Tanimoto), data=TMT0 ) +
    scale_color_gradient(low="white", high=colh, guide=FALSE) +
    scale_y_continuous( position="right", limits=c(0,1) ) +
    ##    facet_wrap( ~"Tanimoto\nSimilarity" ) +
    ggtitle( "Tanimoto\nSimilarity" ) +
    theme(panel.grid.major = ebl, panel.grid.minor = ebl,
          axis.ticks.x = ebl, axis.text.x = ebl,
          axis.title.y = ebl, panel.border = ebl,
          strip.background=ebl, strip.text = etxt(14) )

## Compose the final plot
lgc <- cowplot::plot_grid(NULL, lgnd, NULL, ncol=1,
                          rel_heights=c(.15, 1, .45))
gg <- cowplot::plot_grid(hmp, lgc,
                         rel_widths=c(.9, .1))

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Suppl6-", Sys.Date(), ".pdf") )
ggsave( fnOut, gg, width=9.5, height=9 )
