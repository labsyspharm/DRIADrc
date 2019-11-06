source( "results.R" )
source( "plot.R" )

## All ECDF data for p < 0.05, sorted by p value
ECDF <- syn_csv( "syn21094322" ) %>% filter( p.value < 0.05 ) %>%
    arrange( p.value ) %>% mutate_at( c("Target", "TAS"), as_factor )

TST <- syn_csv( "syn21094322" ) %>% select( Target, TAS, AUC, p.value ) %>% distinct %>%
    filter( !grepl( "10 \\(", TAS) ) %>% group_by( Target ) %>%
    filter( AUC[1] > 0.5 ) %>%
    select( Target, p.value ) %>% distinct %>% ungroup %>%
    mutate( FDR = p.adjust(p.value, "BH") ) %>% arrange(p.value)

## A single p value entry for each facet
ECDFp  <- ECDF %>% mutate_at( "p.value", ~as.character(round(.x,3)) ) %>%
    select( Target, p.value ) %>% distinct

## Plotting elements and palettes
max_rank <- max(ECDF$Rank)
cpal <- set_names( c("#b2182b", "#ef8a62", "#fddbc7", "#d9d9d9"), levels(ECDF$TAS) )
ltpal <- set_names( c("solid", "solid", "solid", "dotted"), levels(ECDF$TAS) )

## Plot each target on a separate facet
gg <- ggplot( ECDF, aes(x=Rank, y=CumProb) ) +
    theme_bw() + theme_bold() +
    facet_wrap( ~Target, nrow=3 ) +
    geom_step( aes(color=TAS, linetype=TAS), size=1.1 ) +
    geom_text( aes(label=p.value), data=ECDFp, fontface="bold",
              x=max_rank, y=0, vjust="inward", hjust="inward" ) +
    xlab( "Drug Rank" ) + ylab( "Cumulative Probability" ) +
    scale_color_manual( values = cpal ) +
    scale_linetype_manual( values = ltpal ) +
    theme( panel.grid.major=element_blank(),
          panel.grid.minor=element_blank() )

gg2 <- lemon::reposition_legend( gg, "center", panel= "panel-5-3" )
ggsave( "test.pdf", gg2, width=10, height=6 )

##gt <- cowplot::plot_to_gtable(gg) %>% gtable_filter("panel") %>%
##    with( set_names(grobs, layout$name) ) %>% keep( ~identical(.x,zeroGrob()) ) %>%
##    names
