## Correspondence of replicates across DGE plates
##
## by Artem Sokolov

library( here )

source( here("figures","plot.R") )

f <- function( fn, colid )
{
    here( "results", fn ) %>%
        read_csv( col_types=cols() ) %>%
        mutate_at( "Drug", str_to_lower ) %>%
        filter( FDR <= 0.05 ) %>%
        select( Drug, Gene, !!rlang::ensym(colid) := logFC )
}

DFX1 <- f( "DGE1-dfx.csv.gz", logFC1 ) 
DFX2 <- f( "DGE2-dfx.csv.gz", logFC2 )
DFX <- inner_join(DFX1, DFX2, by=c("Drug","Gene"))
CR <- DFX %>% group_by( Drug ) %>% summarize( R=cor(logFC1, logFC2, method="sp") ) %>%
    mutate( Lbl=str_c(round(R,3), "  \n") )

gg <- ggplot(DFX, aes(x=logFC1, y=logFC2)) + theme_bw() + theme_bold() +
    geom_point(alpha=0.5, color="steelblue") +
    facet_wrap( ~Drug, nrow=2 ) +
    geom_abline( slope=1, intercept=0, lty="dashed" ) +
    geom_text( aes(label=Lbl), data=CR, x=Inf, y=-Inf, hjust=1, vjust=0.5, size=5 ) +
    xlab( "log(Fold Change) in Experiment 1" ) +
    ylab( "log(Fold Change) in Experiment 2" )

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Suppl2-", Sys.Date(), ".pdf") )
ggsave( fnOut, gg, width=7.25, height=5.5 )
