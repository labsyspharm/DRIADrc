## Drug toxicity assayed by Deep-Dye-Drop
##
## by Artem Sokolov

source( "results.R" )
source( "plot.R" )

main <- function()
{
    ## Download the nuclei count data and compute summary statistics
    TX <- syn_csv( "syn18496657" ) %>% rename(Drug=`Fluid name`)
    TXm <- TX %>% group_by(Drug) %>% summarize_at( "Nuclei counts", median ) %>%
        mutate_at( "Drug", str_to_lower )

    ## Load the counts data from both DGE experiments, and the matching metadata
    XY <- list(DGE1="syn20820769", DGE1Y="syn20820771",
               DGE2="syn20820825", DGE2Y="syn20820826") %>% map( syn_csv )
    
    ## Compute total counts per well and map the values to the corresponding drugs
    DGE <- map( XY[c("DGE1","DGE2")], summarize_at, vars(-HUGO), sum ) %>%
        map( gather, Well, DGEcount ) %>%
        map2( XY[c("DGE1Y","DGE2Y")], inner_join, by="Well" ) %>%
        bind_rows(.id="DGE") %>% mutate_at( "Drug", str_to_lower ) %>%
        group_by( Drug, DGE ) %>% filter( Concentration == max(Concentration) ) %>%
        summarize_at( "DGEcount", median ) %>% ungroup() %>%
        inner_join( TXm, by="Drug" ) %>%
        mutate( Highlight = ifelse(`Nuclei counts` < 2400, "toxic", "no") )
    
    ## Plot the bi-modal distribution
    f <- function(v){rep("    ", length(v))}
    gg1 <- ggplot( TXm, aes(x=`Nuclei counts`) ) + theme_bw() + theme_bold() +
        geom_density() + geom_vline( xintercept=2400, lty="dashed" ) +
        xlim( c(0, 5000) ) + ylab( "Density" ) +
        scale_y_continuous( labels=f ) +
        theme( axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks=element_blank() )

    ## Plot the correspondence between DGE transcript counts and Deep-Dye-Drop assay
    gg2 <- ggplot( DGE, aes(y=DGEcount, x=`Nuclei counts`, color=Highlight) ) +
        geom_point() + theme_bw() + theme_bold() +
        scale_y_log10(labels=numform::f_denom) + 
        xlim( c(0, 5000) ) + ylab( "DGE Transcript Count" ) +
        geom_vline( xintercept=2400, lty="dashed" ) +
        scale_color_manual( values=c(toxic="tomato",no="black"), guide=FALSE )

    gg <- egg::ggarrange( gg1, gg2 )
    ggsave( str_c("FigS4-",Sys.Date(),".pdf"), gg, width=8.5, height=5.5 )
}

main()
