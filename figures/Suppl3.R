## Drug toxicity assayed by Deep-Dye-Drop
##
## by Artem Sokolov

library( here )

source( here("figures","results.R") )
source( here("figures","plot.R") )

theme_noaxes <- function() {
    theme(axis.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank()) }

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
    
## Plot the bi-modal distribution of nuclei counts
gg1 <- ggplot( TXm, aes(x=`Nuclei counts`) ) + theme_bw() + theme_bold() +
    geom_density() + geom_vline( xintercept=2400, lty="dashed" ) +
    xlim( c(0, 5000) ) + theme_noaxes()

## Plot the correspondence between DGE transcript counts and Deep-Dye-Drop assay
gg2 <- ggplot( DGE, aes(y=DGEcount, x=`Nuclei counts`, color=Highlight) ) +
    geom_point() + theme_bw() + theme_bold() +
    scale_y_log10(labels=numform::f_denom, limits=c(1000, 1e6)) + 
    xlim( c(0, 5000) ) + ylab( "DGE Transcript Count" ) +
    geom_vline( xintercept=2400, lty="dashed" ) +
    scale_color_manual( values=c(toxic="tomato",no="black"), guide=FALSE )

## Plot the bi-modal distribution of DGE transcript counts
gg3 <- ggplot( DGE, aes(x=DGEcount) ) + theme_bw() + theme_bold() + theme_noaxes() +
    geom_density() + coord_flip() + scale_x_log10(limits=c(1000, 1e6))

## Put everything together
gg <- egg::ggarrange( gg1, cowplot::ggdraw(), gg2, gg3, ncol=2, heights=c(1,2), widths=c(3,1) )

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Suppl3-", Sys.Date(), ".pdf") )
ggsave( fnOut, gg, width=8.25, height=4.75 )

