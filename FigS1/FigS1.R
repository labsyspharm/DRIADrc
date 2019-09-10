## Performance on cell-type specific gene sets
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

FigS1 <- function()
{
    load(syn( "syn18565498" ))

    v <- c( AB = "A-vs-B", AC = "A-vs-C", BC = "B-vs-C" )
    RS <- RS %>% mutate_at( "Task", recode, !!!v )
    RBK <- RBK %>% mutate_at( "Task", recode, !!!v )
    
    SBK <- RBK %>% group_by( Task ) %>% summarize_at( "AUC", list(BK=list) ) %>%
        inner_join( RS, ., by="Task" ) %>% mutate( pval = map2_dbl(AUC, BK, ~mean(.x < .y)) ) %>%
        mutate( sig = ifelse( pval < 0.05, "Yes", "No" ) ) %>% arrange( pval ) 
   
    ## Additional plotting elements
    xrng <- bind_rows(RBK, SBK) %>% pull( AUC ) %>% range
    
    ## Plot everything together
    ggplot() + theme_bw() + facet_wrap( ~Task, ncol=1 ) + guides( color=FALSE ) +
        ggthemes::scale_fill_few() + theme_bold() + ylab( "Density" ) +
        scale_x_continuous( limits = c(1,1.1) * xrng ) + ylim( c(0,15) ) +
        ggrepel::geom_text_repel( aes(x=AUC, label=Name, color=sig), SBK, y=3,
                                 nudge_y=100, fontface="bold" ) +
        geom_density( aes(x=AUC, fill=Task), RBK, alpha=0.65, lwd=1 ) +
        geom_segment( aes(x=AUC, xend=AUC, color=sig), SBK, y=0, yend=3, lwd=1 ) +
        scale_color_manual( values=c("Yes" = "red", "No" = "black") ) +
        theme( strip.background = element_blank(), strip.text.x = element_blank(),
              legend.position=c(0.98,0.98), legend.justification=c(1,1),
              legend.background=element_rect(fill=NA), legend.title.align=0.75 ) +
        ggsave( str_c("FigS1-",Sys.Date(),".pdf"), width=8, height=5 )
}

