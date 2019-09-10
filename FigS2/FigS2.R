## FigS2: Performance on background gene sets
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

FigS2 <- function()
{
    BK <- indexBackground() %>% bkFromIndex() %>% retag() %>% unnest()
    TXT <- BK %>% select(Task) %>% distinct

    ggplot( BK, aes(x=Size, y=AUC, color=Dataset) ) + theme_bw() + theme_bold() +
        facet_wrap( ~Task ) + geom_smooth( se=FALSE ) +
        geom_text( aes(label=Task), data=TXT, color="black", fontface="bold",
                  size=5, x=log10(1000), y=0.45, hjust=1, vjust=0.5 ) +
        scale_color_manual( values=dsPal() ) + #guides( color=FALSE ) +
        theme( strip.background = element_blank(),
              strip.text = element_blank() ) +
        scale_x_log10(name="Number of genes in set") +
        ggsave( str_c("FigS2-",Sys.Date(),".pdf"), width=9, height=3 ) +
        ggsave( str_c("FigS2-",Sys.Date(),".png"), width=9, height=3 )
}
