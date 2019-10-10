## FigS2: Performance on background gene sets
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

FigS1 <- function()
{
    load( syn("syn20835686") )
    load( syn("syn20928450") )

    v <- c( AB = "A-vs-B", AC = "A-vs-C", BC = "B-vs-C" )
    BK0 <- allRes %>% mutate( nFeats = map_int(Feats, length), Feats=NULL ) %>%
        bind_rows( mayoRes ) %>% mutate_at( "Task", recode, !!!v ) %>%
        select( -Plate, -Drug, -AUC, -pval, Vals=BK ) %>%
        arrange( nFeats ) %>% unnest() %>%
        group_by( Dataset, Task, nFeats ) %>% summarize_at( "Vals", list ) %>% ungroup()
    TXT <- BK0 %>% select(Task) %>% distinct

    BK <- BK0 %>% mutate( AUC=map_dbl(Vals,mean), s=map_dbl(Vals,sd) )

    f <- function( p, k, a, step=0.5 )
    {
        p + geom_ribbon( aes(ymin=AUC+k*s, ymax=AUC+(k+step)*s), alpha=a ) +
            geom_ribbon( aes(ymin=AUC-k*s, ymax=AUC-(k+step)*s), alpha=a )
    }
    
    pl <- ggplot( BK, aes(x=nFeats) ) + theme_bw() + theme_bold() +
        facet_wrap( ~Task ) + geom_line( aes(y=AUC, color=Dataset) ) +
##        geom_text( aes(label=Task), data=TXT, color="black", fontface="bold",
##                  size=5, x=log10(300), y=0.45, hjust=1, vjust=0.5 ) +
        scale_color_manual( values=dsPal() ) +
        theme( strip.background = element_blank(), strip.text = element_blank() ) +
        scale_x_log10(name="Number of genes in set")

    nSD <- seq(0,2.5,by=0.5)
    va <- c(0.4, 0.3, 0.2, 0.1, 0.05, 0.025)
    pl <- reduce2( nSD, va, f, .init=pl )

    Y <- mayoRes %>% select( -AUC, AUC=BK ) %>% unnest()
    ggplot( Y, aes(x=nFeats, y=AUC) ) + geom_smooth(se=FALSE) +
        facet_wrap( ~Task )
    
#        ggsave( str_c("FigS1-",Sys.Date(),".pdf"), width=9, height=3 ) +
#        ggsave( str_c("FigS1-",Sys.Date(),".png"), width=9, height=3 )
}
