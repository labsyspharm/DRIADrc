## Performance on background gene sets
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

panelB <- function()
{
    library( ggridges )

    ## A closer look some of the following drugs (#non-MAYO datasets where p < 0.05):
    ##  lapatinib (3), baricitinib (MAYO+2), nvp-tae684 (4), torin1 (3)

    ## Identify performance for the drugs of interest
    vDrugs <- c("lapatinib", "nvp-tae684")
    X <- indexDGE() %>% filter( Task == "AC" ) %>% resFromIndex() %>%
        retag() %>% unnest() %>% filter( Drug %in% vDrugs ) %>%
        mutate( Label = ifelse(p_value == 0, " p<0.01", str_c(" p=",p_value)),
               Size = round(Size/10)*10, Highlight=ifelse(p_value <= 0.05, "yes", "no") ) %>%
        select( -Task, -Plate, -Approval, -IsToxic, -LINCSID, -p_value )

    ## Get the matching background
    XX <- indexBackground() %>% filter( Task == "AC" ) %>% bkFromIndex() %>%
        unnest() %>% nest( AUC, .key="BK" ) %>% retag() %>% select( -Task ) %>%
        inner_join(X, ., by=c("Dataset","Size")) %>%
        mutate( Name = glue::glue("{Drug} ({Target})") ) %>%
        select( -Drug, -Target, -Size ) %>%
        arrange( Dataset ) %>% mutate_at( "Dataset", as_factor )
    BK <- XX %>% select( Name, Dataset, BK ) %>% unnest

    ## Generate the ridge plots
    ggplot( BK, aes(x=AUC, y=Dataset, fill=Dataset) ) +
        facet_wrap( ~Name, nrow=1 ) +
        theme_ridges(center_axis_labels=TRUE) +
        geom_density_ridges2(scale=1.25, size=1, alpha=0.5) +
        geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Dataset),
                          yend=as.numeric(Dataset)+0.9),
                     data=XX, color="black", lwd=2 ) +
        geom_text( aes(y=as.numeric(Dataset)+0.7, label=Label, color=Highlight),
                  x=0.95, hjust=0, data=XX, fontface="bold", size=4 ) +
        coord_cartesian(clip="off") + xlim(0.6, 1.05) +
        scale_fill_manual( values=dsPal(), guide=FALSE ) +
        scale_color_manual( values=c("yes"="red","no"="black"), guide=FALSE ) +
        theme( strip.text.x = element_text(margin=margin(b=5,t=5), face="bold"),
              strip.background = element_blank() )
##        ggsave( str_c("Fig2B-",Sys.Date(),".pdf"), width=8, height=5 ) +
##        ggsave( str_c("Fig2B-",Sys.Date(),".png"), width=8, height=5 )
}

Fig2 <- function()
{
    ## Plot individual panels
    pA <- pdfGrob("syn20540118")
    pB <- panelB()

    ## Place everything onto the same figure
    ff <- cowplot::plot_grid( pA, pB, ncol=1, rel_heights=c(0.65,1),
                             labels=c("A","B"), label_size=24 )
    ggsave( str_c("Fig2-",Sys.Date(),".pdf"), ff, width=10, height=7.5 )
    ggsave( str_c("Fig2-",Sys.Date(),".png"), ff, width=10, height=7.5 )
}

