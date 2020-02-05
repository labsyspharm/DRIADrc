## DGE gene sets
## A - Schematic explaining DGE experiment
## B - Exemplar on two drugs
##
## by Artem Sokolov

library( here )

source( here("figures","results.R") )
source( here("figures","plot.R") )

panelB <- function()
{
    library( ggridges )

    R <- DGEcomposite() %>% select( LINCSID, Drug, Target ) %>% distinct
    
    ## Load all the results
    load( syn("syn20928450") )
    load( syn("syn20948620") )

    ## Isolate the slice corresponding to:
    ##   1. NVP-TAE684 (DGE1) - HMSL10024
    ##   2. Ruxolitinib (DGE2) - HMSL10138
    X1 <- allRes %>% filter(Task == "AC") %>% mutate(nFeats = map_int(Feats, length), Feats=NULL)
    X2 <- mayoRes %>% filter(Task == "AC")
    XX <- bind_rows( X1, X2 ) %>% rename( LINCSID=Drug ) %>% inner_join(R, by="LINCSID") %>%
        mutate_at( "Target", recode, JAK1 = "JAK1/2" ) %>%
        filter( (Drug == "ruxolitinib" & Plate=="DGE2") | Drug == "nvp-tae684" ) %>%
        select( -Task, -LINCSID, -Plate, -nFeats ) %>%
        mutate( Name = glue::glue("{Drug} ({Target})"), Drug=NULL, Target=NULL,
               Label = ifelse(pval == 0, " p<0.001", str_c(" p=",pval)),
               Highlight=ifelse(pval <= 0.05, "yes", "no"), pval=NULL ) %>%
        arrange( Dataset ) %>% mutate_at( "Dataset", as_factor )

    ## Unroll the background information
    BK <- XX %>% select( Name, Dataset, AUC=BK ) %>% unnest(AUC)

    ## Generate the ridge plots
    ggplot( BK, aes(x=AUC, y=Dataset, fill=Dataset) ) +
        facet_wrap( ~Name, nrow=1 ) +
        theme_ridges(center_axis_labels=TRUE) +
        geom_vline( xintercept=seq(0.5, 0.9, by=0.1), color="gray90" ) +
        geom_density_ridges2(scale=1.25, size=1, alpha=0.5) +
        geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Dataset),
                          yend=as.numeric(Dataset)+0.9),
                     data=XX, color="black", lwd=2 ) +
        geom_text( aes(y=as.numeric(Dataset)+0.7, label=Label, color=Highlight),
                  x=0.94, hjust=0, data=XX, fontface="bold", size=4 ) +
        coord_cartesian(clip="off") + xlim(0.5, 0.99) +
        scale_fill_manual( values=dsPal(), guide=FALSE ) +
        scale_color_manual( values=c("yes"="red","no"="black"), guide=FALSE ) +
        theme( strip.text.x = element_text(margin=margin(b=5,t=5), face="bold"),
              strip.background = element_blank(), panel.grid.major.x=element_blank() )
}

## Plot individual panels
pA <- pdfGrob(here("schematics","Fig2A.pdf"))
pB <- panelB()

## Place everything onto the same figure
ff <- cowplot::plot_grid( pA, pB, ncol=1, rel_heights=c(0.65,1),
                         labels=c("a","b"), label_size=24 )

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Fig2-", Sys.Date(), ".pdf") )
ggsave( fnOut, ff, width=10, height=7.5 )
