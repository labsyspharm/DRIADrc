## Figure 1
##
## by Artem Sokolov

source( "results.R" )
source( "plot.R" )

panelC <- function()
{
    library( ggridges )
    
    ## Load all results associated with previous AMP-AD gene sets
    load( syn("syn20834625") )

    ## Changes set names from PMID-S# to PMID (S#)
    f <- function(id) { str_split(id, "-") %>% map_chr(~str_c(.x[1], " (",.x[2],")")) }
    
    ## Look in the original set definitions for keywords
    V <- syn( "syn19032178" ) %>% read_lines() %>% str_split( "\t" ) %>%
        map( ~set_names(.x[2], .x[1]) ) %>% unlist %>% enframe( "Set", "Description" )

    ## Compose a joint results frame
    R <- litRes %>% select(-Feats) %>% arrange( pval ) %>%
        group_by(Set) %>% slice(1) %>% ungroup() %>%
        inner_join(V, by="Set") %>% mutate_at("Set", f) %>%
        mutate_at( "Description", recode, "Tau Neuropathology" = "PSP Tau Neuropathology" ) %>%
        mutate( Name = glue::glue("{Set}\n{Description}"),
               Label = glue::glue("p={pval}") ) %>%
        mutate_at( "Label", recode, `p=0` = "p<0.01" ) %>%
        select( -Set, -Description ) %>%
        arrange( AUC ) %>% mutate_at( "Name", as_factor )

    ## Unroll the background information
    BK <- R %>% select( Name, Dataset, AUC=BK ) %>% unnest(AUC)

    ## Generate the ridge plots
    ggplot( BK, aes(x=AUC, y=Name, fill=Dataset) ) +
        theme_ridges(center_axis_labels=TRUE) +
        geom_density_ridges2(scale=1.25, size=1) +
        geom_segment( aes(x=AUC, xend=AUC, y=as.numeric(Name), yend=as.numeric(Name)+0.9),
                     data=R, color="darkgray", lwd=2 ) +
        geom_text( aes(y=as.numeric(Name)+0.75, label=Label, hjust=0),
                  x=0.85, hjust=0, data=R, size=4 ) +
        scale_y_discrete( name=NULL ) + coord_cartesian(clip="off") +
        scale_fill_manual( values=dsPal(), guide=FALSE ) 
    ##        theme( axis.text=etxt(12), axis.title=etxt(14) )
}

## Identify individual panels
pA <- pdfGrob("syn21212910")
pB <- pdfGrob("syn20506949")
pC <- panelC()

## Put everything together
fAB <- cowplot::plot_grid( NULL, pA, NULL, pB, ncol=2, labels=c("a","","b",""),
                          rel_widths=c(0.02,1), rel_heights=c(1,0.8), label_size=24 )

ff <- cowplot::plot_grid( fAB, NULL, pC, nrow=1, labels=c("","c",""),
                         rel_widths=c(1.5,0.02,1), label_size=24 )

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Fig1-", Sys.Date(), ".pdf") )
ggsave( fnOut, ff, width=14, height=7 )
