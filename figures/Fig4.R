library( here )

source( here("figures","results.R") )
source( here("figures","plot.R") )

panelB <- function()
{
    ## Load drug ranking combining scores across plates
    R <- DGEcomposite() %>% group_by( Drug, LINCSID, Approval ) %>%
        summarize( HMP=exp(mean(log(HMP))) ) %>% ungroup() %>%
        arrange( HMP ) %>% mutate( Rank=1:n() )
    max_rank = max(R$Rank)

    ## Load TAS information
    TAS <- TASvalues() %>% select( LINCSID=compound_id, Target=entrez_symbol, TAS=tas )

    ## Combine with TAS data and isolate JAK targets
    vJAK <- c("JAK1", "JAK2", "JAK3", "TYK2")
    X <- left_join( R, TAS, by="LINCSID" ) %>% select( -HMP, -LINCSID ) %>%
        spread( Target, TAS ) %>% select( Drug, Approval, Rank, !!!syms(vJAK) ) %>%
        mutate( combined = pmin(JAK1, JAK2, JAK3, TYK2, na.rm=TRUE) ) %>%
        arrange( Rank ) %>% gather( Target, TAS, JAK1, JAK2, JAK3, TYK2, combined ) %>%
        mutate_at( "TAS", as_factor ) %>%
        mutate_at( "Target", factor, c("combined", "TYK2", str_c("JAK",3:1)) )

    ## Compute ECDFs over the combined ranking
    ECDF <- X %>% filter( Target == "combined", TAS != 10 ) %>%
        group_by(TAS) %>% summarize( RR = list(Rank) ) %>%
        mutate(Rank = map(RR, ~unique(c(1L, .x, max_rank))),
               CumProb = map2( RR, Rank, ~ecdf(.x)(.y) )) %>%
        select( -RR ) %>% unnest( c(Rank, CumProb) ) %>%
        mutate_at( "TAS", fct_drop ) %>%
        mutate_at( "TAS", fct_recode, "1 (<100 nM)"="1",
                  "2 (100-999 nM)"="2", "3 (1-10 uM)"="3" )

    ## Approved drugs that bind JAK family
    A <- X %>% filter( Approval == "approved", TAS != 10 ) %>%
        distinct( Drug, Rank )
    
    ## Plotting elements
    cpal <- set_names( c("#b2182b", "#ef8a62", "#fddbc7"), levels(ECDF$TAS) )
    fpal <- set_names( c("#b2182b", "#ef8a62", "#fddbc7", "#d9d9d9"), c(1:3,10) )
    elbl <- element_blank()
    crd <- function() coord_cartesian( xlim=c(1, 80), expand=FALSE, clip="off" )

    ## ECDF plot
    gg1 <- ggplot( ECDF, aes(x=Rank, y=CumProb) ) +
        theme_bw() + theme_bold() + crd() +
        geom_step( aes(color=TAS), size=1.1 ) +
        xlab("Drug Rank") + ylab("Cum. frac. of drugs") +
        scale_color_manual( values=cpal ) +
        scale_x_continuous(
            breaks = c(1,20,40,60,80),
            minor_breaks = c(10,30,50,70) ) +
        theme(  panel.border=elbl )
#            legend.position="bottom", legend.direction="horizontal",

    ## TAS plot
    gg2 <- ggplot( X, aes(x=Rank, y=Target, fill=TAS) ) +
        theme_bw() + theme_bold() + crd() +
        geom_tile() + labs(x=NULL) +
        scale_fill_manual( values=fpal, guide=FALSE ) +
        theme( axis.line.x = elbl, axis.ticks = elbl, axis.text.x = elbl,
              plot.margin = unit(c(0, 7.5, 0, 5.5), "pt"), panel.border=elbl )

    ## Highlight approved drugs
    gg3 <- ggplot( A, aes(x=Rank, y=0.01) ) + theme_void() + crd() +
        geom_col( aes(y=0.01), fill="black" ) +
        ggrepel::geom_text_repel( aes(label=Drug), fontface="bold",
                                 nudge_y=1, direction="x",
                                 angle=90, size = 11 / .pt ) +
        scale_y_continuous( lim=c(0,.1), expand=c(0,0) )

    egg::ggarrange( plots=list(gg3, gg2, gg1), ncol=1, heights = c(1.5,1,3.5) )
}

panelC <- function()
{
    ## All ECDF data for p < 0.05, sorted by p value
    ECDF <- read_csv( here("results","TAS-ecdf.csv"), col_types=cols() ) %>%
        filter( p.value < 0.05 ) %>% arrange( p.value ) %>%
        mutate_at( c("Target", "TAS"), as_factor )

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
        facet_wrap( ~Target, nrow=2 ) +
        geom_step( aes(color=TAS, linetype=TAS), size=1.1 ) +
        geom_text( aes(label=p.value), data=ECDFp, fontface="bold",
                  x=max_rank, y=0, vjust="inward", hjust="inward" ) +
        xlab( "Drug Rank" ) + ylab( "Cumulative Probability" ) +
        scale_color_manual( values = cpal, guide=FALSE ) +
        scale_linetype_manual( values = ltpal, guide=FALSE ) +
        theme( panel.grid.major=element_blank(),
              panel.grid.minor=element_blank() )

    gg
##    lemon::reposition_legend( gg, "center", panel= "panel-5-3" )
}

pA <- pdfGrob(here("schematics","Fig4A.pdf"))
pB <- panelB()
pC <- panelC()

ff <- cowplot::plot_grid( NULL, pA, NULL, pB, NULL, pC, ncol=1,
                         labels=c("a","","b","","c",""),
                         rel_heights=c(0.001,1, 0.05,1.1, 0.1,0.75),
                         label_size=24 )

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Fig4-", Sys.Date(), ".pdf") )
ggsave( fnOut, ff, width=9, height=13 )
