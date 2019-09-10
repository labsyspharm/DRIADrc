## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

## Returns a score legend grob
legendGrobScore <- function( cbr, mar = margin(b=0.5, l=0.5, unit="cm"))
{
    X <- tibble( x=c(0,1), y=c(0,1), HMP = range(cbr) )
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% colorRampPalette
    gcb <- guide_colorbar(title.vjust = .85, barwidth=unit(4,"cm"))
    gg <- ggplot( X, aes(x=x, y=y, color=HMP) ) + geom_bar(stat="identity") +
        scale_color_gradientn( colors=pal(100), guide=gcb, trans="log",
                              breaks=c(0.01,0.05,0.3) ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom", legend.margin=mar )
    cowplot::get_legend( gg )
}

## Helper function that prepares a data frame for heatmap plotting
fprep <- function( .df ) {
    .df %>% select( Drug:HMP ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" ) %>% as.matrix()
}

## Prepares row annotations (toxicity) for the heatmap
aprep <- function( .df ) {
    .df %>% select( Drug, Toxicity = IsToxic, Approval ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" )
}
        
## Produces a set of values to print over the heatmap
nprep <- function( .df ) {
    .n  <- fprep(.df) %>% round(3)
    storage.mode(.n) <- "character"
    .n[,1:5] <- ""
    .n
}

## Additional MoA information for FDA-approved compounds
FDA_MOA <- function()
{
    tribble( ~Drug, ~MoA,
            "baricitinib", "JAK family",
            "lapatinib", "EGFR, HER2",
            "regorafenib", "RET, VEGFRs, a.o.",
            "tofacitinib", "JAK family",
            "ruxolitinib", "JAK1/2",
            "dasatinib", "BCR-ABL, SRCfamily, a.o.",
            "ponatinib", "BCR-ABL, SRCfamily, a.o.",
            "bortezomib", "26S proteazome",
            "cabozantinib", "MET, VEGFRs, AXL, a.o.",
            "vorinostat", "HDAC1/2/3/6",
            "palbociclib", "CDK4/6",
            "nilotinib", "BCR-ABL",
            "ibrutinib", "BTK" )            
}

## Plots the top FDA-approved candidates ranked by composite score
Fig3 <- function()
{
    ## Fetch the composite score matrix and separate drugs in FDA-approved and non-approved
    XX <- DGEcompositePre() %>% select( -LINCSID ) %>%
        mutate( IsApproved = ifelse( Approval %in% c("approved","vet_approved"),
                                    "FDA-Approved", "Non-Approved" ) ) %>%
        mutate_at( "IsToxic", recode, `0` = "Non-Toxic", `1` = "Toxic" ) %>%
        mutate_at( "Approval", str_to_title ) %>% left_join( FDA_MOA(), by="Drug" ) %>%
        mutate( Target = ifelse(is.na(MoA), Target, MoA), MoA=NULL ) %>%
        mutate( Drug = str_c(Drug, " (", Target, ") [", str_sub(Plate, 4, 4), "]"),
               Plate=NULL, Target=NULL ) %>% split( ., .$IsApproved ) %>%
        map( arrange, HMP ) %>% map( head, 15 )

    ## Define palettes
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% colorRampPalette
    palA <- list( Toxicity = c("Toxic"="tomato", "Non-Toxic"="steelblue"),
                 Approval = set_names(ggthemes::few_pal()(3),
                                      c("Approved", "Experimental", "Investigational")) )

    ## Define the value -> color mapping
    xmn <- min( fprep(XX[[1]]) )
    xmx <- max( fprep(XX[[1]]) )
    cbr <- seq( log(xmn), log(xmx), length.out=101 ) %>% exp
        
    ## Set up the plotting mechanism
    fplot <- partial( pheatmap::pheatmap, cluster_rows=FALSE, cluster_cols=FALSE,
                     color=pal(100), breaks=cbr, fontsize=11, fontface="bold", gaps_col=5,
                     number_color="black", width=6.5, height=6, silent=TRUE,
                     annotation_colors=palA, legend=FALSE, annotation_legend=FALSE )
    
    ## Plot approved and non-approved drugs separately
    ggh <- map2( XX, c("FDA-approved\n", "Experimental and Investigational\n"),
                ~fplot(fprep(.x), display_numbers=nprep(.x), annotation_row=aprep(.x), main=.y) ) %>%
        map( pluck, "gtable" )

    ## Faux plot to generate the Score legend
    gl1 <- legendGrobScore( cbr )
    
    ## Faux plot to generate the Toxicity legend
    X2 <- XX[[1]] %>% select( Drug, Toxicity=IsToxic )
    gg2 <- ggplot( X2, aes(x=Drug, fill=Toxicity) ) + geom_bar(aes(y=1), stat="identity") +
        scale_fill_manual( values=palA$Toxicity,
                          guide=guide_legend(title.position="top") ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom", legend.margin=margin(b=0.5, l=0.5, unit="cm") )
    gl2 <- cowplot::get_legend( gg2 )
    
    ## Another faux plot to generate the Approval legend
    X3 <- bind_rows(XX) %>% select( Drug, Approval )
    gg3 <- ggplot( X3, aes(x=Drug, fill=Approval) ) +
        geom_bar( aes(y=1), stat="identity" ) +
        scale_fill_manual( values=palA$Approval,
                          guide=guide_legend(title.position="top") ) +
        theme( legend.title=etxt(12), legend.text=etxt(10),
              legend.position="bottom", legend.margin=margin(l=0, b=0.5, unit="cm"))
    gl3 <- cowplot::get_legend( gg3 )

    ## Put everything together
    gleg <- gridExtra::arrangeGrob( grobs = list(gl3,gl2,gl1), nrow=1,
                                   widths=c(4,2,3) )
    ly <- matrix( c(1,3,2,3), 2, 2 )
    gg <- gridExtra::arrangeGrob( grobs = c(ggh,list(Legend=gleg)), layout_matrix=ly,
                                 heights=c(0.9,0.1), widths=c(1.3,1) )

    ## Add custom annotations
    cowplot::ggdraw(gg) +
        cowplot::draw_text( "Drug (FDA-proposed MoA) [Plate Index]",
                          0.4, 0.93, fontface="bold", size=11 ) +
            cowplot::draw_text( "Drug (Vendor Target) [Plate Idx]",
                               0.87, 0.93, fontface="bold", size=11 )
    cowplot::ggsave( str_c("Fig3-", Sys.Date(), ".pdf"), width=10.4, height=7.1 )
}
