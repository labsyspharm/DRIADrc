## Computes a composite score for each DGE drug
##
## by Artem Sokolov

source( "../results.R" )
source( "../plot.R" )

## Helper function that prepares a data frame for heatmap plotting
fprep <- function( .df ) {
    .df %>% select( Drug, MSBB10, MSBB22, MSBB36, MSBB44, `MSBB ave` = MSBB,
                   ROSMAP, Composite ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" ) %>% as.matrix()
}

## Prepares row annotations (toxicity) for the heatmap
aprep <- function( .df ) {
    .df %>% select( Drug, Toxicity = IsToxic, Approval ) %>%
        as.data.frame() %>% column_to_rownames( "Drug" )
}
        
## Produces a set of values to print over the heatmap
nprep <- function( .df ) {
    .n  <- fprep(.df) %>% round(2)
    storage.mode(.n) <- "character"
    .n[,1:6] <- ""
    .n
}

## Plots the top FDA-approved candidates ranked by composite score
FigS3 <- function()
{
    ## Fetch the composite score matrix and separate drugs in FDA-approved and non-approved
    XX <- DGEcomposite() %>%
        mutate( IsApproved = ifelse( Approval %in% c("approved","vet_approved"),
                                    "FDA-Approved", "Non-Approved" ) ) %>%
        mutate_at( "IsToxic", recode, `0` = "Non-Toxic", `1` = "Toxic" ) %>%
        mutate_at( "Approval", str_to_title ) %>%
        mutate( Drug = str_c(Drug, " (", Target, ") [", str_sub(Plate, 4, 4), "]"),
               Plate=NULL, Target=NULL ) %>% split( ., .$IsApproved ) %>%
        map( arrange, desc(Composite) ) %>% map( head, 15 )

    ## Set up the plotting mechanism
    pal <- RColorBrewer::brewer.pal( 9, "YlOrBr" ) %>% rev %>% colorRampPalette
    palA <- list( Toxicity = c("Toxic"="tomato", "Non-Toxic"="steelblue"),
                 Approval = set_names(ggthemes::few_pal()(3),
                                      c("Approved", "Experimental", "Investigational")) )
    fplot <- partial( pheatmap::pheatmap, cluster_rows=FALSE, cluster_cols=FALSE,
                     color=pal(100), fontsize=11, fontface="bold", gaps_col=c(4,6),
                     number_color="black", width=6.5, height=6, silent=TRUE,
                     annotation_colors=palA, legend=FALSE, annotation_legend=FALSE )
    
    ## Plot approved and non-approved drugs separately
    ggh <- map( XX, ~fplot(fprep(.x), display_numbers=nprep(.x), annotation_row=aprep(.x)) ) %>%
        map( pluck, "gtable" )

    ## Faux plot to generate the legend (Score and Toxicity)
    X1 <- XX[[1]] %>% select( Drug, Score=ROSMAP, Toxicity=IsToxic )
    smx <- max(X1$Score)
    gg1 <- ggplot( X1, aes(x=Drug, color=Score, fill=Toxicity) ) +
        geom_bar(aes(y=1), stat="identity") +
        scale_color_gradientn( colors=pal(100), limits=c(0, smx) ) +
        scale_fill_manual( values=palA$Toxicity ) +
        guides( color=guide_colorbar(title.vjust = .85) ) +
        theme( legend.title = etxt(12), legend.text = etxt(10),
              legend.position="bottom" )
    gl1 <- cowplot::get_legend( gg1 )

    ## Another faux plot to generate the Approval legend
    X2 <- bind_rows(XX) %>% select( Drug, Approval )
    gg2 <- ggplot( X2, aes(x=Drug, fill=Approval) ) +
        geom_bar( aes(y=1), stat="identity" ) +
        scale_fill_manual( values=palA$Approval ) +
        theme( legend.title=etxt(12), legend.text=etxt(10),
              legend.position="bottom", legend.margin=margin(b=0.5, unit="cm") )
    gl2 <- cowplot::get_legend( gg2 )

    ## Put everything together
#    ly <- matrix(c(1,4,2,3),2,2)
    gg <- gridExtra::arrangeGrob( grobs = c(ggh,list(gl2,gl1)), # layout_matrix=ly,
                                 heights=c(0.9,0.1), widths=c(0.96,1) )
    
    ## Write out to file
    ggsave( str_c("FigS3-", Sys.Date(), ".pdf"), gg, width=10.4, height=7 )
    ggsave( str_c("FigS3-", Sys.Date(), ".png"), gg, width=10.4, height=7 )
}
