library( tidyverse )
library( DRIAD )

synapser::synLogin()
syn <- synExtra::synDownloader( "~/data/DRIAD/rebuttal", ifcollision="overwrite.local" )

## Set up the prediction task
XY <- prepareTask( "~/data/ROSMAP/rosmap.tsv.gz", "AC" )
lP <- preparePairs(XY)

## Load relevant data
## Prefilter to the common space of genes
read_csv2 <- function( ... ) { read_csv( ..., col_types=cols() ) %>% filter(Gene %in% colnames(XY)) }
dfx1 <- syn( "syn18145776" ) %>% read_csv2
dfx2 <- syn( "syn17167348" ) %>% read_csv2

## Count the number of differentially expressed genes for each drug
nn1 <- dfx1 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )
nn2 <- dfx2 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% summarize( nn = n() ) %>% arrange( desc(nn) )

## Identify drugs with a large number of differentially expressed genes
d1 <- dfx1 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% filter( n() > 300 ) %>% pull(Drug) %>% unique
d2 <- dfx2 %>% group_by(Drug) %>% filter( FDR < 0.05 ) %>% filter( n() > 300 ) %>% pull(Drug) %>% unique

## Selects the top k differentially expressed genes for drugs d
topDFX <- function( .dfx, d, k )
    .dfx %>% filter(Drug %in% d) %>% group_by(Drug) %>% arrange(FDR) %>% slice(1:k) %>%
        summarize_at("Gene", list) %>% with( set_names(Gene, Drug) )

## Select increasingly large sets of differentially expressed genes for these drugs
v <- seq(100, 900, by=100)
DFX1 <- set_names( v, str_c("L",v) ) %>% map( ~topDFX(dfx1, d1, .x) )
DFX2 <- set_names( v, str_c("L",v) ) %>% map( ~topDFX(dfx2, d2, .x) )
DFX <- map2( DFX1, DFX2, c )

## Evalute all gene sets
##res <- imap( DFX, ~{cat("Evaluating ",.y,"\n"); evalGeneSets(.x, XY, lP)} )
##saveRDS( res,"driad-geq300.rds" )
res <- readRDS("driad-geq300.rds" )

RR <- bind_rows( res, .id="nn" ) %>% mutate( nFeats = map_int(Feats, length) ) %>%
    select( nn, nFeats, Drug=Set, AUC ) %>% filter( nFeats %% 100 == 0 ) %>%
    arrange(Drug) %>% group_by( nn, nFeats ) %>% summarize_at( "AUC", list )

CM <- expand_grid(
    rename_all( RR, str_c, "x" ),
    rename_all( RR, str_c, "y" )) %>%
    mutate(cr = map2_dbl(AUCx, AUCy, cor, method="sp"),
           lbl = ifelse( nFeatsx == 300 | nFeatsy == 300,
                        as.character(round(cr,2)), "" ) )

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))

## Panel A
pA <- ggplot( bind_rows(nn1, nn2), aes(x=nn) ) + theme_bw() +
    geom_histogram(binwidth=100, fill="steelblue") +
    geom_vline(xintercept=300) +
    scale_x_continuous( breaks=c(0, 300, 1000, 2000, 3000) ) +
    xlab( "Number of diff. expressed genes with FDR < 0.05" ) +
    ylab( "Number of drugs" )

## Panel B
pB <- ggplot( CM, aes(x=as.character(nFeatsx),
                      y=as.character(nFeatsy), fill=cr) ) +
    geom_tile() + theme_bw() +
    geom_text( aes(label=lbl) ) +
    scale_fill_gradientn(colors=pal, limits=c(-1,1),
                         name="Spearman\ncorrelation") +
    xlab( "Number of top diff. expressed genes" ) +
    ylab( "Number of top diff. expressed genes" )

## Generate a composite figure
ff <- cowplot::plot_grid(pA, NULL, pB, nrow=1, rel_widths=c(1,0.1,1.25),
                         labels=c("a","b",""), label_size=24 )
ggsave( "R4-6.pdf", ff, width=9.5, height=4 )
ggsave( "R4-6.png", ff, width=9.5, height=4 )
