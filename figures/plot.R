## Plotting elements that are shared across multiple figures
##
## by Artem Sokolov

suppressMessages(library( tidyverse ))

## Composes a grob from an imported .pdf
pdfGrob <- function( fnPdf )
{
    ## Fetch the .pdf source and convert it to .svg format
    fnSvg <- gsub( "pdf", "svg", fnPdf )
    grConvert::convertPicture(fnPdf, fnSvg)

    ## Import the resulting .svg and construct the grob
    grImport2::readPicture(fnSvg) %>%
        grImport2::pictureGrob()
}

## Short-hand for bold element_text of desired size
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Custom ggplot theme that boldifies text elements
theme_bold <- function()
{
    theme( axis.text = etxt(12), axis.title = etxt(14), legend.text = etxt(12),
          legend.title = etxt(14), strip.text = etxt(12),
          strip.background = element_blank() )
}

## Define the dataset palette
dsPal <- function()
{
##    c( ROSMAP = "#882256", MAYO = "#cc6677", `MSBB10` = "#89cced",
##      `MSBB44` = "#332f85", `MSBB36` = "#13783d", `MSBB22` = "#ddcb76" )
    c( ROSMAP = "#872756", MAYO = "#CA6676", `MSBB10` = "#87CAEC",
      `MSBB44` = "#363283", `MSBB36` = "#46A898", `MSBB22` = "#DBCA75" )
}

