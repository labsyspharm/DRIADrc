suppressMessages(library(tidyverse))
library(grid)
library(gridExtra)
library(cowplot)
library(here)

wd <- here("mechanism", "polypharmacology")
source(here("figures", "plot.R"))

gmean <- function(x) {exp(mean(log(x)))}

calculate_ecdf <- function( drug_ranking, ..., max_rank = NULL) {
  max_rank <- if ( is.null(max_rank) ) max( drug_ranking$Rank ) else max_rank
  as_tibble(drug_ranking) %>%
    group_by(...) %>%
    summarize( RR=list(Rank) ) %>%
    ## Always include 1 and max rank for visualization
    mutate(Rank = map(RR, ~unique(c(0L, .x, max_rank))),
           CumProb = map2( RR, Rank, ~ecdf(.x)(.y) )) %>%
    select( ..., Rank, CumProb ) %>%
    unnest( c(Rank, CumProb) ) %>%
    ungroup()
}

hmean <- function(v) {length(v)/sum(1/v)}

drug_set_venns <- function() {
  vp <- viewport(x = .5, y = .5, width = unit(1, "snpc"), height = unit(1, "snpc"), name = "venn_diagram")
  r <- .2

  c1 <- circleGrob(
    unit(.4, "npc"), unit(.5, "npc"), r,
    gp = gpar(col = "black", fill = "red"),
    name = "circle_1",
    vp = vp
  )
  c2 <- circleGrob(
    unit(.6, "npc"), unit(.5, "npc"), r,
    gp = gpar(col = "black", fill = "blue"),
    name = "circle_2",
    vp = vp
  )

  c30 <- .5*sqrt(3)

  x <- unit(c(.5, .4 + r*c30, .6, .4 + r*c30, .5, .6 - r*c30, .4, .6 - r*c30), "npc")
  y <- unit(c(0.5 + r*c30, .6, .5, .4, .5 - r*c30, .4, .5, .6), "npc")
  s <- c(0, -1, -1, -1, 0, -1, -1, -1)
  m <- xsplineGrob(
    x = x,
    y = y,
    shape = s,
    open = FALSE,
    gp = gpar(col = "black", fill = "white"),
    name = "intersection",
    vp = vp
  )
  combined <- gTree(children = gList(c1, c2, m), name = "venn")

  pal <- ggthemes::few_pal()(3)
  
  venn_configs <- list(
    "s1_and_s2" = c(circle_1 = "white", circle_2 = "white", intersection = pal[1]),
    "s1_xor_s2" = c(circle_1 = "red", circle_2 = "red", intersection = "white"),
    "s1_not_s2" = c(circle_1 = pal[2], circle_2 = "white", intersection = "white"),
    "s2_not_s1" = c(circle_1 = "white", circle_2 = pal[3], intersection = "white")
  )

  venns <- imap(
    venn_configs,
    function(config, name) {
      edits <- imap(
        config,
        ~gEdit(gPath(.y), gp = gpar(fill = .x))
      ) %>%
        unname() %>%
        {do.call(gEditList, .)}
      nv <- applyEdits(combined, edits) %>%
        editGrob(name = name)
    }
  )
  venns
}

plot_single_drug_combo_density <- function(target_combos, targets) {
    ## Mapping T1/T2 to actual target names
    tgtmap <- c(T1_AND_T2 = paste0(targets[1], " AND\n", targets[2]),
                T1_AND_NOT_T2 = paste0(targets[1], " AND\nNOT ", targets[2]),
                T2_AND_NOT_T1 = paste0(targets[2], " AND\nNOT ", targets[1]))

    ## Pull out the relevant slice of results data frame
    drug_sets <- target_combos %>%
        filter(Target_1 == targets[[1]], Target_2 == targets[[2]]) %>%
        chuck("data", 1) %>%
        bind_rows( .id="Drug_Set" ) %>%
        filter( Drug_Set != "T1_XOR_T2" ) %>%
        mutate_at( "Drug_Set", recode, !!!tgtmap ) %>%
        mutate_at( "Drug_Set", factor, levels=tgtmap )

    ## Compose text labels for each facet
    TXT <- drug_sets %>% select( Drug_Set ) %>%
        distinct() %>% mutate( Lbl = as.character(Drug_Set) ) %>%
        mutate_at( "Lbl", map_chr, gsub, pattern=" AND ", replacement=" AND\n" )
    
    ggplot(drug_sets, aes(Rank)) + theme_bw() +
        geom_density(aes(y = stat(scaled), color = Drug_Set, fill = Drug_Set), alpha = .8) +
        geom_tile(aes(y = -0.05, x = Rank, fill = Drug_Set), width = 0.80, height = 0.1) +
        scale_y_continuous(breaks = NULL, minor_breaks = NULL) +
        scale_x_continuous(breaks = c(1, 20, 40, 60, 77)) +
        labs(x = "Drug rank", y = "Density estimate") +
        coord_cartesian(ylim = c(-.1, 1.05), xlim = c(0, 77), expand = FALSE) +
        facet_wrap(~Drug_Set) +
        ggthemes::scale_color_few(guide = FALSE) +
        ggthemes::scale_fill_few(name = "Drug Set", labels=c("A AND B", "A AND NOT B", "B AND NOT A")) +
        theme( strip.background=element_blank(), panel.grid=element_blank(),
              legend.key.size=unit(2,"lines") )
}

panelA <- function() {
    gg <- plot_single_drug_combo_density(target_combos, c("RPS6KA1", "TYK2"))
    gt <- ggplotGrob(gg)
    hj <- c(0, 0.5, 0.5, 0.7, 0.7)
    for( j in grep("axis-b", gt$layout$name) )
        pluck( gt, "grobs", j, "children", 2, "grobs", 2, "children", 1, "hjust" ) <- hj

    
    ## Add venn diagrams
    vns <- drug_set_venns()[c(1L, 3L, 4L)]
    gtl <- gtable::gtable_filter( gt$grobs[[24]]$grobs[[1]], "key", invert=TRUE ) %>%
        gtable::gtable_add_grob( vns[[1]], 4, 2 ) %>%
        gtable::gtable_add_grob( vns[[2]], 5, 2 ) %>%
        gtable::gtable_add_grob( vns[[3]], 6, 2 )
    gt$grobs[[24]]$grobs[[1]] <- gtl

    gt
}

panelB <- function( target_combo_significance_aggregated, target_combos,
                   ns = c(synergistic = 10, antagonistic = 10, neutral = 5) ) {
  plot_data <- target_combo_significance_aggregated %>%
    mutate(
      Class = factor(Class, levels = c("synergistic", "NS", "antagonistic"))
    ) %>%
    arrange(Class, Combined) %>%
    group_by(Class) %>%
    group_modify(
      ~switch(
        as.character(.y$Class),
        "synergistic" = head(.x, n = ns[["synergistic"]]),
        "antagonistic" = head(.x, n = ns[["antagonistic"]]),
        "NS" = sample_n(.x, size = ns[["neutral"]])
      )
    ) %>%
    ungroup() %>%
    left_join(target_combos, by = c("Target_1", "Target_2")) %>%
    unnest_longer(data, indices_to = "Drug_Set") %>%
    filter(Drug_Set != "T1_XOR_T2") %>%
      mutate(
          Drug_Set = recode(Drug_Set, T1_AND_T2 = "A AND B",
                            T1_AND_NOT_T2 = "A AND NOT B", T2_AND_NOT_T1 = "B AND NOT A") %>%
              factor(levels = c("A AND B", "A AND NOT B", "B AND NOT A")),
          Combination = paste(Target_1, Target_2, sep = "\n") %>%
              fct_inorder() ) %>%
    unnest(data)

  # browser()

  effect_color_map <- c(
    "synergistic" = "#f6f9fc",
    "antagonistic" = "#fcf9f6",
    "neutral" = "#f9f9f9"
  )

  density_plot <- plot_data %>%
    ggplot(aes(Rank)) +
    facet_wrap(
      vars(Combination),
      strip.position = "left"
    ) +
#    geom_density(aes(y = stat(scaled), color = Drug_Set), fill = "gray90", alpha = .6) +
    geom_density(aes(y = stat(scaled), color = Drug_Set, fill = Drug_Set), alpha = .6) +
    geom_text(
      aes(label = Combined_Text),
      data = function(data) {
        data %>%
          distinct(Combination, Combined) %>%
          mutate(Combined_Text = sprintf("%.2g", Combined))
      },
      x = 76, y = 0.05, hjust = 1, vjust = 0, size = 8 / .pt
    ) +
    geom_text(
      aes(label = Class_Rank),
      data = function(data) {
        data %>%
          filter(Class != "NS") %>%
          group_by(Class) %>%
          distinct(Combination, Combined) %>%
          arrange(Combined) %>%
          mutate(Class_Rank = map_chr(1:n(), toOrdinal::toOrdinal, language = "English"))
      },
      x = 3, y = .05, hjust = 0, vjust = 0, size = 8 / .pt
    ) +
    geom_tile(
      aes(y = -.1, fill = Drug_Set),
      width = 0.80, height = 0.2
    ) +
      ggthemes::scale_color_few(guide = FALSE) +
    ggthemes::scale_fill_few(guide=FALSE) +
    theme(
      strip.text = element_text(color = "black"),
      strip.text.y = element_text(angle = 180),
      strip.background = element_rect(fill = NA),
      strip.switch.pad.wrap = unit(0, "pt"),
      panel.spacing.x = unit(0, "pt"),
      # Put more space between rows with different effects
      panel.spacing.y = unit(c(1, 4, 4, 1), "pt"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing = unit(0, "pt"),
      legend.box = "vertical",
      legend.box.margin = margin(2, 2, 2, 2, "pt"),
      legend.box.spacing = unit(0, "pt")
    ) +
    scale_y_continuous(
      breaks = NULL,
      minor_breaks = NULL,
      limits = c(-0.2, 1.1)
    ) +
    scale_x_continuous(breaks = c(1, 20, 40, 60, 77)) +
    coord_cartesian(xlim = c(1, 77), expand = FALSE) +
    labs(
      x = "Drug rank",
      y = ""
    )

  # Add shaded background
  density_plot_grob <- density_plot %>%
    ggplotGrob() %>%
    # Synergistic
    gtable::gtable_add_grob(
      grid::rectGrob(gp = grid::gpar(fill = effect_color_map[["synergistic"]], col = NA)),
      7, 5, 11, 28, z = 0.5
    ) %>%
    # Neutral
    gtable::gtable_add_grob(
      grid::rectGrob(gp = grid::gpar(fill = effect_color_map[["neutral"]], col = NA)),
      15, 5, 15, 28, z = 0.5
    ) %>%
    # Antagonistic
    gtable::gtable_add_grob(
      grid::rectGrob(gp = grid::gpar(fill = effect_color_map[["antagonistic"]], col = NA)),
      19, 5, 23, 28, z = 0.5
    ) %>%
    # Add effect labels on y-axis
    {
        labels <- list(
            grid::textGrob("synergistic", rot = 90, gp = gpar(fontsize = 12)),
            grid::textGrob("neutral", rot = 90, gp = gpar(fontsize = 12)),
            grid::textGrob("antagonistic", rot = 90, gp = gpar(fontsize = 12))
        ) %>%
            map(
                ~ggplot2:::add_margins(
                               gList(.x),
                               height = grobHeight(.x), width = grobWidth(.x),
                               margin = margin(0, 5, 0, 5, "pt"), margin_x = TRUE, margin_y = FALSE
                           )
            )
        i <- 3
        .$widths[i] <- exec(max, !!!map(labels, grobWidth))
        gtable::gtable_add_grob( ., labels[[1]],  7, i, 11, z = 0.5 ) %>%
        gtable::gtable_add_grob(    labels[[2]], 15, i, 15, z = 0.5 ) %>%
        gtable::gtable_add_grob(    labels[[3]], 19, i, 23, z = 0.5 )
    }

  density_plot_grob
}

Fig5 <- function() {
    pB <- panelB( target_combo_significance_aggregated, target_combos )
    plot_grid( plot_grid( panelA(), cowplot::ggdraw(), rel_widths=c(3,1) ),
              pB, labels = c("a", "b"), ncol = 1, rel_heights = c(1, 2) )
}

## Load significance of target combinations for each evaluated target pair
target_combo_significance <- read_rds(file.path(wd, "target_combo_significance.rds"))

## Load significance aggregated across all drug set tests for each target pair
target_combo_significance_aggregated <- read_rds(
  file.path(wd, "target_combo_significance_aggregated.rds")
)

## Load drug sets for each target pair
target_combos <- read_rds(file.path(wd, "target_combos.rds")) %>%
  mutate_at( "data", map, function(x) {
      x[["T1_XOR_T2"]] <- bind_rows(x[["T1_AND_NOT_T2"]], x[["T2_AND_NOT_T1"]])
      x[["NOT_T1_AND_NOT_T2"]] <- NULL
      x
  } )

set.seed(42)
fig5_plot <- Fig5()

## Compose the filename or extract it from the command line
cmd <- commandArgs( trailingOnly=TRUE )
fnOut <- `if`( length(cmd) > 0, cmd[1], str_c("Fig5-", Sys.Date(), ".pdf") )

ggsave2( fnOut, fig5_plot, width = 9, height = 7 )

