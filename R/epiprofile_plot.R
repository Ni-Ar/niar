#' Plot the histone PTMs fold changes as an MA plot.
#'
#' @param data A tibble processed with `get_hPTM_FC()`.
#' @param sample_name Character. Name of the sample to plot.
#' @param RT_min_sd_thrshld Numeric. A minimum threshold for the retention time standard deviation across all samples in the tibble for an individual PTM. Defatult `2` min.
#' @param fc_thrshld Numeric. Label histone PTMs in the plot with absolute log2 fold change higher than `fc_thrshld`. Default `2`.
#' @param area_max_thrshld Numeric. Label histone PTMs in the plot with with average area across samples higher higher than `area_max_thrshld`. Default `2e12`.
#' @param label_txt_size Numeric. Size of the labels in the plot.
#' @param Pseudocount_value Numeric. Pseudo count value used in `get_hPTM_FC()`, used only for plotting info, not actually added. Default `1`.
#'
#' @return A nice ggplot. 
#' @import ggplot2
#' @import ggrepel
#' @export
#' @details This plotting function is inspired by bulk mRNA-seq MA plots as in DESeq2 analysis. Histone PTMs are coloured following Karolin Luger's histone proteins palette.
#'
#' @examples
#'read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
#'     tidy_hPTMs() |>
#'     get_hPTM_FC(Reference_Sample_Name = 'WT_1_DIA') |>
#'     plotMA_hPTM(data = subset(RT_sd <= 2 & Ratio > 0.000004 & Area != 1),
#'                 sample_name = "KO_2_DIA")
plotMA_hPTM <- function(data, sample_name, RT_min_sd_thrshld = 2, fc_thrshld = 2, 
                        area_max_thrshld = 2e12, label_txt_size = 2, 
                        Pseudocount_value = 1) {
  # subset for only a
  data <- subset(data, Sample == sample_name)
  
  # Define PTMs to plot in MA
  data$label <- F
  data[data$Area_log2FC <= -fc_thrshld, ]$label <- T
  data[data$Area_log2FC >= fc_thrshld, ]$label <- T
  data[data$Area_mean >= area_max_thrshld, ]$label <- T
  
  # MA plot
  ggplot(data = data ) +
    aes(x = Area_mean, y = Area_log2FC, fill = Histone, size = Ratio) +
    facet_wrap(~Sample) +
    geom_hline(yintercept = 0, linetype = 'solid', colour = "black", size = 0.5) +
    
    # # Label on top
    geom_label_repel(data = subset(data, label & Area_log2FC > 0),
                     aes(label = PTM,
                         colour = ifelse(Histone %in% c("H4", "H3", "H3.3"), yes = "A", no = "B") ),
                     family = "Arial", verbose = F, segment.color = 'black',
                     segment.curvature = -1e-20,
                     force = 2, show.legend = F, size = label_txt_size, nudge_y = 2,
                     label.padding = grid::unit(0.5, "mm")) +
    # # Labels at the bottom
    geom_label_repel(data = subset(data, label & Area_log2FC < 0),
                     aes(label = PTM,
                         colour = ifelse(Histone %in% c("H4", "H3", "H3.3"), yes = "A", no = "B")),
                     force = 2, show.legend = F, size = label_txt_size,
                     family = "Arial", segment.color = 'black',
                     segment.curvature = -1e-20,
                     nudge_y = -2, verbose = F,
                     label.padding = grid::unit(0.5, "mm")) +
    
    # # Few very abundant points above midline
    geom_label_repel(data = subset(data, label & Area_mean > area_max_thrshld & Area_log2FC >= 0.5),
                     aes(label = PTM,
                         colour = ifelse(Histone %in% c("H4", "H3", "H3.3"), yes = "A", no = "B")),
                     family = "Arial", verbose = F, segment.color = 'black',
                     force = 2, show.legend = F, size = label_txt_size,
                     nudge_x = 2, nudge_y = 0.75, min.segment.length = 1,
                     segment.curvature = -1e-20, hjust = 0,
                     label.padding = grid::unit(0.75, "mm") ) +
    
    # # Few very abundant points below midline
    geom_label_repel(data = subset(data, label & Area_mean > area_max_thrshld & Area_log2FC < -0.5),
                     aes(label = PTM,
                         colour = ifelse(Histone %in% c("H4", "H3", "H3.3"), yes = "A", no = "B") ),
                     family = "Arial", verbose = F, segment.color = 'black',
                     force = 2, show.legend = F, size = label_txt_size,
                     nudge_y = -0.5, nudge_x = 0.7, min.segment.length = 1,
                     segment.curvature = -1e-20, hjust = 0,
                     label.padding = grid::unit(0.75, "mm") ) +
    geom_point(shape = 21) +
    scale_x_log10(n.breaks = 10) +
    coord_cartesian(xlim = c(1e6, NA)) +
    scale_color_manual(values = c('A' = "white", 'B' = "black")) +
    scale_fill_manual(values = c('H1' = "gray", 'H14' = "gray", 'H12' = "gray", 'H1H15' = "gray",
                                 'H3' = "dodgerblue", 'H3.3' = "dodgerblue",
                                 'H4' = "forestgreen", 
                                 'H2B' = "firebrick2", 'H2B1B' = "firebrick2", 
                                 'H2A' = "goldenrod1",  'H2A1' = "goldenrod1",
                                 'H2A3' = "goldenrod1", 'H2AX' = "goldenrod1", 'H2AJ' = "goldenrod1",
                                 'H2AV' = "goldenrod1", 'H2AZ' = "goldenrod1")) +
    labs(x = "mean chromatogram peak area", 
         y = paste0("log2 Sample / WT ( + ", Pseudocount_value, " psuedocount)"),
         title = "EpiProfile2.1") +
    guides(fill = "none", size = guide_legend(nrow = 1, byrow = TRUE, 
                                              title =  "PTM proportion to total peptide")) +
    annotation_logticks(sides = "b") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.background = element_rect(colour = 'black'),
          panel.grid.major = element_line(colour = 'gray84', size = 0.2),
          axis.text = element_text(colour = 'black'),
          axis.ticks.length.y = grid::unit(2, "mm"), 
          axis.ticks.x = element_blank(),
          strip.background = element_blank())
}


#' Make an alluvial plot with the histone PTMs ratio in each sample
#'
#' @param data 
#' @param peptide_name Name of the peptide to plot
#' @param right_sample_name name of the sample on the RHS that will have some extra geom_labels
#' @param num_col number of colours to use for the plot
#'
#' @return A nice alluvial plot
#' @import ggplot2
#' @import ggrepel
#' @import ggalluvial
#' @import ggfittext
#' @import MetBrewer
#' @export
#'
#' @examples
# read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
#'     tidy_hPTMs() |>
#'     get_hPTM_FC(Reference_Sample_Name = 'WT_1_DIA') |>
#'     plotALLU_hPTM(peptide_name = "KSAPATGGVKKPHR(H3_27_40)",
#'                   right_sample_name = "Last_Sample_DIA", 
#'                   num_col = 12)
plotALLU_hPTM <- function(data, peptide_name, right_sample_name, num_col = 12) {
  ungroup(data) |>
    subset(Peptide_Type == peptide_name) %>%
    select( c(Histone, Peptide_Start, Peptide_End, Peptide_Sequence, Ratio, Sample, 
              Modification) ) -> tmp_data
  
  num_mod <- length(unique(tmp_data$Modification))
  # plt <- brewer.pal(n = num_col, name = "Accent")
  plt <- met.brewer(name = "Hiroshige", n = num_col, type = "continuous")
  rep_plt <- rep(plt, ceiling(num_mod / num_col))
  stopifnot(num_mod > rep_plt)
  min_ratio_label <- 0.015
  pep_plot_title <- unique( paste0(tmp_data$Histone, " ", tmp_data$Peptide_Start,
                                   "-", tmp_data$Peptide_End, " ",
                                   tmp_data$Peptide_Sequence) )
  
  
  ggplot(tmp_data) +
    aes(x = Sample, stratum = Modification, alluvium = Modification, 
        y = Ratio, fill = Modification) +
    geom_alluvium(alpha = 0.4, width = 0.4) +
    geom_stratum(alpha = 1, width = 0.6, size = 0.1) +
    geom_label_repel(data = subset(tmp_data, Sample == right_sample_name ),
                              # Repel Label all the PTMs with Ratio <= min_ratio_label 
                              # All other are NA, not shown
                              aes(label = ifelse(test = Ratio <= min_ratio_label, 
                                                 yes = Modification,
                                                 no = NA) ), 
                              stat = "stratum", size = 2, alpha = 1,
                              seed = 16, family = "Arial",
                              direction = "y",
                              position = position_nudge_repel(x = 1, y = 0),
                              hjust = 0.5,
                              label.padding = grid::unit(0.5, "mm"),
                              min.segment.length = 0.1,
                              box.padding = grid::unit(0.2, "mm"),
                              segment.color = 'black',
                              segment.size = grid::unit(0.1, "mm"),
    ) +
    geom_fit_text(aes(label = Modification ),
                             stat = "stratum", width = 0.6, min.size = 2,
                             size = 7, family = "Arial") +
    scale_fill_manual(values = rep_plt) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.2), add = 0)) +
    scale_y_continuous(expand = expansion(mult = 0, add = 0) ) +
    labs(title = pep_plot_title) +
    theme_classic() +
    theme(plot.title = element_text(size = 8, vjust = -1, hjust = 0),
          plot.background = element_blank(),
          panel.background = element_blank(),
          text = element_text(family = "Arial"),
          axis.title = element_blank(),
          axis.text = element_text(colour = "black", size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(2, "mm"),
          legend.position = "none",
          legend.box.margin = margin(l = -15),
          legend.text = element_text(size = 3),
          legend.title = element_text(size = 3.5),
          legend.box.background = element_blank(),
          legend.key.height =  unit(2.5, "mm"),
          legend.key.width =  unit(3, "mm"),
          legend.background = element_blank(),
          legend.text.align = 0) -> p_temp
  return(p_temp)
} 
