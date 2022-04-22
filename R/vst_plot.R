#' Plot the dPSI as an histogram using ggplot2
#'
#' @param data A data.frame that has at least one column called 'dPSI'
#' @param by_complex Logical, whether to split the histogram by AS event complexity/type. Uses vast-tools tables "COMPLEX" column.
#'
#' @return A ggplot2 plot
#' @import ggplot2
#' 
#' @export
#'
#' @examples
#' read_vst_tbl(path = 'path/to/compare/inclusion/tbl', show_col_types = FALSE) %>%
#'    tidy_vst_psi(verbose = FALSE) %>%
#'    hist_dPSI(by_complex = TRUE)
hist_dPSI <- function(data, by_complex = FALSE) {
    if ( !any(colnames(data) == 'dPSI') ) {
        stop("The input data.frame must have a column called 'dPSI'. I only see:\n",
             colnames(data) )
    } 
    col_to_select <- c("dPSI")
    if (by_complex) {
        if ( !any(colnames(data) == 'COMPLEX') ) {
            stop("The input data.frame must have a column called 'COMPLEX' ", 
                 "in order to use by_complex = TRUE." )
        }
        col_to_select <- c("dPSI", "COMPLEX")
    }
    
    ggplot(unique(data[, col_to_select])) +
        aes(x = dPSI, fill = dPSI > 0 ) +
        geom_histogram(binwidth = 1, show.legend = F) +
        scale_x_continuous(n.breaks = 8) +
        scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0, 0.01))) +
        scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "dodgerblue3")) +
        labs(x = "\u0394PSI") +
        theme_classic() +
        theme(panel.grid.major = element_line(colour = 'gray84', size = 0.1),
              axis.text = element_text(colour = 'black'), 
              strip.background = element_blank()) -> p_dPSI
    
    if (by_complex) { 
        p_dPSI <- p_dPSI + facet_wrap(~ COMPLEX, scales = "free_y", ncol = 4) 
        
    }
    return(p_dPSI)
}

#' Plot the PSI of all samples in a data.frame
#'
#' @param data A dataframe with certain columns. Do NOT use this function for plotting a lot of AS events. 
#' @param simplify_names Logical, Should names be shortened? Helps with shorter X-axis text names.
#'
#' @return A ggplot2 plot
#' @import ggplot2
#' @importFrom dplyr mutate
#' @export
#'
#' @examples
#' read_vst_tbl(path = "file/inclusion/table.tab", show_col_types = FALSE) %>%
#'      tidy_vst_psi(verbose = F) %>%
#'      subset(abs(dPSI) >= 80 ) %>%
#'      plot_samples_PSI()
plot_samples_PSI <- function(data, simplify_names = TRUE) {
    
    required_cols <- c("Sample", "PSI", "Quality_Score_Value", "GENE", "EVENT")    
    if (!all(required_cols %in% colnames(data) )) {
        stop("The input data needs the following columns:\n", 
             paste0(required_cols, collapse = ", "))
    }
    
    if (simplify_names) {
        data <- mutate(data, Sample = longest_common_prefix(names = Sample))
    }
    
    # set colour palette for AS event quality
    quality_score_colors <- c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837')
    quality_score_values <- c("N", "VLOW", "LOW", "OK", "SOK")
    names(quality_score_colors) <- quality_score_values
    
    ggplot(data) +
        aes(x = Sample, y = PSI, fill = Quality_Score_Value) +
        facet_wrap(  ~ GENE + EVENT, scales = "free_x") +
        geom_point(size = 2, shape = 21) +
        scale_fill_manual(values = quality_score_colors, name = "Score") +
        scale_y_continuous(breaks = seq(0, 100, 10), 
                           expand = expansion(mult = 0, add = 0.1)) +
        guides(fill = guide_legend(override.aes = list(shape = 21))) +
        coord_cartesian(ylim = c(0, 100), clip = 'off') +
        theme_classic() +
        theme(axis.text = element_text(colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank(),
              strip.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray84', size = 0.15),
              legend.position = "left") -> p_PSI
    
    return(p_PSI)
}