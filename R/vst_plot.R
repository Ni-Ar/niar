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
#' read_vst_tbl(path = 'path/to/compare/inclusion/tbl', show_col_types = FALSE) |>
#'    tidy_vst_psi(verbose = FALSE) |>
#'    hist_dPSI(by_complex = TRUE)
plot_hist_dPSI <- function(data, by_complex = FALSE) {
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
#' read_vst_tbl(path = "file/inclusion/table.tab", show_col_types = FALSE) |>
#'      tidy_vst_psi(verbose = F) |>
#'      subset(abs(dPSI) >= 80 ) |>
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

#' Plot gene expression vs alternatively spliced event PSI
#'
#' @param data A `data.frame` with these column names: "Quality_Score_Value", "Gene_Expr", "PSI", "GENE", and "Sample". Header is case-sensitive.
#' @param quality_thrshld vast-tools event quantification quality score threshold. Must be one of "N", "VLOW", "LOW", "OK", "SOK". For more info read the official documentation [here](https://github.com/vastgroup/vast-tools#combine-output-format) under "Column 8, score 1".
#' @param external_gene_name An ensembl-gene-name. 
#' @param vst_id vast-tools AS event ID.
#' @param unit Was the vast-tools gene expression quantified in `cRPKMs` or `TPMs`?
#' @param text Logical. Do you want to label the point in the plot. Uses `ggrepel`. Default `FALSE`.
#' @param beautify Remove pre-fixed or post-fixes from sample names. Kinda experimental.
#' @param xzero Logical. Should x-axis start from zero?
#' @param colour Either 'score' or 'PSI' to indicate if the points should be coloured by the AS event score (see `quality_thrshld`) or the PSI level. 
#' @param save_plot Logical. Do you wanna save the plot to pdf? Uses `Cairo` as device.
#' @param out_plot_dir Path specifing location where to save the plot pdf.
#' @param verbose Lofical, print info on correlation.
#' @param return_data Logical. Do not plot the data but just return the data used to plot. Default `FALSE`.
#' @param subttl Character string in case you wanna add some info to the subtitle.
#'
#' @return A ggplot2 plot or a `data.frame`.
#' @import ggplot
#' @import ggrepel
#' @import Cairo
#' @import dplyr
#' @import XICOR
#' @export
#'
#' @examples
#' # 1. Get PSI values
#' grep_psi(inclusion_tbl = psi_tbl_path, vst_id = "HsaEX0000001") |>
#'        tidy_vst_psi() -> psi_tbl
#'        
#' # 2. Get gene expression values
#' grep_gene_expression(vst_expression_tbl = expression_tbl_path,
#'                      ensembl_gene_id = gene_name ) |>
#'     tidy_vst_expr() |>
#' # 3. Join gene expression table to PSI table
#'     left_join(psi_tbl, by = "Sample") |> 
#' # 4. Calculate correlations and plot
#'     plot_corr_gene_expr_psi(external_gene_name = gene_name, unit = "TPM",
#'                             vastid = "HsaEX0000001", text = TRUE,
#'                             beautify = FALSE, xzero = FALSE, verbose = FALSE,
#'                             subttl = "A great subtitle", 
#'                             out_plot_dir = path_out, save_plot = TRUE)  
plot_corr_gene_expr_psi <- function(data, quality_thrshld = "N", 
                                    external_gene_name, vst_id, 
                                    unit, text = FALSE,
                                    beautify = FALSE, xzero = TRUE, 
                                    colour = c('score', 'PSI'), 
                                    save_plot = FALSE, out_plot_dir,
                                    verbose = TRUE, return_data = FALSE,
                                    subttl = NULL) {
  # --- CHECK PARAMS ----
  if ( !any(quality_thrshld == c("N", "VLOW", "LOW", "OK", "SOK")) ) {
    stop("Can't understand the filtering option:\t", quality_thrshld,
         "\nThe parameter quality_thrshld must be one of N, VLOW, LOW, OK, or, SOK")
  }
  
  if ( !any( c('TPM', 'cRPKM') %in% unit) ) {
    stop("unit must be either 'TPM' or 'cRPKM'")
  }
  
  if ( !any( c('score', 'PSI') %in% colour ) ) {
    stop("colour must be either 'score' or 'PSI'")
  }
  
  # if using default specification set colour-coding to 'score'
  if ( all(colour == c('score', 'PSI') )) { colour <- 'score' }
  
  
  data_required_cols <- c("Quality_Score_Value", "Gene_Expr", "PSI", "GENE", "Sample" )
  if ( all(data_required_cols %in% colnames(data) ) ) {
    
    # # If XICOR is not yet installed, try to install it now
    # if (! "XICOR" %in% installed.packages()[,"Package"] ) { 
    #   warning("The R package 'XICOR' is not installed, I'm gonna install it now.")
    #   install.packages("XICOR")
    # }
    
  } else{
    missing_col <- which( ! data_required_cols %in% colnames(data_required_cols) )
    message("Input dataframe is missing the required columns: ", 
            paste0(data_required_cols[missing_col], collapse = " ") )
  }
  
  # Subset input data only for the one vst_id, in case data contains more than one.
  data <- subset(data, EVENT == vst_id)
  
  # --- PLOT LAYOUT ----
  theme_classic() +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.title = element_text(size = 8, vjust = 1),
          legend.key.width = unit(16, "mm"),
          legend.key.height = unit(1.8, "mm"),
          legend.margin = margin(t = -2, b = -3.5, unit = "pt"),
          
          plot.background = element_blank(),
          plot.caption = element_text(size = 5, margin = margin(t = -8, unit = "pt")),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          plot.tag = element_text(size = 5.75, face = "bold", color = "black"),
          
          panel.background = element_blank(),
          panel.grid.major.y = element_line(colour = 'grey73', size = 0.5),
          
          axis.text = element_text(colour = "black", size = 8),
          axis.title = element_text(colour = "black", size = 8),
          axis.ticks.length.y = unit(2, units = "mm"),
          
          strip.background = element_blank() ) -> niar_theme
  
  # --- KEEP SAMPLES WITH AS QUALITY >= quality_thrshld ---
  # Set the Score 1 quality values as factors 
  Quality_Score_Values <- c("N", "VLOW", "LOW", "OK", "SOK")
  data$Quality_Score_Value <- factor(data$Quality_Score_Value,
                                     Quality_Score_Values)
  
  num_quality_thrshld <- as.numeric(factor(quality_thrshld, levels = Quality_Score_Values))
  data <- subset(data, as.numeric(Quality_Score_Value) >= num_quality_thrshld)
  
  # --- CORRELATION EXPRESSION vs PSI  --- 
  # Spearmann
  expr_psi_sprmn <- cor(data$Gene_Expr, data$PSI, use = "complete.obs",
                        method = "spearman") 
  expr_psi_sprmn <- round(expr_psi_sprmn, 3)
  
  # Pearson
  expr_psi_prsn <- cor(data$Gene_Expr, data$PSI, use = "complete.obs",
                       method = "pearson")
  expr_psi_prsn <- round(expr_psi_prsn, 3)
  
  # Chatterjee ( http://arxiv.org/abs/1909.10140 )
  # Remember this is not symmetric 
  expr_psi_chttrj <- calculateXI(xvec = data$Gene_Expr, yvec = data$PSI, 
                                 simple = T, seed = 16)
  expr_psi_chttrj <- round(expr_psi_chttrj, 3)
  
  if ( verbose ) {
    message("Expression vs PSI correlations:\n",
            "GENE: ", external_gene_name, "\tID: ", vst_id, "\n",
            "Spearman: ", expr_psi_sprmn, " Pearson: ", expr_psi_prsn,
            " Chatterjee: ", expr_psi_chttrj)
  }
  
  if ( beautify ) {
    suppressPackageStartupMessages( require(Biostrings) )
    names <- sort(unique(data$Sample))
    random_names <- base::sample(x = names, size = 2)
    # Try to remove the longest common prefix from the sample names
    # This looks 2 random elements in the names and checks what prefix they have 
    # in common and removes if from the plotted sample names.
    # PREFIX
    nchar_common_prefix <- Biostrings::lcprefix(s1 = random_names[1],
                                                s2 = random_names[2])
    
    if ( nchar_common_prefix > 0 ) {
      common_prefix <- substr(random_names[1], start = 1, stop = nchar_common_prefix)
      if (verbose) { message("Found prefix: ", common_prefix) }
      data <- data |>
        dplyr::mutate(Pretty_Sample = gsub(pattern = paste0("^", common_prefix),
                                           replacement = "", x = Sample,
                                           ignore.case = T, perl = F))
    } else {
      if (verbose) { message("Couldn't find common prefix.") }
      data <- data |> dplyr::mutate(Pretty_Sample =  Sample)
    }
    
    # Try to remove the longest AND most abundant common suffix from 
    # sample names with the a somewhat similar approach as for the prefixes
    # SUFFIX
    common_suffix <- longest_most_abundant_common_suffix(x = names, k = 80)
    
    if ( length(nchar(common_suffix)) > 0 ) {
      
      if (verbose) { message("Found suffix: ", common_suffix) }
      data <- data |>
        dplyr::mutate(Pretty_Sample = gsub(pattern = paste0(common_suffix, "$"),
                                           replacement = "", x = Sample,
                                           ignore.case = T, perl = F) )
    } else {
      if (verbose) { message("Couldn't find common prefix.") }
      data <- data |> dplyr::mutate(Pretty_Sample =  Sample)
    }
    # Maybe in the future this could be improved...?
  }
  
  # --- PLOT --- 
  AS_EVENT_GENE <- unique(data$GENE)
  
  info_ttl <- paste0("GENE: ", external_gene_name, "  ~  ID: ", vst_id, " (",
                     AS_EVENT_GENE, ")", "      ", "Spearman: ", expr_psi_sprmn,
                     ", Pearson: ", expr_psi_prsn, 
                     ", Chatterjee: ", expr_psi_chttrj)
  
  ggplot(data) +
    aes(x = Gene_Expr, y = PSI) +
    stat_smooth(method = 'lm', formula = 'y ~ x', se = T, level = 0.95,
                colour = 'black', size = 0.5) +
    coord_cartesian(ylim = c(0, 100), clip = 'off', default = TRUE)  +
    scale_x_continous(n.breaks = 10) +
    labs(title = info_ttl, x = paste0(external_gene_name, " ", unit),
         subtitle = subttl)  +
    ylab( bquote(.(AS_EVENT_GENE)~.(vst_id)~~Psi) ) +
    niar_theme -> cor_plot
  
  # Decide if gene expression should start from zero or not
  if( xzero ) { 
    cor_plot <- cor_plot + coord_cartesian(xlim = c(0, NA), 
                                           ylim = c(0, 100), clip = 'off' )
  }
  
  # Decide if showing sample names as text
  if ( text ) {
    suppressPackageStartupMessages( require( ggrepel ) )
    # Try to plot at best half of the sample names
    n_samples <- length(unique(data$Sample))
    k <- 2
    max_n_text <- ( floor(n_samples/k) + k)
    
    if ( beautify ) { # Plot pretty names
      cor_plot <- cor_plot + 
        geom_text_repel(aes(label = Pretty_Sample), size = 2.25, 
                        max.overlaps = max_n_text, min.segment.length = 0.25,
                        segment.size = 0.25, segment.alpha = 0.75, max.time = 2,
                        seed = 16, na.rm = T, show.legend = F, verbose = verbose)
    } else { # or plot regular names
      cor_plot <- cor_plot + 
        geom_text_repel(aes(label = Sample), size = 2.25, 
                        max.overlaps = max_n_text, min.segment.length = 0.25,
                        segment.size = 0.25, segment.alpha = 0.75, max.time = 2,
                        seed = 16, na.rm = T, show.legend = F, verbose = verbose)
    }
  }
  
  if ( colour == 'score' ) {
    # Palette of greens
    quality_score_colors <- c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837')
    names(quality_score_colors) <- Quality_Score_Values
    
    cor_plot <-  cor_plot +
      geom_point(aes(fill = Quality_Score_Value), shape = 21, size = 3) +
      scale_fill_manual(values = quality_score_colors) 
    
  } else if ( colour == 'PSI' ) {
    
    cor_plot <- cor_plot +
      geom_point(aes(fill = PSI), shape = 21, size = 3) +
      scale_fill_gradient2(low = "#E317BF", mid = "#E3C22D", high = "#11E3DA",
                           midpoint = 50, limits = c(0, 100) ) + 
      guides(fill = guide_colorbar(barwidth = grid::unit(x = 11, "cm") ) ) 
    
  } else {
    stop("Colour option is not defined correctly! colour = must be either 'score' or 'PSI'.")
  }
  
  # --- RETURN OR EXPORT PLOT AS PDF
  if (save_plot) {
    # Plot name contains filtering decisions
    if ( quality_thrshld == "N") {
      filtering_name <- "unfltrd"
    } else {
      filtering_name <- paste0("fltrd", quality_thrshld)
    }
    
    sub_ttl_name <- gsub(pattern = " ", replacement = "_", subttl)
    
    snazzy_plt_name <- paste(external_gene_name, unit, 
                             vst_id, "PSI", filtering_name, sub_ttl_name,
                             "expr_psi_correlations.pdf", sep = "_")
    
    # If output plot dir doesn't exists create it.
    if (!dir.exists(out_plot_dir)) { dir.create(out_plot_dir, recursive = T) }
    # Warn if file already exists
    if ( file.exists(file.path(out_plot_dir, snazzy_plt_name) ) ) {
      warning("The output plot already exists, I'm gonna overwrite it!")
    }
    ggsave(filename = snazzy_plt_name, 
           plot = cor_plot,
           path = out_plot_dir,
           units = "cm",
           width = 18, height = 10,
           device = cairo_pdf)
    
    if (verbose) { message("Plot: ", snazzy_plt_name, " saved in:\n", out_plot_dir) }
    
  }
  
  # Return the tidy data or the plot. Can return data independently from save_plot 
  if (return_data) {  return(data) } else { return(cor_plot) }
}




#' Nice plot to check expression and PSI across ENCODE mouse developmental stages of different tissues.
#'
#' @param data_tbl A tibble generated with `get_mouse_tissue_devel_tbl`.
#' @param title A title for the plot. By default (`NULL`) is "gene ~ vast-id".
#' @param legend Either `inside` or `outside` for specifying the position of the legend relative to the plot.
#' @param save_plot Logical, whether or not to save the plot as pdf.
#' @param plot_name Name of the pdf plot file. By default (`NULL`) is a meaningful name.
#' @param out_plot_dir Path specifying where to save the plot. By default (`NULL`) is dated subfolder.
#' @param width Width of the plot in cm (numeric). Defaults are already very good.
#' @param height Height of the plot in cm (numeric). Defaults are already very good.
#'
#' @return Either a plot of a pdf.
#' @import ggplot2 
#' @import scales
#' @import Cairo
#' @import MetBrewer
#' @export
#'
#' @examples
#' # Save plot of Pax6 expression and its exon 6 PSI 
#' get_mouse_tissue_devel_tbl(vst_id = "MmuEX0033804", ensembl_gene_id = "ENSMUSG00000027168")  |>
#'          plot_mouse_tissue_devel(legend = "outside", save_plot = T)
plot_mouse_tissue_devel <- function(data_tbl, title = NULL, legend = c('inside', 'outside'), 
                                    save_plot = FALSE, plot_name = NULL, 
                                    out_plot_dir = NULL, width = 7, height = 7) {

    # Check params 
    if ( missing(data_tbl) ) { stop("You didn't specified an table generated with `get_mouse_tissue_devel_tbl`!") } 
    
    if( is.null(title) ) {
        title <- unique(paste0(data_tbl$GENE, " ~ ", data_tbl$EVENT))
    }
    
    ggplot(data_tbl) +
        aes(x = Stage, y = Tissue, fill = mean_PSI, size = log2(mean_Gene_Expr + 1) ) +
        geom_point(shape = 21, stroke = 0.2) +
        geom_text(aes(label = round(mean_PSI, 0)), size = 1.75, family = "Arial") + 
        scale_size(range = c(2, 6), breaks = c(1:10), name = 'log2(TPMs)') +
        scale_x_discrete(expand = expansion(mult = 0.05, add = 0.1) ) +
        coord_cartesian(clip = 'off') +
        labs(title = title ) +
        theme_classic(base_family = "Arial") +
        theme(axis.text = element_text(colour = 'black', size = 5),
              axis.text.x = element_text(margin = margin(t = 0, unit = "mm")),
              axis.text.y = element_text(hjust = 1, margin = margin(r = 0, unit = "mm")),
              axis.title = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(size = 0.25, colour = "black"),
              axis.ticks.length.y = unit(x = 1, units = "mm"),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 5, vjust = 0, hjust = 0.05, margin = margin(b = -1, unit = "mm")),
              plot.background = element_blank(),
              panel.background = element_blank(),
              panel.grid.major.y = element_line(colour = "gray84", size = 0.2) ) -> core_plot
    
    
    if (legend == 'inside') {
        core_plot +
            scale_fill_gradientn(colours = met.brewer("Hiroshige", 9, direction = -1),
                                 breaks = c(0, 25, 50, 75, 100), name = "Mean PSI", 
                                 limits = c(0, 100), na.value = "gray84",
                                 labels = c("0\nSkipping", 25, 50, 75, "100\nInclusion") ) +
            guides(
                fill = guide_colourbar(
                    barwidth = unit(16, units = "mm"),
                    barheight = unit(1, units = "mm"),
                    title.position = "top",
                    title.hjust = 0,
                    title.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5,
                        vjust = 0,
                        margin = margin(t = -1, b = -0.5, unit = "mm")
                    ),
                    label.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5,
                        margin = margin(t = 0, unit = "mm")
                    )
                ),
                size = guide_legend(
                    reverse = F,
                    override.aes = list(fill = "gray84"),
                    keyheight = grid::unit(1, 'mm'),
                    keywidth = grid::unit(1, 'mm'),
                    nrow = 1,
                    title.vjust = 0.55,
                    title.position = "top",
                    title.hjust = 0,
                    title.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5,
                        vjust = 0,
                        margin = margin(t = -1, b = -1.5, unit = "mm")
                    ),
                    label.position = "bottom",
                    label.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5,
                        hjust = 0.5,
                        vjust = 0,
                        margin = margin(t = -2, unit = "mm")
                    )
                )
            )  +
    theme(legend.position = c(0.225, 0.80),
          legend.direction = "horizontal",
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.background = element_rect(colour = "white"),
          legend.margin = margin(t = -1, b = -1, l = 0, r = -1, unit = "mm"),
          legend.box.spacing = grid::unit(0, 'mm'),
          legend.key = element_blank() ) -> final_plot
        
    } else if ( legend == "outside") {
        core_plot +
            scale_fill_gradientn(colors = MetBrewer::met.brewer("Hiroshige", 9, direction = -1),
                                 breaks = c(0, 25, 50, 75, 100), name = "Mean PSI", 
                                 limits = c(0, 100), na.value = "gray84",
                                 labels = c("0 (Skipping)", 25, 50, 75, "100 (Inclusion)") ) +
            guides(
                fill = guide_colourbar(
                    barheight = unit(25, units = "mm"),
                    barwidth = unit(1, units = "mm"),
                    title.position = "top",
                    title.hjust = 0,
                    title.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5),
                    label.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5) ),
                size = guide_legend(
                    override.aes = list(fill = "gray90"),
                    keyheight = grid::unit(1, 'mm'),
                    keywidth = grid::unit(1, 'mm'),
                    title.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5 ),
                    label.theme = element_text(
                        family = "Arial",
                        colour = "black",
                        size = 5
                    )
                )
            )  +
            theme(legend.title = element_text(size = 5),
                  legend.text = element_text(size = 5),
                  legend.background = element_blank(),
                  legend.margin = margin(t = -1, b = -1, l = 2, r = -1, unit = "mm"),
                  legend.box = "vertical",
                  legend.box.spacing = grid::unit(0, 'mm'),
                  legend.box.background = element_blank(),
                  legend.key = element_blank() ) -> final_plot
        
        width <- width + 0.5
        
    } else{
        stop("'legend' must specify the position of the legend relative to the plot.")
    }

    if (is.null(plot_name) ) {
        plot_name <- paste0("Mouse_Devel_", unique(data_tbl$GENE), "_Expr_", 
                            unique(data_tbl$EVENT), "_PSI_", 
                            width, "x", height, "cm.pdf" ) 
    }
    
    if (save_plot) {
        if( is.null (out_plot_dir) ) { 
            anal_dir <- file.path('/users/mirimia/narecco/projects/07_Suz12AS/analysis')
            out_plot_dir <- file.path(anal_dir, 'tools_output/ENCODE_Mouse_Development', format(Sys.Date(), "%Y_%m_%d"))
        }
    }

    if ( save_plot == TRUE ) {
        # If output plot dir doesn't exists create it.
        if (!dir.exists(out_plot_dir)) { dir.create(out_plot_dir, recursive = T) }
        # Warn if file already exists
        if ( file.exists(file.path(out_plot_dir, plot_name) ) ) {
            warning("The output plot already exists, I'm gonna overwrite it!")
        }
        
        ggsave(filename = plot_name, plot = final_plot, device = cairo_pdf, 
               path = out_plot_dir, units = "cm", width = width, height = height)
        
    } else if ( save_plot == FALSE ) {
        return(final_plot)
        
    } else {
        stop("`save_plot` must be logical! Either TRUE or FALSE")
    }
       
}