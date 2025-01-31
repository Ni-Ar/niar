#' Plot a Histogram of ΔPSI Values
#'
#' @description 
#' This function creates a histogram of ΔPSI values from a given data frame, with optional faceting by splicing type or splicing complexity. It is particularly useful for visualising the distribution of ΔPSI values in alternative splicing analysis.
#'
#' @param data A data frame containing the ΔPSI values and optionally the splicing type or complexity. Must include a column named 'dPSI'. If `by` is "complex" or "type", the data frame must also include a 'COMPLEX' column.
#' @param by A character string specifying the grouping variable for faceting. Accepts "complex" (splicing complexity defined by vast-tools), "type" (splicing type: Exon, Intron, Alt. ss), or `NULL` (no faceting). Default is `NULL`.
#' @param num_col An integer specifying the number of columns for faceted plots. Default is 1.
#' @param x_lims A numeric vector of length 2 specifying the x-axis limits for the histogram. Default is |∆PSI| 95% (`c(-95, 95)`).
#' @param ... Additional arguments passed to `facet_wrap()`, such as `strip.position` or `scales = "fixed"`.
#'
#' @details 
#' If `by = "type"`, the function automatically categorises splicing events in the 'COMPLEX' column into one of three types: Exon, Intron, or Alternative splice site (Alt. ss). Valid splicing complexity values are:
#' - Exon types: S, C1, C2, C3, ANN, MIC
#' - Intron types: IR
#' - Alternative splice sites: Alt3, Alt5
#'
#' The function also validates the input parameters and ensures that necessary columns are present in the data frame. It throws an informative error if any required input is missing or invalid.
#' The `x_lims` parameter allows users to define custom x-axis limits for the histogram, providing flexibility in visualising the ΔPSI distribution over a specific range. The default is set to `c(-95, 95)`, which is suitable for most ΔPSI datasets.
#'
#' The `...` parameter enables the user to pass additional arguments to the `facet_wrap()` function when faceting by splicing complexity or type. This can include options such as `scales = "free"` for adjusting individual facet scales or `strip.position = "bottom"` to modify the placement of facet strip labels.
#'
#' @return 
#' A `ggplot` object showing a histogram of ΔPSI values, optionally faceted by splicing complexity or type.
#'
#' @importFrom dplyr mutate case_when
#' @import ggplot2
#' @export
#'
#' @examples 
#' read_vst_tbl(path = 'path/to/compare/inclusion/tbl', show_col_types = FALSE) |>
#'    tidy_vst_psi(verbose = FALSE) |>
#'    plot_hist_dPSI(data, by = "type", num_col = 2)
plot_hist_dPSI <- function(data, by = NULL, num_col = 1, x_lims = c(-95, 95), ... ) {
    if ( !any(colnames(data) == 'dPSI') ) {
        stop("The input data.frame must have a column called 'dPSI'. I only see:\n",
             colnames(data) )
    } 
    
    # Sanity check for the 'by' parameter
    valid_by_values <- c(NULL, "complex", "type")
    if (!by %in% valid_by_values) {
      stop("Invalid value for 'by'. It must be one of: NULL, 'complex', or 'type'. ",
           "You provided: ", deparse(by))
    }
  
    
    if (by == "complex" || by == "type") {
        if ( !any(colnames(data) == 'COMPLEX') ) {
            stop("The input data.frame must have a column called 'COMPLEX' ", 
                 "in order to use by_complex = TRUE." ) }
      
        col_to_select <- c("dPSI", "COMPLEX")
        
    } else if ( is.null(by) ) {
      col_to_select <- c("dPSI")
      
    } else {
      # redundant by okay.
      stop("The parameter 'by' must be either 'complex', 'type', or NULL" )
    }
    
    fltrd_data <- unique(data[, col_to_select])
    
    if ( by == 'type'){
  
      # define the AS type from the "COMPLEX" vast-tools column
      exon_type <- c("S", "C1", "C2", "C3", "ANN", "MIC")
      intron_type <- c("IR")
      alt_ss <- c("Alt3", "Alt5")
      
      fltrd_data |>
        mutate(COMPLEX = case_when(COMPLEX %in% exon_type ~ 'Exon',
                                  COMPLEX %in% intron_type ~ 'Intron',
                                  COMPLEX %in% alt_ss ~ 'Alt. ss') ) -> fltrd_data
      
    }
    
    ggplot(fltrd_data) +
        aes(x = dPSI, fill = dPSI >= 0 ) +
        geom_histogram(binwidth = 1, show.legend = F) +
        scale_x_continuous(n.breaks = 8) +
        scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0, 0.01))) +
        scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "dodgerblue3")) +
        labs(x = "\u0394PSI") +
        coord_cartesian(xlim = x_lims) +
        theme_classic() +
        theme(panel.grid.major = element_line(colour = 'gray84', linewidth = 0.1),
              axis.text = element_text(colour = 'black'), 
              strip.background = element_blank()) -> p_dPSI
    
    if ( !is.null(by) ) {
        p_dPSI <- p_dPSI + facet_wrap(~ COMPLEX, ncol = num_col, ...) 
        
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
              panel.grid.major = element_line(colour = 'gray84', linewidth = 0.15),
              legend.position = "left") -> p_PSI
    
    return(p_PSI)
}

#' Plot a gene expression vs an alternatively spliced event PSI ( 1 vs 1)
#'
#' @description Plot gene expression on X-axis vs AS event PSI on Y-axis and calculate correlation between samples.
#'
#' @param data A `data.frame` with these column names: "Quality_Score_Value", "Gene_Expr", "PSI", "GENE", and "Sample". Header is case-sensitive.
#' @param quality_thrshld vast-tools event quantification quality score threshold. Must be one of "N", "VLOW", "LOW", "OK", "SOK". For more info read the official documentation [here](https://github.com/vastgroup/vast-tools#combine-output-format) under "Column 8, score 1".
#' @param external_gene_name An ensembl-gene-name. 
#' @param vst_id vast-tools AS event ID.
#' @param unit Was the vast-tools gene expression quantified in `cRPKMs` or `TPMs`?
#' @param log_expr Logical. Should gene expression values on the X-axis be log2 transformed before calculating the correlations and plotting the data? Default `FALSE`. This influences only the Pearson (`r`) correlation, the other 3 are unaffected. 
#' @param text Logical. Do you want to label the point in the plot. Uses `ggrepel`. Default `FALSE`.
#' @param text_size Number How big should the text be? Default is `8`.
#' @param label_size Number How big should the label text be? Default is `2.5`.
#' @param beautify Remove pre-fixed or post-fixes from sample names. Kinda experimental.
#' @param xzero Logical. Should x-axis start from zero?
#' @param colour Either 'score' or 'PSI' to indicate if the points should be coloured by the AS event score (see `quality_thrshld`) or the PSI level. 
#' @param save_plot Logical. Do you wanna save the plot to pdf? Uses `Cairo` as device.
#' @param out_plot_dir Path specifing location where to save the plot pdf.
#' @param verbose Lofical, print info on correlation.
#' @param return_data Logical. Do not plot the data but just return the data used to plot. Default `FALSE`.
#' @param subttl Character string in case you wanna add some info to the subtitle. For example one could specify the type of cell line or tissue.
#'
#' @return A ggplot2 plot or a `data.frame`.
#' @import ggplot2
#' @import ggrepel
#' @import Cairo
#' @import dplyr
#' @import XICOR
#' @importFrom Biostrings lcprefix
#' @export
#' 
#' @details The plot title reports the correlation values of the:
#' \itemize{
#' \item{ Spearman's rho }
#' \item{ Pearson's r }
#' \item{ Kendall's tau}
#' \item{ Chatterjee's xi}
#' }   
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
#' left_join(psi_tbl, by = "Sample") |> 
#' # 4. Calculate correlations and plot
#' plot_corr_gene_expr_psi(external_gene_name = gene_name, unit = "TPM",
#'                         vst_id = "HsaEX0000001", text = TRUE,
#'                         beautify = TRUE, xzero = FALSE, verbose = FALSE,
#'                         subttl = "A great subtitle", 
#'                         out_plot_dir = path_out, save_plot = TRUE)  
#'
plot_corr_gene_expr_psi <- function(data, quality_thrshld = "N", 
                                    external_gene_name, vst_id, 
                                    unit, log_expr = FALSE, text = FALSE,
                                    text_size = 10, label_size = 3,
                                    beautify = FALSE, xzero = TRUE, 
                                    colour = c('score', 'PSI'), 
                                    save_plot = FALSE, out_plot_dir,
                                    verbose = TRUE, return_data = FALSE,
                                    subttl = NULL) {
  # ---- CHECK PARAMS ----- 
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
  if ( !all(data_required_cols %in% colnames(data) ) ) {
    
    missing_col <- which( ! data_required_cols %in% colnames(data) )
    stop("Input dataframe is missing the required columns: ", 
            paste0(data_required_cols[missing_col], collapse = " ") )
  }
  
  if( !is.logical(log_expr)) {
    stop("Parameter 'log_expr' must be TRUE or FALSE not: ", log_expr)
  }

  # Reduce text size.
  if (save_plot == TRUE & text_size > 8 & label_size > 3) {
    text_size <- 5
    label_size <- 2
  }
    
  # Subset input data only for the one vst_id, in case data contains more than one.
  data <- subset(data, EVENT == vst_id)
  
  # ---- PLOT LAYOUT ----
  theme_classic(base_size = text_size, base_family = "Arial") +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.title = element_text(margin = margin(r = 0, unit = "mm")),
          legend.text = element_text(margin = margin(l = -4, r = -4, unit = 'mm')),
          legend.key = element_blank(),
          legend.margin = margin(t = -1, b = -1, unit = "mm"),
          legend.spacing.x = unit(5, "mm"),
          legend.spacing.y = unit(0, "mm"),
          
          plot.background = element_blank(),
          plot.title = element_text(colour = "black"),
          plot.subtitle = element_text(colour = "black"),
          
          panel.background = element_blank(),
          panel.grid.major.y = element_line(colour = 'grey73', linewidth = 0.5),
          
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          axis.title.x = element_text(margin = margin(t = 0, unit = 'mm')),
          axis.title.y = element_text(margin = margin(r = -0.5, unit = 'mm')),
          axis.ticks.length.y = unit(1, units = "mm"),
          
          strip.background = element_blank() ) -> niar_theme
  
  # ---- KEEP SAMPLES WITH AS QUALITY >= quality_thrshld ----
  Quality_Score_Values <- c("N", "VLOW", "LOW", "OK", "SOK")
  data$Quality_Score_Value <- factor(data$Quality_Score_Value,
                                     Quality_Score_Values)
  
  num_quality_thrshld <- as.numeric(factor(quality_thrshld, levels = Quality_Score_Values))
  data <- subset(data, as.numeric(Quality_Score_Value) >= num_quality_thrshld)
  
  # ---- LOG TRANSFORM GENE EXPRESSION ---- 
  if ( log_expr == TRUE) {
    data$Gene_Expr <- log2(data$Gene_Expr)
  }
  
  
  # ---- CORRELATION EXPRESSION vs PSI ----
  # Spearmann
  expr_psi_sprmn <- cor(data$Gene_Expr, data$PSI, use = "complete.obs",
                        method = "spearman") 
  expr_psi_sprmn <- round(expr_psi_sprmn, 3)
  
  # Pearson
  expr_psi_prsn <- cor(data$Gene_Expr, data$PSI, use = "complete.obs",
                       method = "pearson")
  expr_psi_prsn <- round(expr_psi_prsn, 3)
  
  # Kendall
  expr_psi_kndl <- cor(data$Gene_Expr, data$PSI, use = "complete.obs",
                       method = "kendall")
  expr_psi_kndl <- round(expr_psi_kndl, 3)
  
  # Chatterjee ( http://arxiv.org/abs/1909.10140 )
  # Remember this is not symmetric 
  expr_psi_chttrj <- calculateXI(xvec = data$Gene_Expr, yvec = data$PSI, 
                                 simple = T) #, seed = 16)
  expr_psi_chttrj <- round(expr_psi_chttrj, 3)
  
  if ( verbose ) {
    message("Expression vs PSI correlations:\n",
            "GENE: ", external_gene_name, "\tID: ", vst_id, "\n",
            "Spearman: ", expr_psi_sprmn, " Pearson: ", expr_psi_prsn,
            " Kendall: ", expr_psi_kndl, " Chatterjee: ", expr_psi_chttrj)
  }
  
  # ---- PREFIX & SUFFIX REMOVAL ----
  if ( beautify == TRUE) {
    # FIND PREFIX
    names <- sort(unique(data$Sample))
    random_names <- base::sample(x = names, size = 2)
    # Try to remove the longest common prefix from the sample names
    # This looks 2 random elements in the names and checks what prefix they have 
    # in common and removes if from the plotted sample names.
    nchar_common_prefix <- Biostrings::lcprefix(s1 = random_names[1],
                                                s2 = random_names[2])
    if( nchar_common_prefix > 0 ) {
      common_prefix <- substr(random_names[1], start = 1, stop = nchar_common_prefix)
    } else {
      common_prefix <- substr(random_names[1], start = 1, stop = nchar_common_prefix)
    }
    
    # FIND SUFFIX
    # Try to remove the longest AND most abundant common suffix from 
    # sample names with the a somewhat similar approach as for the prefixes
    common_suffix <- longest_most_abundant_common_suffix(x = names, k = 80)
    nchar_common_suffix <- nchar(common_suffix)
    # Catch when there are no common suffix and the character are zero
    if( identical(nchar_common_suffix, integer(0)) ) { nchar_common_suffix <- 0 }
    
    
    if ( nchar_common_prefix > 0 & nchar_common_suffix > 0 ) {
      
      if (verbose) { message("Found prefix: ", common_prefix) }
      if (verbose) { message("Found suffix: ", common_suffix) }
      
      data <- mutate(data, 
                     Pretty_Sample = gsub(pattern = paste0("^", common_prefix),
                                          replacement = "", x = Sample,
                                          ignore.case = T, perl = F) ) |>
        mutate(Pretty_Sample = gsub(pattern = paste0(common_suffix, "$"),
                                    replacement = "", x = Pretty_Sample,
                                    ignore.case = T, perl = F) )
      
    } else if ( nchar_common_prefix > 0 & nchar_common_suffix <= 0 ) {
      if (verbose) { message("Found prefix: ", common_prefix) }
      if (verbose) { message("Couldn't find common suffix") }
      
      data <- mutate(data, 
                     Pretty_Sample = gsub(pattern = paste0("^", common_prefix),
                                          replacement = "", x = Sample,
                                          ignore.case = T, perl = F) )
      
    } else if ( nchar_common_prefix <= 0 & nchar_common_suffix > 0 ) { 
      if (verbose) {  message("Couldn't find common prefix.") }
      if (verbose) { message("Found suffix: ", common_suffix) }
      
      data <- mutate(data, 
                     Pretty_Sample = gsub(pattern = paste0(common_suffix, "$"),
                                          replacement = "", x = Sample,
                                          ignore.case = T, perl = F) )
      
    } else if ( nchar_common_prefix <= 0 & nchar_common_suffix <= 0 ) {
      if (verbose) {  message("Couldn't find common prefix.") }
      if (verbose) { message("Couldn't find common suffix") }
      data <- mutate(data, Pretty_Sample = Sample)
    } else {
      warning("Coulnd't figure out the prefix and suffix of sample names")
    }
  }
  
  # ---- PLOT ----
  AS_EVENT_GENE <- unique(data$GENE)
  
  info_ttl <- paste0(external_gene_name, " ~ ", vst_id, " (",
                     AS_EVENT_GENE, ")", " ", "\u03c1 ", expr_psi_sprmn,
                     ", r ", expr_psi_prsn, ", \u03c4 ", expr_psi_kndl,
                     ", \u03be ", expr_psi_chttrj)
  
  if ( log_expr == TRUE) {
    X_axis_ttl <- paste0(external_gene_name, " ", unit, " (log2)")
    
  } else if ( log_expr == FALSE ) {
    X_axis_ttl <- paste0(external_gene_name, " ", unit)
  } else {
    warning(' Not sure whether the X-axis was log2 transformed')
    X_axis_ttl <- paste0(external_gene_name, " ", unit)
  }
  
  
  ggplot(data) +
    aes(x = Gene_Expr, y = PSI) +
    stat_smooth(method = 'lm', formula = 'y ~ x', se = T, level = 0.95,
                colour = 'black', size = 0.5) +
    coord_cartesian(ylim = c(0, 100), clip = 'off', default = TRUE)  +
    scale_x_continuous(n.breaks = 10) +
    labs(title = info_ttl, x = X_axis_ttl,
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
    # Try to plot at best half of the sample names
    n_samples <- length(unique(data$Sample))
    k <- 2
    max_n_text <- ( floor(n_samples/k) + k)
    
    if ( beautify == TRUE ) { # Plot pretty names
      cor_plot <- cor_plot + 
        geom_text_repel(aes(label = Pretty_Sample), size = label_size,
                        max.overlaps = max_n_text, min.segment.length = 0.25,
                        segment.size = 0.25, segment.alpha = 0.75, max.time = 2,
                        seed = 16, na.rm = T, show.legend = F, verbose = verbose,
                        family = "Arial")
    } else { # or plot regular names
      cor_plot <- cor_plot + 
        geom_text_repel(aes(label = Sample), size = label_size, 
                        max.overlaps = max_n_text, min.segment.length = 0.25,
                        segment.size = 0.25, segment.alpha = 0.75, max.time = 2,
                        seed = 16, na.rm = T, show.legend = F, verbose = verbose,
                        family = "Arial")
    }
  }
  
  if ( colour == 'score' ) {
    # Palette of greens
    quality_score_colors <- c('#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837')
    names(quality_score_colors) <- Quality_Score_Values
    
    cor_plot <-  cor_plot +
      geom_point(aes(fill = Quality_Score_Value), shape = 21, size = 1.25,stroke = 0.25) +
      scale_fill_manual(values = quality_score_colors, name = "PSI Quality") 
    
  } else if ( colour == 'PSI' ) {
    
    cor_plot <- cor_plot +
      geom_point(aes(fill = PSI), shape = 21, size = 1.25, stroke = 0.25) +
      scale_fill_gradient2(low = "#E317BF", mid = "#E3C22D", high = "#11E3DA",
                           midpoint = 50, limits = c(0, 100), name = "PSI" ) + 
      guides(fill = guide_colorbar(barwidth = grid::unit(x = 6, "cm"),
                                   barheight = grid::unit(x = 2, "mm") ) )
      
    
  } else {
    stop("Colour option is not defined correctly! colour = must be either 'score' or 'PSI'.")
  }
  
  # ---- RETURN OR EXPORT PLOT AS PDF ----
  if (save_plot) {
    # Plot name contains filtering decisions
    if ( quality_thrshld == "N") {
      filtering_name <- "unfltrd"
    } else {
      filtering_name <- paste0("fltrd", quality_thrshld)
    }
    
    if ( !is.null(subttl)) {
      sub_ttl_name <- gsub(pattern = " ", replacement = "_", subttl)
    } else {
      sub_ttl_name <- "samples"
    }
    
    if (log_expr == TRUE) {
      transform_name <- "log2GeneExpr"
    } else if( log_expr == FALSE) {
      transform_name <- ""
    } else {
      transform_name <- ""
      warning("param log_expr is not clear to me.")
    }
    
    snazzy_plt_name <- paste(external_gene_name, unit, vst_id, "PSI", 
                             filtering_name, transform_name, sub_ttl_name,
                             "corr.pdf", sep = "_")
    
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
           width = 10, height = 9,
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
#' @param colour_bar Colour palette to represent the PSI values. Select between `BlueRed` or `viridis`.
#' @param PSI_limits A vector of length 2, specifying the minimum and maximum PSI colour limits in the legend.
#' @param blck_wht_PSI_col_thshld The PSI text is shown in black, PSI values above this threshold will be displayed in white instead. Gives better contrast with dark colours of high PSI.
#'
#' @return A ggplot2 plot or a pdf.
#' @import ggplot2 
#' @import Cairo
#' @import MetBrewer
#' @export
#' 
#' @seealso [get_mouse_tissue_devel_tbl()]
#'
#' @examples
#' # Plot Pax6 expression and its exon 6 PSI 
#' get_mouse_tissue_devel_tbl(vst_id = "MmuEX0033804", ensembl_gene_id = "ENSMUSG00000027168")  |>
#'          plot_mouse_tissue_devel(legend = "inside", colour_bar = 'BlueRed')
plot_mouse_tissue_devel <- function(data_tbl, title = NULL, legend = c('inside', 'outside'), 
                                    save_plot = FALSE, plot_name = NULL, 
                                    out_plot_dir = NULL, width = 7, height = 7,
                                    colour_bar = c('BlueRed', 'viridis'), 
                                    PSI_limits = c(0, 100),
                                    blck_wht_PSI_col_thshld = 70 ) {

    # Check params 
    if ( missing(data_tbl) ) { 
      stop("You didn't specified an table generated with `get_mouse_tissue_devel_tbl`!") 
      } 
    
    if( is.null(title) ) {
        title <- unique(paste0(data_tbl$GENE, " ~ ", data_tbl$EVENT))
    }
  
  PSI_min <- PSI_limits[1]
  PSI_max <- PSI_limits[2]
  
  # this should be tested more and improved or removed as it might be an overcomplication
  if ( legend == 'inside' ) {
    if( PSI_min == 0 ) { zero_label <- "0\nSkipping" }
    if( PSI_max == 100 ) { hundred_label <- "100\nInclusion" }
  } else if ( legend == 'outside' ) {
    if( PSI_min == 0 ) { zero_label <- 0 }
    if( PSI_max == 100 ) { hundred_label <- 100 }
  }
  
  # Create a rich colourbar for the PSI
  num_breaks_colour_bar <- 5
  my_breaks <- round(seq(from = PSI_min, to = PSI_max, length.out = num_breaks_colour_bar), 0)
  my_labels <- as.character(my_breaks)
  
  if ( my_breaks[1] == 0) { my_labels[1] <- zero_label }
  if ( my_breaks[length(my_breaks)] == 100) { my_labels[length(my_labels)] <- hundred_label }
  
  # change PSI colour
  data_tbl <- data_tbl |>
  mutate(txt_colour = case_when(mean_PSI >= blck_wht_PSI_col_thshld ~ "above_white",
                                mean_PSI < blck_wht_PSI_col_thshld ~ "below_black")) 

  ggplot(data_tbl) +
      aes(x = Stage, y = Tissue, fill = mean_PSI, size = log2(mean_Gene_Expr + 1) ) +
      geom_point(shape = 21, stroke = 0.2) +
      geom_text(aes(label = round(mean_PSI, 0), color = txt_colour), size = 1.75, family = "Arial") + 
      scale_size(range = c(2, 6), breaks = c(1:10), name = 'log2(TPMs)') +
      scale_x_discrete(expand = expansion(mult = 0.05, add = 0.1) ) +
      scale_colour_manual(values = c("below_black" = "black", "above_white" = "white"), guide = 'none') +
      coord_cartesian(clip = 'off') +
      labs(title = title ) +
      theme_classic(base_family = "Arial") +
      theme(axis.text = element_text(colour = 'black', size = 5),
            axis.text.x = element_text(margin = margin(t = 0, unit = "mm")),
            axis.text.y = element_text(hjust = 1, margin = margin(r = 0, unit = "mm")),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(linewidth = 0.25, colour = "black"),
            axis.ticks.length.y = unit(x = 1, units = "mm"),
            axis.line = element_line(linewidth = 0.25, colour = "black"),
            plot.title = element_text(size = 5, vjust = 0, hjust = 0.05, margin = margin(b = -1, unit = "mm")),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.grid.major.y = element_line(colour = "gray84", linewidth = 0.2) ) -> core_plot
    
  if ( colour_bar == 'BlueRed' ) {
    core_plot <- core_plot +
      scale_fill_gradientn(colours = met.brewer("Hiroshige", 9, direction = -1),
                           breaks = my_breaks, name = "Mean PSI", 
                           limits = c(PSI_min, PSI_max), na.value = "gray84",
                           labels = my_labels )
  } else if ( colour_bar == 'viridis') {
    library(viridis)
    core_plot <- core_plot +
      scale_fill_continuous(type = "viridis", direction = -1,
                            limits = c(PSI_min, PSI_max),
                            breaks = my_breaks, na.value = "gray84",
                            name = "Mean PSI", labels = my_labels)
  } else {
    stop("The parameter 'colour_bar' must be either: 'BlueRed' or 'viridis'")
  }
  
  if (legend == 'inside') {
      core_plot +
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
      stop("The parameter 'legend' must specify the position of the legend relative to the plot.")
  }

  if ( is.null(plot_name) ) {
      plot_name <- paste0("Mouse_Devel_", unique(data_tbl$GENE), "_Expr_", 
                          unique(data_tbl$EVENT), "_PSI_", 
                          width, "x", height, "cm.pdf" ) 
  }
  ## This needs to be double checked... what was I doing here? Why 2 ifs for save_plot?
  if (save_plot) {
      if( is.null (out_plot_dir) ) { 
          anal_dir <- file.path('/users/mirimia/narecco/projects/07_Suz12AS/analysis')
          out_plot_dir <- file.path(anal_dir, 'tools_output/ENCODE_Mouse_Development', 
                                    format(Sys.Date(), "%Y_%m_%d"))
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
      stop("The parameter `save_plot` must be logical, either TRUE or FALSE!")
  }
       
}

#' Nice plot to check the density distribution of the correlation of one AS event PSI vs all genes expression levels.
#'
#' @param data A tibble generated with `gimme_PSI_expr_corr`.
#' @param binwidth Correlation bins width. Default 0.05. 
#'
#' @return A ggplot2 plot.
#' @import ggplot2 
#' @export
#'
#' @examples
#' gimme_PSI_expr_corr(inclusion_tbl = psi_path, vst_id = "HsaEX0000001", 
#'                     vst_expression_tbl = expr_path, corr_method = "spearman", 
#'                     verbose = F, use = "complete.obs") |>
#'                     plot_corr_dist()
plot_corr_dist <- function(data, binwidth = 0.05 ){
  num_genes <- nrow(data)
  ggplot(data) +
    aes(x = Correlation) +
    geom_histogram(aes(fill = after_stat(x)), binwidth = binwidth) +
    geom_density(aes( y = binwidth * ..count..) ) + 
    geom_vline(xintercept = 0, linetype = 'solid') +
    scale_fill_gradient2(low = "dodgerblue", mid = "white",
                         midpoint = 0, high = "firebrick3",
                         guide = "none") +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    coord_cartesian(expand = F) +
    labs(x = "Correlation", y = "Number of genes",
         title = paste0("Distribution of the correlation of one AS event PSI ", 
                        "against ", num_genes, " genes") ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank())
}

#' Plot to check the top and bottom correlating genes of one AS event PSI vs all genes expression levels. (1 vs many)
#'
#' @param data A tibble generated with `all_gene_expr_corr(num_genes = <NUM>, map_ID_2_names = T)`
#' @param num_genes How many top and bottom genes to plot? A vector of length 2 specifying the number of genes to plot for the bottom and for the top correlating genes. Must be above 0.
#' @param ... Extra parameters passed to `geom_point` 
#'
#' @return A ggplot2 plot.
#' @import ggplot2 
#' @export
#'
#' @examples
#' #' gimme_PSI_expr_corr(inclusion_tbl = psi_path, vst_id = "HsaEX0000001", 
#'                        vst_expression_tbl = expr_path, corr_method = "spearman", 
#'                        verbose = F, use = "complete.obs", quality_thrshld = "LOW",
#'                        num_genes = 10, expr_min_mean_fltr = 10,
#'                        map_ID_2_names = T, species = "hsapiens") |>
#'                    plot_best_corr_genes(size = 5, stroke = 0.25)  
plot_best_corr_genes <- function(data, num_genes, ...) {
  
  required_cols <- c("Correlation", "external_gene_name")    
  if (!all(required_cols %in% colnames(data) )) {
    stop("The input data needs the following columns:\n", 
         paste0(required_cols, collapse = ", "))
  }
  
  if ( !missing(num_genes) ) {
      if ( length(num_genes) == 1 ) {
          num_genes <- rep(x = num_genes, times = 2)
      } else if ( length(num_genes) >= 3 ){
          num_genes <- num_genes[1:2]
      } else if ( any(num_genes < 1) ) {
          num_genes <- c(1,1)
      }
  }
  
  if( missing(num_genes) ) {
      data_plot <- data
  } else {
      data |> (\(x) {
          rbind( head(x, num_genes[2]), tail(x, num_genes[1]) )
      })() -> data_plot
  }
 
   left_nudge <- round((min(data_plot$Correlation) - 0) * 0.01, 3)
  right_nudge <-  round((max(data_plot$Correlation) - 0) * 0.01, 3)
  
  ggplot(data_plot) +
    aes(x = Correlation, y = external_gene_name, fill = Correlation) +
    geom_segment(aes(x = 0, xend = Correlation, yend = external_gene_name)) +
    geom_point(shape = 21, ...) +
    
    # Positively correlating genes in the plot
    geom_text(inherit.aes = F, 
              data = subset(data_plot, Correlation >= 0),
              aes(label = external_gene_name, x = 0, y = external_gene_name),
              nudge_x = left_nudge, hjust = 1 ) +
    # Negatively correlating genes in the plot
    geom_text(inherit.aes = F, 
              data = subset(data_plot, Correlation < 0),
              aes(label = external_gene_name, x = 0, y = external_gene_name),
              nudge_x = right_nudge, hjust = 0 ) +
    
    scale_x_continuous(n.breaks = 10) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white",
                         midpoint = 0, high = "firebrick3",
                         guide = "none") +
    labs(x = "Correlation" ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank() )
}

#' Plot a heatmap with all correlations (many vs many) between an AS event PSI and a set of genes. 
#'
#' @param inclusion_tbl path to vast-tools inclusion table that contains a vst_id event.
#' @param vst_id vast-tools alternative splicing event to grep in the `inclusion_tbl`.
#' @param quality_thrshld vast-tools event quantification quality score threshold. Must be one of "N", "VLOW", "LOW", "OK", "SOK". For more info read the official documentation [here](https://github.com/vastgroup/vast-tools#combine-output-format) under "Column 8, score 1".
#' @param vst_expression_tbl Path to a vast-tools expression table (either cRPKM or TPM).
#' @param min_mean_count Filter out low expressed genes in the table read from `inclusion_tbl`. Defines the minimum row mean expression value across all samples that a gene must have to be selected. 
#' @param gene_ids Character vector containing valid ENSEMBL gene IDs to be used for calculating the correlations.
#' @param map_ID_2_names Logical. Whether or not to map the ENSEMBL gene IDs to gene names. Can be used only if `num_genes` is specified and the table contains ENSEMBL gene ID (check automatically). Default `TRUE`.
#' @param species Species character to use to map the ENSEMBL gene ID. Used by `gimme_mart()` to built a bioMaRt object. Default is guessed from `vst_id`.
#' @param event_position Where the place the `vst_id` in the heatmap. Either `first` for placing it as a first plotted element, or `middle` to plot it in the middle.
#' @param corr_method Either `spearman`, `pearson`, or `kendall` passed to the function `cor()`. 
#' @param return_data Logical. Return the heatmap or the dataframe. Default is `FALSE` returning the heatmap.
#' @param verbose 
#' @param ... 
#'
#' @return A ggplot2 plot or a tibble.
#' @import tibble
#' @import dplyr
#' @importFrom forcats fct_inorder
#' @import ggplot2
#' @export
#' @description The order of the genes in the lower triangle heatmap are defined as the input order of `gene_ids`.
#'
#' @examples
#' # Create a gene ID mapping object.
#' ensembl <- gimme_mart()
#' 
#' # Map SRSF1-12 gene names to ENSEMBL GENE IDs.
#' SRSF_IDs <- gene_name_2_ensembl_id(gene_name = paste0("SRSF", 1:12))
#' 
#' # Plot the correlations between all SRSFs genes expressions and one AS event
#' plot_corr_heatmap(inclusion_tbl = psi_path, vst_id = "HsaEX0000001",
#'                   quality_thrshld = "VLOW",
#'                   vst_expression_tbl = expr_path,
#'                   corr_method = "spearman", use = "complete.obs",
#'                   gene_ids = SRSF_IDs, verbose = F,
#'                   event_position = "middle")
plot_corr_heatmap <- function(inclusion_tbl, vst_id, quality_thrshld,
                              vst_expression_tbl, min_mean_count = 1,
                              gene_ids, 
                              map_ID_2_names = TRUE, species,
                              event_position = "first",
                              corr_method = c("spearman", "pearson", "kendall"),
                              return_data = FALSE,
                              verbose, ...
) {
  # 1 ---- Check input parameters ----
  if ( missing(vst_id) ) { stop("You didn't specified a vst_id!") } 
  if ( missing(gene_ids) ) { stop("You didn't specified a vector of ENSEMBL gene IDs!") } 
  if( length(gene_ids) < 3 ) { stop("Please provide at least 3 ENSEMBL gene IDs") }
  
  # if ( ! grepl(pattern = "^ENS", x = gene_ids[1] ) ) {
  #   stop("Are you sure that the first element of 'gene_ids' is a valid ENSEMBL",
  #        "gene ID?")
  # }
  
  if ( missing(corr_method) ) { 
    stop("You didn't specified an correlation method!",
         "Use '?cor' read more about what method to use.") 
  } 
  if ( ! any( corr_method %in%  c("spearman", "pearson", "kendall") ) ) {
    message("corr_method must be either:",
            "'spearman', 'pearson' or, 'kendall'.",
            "Use '?cor' read more about what method to use.")
  }
  
  if ( ! any( event_position %in%  c("first", "middle") ) ) {
    message("event_position must be either: 'first' or 'middle'")
  }
  
  if ( map_ID_2_names == TRUE) {
    # Guess the species
    if ( missing(species) ) { 
      species <- guess_species(vst_id = vst_id, latin_name = TRUE)
      if (verbose) { message("I believe ", vst_id, " is a ", species, " ID.") }
    }
    
    if ( is.null(species) ) {
      stop("You need to specify the species to use for mapping the ENSEMBL",
           "gene IDs! Use '?gimme_mart()' to check which species are supported")
    }
  }
  
  # 2 ---- GET PSI MATRIX ----
  gimme_psi_mat(inclusion_tbl = psi_path, vst_id = vst_id,
                quality_thrshld = quality_thrshld, 
                verbose = verbose ) -> psi_mat
  
  # 3 ---- GET EXPRESSION MATRIX  ----
  gimme_expr_mat(vst_expression_tbl = vst_expression_tbl, 
                 min_mean_count = min_mean_count, 
                 verbose = verbose) -> gene_expr_mat
  
  # 4 ---- FILTER OUT GENES THAT ARE NOT REQUIRED ---- 
  gene_expr_mat <- gene_expr_mat[, colnames(gene_expr_mat) %in% gene_ids]
  
  # 5 ---- FILTER OUT SAMPLES FOR WHICH THE PSI WAS DISCARDED ----
  gene_expr_mat <- gene_expr_mat[rownames(gene_expr_mat) %in%  rownames(psi_mat), ]
  
  # 6 ---- COMBINE THE 2 MATRIXES ----
  p_e_mat <- cbind(psi_mat, gene_expr_mat)
  
  # 7 ---- CALCULATE FULL CORRELATION MATRIX ----
  M <- cor(x = p_e_mat, method = corr_method, ...)
  # TO DO: check that the matrix is numeric and not empty
  
  # 8 ---- MAP ENSEMBL GENE IDs TO GENE NAMES ----
  if ( map_ID_2_names == TRUE ) {
    # Create an ENSEMBL's BioMaRt ojbect to map ensembl gene IDs to gene names
    ensembl <- gimme_mart(species = species)
    gene_colnames <- colnames(M)[-1]  
    ensembl_id_2_gene_name(ensembl_gene_id = gene_colnames,
                           only_gene_name = F,
                           mRt_objct = ensembl) -> map_df
    external_gene_names <- map_df$external_gene_name
    
    # Test if external gene names are unique
    if( length(unique(external_gene_names)) != length(external_gene_names) ) {
      warning("It could be that the mapped gene names from the ENSEMBL IDs",
              " are not unique, and some are duplicated. This is likely going ",
              " to give you an error at some point.")
    }
    
    # Assign mapped gene names to the matrix
    colnames(M)[-1] <- external_gene_names
    rownames(M)[-1] <- external_gene_names
    
    # Set order as genes appear in the input parameter "gene_ids"
    map_df |> 
      mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = gene_ids)) |>
      arrange(ensembl_gene_id) |>
      mutate(external_gene_name = fct_inorder(external_gene_name) ) |>
      as_tibble() -> map_df_fct
    
    genes_order <- as.character(map_df_fct$external_gene_name)
    
  } else if (map_ID_2_names == FALSE ){
    genes_order <- gene_ids
  } else { 
    stop("Are you sure map_ID_2_names is a logical?")
  }
  
  # ---- DEFINE GENE ORDER WITH AS EVENT FIRST OR IN THE MIDDLE ----
  if ( event_position == "first" ) {
    elements_order <- c(colnames(M)[1], genes_order)
  } else if ( event_position == "middle" ) {
    mid_point <- round(length(genes_order) / 2, 0)
    elements_order <- c(
      # first part OF GENES
      genes_order[1:mid_point],
      # AS event
      colnames(M)[1], 
      # second part
      genes_order[seq(from = mid_point+1, to = mid_point*2, by = 1)]
    ) 
  } else {
    stop("event_position must be either 'first' or 'middle', not: ", 
         event_position, "!")
  }
  
  # ---- TURN CORRELATION MATRIX TO A DATAFRAME KEEPING ONLY THE LOWER TRIANGLE ----
  as.data.frame(M) |>
    rownames_to_column("FROM") |>
    as_tibble() |>
    pivot_longer(cols = !starts_with("FROM"), names_to = "TO",
                 values_to = "Correlation") |>
    subset(!is.na(Correlation)) |>
    mutate(FROM = factor(FROM, levels = elements_order ) ) |>
    mutate(TO = factor(TO, levels = rev(elements_order) ) ) |>
    arrange(FROM, desc(TO)) |> 
    # FROM and TO are sorted in opposite direction as their levels are reversed
    # Filter all the rows below or equal to when TO is equal to FROM.
    filter(row_number() <= which(TO == FROM)) -> df
  # ---- PLOT LOWER TRIANGLE CORRELATION HEATMAP WITH GGPLOT2 ----   
  if ( return_data == TRUE ) {
    return(df)
  } else if ( return_data == FALSE ) {
    
    right_nudge <- length(unique(df$FROM)) * 0.042
    
    ggplot(df, aes(y = FROM, x = TO, fill = Correlation)) +
      geom_tile() +
      geom_text(data = subset(df, FROM == TO),
                aes(label = FROM), hjust = 0,
                nudge_x = right_nudge, size = 3, check_overlap = T, 
                family = "Arial")  +
      scale_fill_gradient2(low = "dodgerblue", mid = "white", 
                           midpoint = 0, high = "firebrick4", 
                           n.breaks = 10, limits = c(-1, 1)) +
      guides(fill = guide_colourbar(title = paste0(corr_method, "\ncorrelation"), 
                                    barheight = unit(5, "cm"), 
                                    barwidth = unit(2.5, "mm"),
                                    title.vjust = 0.85 )) +
      coord_fixed(clip = 'off', ratio = 1) +
      theme_classic(base_size = 8, base_family = "Arial") +
      theme(axis.text.x = element_text(hjust = 1, angle = 45, 
                                       margin = margin(r = -1), 
                                       colour = "black"),
            axis.text.y = element_text(hjust = 1, margin = margin(r = -1), 
                                       colour = "black"),
            axis.title = element_blank(), 
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.position = c(0.9, 0.75),
            legend.text = element_text(family = "Arial"),
            legend.title = element_text(family = "Arial", hjust = 0),
            legend.margin = margin(t = 2, b = 0, unit = "mm"),
            legend.background = element_blank(),
            plot.caption = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            plot.margin = margin(l = -2, r =  -2, unit = "cm"),
      ) -> p_corr_heatmap
    
    return(p_corr_heatmap)
    
  } else {
    warning("The parameter return_data must be a logical, either TRUE or ",
            "FALSE, not: ", return_data, ". I'm returnig just the data ",
            "as it is not clear what you want!")
    return(df)
  }
}
