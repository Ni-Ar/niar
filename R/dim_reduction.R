#' Given a matrix show me the principal components and other neat things.
#'
#' @param mat A matrix 
#' @param x Integer indicating principal component on X-axis. Default 1.
#' @param y Integer indicating principal component on Y-axis. Default `x + 1`.
#' @param mt `data.frame` containing `mat` additional information (i.e. metadata). Not required but necessary if one wants label or colour according to specific parameters. When not specified the function will create a simple metadata based on the `mat` column names.
#' @param mcol Character specifying one column of `mt` that contains the matrix column names.
#' @param m_fill Character specifying one column name of `mt` to use for colouring the samples.
#' @param m_label Character or logical specifying one column name of `mt` to use for labelling the points in the PCA. Set to `FALSE` (default) for omitting labels in plot. Use `'mat_col'` if `mt` is not defined but you want to use the `mat` column names for labelling.
#' @param n_top_var Integer: Number of most variable matrix rows to use for `prcomp`.
#' @param filt_mat Logical. Whether or not to remove certain rows from `mat` that contain to many `NA`. See `NA_filt_thrshld` to specify the threshold for removing.
#' @param NA_filt_thrshld Integer between 0 and 1: maximum % of NA accepted on `mat` rows. 
#' @param show_variance Logical. Show an extra plot with the variance on each component.
#' @param show_stats Logical. Show an extra plot with the summary statistics for the data in `mat`.
#' @param n_loadings Integer indicating how many top and bottom loadings to plot.
#' @param return_data Logical. If `TRUE` returns rotated data used for plotting instead of the actual plot. Can be used with `n_loadings` equal to any positive integer to return all components loadings.
#' @param real_aspect_raio Logical. If `TRUE` (deafault) represent the distances between samples as faithfully as possible. Take into account that the second component is always smaller than the first, sometimes considerably so, thus `TRUE` normalize the axis aspect ratio to the relevant ratio for the PCA plot. Adapted from by: https://f1000research.com/articles/5-1492/v2 .
#' @param ... Set extra parameter for the `prcomp` function like `scale.` (default `FALSE`) and `center` (default `TRUE`).
#'
#' @return Either a plot (created with `ggplot2`), a combination of plots ( created with `patchwork`) or a `data.frame`. 
#' 
#' @import ggplot2
#' @import magrittr
#' @import patchwork
#' @importFrom stats prcomp
#' @importFrom matrixStats rowVars
#' @importFrom dplyr left_join mutate arrange group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom forcats fct_reorder
#' @importFrom scales number_format
#' @importFrom ggrepel geom_text_repel
#' 
#' @export
#'
#' @examples
#' showme_PCA2D(mat)
#' 
#' showme_PCA2D(mat = mat,  mt = mt, mcol = "sample_name", m_fill = "replicate",
#'              x = 3, show_stats = T, m_label = F)
#'              
#' showme_PCA2D(mat = mat, mt = mt, mcol = "sample_name", n_loadings = 12)              
showme_PCA2D <- function(mat, x = 1, y = x + 1, mt, mcol, 
                         m_fill = mcol, m_label = FALSE, 
                         n_top_var = 500, filt_mat = FALSE, 
                         NA_filt_thrshld = 0.95, show_variance = FALSE, 
                         show_stats = FALSE, n_loadings = NULL, 
                         return_data = FALSE, real_aspect_ratio = TRUE, ... ) {
  
  # ----- Check input parameters
  if(missing(mat)){ stop ("Specify a matrix with 'mat = ...'") }
  
  if ( !any(class(mat) %in% "matrix") ) { stop("Input matrix is not a matrix") }
  
  if( missing(mt) ){
    # warning("You didn't specify a metadata!")
    if (is.null(colnames(mat))) {
      # If matrix doesn't have colnames create an increasing alphanumeric string
      colnames(mat) <- paste0(LETTERS, 1:ncol(mat))
    }
    # Create a simple one with just the matrix column names
    mt <- data.frame(mat_col = colnames(mat), stringsAsFactors = F)
    
    # If trying to specify mcol without a metadata return a warning
    if (!missing(mcol)) {
      warning("Dropping 'mcol' specification cause there's not metadata!")
    }
    # Set new parameters when metadata is not provided
    mcol <- colnames(mt)[1] 
    m_fill <- mcol
    if ( missing(m_label) & m_label != FALSE) {
      m_label <- mcol  
    }
  }
  
  # check if the metadata is a tibble and convert to dataframe
  if ( any(class(mt) != "data.frame") ) {
    if ( any(class(mt) != "tbl_df") )  {
      mt <- as.data.frame(mt)
    }
  }
  # stop if metadata is not a tibble
  if ( class(mt) != "data.frame")  { stop("Metadata is not a dataframe") }
  
  if ( !all(colnames(mat) %in% mt[, mcol])) {
    stop("The matrix column names are not all present in the 'mcol' column ", 
         " of the mt dataframe.\nSee for yourself:\n",
         "Matrix columns: ",  paste(colnames(mat), collapse = " " ), "\n",
         "mcol values: ", paste(mt[, mcol], collapse = " " ))
  }
  
  if ( !all(m_fill %in% colnames(mt)) ) {
    stop("The metadata data frame does NOT contain a column called ", m_fill, 
         "\nMetadata columns: ",  paste(colnames(mt), collapse = " " ))
  }
  
  if ( !all(m_label %in% c(colnames(mt), FALSE)) ) {
    stop("The metadata data frame does NOT contain a column called ", m_label, 
         "\nMetadata columns: ",  paste(colnames(mt), collapse = " " ))
  }
  
  # ----- Filter matrix
  if ( filt_mat ) {
    apply(mat, 2, function(x) {
      length(which(is.na(x)))
    } ) / nrow(mat) <= NA_filt_thrshld -> good_cols
    mat <- mat[, good_cols]
  }
  
  # ---- Perform Principal Component Analysis
  # calculate the variance for each row
  row_vars <- matrixStats::rowVars(mat)
  
  # select the n_top_var rows by variance
  indx <- order(row_vars, decreasing = TRUE)[seq_len(min(n_top_var, length(row_vars)))]
  
  # Perform PCA
  pca_data <- stats::prcomp(x = t( mat[indx, ] ), retx = TRUE, ...)
  per100Var <- round( 100*pca_data$sdev^2 / sum(pca_data$sdev^2), 1)
  
  # Check the information content on each component
  if (any(per100Var == 0)) {
    
    # If the number of components with zero variance are more than 30% of all
    # computed components print a message
    if(length(which(per100Var == 0)) >= round(length(per100Var) * 0.3, 0)) {
      message("Components: ", paste(which(per100Var == 0), collapse = ", "), 
              " explain zero variance in the data!")
    }
  }
  
  if (show_variance) {
    ggplot2::qplot(data = as.data.frame(per100Var), geom = "col", 
                   x = as.numeric(rownames(as.data.frame(per100Var))),
                   y = per100Var, show.legend = F,
                   xlab = paste0("Components n = ", length(per100Var) ), 
                   ylab = "Variance explained (%)") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0, NA), 
                                                              mult = c(0, NA))) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = element_text(colour = "black"), 
            axis.line = element_line(color = 'black'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            plot.background = element_blank()) -> p_variance
  }
  
  # ---- Reshape PCA data
  num_components <- max(as.numeric(gsub("PC", "", colnames(pca_data$x))))
  
  # Check for low number of components
  if (x > num_components) { 
    stop("Lower x: there are only ", num_components, " components!") 
  }
  if (y > num_components) { 
    y <- x
    x <- x - 1
    warning("Lowering y: there are only ", num_components, " components!\n",
            "Now y = ", y, " and x = ", x) 
  }
  # Quick check
  if(num_components != length(per100Var)) { stop('Check num components found') }
  
  # Add info column for joining
  pca_df <- as.data.frame(pca_data$x, stringsAsFactors = F)
  colnames(pca_df) <- paste0("PC", 1:num_components)
  pca_df[, mcol] <- rownames(pca_df)
  rownames(pca_df) <- NULL
  # Add metadata to the rotated data
  pca_df <- dplyr::left_join(x = pca_df, y = mt, by = mcol) 
  
  # Check that mcol for colouring plot is in the metadata
  if( !any(colnames(pca_df) == mcol) ){
    stop("The 'mcol' character is not a 'mt' column.\n",
         "Check for yourself: ", mcol, " is not in ",
         paste(colnames(mt), collapse = ", " ) )
  }
  
  # ---- Return data of plot the PCA
  if (return_data) {
    if ( is.null(n_loadings) ) {
      # if not plotting loadings return the PCA points coordinates
      return(pca_df) 
    } else if ( !is.null(n_loadings) ) {
      # if number of loadings it not NULL return the loadings dataframe
      all_loadings_df <- as.data.frame(pca_data$rotation, stringsAsFactors = F)
      return(all_loadings_df)
    } else {
      stop("There's something wrong with returning the loadings data...")
    }
  } else {
    # ---- Plot PCA
    ggplot2::ggplot(data = pca_df) +
      aes(x = pca_df[, x], y = pca_df[, y], fill = pca_df[, m_fill] ) +
      geom_point(size = 4, shape = 21, stroke = 0.2, show.legend = F) +
      labs(x = paste0(colnames(pca_df)[x]," Variance: ", per100Var[x], "%"),
           y = paste0(colnames(pca_df)[y]," Variance: ", per100Var[y], "%") ) +
      theme_bw() + 
      theme(axis.text = element_text(colour = "black"), 
            axis.line = element_line(color = 'black'),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = 0.1, colour = "gray84"),
            panel.background = element_blank(),
            panel.border = element_blank(),
            plot.background = element_blank()) -> p_pca
    
    # ---- Plot PCA: add labels to points
    if (m_label != FALSE) {
      p_pca <- p_pca +
        ggrepel::geom_text_repel(aes(label = pca_df[, m_label] ) )
    }
    # ---- Plot PCA: plot with a realistic aspect ratio 
    # https://figshare.com/articles/figure/Aspect_ratio_for_PCA_plots_/8301197/1
    if (real_aspect_ratio == TRUE) {
      plot_aspect_ratio <- round(per100Var[y] / per100Var[x], 1)
      if (plot_aspect_ratio == 0) {
        warning("There's very little variance in your data! ",
                "The ratio between the variance on the plotted components is: ",
                plot_aspect_ratio,
                "\nSetting PCA plot aspect ratio to 1.")
        plot_aspect_ratio <- 1
      }
      p_pca <- p_pca + coord_fixed(clip = "off", ratio = plot_aspect_ratio )  
    } else if (real_aspect_ratio == FALSE ) {
      p_pca <- p_pca + coord_cartesian(clip = "off")  
    } else { 
      stop("real_aspect_ratio must be a logical! Either TRUE or FALSE)")
    }
    # p_pca <- p_pca + coord_fixed(clip = "off", ratio = plot_aspect_ratio )     
    
    # --- Plot PCA Loadings
    if ( !is.null(n_loadings) ) {
      if (n_loadings > length(pca_data$rotation[, x]) ) {
        warning("Don't ask for more loading to plot than what there ", 
                "actually are.\n", "Loadings in input matrix: ", 
                length(pca_data$rotation[, x]), ". Demanded: ", n_loadings)
        
      }
      pca_loadings <- data.frame( PCx = pca_data$rotation[, x], 
                                  stringsAsFactors = F) 
      pca_loadings$Observations <- rownames(pca_loadings)
      pca_loadings %>%
        dplyr::mutate(Observations = forcats::fct_reorder(Observations, PCx, .desc = T)) %>%
        dplyr::arrange(.data = .,desc(PCx) ) %>% {
          rbind(head(., n_loadings), tail(., n_loadings))
        } -> loadings_data
      
      mid_line <- round(median(loadings_data$PCx), 2)
      
      ggplot2::ggplot(loadings_data) +
        aes(x = Observations, y = PCx) +
        geom_hline(yintercept = mid_line, linetype = 2, colour = 'grey61') +
        geom_point(col = 'black', size = 1.2) +
        geom_segment(aes(x = Observations, xend = Observations,
                         y = mid_line, yend = PCx) ) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
        coord_cartesian(clip = "off") +
        labs(title = paste("Loadings on principal component", x), x = "",
             y = "Loadings eigenvectors")+
        theme_bw() + 
        theme(axis.text = element_text(colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1), 
              axis.line = element_line(color = 'black'),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              plot.background = element_blank()) -> p_loadings
    }
    # --- Plot Summary Statistics
    if ( show_stats ) {
      # Suppress summarise info
      options(dplyr.summarise.inform = FALSE)
      as.data.frame(mat) %>%
        tidyr::pivot_longer(cols = everything()) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(Minimum = min(value), 
                         Median = median(value),
                         Mean = mean(value), 
                         Maximum = max(value),
                         StDev = sd(value)) %>%
        tidyr::pivot_longer( cols = !matches("name"), 
                             names_to = "summary_stats") %>%
        ggplot2::ggplot(aes(y = summary_stats, x = value, fill = summary_stats)) +
        geom_boxplot(show.legend = F) + 
        scale_x_log10() +
        labs(y = "Stats of columns in mat") +
        theme_bw() + 
        theme(axis.text = element_text(colour = "black"),
              axis.line = element_line(color = 'black'),
              axis.title.x = element_blank(),
              panel.border = element_blank(), 
              panel.background = element_blank(),
              plot.background = element_blank()) -> p_stats
    }
    # --- Decide what plots to return
    if ( all( show_variance & show_stats ) ) {
      (p_variance + p_stats ) / p_pca
    } else if ( show_variance ) {
      p_variance + p_pca
    } else if ( show_stats ) {
      p_stats + p_pca
    } else if ( all(!is.null(n_loadings) & show_stats == F & show_variance == F)) {
      p_loadings
    } else {
      p_pca
    }
  }
}
