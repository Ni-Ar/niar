#' Work in PROGRESS. Convert a general DEseq2 result object into a data frame.
#'
#' @param res A DESeq2 object created with `results()` or `lfcShrink()`
#' @param map_gene_IDs Logical. Do you want to convert ensembl gene IDs to gene names?
#' @param histone_PTMs Logical. Are you analysing histone PTMs? Used to name the first column `First_Col` to be compatible with other functions.
#' @param pval_cutoff Numeric threshold to consider observations statistically significant. Default `0.05`.
#' @param log2FC_cutoff Numeric threshold to consider observations moving in a given direction (UP or DOWN). Default `0`.
#'
#' @return A tibble data frame
#' @importFrom dplyr left_join relocate between arrange mutate case_when desc
#' @importFrom tibble rownames_to_column as_tibble
#' @export 
#'
#' @examples
#' 
#' # Somewhere in your DESeq2 analysis
#' resLFC <- lfcShrink(dds, coef = "condition_treated_vs_untreated", type = "apeglm")
#' res_as_df <- res2df(resLFC)
res2df <- function(res, map_gene_IDs = F, histone_PTMs = F, pval_cutoff = 0.05,
                   log2FC_cutoff = 0) {
  
  if (histone_PTMs == TRUE) {
    rownames_colname <- "First_Col"
  } else if ( histone_PTMs == FALSE ) {
    rownames_colname <- "ensembl_gene_id"
  } else{
    stop("Parameter 'histone_PTMs' must be either TRUE or FALSE")
  }
  
  if ( ! between(x = pval_cutoff, left = 0, right = 1) ){
    stop("Parameter 'pval_cutoff' must be a number between 0 and 1.")
  }
  
  res |>
    as.data.frame() |>
    rownames_to_column(rownames_colname) |> 
    as_tibble() -> tbl
  
  if (map_gene_IDs == T) {
    tbl |>
      left_join(ensembl_symbol_df, by = "ensembl_gene_id") |>
      relocate(external_gene_name, .before = baseMean) -> tbl
  }
  
  tbl |>
    arrange(desc(log2FoldChange)) |>
    mutate(Direction = case_when(padj <= pval_cutoff & log2FoldChange >= log2FC_cutoff ~ "UP",
                                 padj <= pval_cutoff & log2FoldChange <= -log2FC_cutoff ~ "DOWN",
                                 padj > pval_cutoff | is.na(padj) ~ "None") ) -> df
  return(df)
}


#' Convert a DESeq2 dds object into a tibble
#'
#' @param deseq_dataset A DESeq2 object.
#' @param tidy Logical, whether or not to pivot the data into a long format
#' @param counts_are_genes Logical, by default this function assumes the counts are from gene expression experiments where the IDs are ensembl gene IDs.
#'
#' @return A tibble
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect contains
#' @export
#'
#' @examples
#' tidy_counts <- dds2counts(dds)
dds2counts <- function(deseq_dataset, tidy = T, counts_are_genes = T, norm_counts = T){
  
  if(missing(x = deseq_dataset) ){
    stop('Please specify a DESeq2 object with dds = ')
  }
  
  if (counts_are_genes == TRUE) {
    rownames_ID <- "ensembl_gene_id"
  } else if( counts_are_genes == FALSE ){
    rownames_ID <- "generic_counts_id"
  } else {
    stop("'counts_are_genes' must be a logical!")
  }
  
  if (norm_counts == TRUE) {
    val_col_name <- "Norm_counts"
  } else if( norm_counts == FALSE) {
    val_col_name <- "Counts"
  } else {
    stop("'norm_counts' must be a logical!")
  }
  
  counts(deseq_dataset, normalized = norm_counts) |>
    as.data.frame() |>
    rownames_to_column(rownames_ID) |>
    as_tibble() -> counts
  
  if (counts_are_genes == TRUE) {
    # remove the version (.X) from the ensembl gene id
    counts <- mutate(counts, 
                     ensembl_gene_id = str_remove(string = ensembl_gene_id, 
                                                  pattern = "\\.[0-9]*$") ) 
  }
  
  if (tidy == TRUE) {
    counts <- pivot_longer(data = counts, cols = !contains(rownames_ID), 
                           names_to = "Sample", values_to = val_col_name)
  }
  
  return(counts)
}

#' Helper function to check distribution of reads in DESeq2 object
#'
#' @param deseq_dataset 
#' @param xlim 
#' @param title 
#' @param min_counts Small number of pseudocounts added to the normalised counts before the `log2()`.
#'
#' @return A nice ggplot2 density plot
#' @import ggplot2
#' @import BiocGenerics
#' @export
#'
#' @examples
#' show_dds_counts_freq(dds)
show_dds_counts_freq <- function(deseq_dataset, xlim = c(0, 20),
                                 title = "Here goes the title", min_counts = 5) {
  # To return normalised counts first estimate size factors.
  # If the dds is coming from tximeta import the average transcript lengths 
  # are used as offsets which are gene- and sample-specific normalization factors
  if( is.null(normalizationFactors(deseq_dataset) ) ) {
    deseq_dataset <- estimateSizeFactors(deseq_dataset, quiet = T)
  } 
  # if ( is.null(sizeFactors(deseq_dataset) ) ) {
  #   
  #     
  #   deseq_dataset <- estimateSizeFactors(deseq_dataset, quiet = T)
  # }

  
  dds2counts(deseq_dataset, tidy = T, norm_counts = T) |>
    ggplot(aes(x = log2(Norm_counts + 0.5), y = after_stat(density), fill = Sample)) +
    geom_vline(xintercept = log2(min_counts + 0.5) , linetype = 'solid',
               color = 'firebrick1') +
    geom_density(show.legend = F, alpha = 0.25) +
    scale_x_continuous(limits = xlim, expand = expansion(add = 0, mult = c(0, 0.05) )) +
    scale_y_continuous(expand = expansion(add = 0, mult = c(0, 0.05))) +
    labs(x = "log2(DESeq2 Normalised counts + 0.5)", y = "Density",
         title = title) +
    theme_bw(base_family = "Arial") + 
    theme(axis.text = element_text(colour = "black"),
          panel.background = element_blank(),
          panel.border = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.25, colour = "black"),
          plot.background = element_blank()) -> density_plot
  return(density_plot)
}


#' Filter lowly abundant genes with the possibility of doing a normalisation on the genes to correct for trended biases. 
#' This function will not estimate size factors of a DESeq2 object, but it will add normalization factors
#' to the DESeq2 object if the cyclic_loess is TRUE.
#' Plotting options allow to see how the read distribution changes using different filtering and normalisation parameters
#'
#' @param deseq_dataset 
#' @param filt_method 
#' @param min_counts 
#' @param verbose 
#' @param show_mltdnsty 
#' @param cyclic_loess 
#'
#' @return a DESeq2 object and a plot
#' @importFrom csaw normOffsets
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @import DESeq2
#' @import patchwork
#' @export
#'
#' @examples
#' dds_fltrd <- filter_dds(dds, filt_method = "mean")
filter_dds <- function(deseq_dataset, filt_method = c("sum", "mean", "max"),
                       min_counts = 5, verbose = T, show_mltdnsty = T, 
                       cyclic_loess = F) {
  
  stopifnot( any( filt_method %in% c("sum", "mean", "max") ) )
  
  prefltr <- nrow(deseq_dataset)
  deseq_dataset <- deseq_dataset[ rowSums(counts(deseq_dataset)) > 1, ]
  minfltr <- nrow(deseq_dataset) 
  
  if( show_mltdnsty ){
    show_dds_counts_freq(deseq_dataset, min_counts = min_counts,
                         title = "Unfiltered gene counts frequencies") -> p_unfltrd
  }
  
  # Filter Lowly Abundant Genes: AT LEAST "min_counts" READS PER GENE ACROSS ALL SAMPELS
  indx <- apply( counts(deseq_dataset, norm = FALSE), 1, function(x) {
    eval( expr = parse(text = paste0(filt_method, "(x) > ", min_counts ) ) )
  })
  deseq_dataset <- deseq_dataset[indx, ]
  postfltr <- nrow(deseq_dataset) 
  
  if( verbose ){
    # summaries filtering steps
    message("Total number of genes identifieds: ", prefltr, "\n",
            "Genes with least one count in one sample: ", minfltr, "\n",
            "Filter using: ", filt_method, " > ", min_counts, " counts: ", postfltr)
  }
  
  rm(prefltr, minfltr, postfltr)
  
  if( show_mltdnsty ){
    show_dds_counts_freq(deseq_dataset, min_counts = min_counts,
                         title = paste0("Filter using: ", filt_method, " > ",
                                        min_counts, " counts") ) -> p_fltr
  }
  
  if( cyclic_loess ){
    # Normalize trended biases across libraries using cyclic loess
    
    # create a SE from the DESeq2 object for csaw function
    dds_se <- SummarizedExperiment(
      assays = list(counts = DESeq2::counts(deseq_dataset, normalized = FALSE)),
      colData = colData(deseq_dataset)
    )
    
    dds_se$totals <- colSums( DESeq2::counts(deseq_dataset) )
    normFacs <- normOffsets( object = dds_se, assay.id = "counts", 
                             iterations = 5L, se.out = F )
    
    # the normalization factors matrix should contain no zeros
    # and have a geometric mean near 1 for each row
    normFacs <- normFacs / exp( rowMeans( log(normFacs) ) )
    
    message("Sanity Check:",
            "Is the geometric mean across samples equal to one? ")
    
    geom_mean_cl_norm_factors <- apply(normFacs, 1, function(x) { exp(mean(log(x)))})
    names(geom_mean_cl_norm_factors) <- NULL
    
    message( all.equal( geom_mean_cl_norm_factors,
                        rep(1, dim(deseq_dataset)[1]) ) )
    
    # Apply these new normalisation factors to the DESeq2 object
    normalizationFactors(deseq_dataset) <- normFacs
    
    if( show_mltdnsty ){
      show_dds_counts_freq(deseq_dataset, min_counts = min_counts,
                           title = "Renormalized with cyclic loess") -> p_cyclc_lss
    }
    
  }
  
  if ( show_mltdnsty ) {
    if ( cyclic_loess ) {
      print( p_unfltrd + p_fltr + p_cyclc_lss )
    } else {
      print( p_unfltrd + p_fltr )
    }
    
  }
  # Return Filtered dds
  return(deseq_dataset)
}



#' Save the counts of a DESeq2 object.
#'
#' @param deseq_dataset A dds
#' @param name A name for the output file
#' @param out_dir Path to the output directory
#' @param tidy Logical, whethere to reshape the data into a long (`TRUE`) or keep it as wide ('FALSE`) format.
#' @param ... extra arguments passed to `dds2counts()` function, like `counts_are_genes = T|F` or `norm_counts = T|F`
#'
#' @return Nothing, this function writes a file
#' @import readr
#' @export
#'
#' @examples
#' save_counts(dds, name = "Your_comparison_long", out_dir = "/path/to/dir")
save_counts <- function(deseq_dataset, name, out_dir, tidy = TRUE, ... ) {
  
  counts <- dds2counts(deseq_dataset, tidy = tidy, ...)
  
  # create output subfolders
  out_dir_counts <- file.path(out_dir, "gene_counts", format(Sys.Date(), "%Y_%m_%d") )
  if (!dir.exists(out_dir_counts)) { dir.create(path = out_dir_counts, recursive = T) }
  
  num_genes <- length(unique(counts$ensembl_gene_id))
  
  if (tidy == TRUE) {
    num_samples <- length(unique(counts$Sample))
    suffix <- "long"
  } else if ( tidy == FALSE ) {
    num_samples <- ncol(counts) - 1 # num of samples - enesembl gene IDs
    suffix <- "wide"
  } else {
    stop('tidy must be a logical!')
  }
  
  out_tab <- file.path(out_dir_counts, 
                       paste0(name, '_n', num_genes, 'xSamples', num_samples, 
                              '_', suffix, '.tab'))
  
  write_delim(x = counts, file = out_tab, delim = '\t', append = F, 
              col_names = T, quote = 'none')
}

#' This function is still experimental as it should be broken down into simpler sub-functions
#'
#' @param res 
#' @param name 
#' @param map_df 
#' @param baseMean_thrshold 
#' @param pval_adj_thrshold 
#' @param out_dir 
#'
#' @return Nothing, this function writes to file.
#' @export
#'
save_df_gsea_list <- function(res, name, map_df, baseMean_thrshold = 300, 
                              pval_adj_thrshold = 0.001, out_dir) {
  
  # the mapper df must have a column with the gene biotype to filter for protein coding genes.
  stopifnot(any(colnames(map_df) == "ensembl_gene_id"))
  stopifnot(any(colnames(map_df) == "gene_biotype"))
  
  # create output subfolders
  out_dir_res <- file.path(out_dir, "gene_tables", format(Sys.Date(), "%Y_%m_%d") )
  if (!dir.exists(out_dir_res)) { dir.create(path = out_dir_res, recursive = T) }
  
  out_dir_gsea <- file.path(out_dir, "gene_lists", format(Sys.Date(), "%Y_%m_%d") )
  if (!dir.exists(out_dir_gsea)) { dir.create(path = out_dir_gsea, recursive = T) }
  
  res <- res2df(res) |>
    mutate(ensembl_gene_id = str_remove(string = ensembl_gene_id, pattern = "\\.[0-9]*$")) 
  
  out_df <- left_join(x = res, y = map_df, by = "ensembl_gene_id") 
  
  out_tab <- file.path(out_dir_res, paste0(name, '_n', nrow(out_df), '.tab'))
  
  write_delim(x = out_df, file = out_tab, delim = '\t', append = F, 
              col_names = T, quote = 'none')
  
  out_xls <- file.path(out_dir_res, paste0(name, '_n', nrow(out_df), '.xls'))
  
  write_excel_csv2(x = out_df, file = out_xls, na = "NA", append = F, 
                   col_names = T, delim = ";", quote = "all", eol = "\n")
  
  subset(out_df, baseMean >= baseMean_thrshold) |>
    arrange(desc(log2FoldChange)) |>
    subset(gene_biotype == "protein_coding") |>
    subset(padj <= pval_adj_thrshold) |>
    subset(Direction == "UP") |>
    select(ensembl_gene_id) -> best_UP_genes
  
  out_UP_GeneList <- file.path(out_dir_gsea, paste0(name, '_UP_genes_n', 
                                                    nrow(best_UP_genes), '.txt'))
  write_delim(x = best_UP_genes, file = out_UP_GeneList, append = F, 
              col_names = F, quote = "none", delim = '\t')
  
  subset(out_df, baseMean >= baseMean_thrshold) |>
    arrange(log2FoldChange) |>
    subset(gene_biotype == "protein_coding") |>
    subset(padj <= pval_adj_thrshold) |>
    subset(Direction == "DOWN") |>
    select(ensembl_gene_id) -> best_DOWN_genes
  
  out_DOWN_GeneList <- file.path(out_dir_gsea, paste0(name, '_DOWN_genes_n', 
                                                      nrow(best_DOWN_genes), '.txt'))
  write_delim(x = best_DOWN_genes, file = out_DOWN_GeneList, append = F, 
              col_names = F, quote = "none", delim = '\t')
  
}

