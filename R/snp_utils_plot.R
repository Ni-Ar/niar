#' Return the bin location of a num given a range from start until end and a set number of bins.
#'
#' @param start 
#' @param end 
#' @param num 
#' @param num_bins 
#' @param return_full_df 
#'
#' @return A number or a dataframe
#' @importFrom dplyr mutate
#' @export
#' 
#' @description
#' A short description...
#' @details
#' If the number (`num`) is not between the start and end (`num > start & num <= end`), the function returns `NA`.
#' 
#' @examples
#' # Example 1
#' find_containing_bin(start = 5515473, end = 5515517, pos = 5515489, return_full_df = F)
#' # Example 2: bin changes when pos increase
#' a <- 10 ; b <- 95
#' sapply(c(1:30), function(x) { find_containing_bin(start = a, end = b, pos = a+x) })
#' # Example 3: returning NA when num not in the range
#' find_containing_bin(start = 50, end = 150, num = 40, return_full_df = F)
find_containing_bin <- function(start, end, num, num_bins = 50, return_full_df = F ) {
  
  # -- 0 -- if num not between start and end return NA
  if (num <= start) { return(NA) }
  if (num > end) { return(NA) }
  
  # -- 1 -- generate bins right-side closed 
  bins <- levels(cut(x = c(start:end), breaks = num_bins, dig.lab = nchar(end)) )
  
  # -- 2 -- Calculate number position inside the bins
  data.frame(interval = bins) |>
    mutate(BinNum = row_number()) |>
    mutate(left_num = as.numeric(str_extract(string = interval, pattern = "(?<=^\\()[0-9].+(?=,)"))) |>
    mutate(right_num = as.numeric(str_extract(string = interval, pattern = "(?<=,)[0-9].+(?=]$)") )) |>
    mutate(Inside = (num > left_num & num <= right_num),
           Dist = round(num - right_num, 4) ) -> binned_df
  
  containing_bin <- subset(binned_df, Inside) 
  if (return_full_df == TRUE) {
    return(binned_df)
  } else if (return_full_df == FALSE) {
    
    # if binning fails because rounding big numbers, pick the one with distance closer to zero
    if ( (nrow(containing_bin) >= 2) ) {
      containing_bin <- containing_bin[which.min(abs(containing_bin$Dist)), ]
    }
    containing_bin_num <- unique(containing_bin$BinNum)
    return(containing_bin_num)
  }else {
    stop('return_full_df must be a logical!')
  }
}

# should add the option to skip certain columns
read_bed_vcf_isec <- function(path, 
                              bed_cols = c('chr', 'start', 'end', 'EVENT', 'GENE', 'strand'),
                              vcf_cols = c('CHROM', 'POS', 'RefSNP_ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'),
                              ...) {
  
  col_names_bed_vcf <-  c(bed_cols, vcf_cols)
  
  # # Define the regex pattern to extract the fist 2 informations from the format column
  # pattern <- "^[0-9]+:(.*?)(?=:)"
  
  imported_df <- read_delim(file = path, col_names = col_names_bed_vcf,
                            delim = "\t", escape_double = FALSE, progress = FALSE,
                            col_types = cols(chr = col_character(), 
                                             start = col_integer(), 
                                             end = col_integer(), 
                                             CHROM = col_skip(), 
                                             POS = col_integer(),
                                             QUAL = col_skip(),
                                             FILTER = col_skip() ), 
                            trim_ws = TRUE, ...) 
  
  return(imported_df)
}


read_bed_vcf_isec1 <- function(path, ...) {
  
  imported_df <- read_delim(file = path, col_names = FALSE, delim = "\t",
                            escape_double = FALSE, progress = FALSE,
                            col_types = cols(X1 = col_character(), 
                                             X2 = col_integer(), 
                                             X3 = col_integer(), 
                                             X7 = col_character(), 
                                             X8 = col_integer(),
                                             X9 = col_integer() ), 
                            trim_ws = TRUE, ...)
  
  col_names <-  c('chr', 'start', 'end', 'EVENT', 'GENE', 'strand', 'CHR',
                  'POS', 'ClinVar_ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
  colnames(imported_df) <- col_names 
  return(imported_df)
}


# imports files generated with:
# `bedtools intersect -a` file.bed -b <(zcat freq.vcf.gz) -wa -wb -sorted > isect.tab`
# Where freq.vcf.gz is the ALFA sorted and chromosome filtered vcf file and 
# bed is a bed file generated with niar::inclusion_tbl2bed(). 
# Both bed and vcf files should have the same sorting order
read_bed_vcf_isec2 <- function(path, min_tot_alle_freq = 0.5, bed_cols, vcf_cols, ...) {
  
  col_names_ALFA <-  c(
    # from bed file:
    'chr', 'start', 'end', 'EVENT', 'GENE', 'strand', 
    # from vcf file:
    'CHR', 'POS', 'RefSNP_ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 
    # Population codes
    'EUR', 'AFO', 'EAS', 'AFA', 'LAC', 'LEN', 'OAS', 'SAS', 'OTR', 'AFR', 'ASN', 'TOT'
  )
  
  # Define the regex pattern to extract the fist 2 informations from the format column
  pattern <- "^[0-9]+:(.*?)(?=:)"
  
  imported_df <- read_delim(file = path, col_names = col_names_ALFA,
                            delim = "\t", escape_double = FALSE, progress = FALSE,
                            col_types = cols(chr = col_character(), 
                                             start = col_integer(), 
                                             end = col_integer(), 
                                             CHR = col_skip(), 
                                             POS = col_integer(),
                                             QUAL = col_skip(),
                                             FILTER = col_skip(), 
                                             INFO = col_skip()), 
                            trim_ws = TRUE, ...) |>
    # extract allele counts from TOTAL population
    mutate(AN_AC = str_extract(string = TOT, pattern = pattern),
           AN = as.integer(str_split_fixed(string = AN_AC, pattern = ":", n = 2)[,1]),
           AC = str_split_fixed(string = AN_AC, pattern = ":", n = 2)[,2]) |>
    mutate(AC = as.integer(str_extract(string = AC, pattern = '^[0-9]+'))) |>
    # calculate allele frequency in total population
    mutate(TOT_ALLELE_FREQ = round(AC/AN, 4)*100, .before = FORMAT ) |>
    # filter out alleles lower than
    subset(TOT_ALLELE_FREQ >= min_tot_alle_freq ) |>
    # arrange(desc(TOT_ALLELE_FREQ)) |>
    select(!c(AN_AC, AN, AC))
  
  return(imported_df)
}


# example 1
# tbv <- read_bed_vcf_isec1(tisas_CV) |>
#        parse_bed_vcf_isec(left_slop = 200, right_slop = 150, 
#                           snp_ID_col = 'ClinVar_ID', left_bins_num = 40, 
#                           centre_bins_num = 50, right_bins_num = 30, 
#                           check_bins = T)
parse_bed_vcf_isec <- function(df, snp_ID_col = 'RefSNP_ID', left_slop, 
                               right_slop, left_bins_num = 50, centre_bins_num = 50, 
                               right_bins_num = 25, check_bins = T) {
  # -- 1 -- check input
  stopifnot(any(snp_ID_col %in% colnames(df)))
  required_cols <- c('GENE', 'EVENT', 'chr', 'start', 'end', 'strand', 'POS')
  
  stopifnot(all(required_cols %in% colnames(df)))
  
  # Slopping had to be applied in a strand-aware fashion.
  if (missing(left_slop)) { stop('Specify the upstream extension you applied to the bed file') }
  if (missing(right_slop)) { stop('Specify the dowstream extension you applied to the bed file') }
  
  # -- 2 -- Parse and de-construct the slopped bed regions
  df |> 
    group_by(GENE, EVENT) |>
    mutate(Aacceptor = ifelse(strand == "+", yes = start + left_slop, no =  end - left_slop), .after = start ) |>
    mutate(Adonor = ifelse(strand == "+", yes = end - right_slop, no = start + right_slop), .after = end ) |>
    mutate(ex_start = ifelse(strand == "+", yes = Aacceptor, no = Adonor), .after = end )  |>
    mutate(ex_end = ifelse(strand == "+", yes = Adonor, no = Aacceptor), .after = ex_start )  |>
    # pass the variable of snp_ID_col to group_by as a column name with: !!sym(snp_ID_col) 
    group_by(EVENT, POS, !!sym(snp_ID_col) ) -> parsed_df
  
  # -- 3 -- Split regions by strand
  pos_df <- subset(parsed_df, strand == "+")
  neg_df <- subset(parsed_df, strand == "-")
  
  # -- 4 -- Map the bin in which the SNP (POS) falls in
  pos_binned <- pos_df |>
    mutate(up_bin = find_containing_bin(start = start, end = ex_start, num = POS, num_bins = left_bins_num) ) |>
    mutate(ce_bin = find_containing_bin(start = ex_start, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
    mutate(do_bin = find_containing_bin(start = ex_end, end = end, num = POS, num_bins = right_bins_num) ) 
  
  # NOTE: reverse bin location on the reverse strand
  # e.g. for the upstream region if a SNV is in the third bin in the plot it will need to be shown at (left_bins_num+1) - 3
  neg_binned <- neg_df |>
    mutate(up_bin = (left_bins_num+1) - find_containing_bin(start = ex_end, end = end, num = POS, num_bins = left_bins_num) ) |>
    mutate(ce_bin = (centre_bins_num+1) - find_containing_bin(start = ex_start, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
    mutate(do_bin = (right_bins_num+1) - find_containing_bin(start = start, end = ex_start, num = POS, num_bins = right_bins_num) )
  
  # -- 5 -- Combined strands SNP binning counts in one df and sort it
  rbind(pos_binned, neg_binned) |>
    ungroup() |>
    arrange(chr, start) -> processed_df

  # -- 6 -- Make a dataframe where the bin numbers of each SNP is reshaped in a long format and ready for plotting
  processed_df |>
    select(GENE, EVENT, strand, POS, !!sym(snp_ID_col), ends_with('_bin')) |>
    pivot_longer(cols = ends_with('_bin'), names_to = 'region', values_to = 'Bin_Number') |>
    subset( !is.na(Bin_Number) ) |> 
    ungroup() -> processed_df
  
  # -- 7 -- Make sure all bins are in the dataframe
  # This is then needed at step 9
  missing_up_bin <- c(left_bins_num:1)[which(!left_bins_num:1 %in% sort(unique(subset(processed_df, region == 'up_bin')$Bin_Number)))]
  missing_ce_bin <- c(centre_bins_num:1)[which(!centre_bins_num:1 %in% sort(unique(subset(processed_df, region == 'ce_bin')$Bin_Number)))]
  missing_do_bin <- c(right_bins_num:1)[which(!right_bins_num:1 %in% sort(unique(subset(processed_df, region == 'do_bin')$Bin_Number)))]
  
  if (length(missing_up_bin) >= 1 ) { message('One bin in upstream region has zero counts')}
  
  if (length(missing_ce_bin) >= 1 ) { message('One bin in centre has zero counts')}
  
  if (length(missing_do_bin) >= 1 ) { 
    message('Bin ', missing_do_bin, ' in dowstream region has zero counts')
  }
  
  
  # -- 8 -- Count how many SNPs are in each bin of each region
  processed_df |>
    group_by(region, Bin_Number) |> mutate(Var_counts = n()) |> ungroup() |>
    group_by(EVENT, region, Bin_Number) |> 
    mutate(Num_SNPs_per_events_per_bin = n()) |> 
    # arrange(desc(Num_SNPs_per_events_per_bin)) |> 
    ungroup() |>
    # make dataframe smaller by dropping some info
    select(!c(POS, !!sym(snp_ID_col) ) ) |> unique() -> processed_df
  
  # -- 9 -- Add empty row for missing bins for consistency in data format
  # TO DO: implement this fix also for upstream and center region.
  
  if (length(missing_do_bin) >= 1 ) { 
   empty_bin_rows_down <- data.frame(GENE = NA, EVENT = NA, strand = '+',
                                      region = 'do_bin', Bin_Number = missing_do_bin,
                                      Var_counts = 0, Num_SNPs_per_events_per_bin = 0)
    processed_df <- rbind(processed_df, empty_bin_rows_down)
  }
  
  # -- 10 -- Count how many events (aka exons) are in a specific bin, 
  # meaning the SNP count in such bin originates out of how many events.
  processed_df |>
    group_by(region, Bin_Number) |>
    mutate(Num_events_in_bin = n(), .after = Var_counts) |>
    ungroup() |>
    # change to prettier names
    mutate(region = gsub('_bin', '', region), region = gsub('ce', 'centre', region),
           region = gsub('up', 'upstream', region), region = gsub('do', 'downstream', region)) |>
    mutate(region = factor(region, levels = c('upstream', 'centre', 'downstream'))) |>
    # make to artificial bins for plotting
    mutate(ReOrdered_Bin_Number = case_when(
      region == 'upstream' ~ Bin_Number - left_bins_num,
      region == 'centre' ~ Bin_Number,
      region == 'downstream' ~ Bin_Number + centre_bins_num)
    ) |> arrange(region, Bin_Number) -> tidy_counts
  
  # -- 11 -- Optionally check stats
  if (check_bins) {
    tidy_counts |> select(region, Bin_Number) |>  unique() |> 
      group_by(region) |> summarise(number_of_bins = n() ) |> print()
  }
  return(tidy_counts)
}


#' Get how many SNPs variation fall within each bin in the different regions of the bed file.
#'
#' @param df 
#' @param snp_ID_col 
#' @param left_slop 
#' @param right_slop 
#' @param left_bins_num 
#' @param centre_bins_num number of bins to divide the centre region (before slopping).
#' @param right_bins_num 
#' @param start_centre_len Do not bin all centre region at once. Just consider the first N part of the centre region.
#' @param end_centre_len Just consider the last N part of the centre region.
#' @param check_bins
#'
#' @return A data.frame
#' @export
#'
count_snps_in_bins <- function(df, snp_ID_col = 'RefSNP_ID', 
                               left_slop, 
                               right_slop, 
                               left_bins_num = 50, 
                               centre_bins_num = 50, 
                               right_bins_num = 25, 
                               start_centre_len = 39, 
                               end_centre_len = 39,
                               check_bins = T) {
  # -- 1 -- check input
  stopifnot(any(snp_ID_col %in% colnames(df)))
  required_cols <- c('GENE', 'EVENT', 'chr', 'start', 'end', 'strand', 'POS')
  
  stopifnot(all(required_cols %in% colnames(df)))
  
  # Slopping had to be applied in a strand-aware fashion.
  if (missing(left_slop)) { stop('Specify the upstream extension you applied to the bed file') }
  if (missing(right_slop)) { stop('Specify the dowstream extension you applied to the bed file') }
  
  # -- 2 -- Parse and de-construct the slopped bed regions
  df |> 
    group_by(GENE, EVENT) |>
    mutate(Aacceptor = ifelse(strand == "+", yes = start + left_slop, no =  end - left_slop), .after = start ) |>
    mutate(Adonor = ifelse(strand == "+", yes = end - right_slop, no = start + right_slop), .after = end ) |>
    mutate(ex_start = ifelse(strand == "+", yes = Aacceptor, no = Adonor), .after = end )  |>
    mutate(ex_end = ifelse(strand == "+", yes = Adonor, no = Aacceptor), .after = ex_start ) |>
    mutate(centre_first = ex_start + start_centre_len) |>
    mutate(centre_last = ex_end - end_centre_len) |>
    # pass the variable of snp_ID_col to group_by as a column name with: !!sym(snp_ID_col) 
    group_by(EVENT, POS, !!sym(snp_ID_col) ) -> parsed_df
  
  # -- 3 -- Split regions by strand
  pos_df <- subset(parsed_df, strand == "+")
  neg_df <- subset(parsed_df, strand == "-")
  
  
  # -- 4 -- Map the bin in which the SNP (POS) falls in
  # cf: centre first
  # cl: centre last
  if (nrow(pos_df) >= 1 ) {
    pos_binned <- pos_df |>
      mutate(up_bin = find_containing_bin(start = start, end = ex_start, num = POS, num_bins = left_bins_num) ) |>
      # mutate(ce_bin = find_containing_bin(start = ex_start, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
      mutate(cf_bin = find_containing_bin(start = ex_start, end = centre_first, num = POS, num_bins = centre_bins_num) ) |>
      mutate(cl_bin = find_containing_bin(start = centre_last, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
      mutate(do_bin = find_containing_bin(start = ex_end, end = end, num = POS, num_bins = right_bins_num) ) 
  } else if ( nrow(pos_df) == 0) {
    pos_binned <- pos_df |>
      mutate(up_bin = 0, cf_bin = 0, cl_bin = 0, do_bin = 0)
  } else {
    stop('Not sure how many rows are in pos_df')
  }
  
  # NOTE: reverse bin location on the reverse strand
  # e.g. for the upstream region if a SNV is in the third bin in the plot it will need to be shown at (left_bins_num+1) - 3
  if ( nrow(neg_df) >= 1 ) {
    neg_binned <- neg_df |>
      mutate(up_bin = (left_bins_num+1) - find_containing_bin(start = ex_end, end = end, num = POS, num_bins = left_bins_num) ) |>
      # mutate(ce_bin = (centre_bins_num+1) - find_containing_bin(start = ex_start, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
      mutate(cf_bin = (centre_bins_num+1) - find_containing_bin(start = centre_last, end = ex_end, num = POS, num_bins = centre_bins_num) ) |>
      mutate(cl_bin = (centre_bins_num+1) - find_containing_bin(start = ex_start, end = centre_first, num = POS, num_bins = centre_bins_num) ) |>
      mutate(do_bin = (right_bins_num+1) - find_containing_bin(start = start, end = ex_start, num = POS, num_bins = right_bins_num) )
  
  } else if ( nrow(neg_df) == 0 ) {
    neg_binned <- neg_df |>
      mutate(up_bin = 0, cf_bin = 0, cl_bin = 0, do_bin = 0)
  } else {
    stop('Not sure how many rows are in neg_df')
  }
  
  # -- 5 -- Combined strands SNP binning counts in one df and sort it
  rbind(pos_binned, neg_binned) |>
    ungroup() |>
    arrange(chr, start) -> processed_df
  
  # -- 6 -- Make a dataframe where the bin numbers of each SNP is reshaped in a long format and ready for plotting
  processed_df |>
    select(GENE, EVENT, strand, POS, !!sym(snp_ID_col), ends_with('_bin')) |>
    pivot_longer(cols = ends_with('_bin'), names_to = 'region', values_to = 'Bin_Number') |>
    subset( !is.na(Bin_Number) ) |> 
    ungroup() -> processed_df
  
  # -- 7 -- Make sure all bins are in the dataframe
  # TO DO: implement this fix also for cf and cf
  # This is then needed at step 9
  missing_up_bin <- c(left_bins_num:1)[which(!left_bins_num:1 %in% sort(unique(subset(processed_df, region == 'up_bin')$Bin_Number)))]
  # missing_ce_bin <- c(centre_bins_num:1)[which(!centre_bins_num:1 %in% sort(unique(subset(processed_df, region == 'ce_bin')$Bin_Number)))]
  missing_do_bin <- c(right_bins_num:1)[which(!right_bins_num:1 %in% sort(unique(subset(processed_df, region == 'do_bin')$Bin_Number)))]
  
  if (length(missing_up_bin) >= 1 ) { 
    message('Bin ', missing_up_bin, ' in upstream region has zero counts')
    }
  
  # if (length(missing_ce_bin) >= 1 ) { message('One bin in centre has zero counts')}
  
  if (length(missing_do_bin) >= 1 ) { 
    message('Bin ', missing_do_bin, ' in dowstream region has zero counts')
  }
  
  # -- 8 -- Count how many SNPs are in each bin of each region
  processed_df |>
    group_by(region, Bin_Number) |> mutate(Var_counts = n()) |> ungroup() |>
    group_by(EVENT, region, Bin_Number) |> 
    mutate(Num_SNPs_per_events_per_bin = n()) |> 
    # arrange(desc(Num_SNPs_per_events_per_bin)) |> 
    ungroup() |>
    # make dataframe smaller by dropping some info
    select(!c(POS, !!sym(snp_ID_col) ) ) |> unique() -> processed_df
  
  # -- 9 -- Add empty row for missing bins for consistency in data format
  if (length(missing_up_bin) >= 1 ) { 
    empty_bin_rows_up <- data.frame(GENE = NA, EVENT = NA, strand = '+',
                                    region = 'up_bin', Bin_Number = missing_up_bin,
                                    Var_counts = 0, Num_SNPs_per_events_per_bin = 0)
    processed_df <- rbind(processed_df, empty_bin_rows_up)
  }
  
  # TO DO: implement this fix also for cf and cf
  
  if (length(missing_do_bin) >= 1 ) { 
    empty_bin_rows_down <- data.frame(GENE = NA, EVENT = NA, strand = '+',
                                      region = 'do_bin', Bin_Number = missing_do_bin,
                                      Var_counts = 0, Num_SNPs_per_events_per_bin = 0)
    processed_df <- rbind(processed_df, empty_bin_rows_down)
  }
  
  # -- 10 -- Count how many events (aka exons) are in a specific bin, 
  # meaning the SNP count in such bin originates out of how many events.
  processed_df |>
    group_by(region, Bin_Number) |>
    mutate(Num_events_in_bin = n(), .after = Var_counts) |>
    ungroup() |>
    # change to prettier names
    mutate(region = gsub('_bin', '', region), 
           # region = gsub('ce', 'centre', region),
           region = gsub('cf', 'first', region),
           region = gsub('cl', 'last', region),
           region = gsub('up', 'upstream', region), region = gsub('do', 'downstream', region)) |>
    mutate(region = factor(region, levels = c('upstream', 'first', 'last', 'downstream'))) |>
    # make to artificial bins for plotting
    mutate(ReOrdered_Bin_Number = case_when(
      region == 'upstream' ~ Bin_Number - left_bins_num,
      region == 'first' ~ Bin_Number,
      region == 'last' ~ Bin_Number,
      region == 'downstream' ~ Bin_Number + centre_bins_num)
    ) |> arrange(region, Bin_Number) -> tidy_counts
  
  # -- 11 -- Optionally check stats
  if (check_bins) {
    tidy_counts |> select(region, Bin_Number) |>  unique() |> 
      group_by(region) |> summarise(number_of_bins = n() ) |> print()
  }
  return(tidy_counts)
}


#' Title
#'
#' @param df 
#' @param left_slop A number defining how much you slopped the bed file upstream (strand aware).
#' @param right_slop A number defining how much you slopped the bed file downstream (strand aware).
#' @param left_bins_num Number of bins you used for binning the upstream region.
#' @param centre_bins_num Number of bins you used for binning the central region(s).
#' @param right_bins_num Number of bins you used for binning the downstream region.
#' @param facet_regions Logical defining whether to slightly separate (facet) the individual regions. Can be either `TRUE` or `FALSE`.
#' @param legend_position Where to put the legend. Goes to `legend.postion` of the `ggplot2::theme()`.
#' @param fill_by How to colour the histogram, one of `region`, `event`, `bin`, `count`. See details for explanation.
#' @param x_labels How to label the X-axis text, still not perfect.
#' @param event_type Text to put in the title.
#'
#' @return A ggplot
#' @import ggplot2
#' @export
#' @note save the plots with ggsave( device = cairo_pdf) and not with device = 'pdf'.
#' @details
#' The parameter `fill_by` defines how to colour the plot: `region` colours by the type of region, meaning upstream, centre or downstream region. 
#' `event` colours by the number of events (i.e. exons) with an SNP in that specific bin. It serves to check that high bars are due to many different event regions and not one single hot spot in few events.
#' `bin` colours by the bin number in each region and `count` fills based on how tall is the bar in the histogram.
#'
hallgrimskirkja_plot <- function(df, left_slop, right_slop, left_bins_num = 50, 
                                 centre_bins_num = 50, right_bins_num = 25, 
                                 facet_regions = F, legend_position = c(0.94, 0.7),
                                 fill_by = c('region', 'event', 'bin', 'count'),
                                 x_labels = c('exon', 'generic', 'bins'),
                                 event_type = NULL ){
  # -- 1 -- check input
  if( missing(df) ) { stop('Specify a dataframe created with: parse_bed_vcf_isec()')}
  # Slopping had to be applied in a strand-aware fashion.
  if (missing(left_slop)) { stop('Specify the upstream extension you applied to the bed file') }
  if (missing(right_slop)) { stop('Specify the dowstream extension you applied to the bed file') }
  
  required_cols <- c('GENE', 'EVENT', 'strand', 'region', 'Bin_Number', 
                     'Var_counts', 'Num_events_in_bin', 'Num_SNPs_per_events_per_bin',
                     'ReOrdered_Bin_Number')
  stopifnot(all(required_cols %in% colnames(df)))
  
  # make sure the 'region' column is a factor
  if ( is.factor(df$region) == FALSE ) {
    df$region <- factor(df$region, c("upstream", "centre", "downstream"))
    
    
  }
  
  # -- 2 -- Define breaks and labels for the plot. (Probably can be simplified)
  extensive_bins_breaks <- unique(sort( c( seq(-left_bins_num, 0, 5), 1, seq(0, centre_bins_num, 5), seq(centre_bins_num, centre_bins_num+right_bins_num,  5) ) ))
  reduced_bins_breaks <- unique(sort( c( -left_bins_num, 1, centre_bins_num, centre_bins_num+right_bins_num)))
  reduced_bins_labels <- c(paste0(-left_slop, 'bp'), 'start', 'end', paste0(right_slop, 'bp') )
  reduced_exons_bins_labels <- c(paste0(-left_slop, 'bp'), "3'", "5'", paste0(right_slop, 'bp') )
  
  if (x_labels == 'exon') {
    plot_breaks <- reduced_bins_breaks
    plot_labels <- reduced_exons_bins_labels
  } else if (x_labels == 'generic' ) {
    plot_breaks <- reduced_bins_breaks
    plot_labels <- reduced_bins_labels
  } else if (x_labels == 'bins' ) {
    plot_breaks <- extensive_bins_breaks
    plot_labels <- extensive_bins_breaks
  } else {
    stop('Parameter "x_labels" must be one of: "exon", "generic", "bins"')
  }
  
  # -- 3 -- Get some statistics for the plot
  Num_Events <- length(unique(df$EVENT))
  Num_Genes <- length(unique(df$GENE))
  Num_Neg_Strand <- table(unique(select(df, GENE, EVENT, strand))$strand)[1]
  Num_Pos_Strand <- table(unique(select(df, GENE, EVENT, strand))$strand)[2]
  
  # -- 4 -- Set bins filling variable
  if (fill_by == 'region') {
    fill_variable <- 'region'
    SFL <- F
  } else if (fill_by == 'event') {
    fill_variable <- 'Num_events_in_bin'
    SFL <- T
  } else if (fill_by == 'count') {
    fill_variable <- 'Var_counts'
    SFL <- F
  } else if (fill_by == 'bin') {
    fill_variable <- 'Bin_Number'
    SFL <- F
  } else {
    stop('Parameter "fill_by" must be either: "region", "event", "count", or "bin"')
  }
  
  # -- 5 -- Remove redundant info not needed for the plot
  df <- unique(select(df, region, Bin_Number, Var_counts, 
                      ReOrdered_Bin_Number, Num_events_in_bin) )
  
  # -- 6 -- Plot
  ggplot(df) +
    aes(x = ReOrdered_Bin_Number, y = Var_counts, fill = !!sym(fill_variable))  +
    geom_bar(stat = 'identity', colour = 'black', width = 0.75, 
             linewidth = 0.1, show.legend = SFL) + 
    scale_y_continuous(expand = expansion(add = c(0, 0.5), mult = c(0, 0.01)),
                       n.breaks = 6) +
    scale_x_continuous(expand = expansion(add = 0.5, mult = 0.01),
                       breaks = plot_breaks, labels = plot_labels) +
    labs(x = 'bins',  y = 'Counts') +
    coord_cartesian(clip = 'off')+
    theme_bw(base_family = 'Arial', base_size = 10) +
    theme(axis.line = element_line(linewidth = 0.2),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.2),
          panel.grid.minor.y = element_blank(),
          plot.background = element_blank(),
          plot.title = element_text(size = 10-1.5),
          axis.ticks.x = element_blank(),
          axis.ticks.length.y = unit(1, units = 'mm'),
          axis.text = element_text(colour = 'black'),
          axis.text.x = element_text(margin = margin(t = -0.85, b = -1)),
          axis.text.y = element_text(vjust = 0.25, hjust = 0.75, margin = margin(l = 0)),
          axis.title.x = element_blank()) -> Hallgrimur_plot
  
  # -- 7 -- Additional features to the plot: Fill colour representation
  if (fill_by == 'region') {
    Hallgrimur_plot <- Hallgrimur_plot +
      scale_fill_manual(values = c('lightblue', 'firebrick1', 'dodgerblue')) 
  } else if (fill_by == 'event') {
    Hallgrimur_plot <- Hallgrimur_plot +
      scale_fill_continuous(type = 'viridis', name = 'Num. exons w/\nSNP per bin') +
      guides(fill = guide_colorbar(barwidth = unit(2, units = 'mm') ) ) +
      theme(legend.position = legend_position,
            legend.justification = 'centre',
            legend.title = element_text(hjust = 0.5, margin = margin(l = -15), size = 8.5),
            legend.text = element_text(colour = 'black', margin = margin(l = -1.2) ),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.box.background = element_blank(),
            legend.text.align = 0 )
  } else if (fill_by == 'count') {
    Hallgrimur_plot <- Hallgrimur_plot +
      scale_fill_gradientn(colours = c(
        "#574571", "#90719F", "#B695BC", "#DEC5DA", 
        "#C1D1AA", "#7FA074", "#466C4B", "#2C4B27") ) 
  } else if (fill_by == 'bin') {
    Hallgrimur_plot <- Hallgrimur_plot +
      scale_fill_gradient(low = 'gray85', high = 'gray15') 
  } else {
    stop('Parameter "fill_by" must be either: "region", "event", "count",  or "bin"')
  }
  
  # -- 8 -- Additional features to the plot: title with number
  if (is.null(event_type)) {
    # show no title
  } else if (is.character(event_type) ) {
    plt_ttl <- paste0(event_type, ' n = ', Num_Events, ' in ', Num_Genes, ' genes')
    Hallgrimur_plot <- Hallgrimur_plot + labs(title = plt_ttl)
  } else {
    stop('Parameter: ', event_type, ' is either NULL or a character!')
  }
  # -- 9 -- Additional features to the plot: Faceting
  if( facet_regions == TRUE ) {
    Hallgrimur_plot <- Hallgrimur_plot + 
      facet_grid(~ region, scales = 'free_x', space = 'free_x') +
      guides(fill = guide_colorbar(barwidth = unit(2, units = 'mm') ),
             x = guide_axis(n.dodge = 2) )+
      theme(strip.background = element_blank(), 
            strip.text = element_text(vjust = 0, margin = margin(b = 0)), 
            strip.clip = 'off', 
            panel.spacing = unit(1, units = 'mm'))
  } else if( facet_regions == FALSE ) {
    Hallgrimur_plot <- Hallgrimur_plot
  } else {
    stop('Parameter: "facet_regions" must either TRUE or FALSE')
  }
  return(Hallgrimur_plot)
}
