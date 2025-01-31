#' Helper function to check a file that is going to be imported into `R`.
#'
#' @param path A character indicating a `File/path/to/a/file`
#'
#' @return Nothing. It returns an error message only if the file doesn't exist or is empty.
#' @export
#'
check_file <- function(path) {
  # -- 1 -- Check that file actually exists
  if ( !file.exists(path) ) { stop('The file does not exists!') }
  
  # -- 2 -- Check that file is not empty
  fasta_info <- file.info(path)
  if ( fasta_info$size == 0) { stop('The file is empty!') }
}

## ------ BED files functions ------

#' Import a bed file into R as a tibble
#'
#' @param path Path to bed file in your system specified as character
#' @param header Is the bed file first line specifying the column names?
#'
#' @return A tibble
#' @importFrom readr read_delim
#' @export
#' 
#' @details
#' If header = `F` than the UCSC official bed files column names are added.
#'
#' @examples
#' read_bed(path = 'path/to/bed/file', header = F)
read_bed <- function(path, header) {
  
  check_file(path)
  if (missing(header) ) { stop('Specify if input bed file contains a header!') }
  
  # Official columns names: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
  bed_cols <- c('chr', 'start', 'end', 'name', 'score', 'strand', 
                'thickStart', 'thickEnd', 'itemRGB', 
                'blockCount', 'blockSize', 'blockStarts')
  
  bed_df <- read_delim(file = path, delim = '\t', quote = '', trim_ws = T, 
                       progress = F, show_col_types = FALSE, col_names = header)
  
  # add official column names if missing
  if (header == FALSE) {
    num_cols_in <- ncol(bed_df)
    colnames(bed_df) <- bed_cols[1:num_cols_in]
  }
  return(bed_df)
}

#' Extend upstream or downstream the coordinates of a bed file
#'
#' @param bed a tibble imported with `read_bed()`
#' @param upstream How many nucleotides to add upstream to the start coordinate
#' @param downstream How many nucleotides to add upstream to the end coordinate
#' @param strand_aware Logical, whether or not to slop in a strand aware fashion.
#'
#' @return A tibble
#' @export
#'
#' @description
#' If `strand_aware` then for the negative strands start position is slopped by 
#' the negative `downstream` value and
#' the end position is slopped by adding the upstream value.
#' 
#' @details
#' strand can only be specified with `+` or `-`,
#' 
#'
#' @examples
#' read_bed(path = test_coord, header = F)
#' 
#' # A tibble: 3 × 6
#' # chr   start   end name    score strand
#' # <chr> <dbl> <dbl> <chr>   <dbl> <chr> 
#' # chr10   100   200 test1     0 +     
#' # chr14   500   550 test2     0 -     
#' # chr1    150   190 test3     0 -   
#' 
#' read_bed(path = exons_coord, header = F) |> slop_bed(upstream = 10, downstream = 5, strand_aware = T)
#' 
#' # A tibble: 3 × 6
#' # chr   start   end name    score strand
#' # <chr> <dbl> <dbl> <chr>   <dbl> <chr> 
#' # chr10    90   205 test1     0 +     
#' # chr14   495   560 test2     0 -     
#' # chr1    145   200 test3     0 -     
#' 
slop_bed <- function(bed, upstream = 10, downstream = 5, strand_aware = T) {
  
  required_cols <- c('start', 'end')
  
  if (strand_aware == T) { required_cols <- c(required_cols, 'strand') }
  
  if (!any(grepl(pattern = '\\+|\\-', x = unique(bed$strand), perl = T)) ) {
    stop('Strand must be specified as "+" or "-", not as: ', table(bed$strand) )
  }
  
  if (strand_aware == TRUE ) {
    slopped_bed <- bed |>
      mutate(start = ifelse(strand == '+', start - upstream, start - downstream),
             end = ifelse(strand == '+', end + downstream, end + upstream + 1) ) # add + 1 for negative strand.
  } else if (strand_aware == FALSE ) {
    stop('Strand unaware slopping is not yet supported!')
  } else {
    stop('strand_aware must be either TRUE or FALSE')
  }
  
  return(slopped_bed)
}

## ------ GTF files functions ------

#' Helper function to handle the last column of a GTF file.
#'
#' @param data A data frame of a GTF file with the 9th column labelled 'attribute' 
#'
#' @return A tibble
#' @importFrom stringr str_split str_extract
#' @export
#' 
#' @details
#' This function helps convert the "attributes" column of a GTF file from character to a named list.
#' 
#' @examples
#' gtf <- read_delim(file = path, delim = '\t', comment = '#!',
#'                   col_names = c('seqname', 'source', 'feature', 'start', 
#'                   'end', 'score', 'strand', 'frame', 'attribute'), 
#'                   col_types = c('i', 'c', 'c', 'i', 'i', 'c', 'c', 'c', 'c'),
#'                   n_max = 100) |>
#'                   gtf_attributes_as_list()
gtf_attributes_as_list <- function(data) {
  
  list_attributes <- list()
  for (n in 1:nrow(data)) {
    # get key/tag/names of the attributes
    str_split(string = data$attribute[n], pattern = "; ", simplify = T ) |>
      str_extract( pattern = "[a-z].+(?= \")") -> attributes_names
    
    # get values of the attributes
    str_split(string = data$attribute[n], pattern = "; ", simplify = T ) |>
      str_extract( pattern = "(?<=\")[aA0-zZ9].*(?=\")") -> attributes_values
    
    stopifnot(length(attributes_values) == length(attributes_names))
    # turn attributes into a named character vector with unique names
    names(attributes_values) <- make.unique(attributes_names, sep = "_")
    # populate an empty list
    list_attributes[[n]] <- attributes_values
  }
  stopifnot(nrow(data) == length(list_attributes))
  # create a new 'attributeS' column
  data$attributes <- NULL
  data$attributes <- list_attributes
  # remove imported 'attribute' column
  data$attribute <- NULL
  return(data)
}


#' Helper function to convert the last column list of a GTF file to a different separate columns
#'
#' @param data A data frame of a GTF file with the 9th column labelled 'attribute' 
#'
#' @return A tibble
#' @importFrom dplyr bind_rows mutate
#' @importFrom tibble as_tibble
#' @export
#'
#' @details
#' This function works only if the "attributes" column was first converted to a list with `gtf_attributes_as_list()`
#' 
#' @examples
#' gtf <- read_delim(file = path, delim = '\t', comment = '#!',
#'                   col_names = c('seqname', 'source', 'feature', 'start', 
#'                   'end', 'score', 'strand', 'frame', 'attribute'), 
#'                   col_types = c('i', 'c', 'c', 'i', 'i', 'c', 'c', 'c', 'c'),
#'                   n_max = 100) |>
#'                   gtf_attributes_as_list() |>
#'                   gtf_attributes_as_cols()
gtf_attributes_as_cols <- function(data) {
  data_attr <- bind_rows(data$attributes) |>
    relocate(starts_with('tag'), .after = last_col())
  stopifnot( nrow(data) == nrow(data_attr) )
  data$attributes <- NULL
  data_attr$attributes <- NULL
  # coerce coordinates to integer
  data <- cbind(data, data_attr) |> mutate(across(c(start, end), as.integer)) |> as_tibble()
  return(data)
}

#' Import a GTF as tibble
#'
#' @param A valid path to a GTF file. 
#' @param max_rows The maximum rows to import. Helps avoiding memory issues. Default 10000.
#' @param col_attribute_type Specify how to handle the last column of a GTF file. See details.
#' #' \itemize{
#' \item{\code{character} - Default. Basically do nothing. }
#' \item{\code{list} - Convert the column to a list }
#' \item{\code{split} - Split the column into different separate columns}
#' }  
#' @param ... Additional parameters to pass to `read_delim()`
#'
#' @return A GTF file as a tibble
#' @importFrom readr read_delim
#' @export
#' 
#' @description
#' Import a genomics GTF file and offer different way to handle the last 9th column (with a bunch of extra info).
#' 
#' @details
#' The parameter `col_attribute_type` can be either: "character", "list", or "split". 
#' The first option (`character`) imports the data as is.
#' The option `list` coerce the column into a list column in the data frame.
#' The last option `split` divides all the attributes into separate columns with a header.
#'
#' @examples
#' hg38gtf <- read_gtf(path = hg38_gtf_path)
read_gtf <- function(path, max_rows = 1e4, col_attribute_type = 'character', ...) {
  check_file(path)
  gtf <- read_delim(file = path, delim = '\t', comment = '#!',
                    col_names = c('seqname', 'source', 'feature', 'start', 
                                  'end', 'score', 'strand', 'frame', 'attribute'), 
                    col_types = c('c', 'c', 'c', 'i', 'i', 'c', 'c', 'c', 'c'),
                    n_max = max_rows) 
  if ( col_attribute_type == 'character' ) {
    # do nothing
  } else if ( col_attribute_type == 'list') {
    gtf <- gtf |> gtf_attributes_as_list()
  } else if ( col_attribute_type == 'split' ) {
    gtf <- gtf |> gtf_attributes_as_list() |> gtf_attributes_as_cols()
  } else {
    stop('The parameter "col_attribute_type" must be one of:\n', 
         ' - character\n', ' - list\n', ' - split')
  }
  return(gtf)
}

