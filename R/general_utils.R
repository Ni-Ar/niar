#' Find longest and most abundant common suffix
#'
#' @param x A character vector
#' @param k An integer between 1 and 99. How much more abundant should the first suffix be relative to the second most abundant suffix?
#' @param max_suffixes An integer. If there are more than `max_suffixes` in `x` don't check for relative percentage abundance.
#' @param ignore.case Logical. Convert all `x` to upper case.
#' @param verbose Logical. Print out info.
#'
#' @return A character corresponding to the longest common suffix.
#' @export
#'
#' @examples
#' sample_names <- c("A_sample1", "B_sample1", "C_sample1", "D_sample2", "E_sample2")
#' longest_most_abundant_common_suffix(x = sample_names)
#' # Returns: "_sample1"
longest_most_abundant_common_suffix <- function(x, k = 80, max_suffixes = 4, 
                                                ignore.case = FALSE, 
                                                verbose = FALSE) {
    # k must be between 1 and 99
    stopifnot(k > 0 | k < 100)
    
    x <- as.character(x)
    if ( ignore.case ) {  x <- toupper(x) }
    
    nc <- nchar(x, type = "char")
    
    i <- 1
    
    suffixes <- c()
    for ( i in 1:min( nc ) ) {
        ss <- substr(x, start = nchar(x) - i + 1, stop = nchar(x) )
        
        # If there's more than one suffix
        if ( length(unique(ss)) >= 2 ) {
            # message(i, ") Found more than one suffix. Most common: ", names( which.max( table(ss) ) ) )
            
            # If there are less than max_suffixes common suffix check the frequencies
            if ( length(unique(ss)) <= max_suffixes ) {
                
                
                if (  ((table(ss)[1] * k) / 100) > table(ss)[2] ) {
                    if (verbose) {
                        message("Suffix 1 ('", names(table(ss)[1]), "') occurs at least ",
                                k, "% more abundant than suffix 2('", names(table(ss)[2]), "')" )
                    }
                    
                } else{
                    if (verbose) {
                        message("Suffix 1 ('", names(table(ss)[1]), "') occurs at least ",
                                k, "% more abundant than suffix 2('", names(table(ss)[2]), "')" )
                    }
                    # skip this longest suffix
                    next
                }
            } else {
                # if there are more than max_suffixes don't even bother checking.
                next
            }
        }
        
        suffixes[i] <- names(table(ss)[1])
    }
    
    return( suffixes[which.max(nchar(suffixes, type = "char"))] )
}

#' Find the longest common prefix between 2 random names. If no common prefix is found the same input values are returned.
#'
#' @param names A character vector with length higher than 1.
#' @param uniquify Logical, should `names` become a unique vector before sampling 2 random `names`? Default `FALSE`, setting to `TRUE`, could help finding a better longest common prefix.
#' @param verbose Logical, print common prefix found.
#'
#' @return A character vector without common prefix that is the same length as the input.
#' @importFrom Biostrings lcprefix
#' @export
#'
#' @examples
#' new_names <- longest_common_prefix(names)
longest_common_prefix <- function(names, uniquify = TRUE, verbose = TRUE) {
    
    if (length(names) <= 1) { stop("I need 2 or more names!")}
    
    # Try to remove the longest common prefix from the sample names
    # This picks 2 random elements in the names and checks what prefix they have 
    # in common and removes it from names.
    sampled_names <- names
    if (uniquify) { sampled_names <- unique(sampled_names) }
    sampled_names <- sort(sampled_names)
    
    set.seed(16)
    random_names <- base::sample(x = sampled_names, size = 2)
    # PREFIX
    nchar_common_prefix <- lcprefix(s1 = random_names[1],
                                    s2 = random_names[2])
    
    if ( nchar_common_prefix > 0 ) {
        common_prefix <- substr(random_names[1], start = 1, stop = nchar_common_prefix)
        if (verbose) { message("Found prefix: ", common_prefix) }
        new_names <- gsub(pattern = paste0("^", common_prefix),
                          replacement = "", x = names, 
                          ignore.case = T, perl = F)
    } else {
        if (verbose) { message("Couldn't find common prefix.") }
        new_names <- names
    }
    return(new_names)
}

## Write a function that combines both removal of longest prefix and suffix

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
  
  if (!file.exists(path) ) { stop('File path does not exist!') }
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

#' Perform one-hot encoding of any character string
#'
#' @param string A string of letters to be one-hot encoded. Case-sensitive!
#' @param alphabet_order A one-letter character vector specifying the alphabet order. 
#' It must contain every letter present in `string`. 
#' If missing the alphabetical order will be used. Default is DNA alphabet: `c('A', 'T', 'C', 'G')`.
#'
#' @return a data.frame
#' @export
#' @description
#' Given any sequence return a data.frame where every column corresponds to one
#' letter of the input `string` and each row corresponds to a letter as set in the
#' `alphabet_order`. Each cell in the data.frame will be a zero, only when the column 
#' equals the row it will be a `1`. 
#' For example, in the case of DNA each nucleotide will be represented as a vector of
#'  length 4, where 3 positions are 0 and only one position is 1, depending on the nucleotide.
#' 
#' @details
#' If `alphabet_order` is not specified the function will create one on the fly 
#' sorting the input `string` by alphabetical order. If a letter in the `string` 
#' character vector is not present in the `alphabet_order` the function will 
#' return an error. Special characters and numbers are allowed (see examples).
#'
#' @author Niccolò Arecco
#' @keywords onehot
#' 
#' @examples
#' # RNA example
#' encode_onehot(string = "UUUAAACCCGG", alphabet_order = c('A', 'U', 'C', 'G'))
#' # Returns the following data.frame where the rows are ordered as in the 
#' # alphabet_order and the columns are ordered as the input string.
#'   U U U A A A C C C G G G
#' A 0 0 0 1 1 1 0 0 0 0 0 0
#' U 1 1 1 0 0 0 0 0 0 0 0 0
#' G 0 0 0 0 0 0 0 0 0 1 1 1
#' C 0 0 0 0 0 0 1 1 1 0 0 0
#' 
#' # Alphabet is optional 
#' encode_onehot(string = 'ATHCAY')
#' # Returns the following data.frame where the input string was sorted 
#' # alphabetically to generate the order on the rows
#'   A T H C A Y
#' A 1 0 0 0 1 0
#' C 0 0 0 1 0 0
#' H 0 0 1 0 0 0
#' T 0 1 0 0 0 0
#' Y 0 0 0 0 0 1
#' 
#' # Case sensitive input
#' encode_onehot(string = 'acACCnN')
#' # Returns the following data.frame where lower case letters appear before upper case one
#'   a c A C C n N
#' a 1 0 0 0 0 0 0
#' A 0 0 1 0 0 0 0
#' c 0 1 0 0 0 0 0
#' C 0 0 0 1 1 0 0
#' n 0 0 0 0 0 1 0
#' N 0 0 0 0 0 0 1
#' 
#' # Special characters and numbers are encoded just fine
#' encode_onehot(string = 'MaQ8T!S-K C2C*')
#' # Returns a data.frame where symbols are sorted first. 
#' # Note how the space (' ') is both a row name and column name
#'   M a Q 8 T ! S - K   C 2 C *
#'   0 0 0 0 0 0 0 0 0 1 0 0 0 0
#' - 0 0 0 0 0 0 0 1 0 0 0 0 0 0
#' ! 0 0 0 0 0 1 0 0 0 0 0 0 0 0
#' * 0 0 0 0 0 0 0 0 0 0 0 0 0 1
#' 2 0 0 0 0 0 0 0 0 0 0 0 1 0 0
#' 8 0 0 0 1 0 0 0 0 0 0 0 0 0 0
#' a 0 1 0 0 0 0 0 0 0 0 0 0 0 0
#' C 0 0 0 0 0 0 0 0 0 0 1 0 1 0
#' K 0 0 0 0 0 0 0 0 1 0 0 0 0 0
#' M 1 0 0 0 0 0 0 0 0 0 0 0 0 0
#' Q 0 0 1 0 0 0 0 0 0 0 0 0 0 0
#' S 0 0 0 0 0 0 1 0 0 0 0 0 0 0
#' T 0 0 0 0 1 0 0 0 0 0 0 0 0 0
#' 
encode_onehot <- function(string, alphabet_order = c('A', 'T', 'C', 'G') ) {
  # 1 - Check input ------------
  if (missing(string)) { stop('Input string is missing!') }
  if (missing(alphabet_order)) { alphabet_order <- sort(input_letters) }
  
  # Check that all input letters are in the ordered alphabet
  input_letters <- unique(unlist(strsplit(x = string, split = '')))
  if( !all(input_letters %in% alphabet_order) ) {
    extra_letter <- input_letters[which(!input_letters %in% alphabet_order)]
    stop('Input letter ', extra_letter, ' is not part of the alphabet')
  }
  
  # 2 - Construct a one hot decoder based on the alphabet order ------------
  num_letters <- length(unique(alphabet_order))
  decoder <- data.frame(letter = alphabet_order, onehot = c(1:num_letters))
  
  # 3 - Create matrix of zeros: cols = string x rows = alphabet ------------
  mat <- matrix( rep(0, num_letters * nchar(string) ), 
                nrow = num_letters, ncol = nchar(string), byrow = T )
  colnames(mat) <- unlist(strsplit(x = string, split = '') )
  rownames(mat) <- alphabet_order
  df <- as.data.frame(mat) # coerce matrix to data.frame
  
  # 4 - Start a double for loop to encode ------------
  # For every column in data.frame (that is a letter in the sequence)
  for (c in 1:ncol(df) ) { 
    # For every letter in the decoder (that is a letter in the alphabet)
    for (l in 1:num_letters) { 
      # If column and row letters are equal
      if ( colnames(df)[c] ==  decoder[l, ]$letter) {
        # Create an empty string of zeros of the same length of the alphabet
        string <- rep(0, num_letters)
        # Convert to 1 only the position that matches
        string[decoder[l, ]$onehot] <- 1
        df[, c] <- string
      }
    }
  }
  return(df)
}
