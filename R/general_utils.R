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
  
  # Check that all input letters are in the ordered alphabet
  input_letters <- unique(unlist(strsplit(x = string, split = '')))
  if (missing(alphabet_order)) { alphabet_order <- sort(input_letters) }
  
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

#' My local function to enumerate the Permutations of the Elements of a Vector
#' 
#' @description
#' This is a copy of the `gtools::permutations()` function. So I have one less dependency to handle since I only need this great function from that package.
#' 
#' @param n Size of the source vector
#' @param r Size of the target vectors
#' @param v Source vector. Defaults to 1:n
#' @param set Logical flag indicating whether duplicates should be removed from the source vector v. Defaults to TRUE.
#' @param repeats.allowed Logical flag indicating whether the constructed vectors may include duplicated values. Defaults to FALSE.
#'
#' @return Returns a matrix where each row contains a vector of length r.
#' @export
#' 
#' @details
#' Caution: The number of combinations and permutations increases rapidly with n and r!.
#' To use values of n above about 45, you will need to increase R's recursion limit. See the expression argument to the options command for details on how to do this.
#' Taken from an email by Brian D Ripley to r-help dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to Alex Ahgarin. Original version was named "subsets" and was Written by Bill Venables.
#'  
#' @author See the authors in `gtools::permutations()`
#'
#' @examples
#' # permutationz(n = length(c('A', 'T', 'C', 'G')), r = 6, v = c('A', 'T', 'C', 'G'), repeats.allowed = TRUE)
permutationz <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) != 
      0) 
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) != 
      0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                                                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                                                  byrow = TRUE))
      }
    }
  else sub <- function(n, r, v) {
    if (r == 1) 
      matrix(v, n, 1)
    else if (n == 1) 
      matrix(v, 1, r)
    else {
      X <- NULL
      for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                        1, r - 1, v[-i])))
      X
    }
  }
  sub(n, r, v[1:n])
}

#' Create same length sequence permutations for each letter 
#'
#' @description
#' Given a set of letters in an alphabet return a dataframe with all the possible permutations of those letters.#' 
#'
#' @param sequence_length How long the sequences should be (keep this below 40)
#' @param alphabet Letters to permutate in order to build different sequences
#' @param verbose Show a pretty print summary data frame. Defatult `FALSE`.
#' @param k Pretty print option: show the first k rows per letter of the alphabet
#' @param scramble Randomize the order of the rows.
#' @param seed Integer passed to `set.seed()` just before data is scrambled to get reproducible results. Default is `NULL`.
#' 
#' @return A data frame
#' @importFrom gtools permutations
#' @export
#' 
#' @details
#' This is a handy wrapper of `gtools::permutation()` function to generate saturated permutation sequences of DNA.
#' The number of permutations increases rapidly with the length of the `alphabet` and `sequence_length`!, I made a pretty print function
#' that can be used to slice the permutation sequences and check the results.
#' To show the pretty print data frame use `verbose = TRUE`.
#' 
#' By default the sequences are returned in alphabetical order. With `scramble = TRUE` one can reshuffle the sequences in a random order.
#' 
#' @seealso [print_seq_perm]
#' 
#' @examples
#' dat <- permutate_seq(sequence_length = 5, k = 3, verbose = T)
#' dat2 <- permutate_seq(sequence_length = 3, alphabet = c("W", "*", "X", "!", "%", "7"), k = 7, verbose = T)
#' 
#' # To make many sequencing barcodes that follow this patter: 
#' # 'NNNN', 'AGCT', 'NNNN', 'TCAG', 'NNNN', 'TAGC', 'NNN', 'CAGT', 'NNN'
#' 
#' barcodes <- list()
#' for (i in 1:100) {
#'   tmp <- cbind( permutate_seq(sequence_length = 4, scramble = T), 'AGCT', 
#'                 permutate_seq(sequence_length = 4, scramble = T), 'TCAG', 
#'                 permutate_seq(sequence_length = 4, scramble = T), 'TAGC', 
#'                 permutate_seq(sequence_length = 3, scramble = T), 'CAGT', 
#'                 permutate_seq(sequence_length = 3, scramble = T) ) 
#'   # concatenate into one single sequence
#'   tmp <- apply(tmp, 1, paste0, collapse = "") |> data.frame() |> setNames('BC')
#'   
#'   barcodes[[i]] <- tmp
#' }
#' 
#' do.call('rbind', barcodes) |> unique() |> nrow()
#' # 25600 ( (N^4) * 100, where N = 4 )
#' 
permutate_seq <- function(sequence_length, alphabet = c('A', 'C', 'G', 'T'), k = 5, verbose = FALSE, 
                          scramble = FALSE, seed = NULL) {
  stopifnot(sequence_length >= 1)
  alphabet <- sort(unique(alphabet))
  
  # this is like gtools::permutations
  permutationz(n = length(alphabet), r = sequence_length, v = alphabet, repeats.allowed = TRUE) |>
    apply( 1, paste0, collapse = "") |> data.frame() -> dat
  colnames(dat) <- 'Sequence'
  
  # scramble order
  if (scramble == TRUE) {
    set.seed(seed) # for deterministic results if different from NULL
    dat$Sequence <- dat$Sequence[sample(nrow(dat))]
  }
  
  if (verbose == TRUE) {
    print(print_seq_perm(perm_dat = dat, alphabet = alphabet, k = k) )
  }
  invisible(dat)
}

#' Pretty print a sequences permutation data frame
#' 
#' @description
#' Show the first k rows of all the words starting with each letter of the `alphabet`.
#'
#' @param perm_dat the data frame generated with `permutate_seq()`
#' @param alphabet Letters to permutate in order to build different sequences
#' @param k Pretty print option: show k rows per letter of the alphabet
#'
#' @return A data frame slicing the permutation data by letter
#' @export
#'
#' @details
#' In addition to the first k words of each letter also the last k rows are always returned.
#'
#' @seealso [permutate_seq]
print_seq_perm <- function(perm_dat, alphabet = c('A', 'C', 'G', 'T'), k = 6) {
  # Get first word (in alphabetical order) for each letter in the alphabet
  lapply(alphabet, function(x) {
    first_letter_regex <- paste0('^', x)
    head( which( grepl(pattern = first_letter_regex, x = perm_dat$Sequence) ), k)  
  }) -> indx
  
  # add last rows as well 
  srtd_indx <- sort(c(unlist(indx), c(
    ( nrow(perm_dat) - (k - 1)  ):nrow(perm_dat) ) )
  )
  
  # divide indexes numbers in a list with sub-vectors of length k elements
  srtd_indx_li <- split(srtd_indx, ceiling(seq_along(srtd_indx) / k))
  
  # get sequences for the selected indexes
  sapply(srtd_indx_li, FUN = function(x){
    perm_dat[unlist(x), ]
  }, simplify = T ) |> as.data.frame() -> sampled_seq
  # add dot dot dot at the end of each group
  sampled_seq[nrow(sampled_seq) + 1, ] <- '...'
  
  # turn into vectors
  seq <- unlist(as.vector(sampled_seq))
  
  lapply(srtd_indx_li, function(x) {
    as.character(c(x, " ..."))
  }) |> unlist() -> indx_dot
  
  # sanity check
  length(seq) == length(indx_dot)
  
  # make a pretty print data frame
  pp_dat <- data.frame(Number = indx_dot, Sequence = seq, row.names = NULL)
  # remove last row
  pp_dat <- pp_dat[ 1:(nrow(pp_dat) - 1), ]
  return(pp_dat)
}

