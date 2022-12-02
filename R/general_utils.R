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
