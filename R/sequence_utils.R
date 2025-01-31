#' Generate Single Nucleotide Variants
#'
#' @title Generate all possible single nucleotide substitutions and deletions for a given sequence
#'
#' @description
#' This function takes a nucleotide sequence and generates all possible single nucleotide
#' substitutions and deletions. It can work with both DNA and RNA sequences, depending on
#' the provided alphabet. The function can either return the variants as a data frame or
#' save them to a VCF (Variant Call Format) file.
#'
#' @param sequence A character string representing the nucleotide sequence to be mutated.
#' @param alphabet A character vector specifying the nucleotide alphabet to use.
#'   Default is c('A', 'C', 'G', 'T') for DNA. Use c('A', 'C', 'G', 'U') for RNA.
#' @param save_to_vcf A logical value indicating whether to save the results to a VCF file.
#'   Default is FALSE.
#' @param out_basename A character string specifying the base name for the output VCF file.
#'   Only used if save_to_vcf is TRUE. Default is 'Temp'.
#' @param out_path A character string specifying the directory path for the output VCF file.
#'   Only used if save_to_vcf is TRUE. Default is tempdir().
#'
#' @details
#' The function generates all possible single nucleotide substitutions by replacing each
#' nucleotide in the sequence with every other nucleotide in the provided alphabet. It also
#' generates deletion variants by removing each nucleotide one at a time. The resulting
#' variants are represented in a data frame or VCF file with the following columns:
#' CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and SEQ.
#'
#' @return
#' If save_to_vcf is FALSE, the function returns a data frame containing all generated
#' variants. If save_to_vcf is TRUE, the function saves the variants to a VCF file and
#' returns nothing, but prints a message indicating where the file was saved.
#' @importFrom utils write.table
#' 
#' @export
#'
#' @examples
#' # Generate variants for a DNA sequence
#' seq <- "ATGCATCA"
#' variants <- generate_single_nt_variants(seq)
#' head(variants)
#'
#' # Generate variants for an RNA sequence
#' rna_seq <- "AUGCAUCA"
#' rna_variants <- generate_single_nt_variants(rna_seq, alphabet = c('A', 'C', 'G', 'U'))
#' head(rna_variants)
#'
#' # Save variants to a VCF file
#' generate_single_nt_variants(seq, save_to_vcf = TRUE, out_basename = "my_variants")
#'
generate_single_nt_variants <- function(sequence, alphabet = c('A', 'C', 'G', 'T'), 
                                        save_to_vcf = FALSE, out_basename = 'Temp',
                                        out_path = tempdir()) {
  # Define the possible DNA/RNA nucleotides
  nucleotides <- alphabet
  
  # Initialize a list to store the mutations
  mutations <- list()
  
  # Iterate through each position in the sequence
  for (i in 1:nchar(sequence)) {
    original_base <- substr(sequence, i, i)
    
    # Generate substitution mutations
    for (new_base in setdiff(nucleotides, original_base)) {
      mutated_seq <- paste0(
        substr(sequence, 1, i-1),
        new_base,
        substr(sequence, i+1, nchar(sequence))
      )
      mutations <- c(mutations, list(c(CHROM = "chrN", 
                                       POS = i, 
                                       ID = "TBA", 
                                       REF = original_base, 
                                       ALT = new_base, 
                                       QUAL = ".", 
                                       FILTER = "PASS", 
                                       INFO = "TYPE=SNP",
                                       SEQ = mutated_seq)))
    }
    
    # Generate deletion mutation
    mutated_seq <- paste0(
      substr(sequence, 1, i-1),
      "-",
      substr(sequence, i+1, nchar(sequence))
    )
    mutations <- c(mutations, list(c(CHROM = "chrN", 
                                     POS = i, 
                                     ID = "TBA", 
                                     REF = original_base, 
                                     ALT = "-", 
                                     QUAL = ".", 
                                     FILTER = "PASS", 
                                     INFO = "TYPE=DEL",
                                     SEQ = mutated_seq)))
  }
  
  # Convert the list to a data frame
  result <- do.call(rbind, mutations)
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  
  # save to file
  if (save_to_vcf == TRUE) {
    # Add VCF header
    header <- c("##fileformat=VCFv4.5",
                "##source=niar::generate_single_nt_variants",
                paste0("##fileDate=", format(Sys.Date(), "%Y%m%d") ),
                paste0("##reference=", sequence),
                "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Type of variant\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tSEQ")
    
    # Write the VCF file
    output_file <- file.path(out_path, paste0(out_basename, ".vcf"))
    writeLines(header, output_file)
    write.table(result, output_file, append = TRUE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    message("VCF file saved at\n", output_file)
  } else if (save_to_vcf == FALSE) {
    return(result)
  } else {
    stop('The parameter "save_to_vcf" must be either TRUE or FALSE !!')
  }
}


# get_kmers(sequence = 'ACTGCAGCATCAGACGCA', k = 25)
# Helper function to generate k-mers from a sequence
get_kmers <- function(sequence, k = 2) {
  n <- nchar(sequence)
  if (k > n) {
    stop('The k-mer length (', k,') is longer than the sequence length (', n,')!')
  } else if (k <= n) {
    kmers <- sapply(1:(n - k + 1), function(i) {
      substr(sequence, i, i + k - 1)
    })
  } else {
    stop("I can't understand the k-mer length and sequence length...?")
  }
  return( unique(kmers) )
}

# Function to calculate the Jaccard distance between two sequences
# seq1, seq2: character strings representing the sequences
# k: length of k-mers to consider (default is 2)
# calculate_jaccard(seq1 = 'ATGCATATGCAT', seq2 = 'CTGCATATCGAACAT', k = 6, return = 'index')
calculate_jaccard <- function(seq1 = '', seq2 = '', k = 3, return = c('index', 'distance')) {
  if ( missing(return) ) {
    stop('Please use "return = index" or "return = distance"')
  }
  # Generate k-mers for each sequence
  kmers1 <- get_kmers(seq1, k)
  kmers2 <- get_kmers(seq2, k)
  
  # Calculate the intersection and union of the k-mers
  intersection <- intersect(kmers1, kmers2)
  union <- union(kmers1, kmers2)
  
  # Handle the case where the union is empty
  if (length(union) == 0) {
    jaccard_index <- 1  # If both sequences are empty or k > length of sequences
  } else {
    # Compute the Jaccard index
    # Measures the similarity between two sets (A and B) 
    jaccard_index <- length(intersection) / length(union)
  }
  
  # Compute the Jaccard distance
  # Measures the dissimilarity between two sets
  jaccard_distance <- 1 - jaccard_index
  
  if (return == 'index') {
    return(jaccard_index)
    
  } else if (return == 'distance') {
    return(jaccard_distance)  
  } else {
    stop('The argument return must be either "index" or "distance, not: ', return)
  }
}

