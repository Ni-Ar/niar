
#' Helper function to check a fasta file that is going to be imported into `R`.
#'
#' @param path A character indicating a `File/path/to/a/fasta/file.fasta`
#'
#' @return Nothing. It returns an error message only if the file doesn't exist or is empty.
#' @export
#'
check_fasta <- function(path) {
  # -- 1 -- Check that fasta file actually exists
  if ( !file.exists(path) ) { stop('The fasta file does not exists!') }
  
  # -- 2 -- Check that fasta file is not empty
  fasta_info <- file.info(path)
  if ( fasta_info$size == 0) { stop('The fasta file is empty!') }
}

#' Handy wrapper to import fasta files into a data.frame
#'
#' @param path A character indicating a `File/path/to/a/fasta/file.fasta`
#' @param ID_col A character specifying a column name for the info contained in the fasta header. Default is `'fasta_header'`.
#' @param ... Additional parameters passed to `seqinr::read.fasta` like `seqtype = "DNA"`
#'
#' @return A data.frame
#' @importFrom seqinr read.fasta
#' @export
#'
#' @examples
#' df1 <- fasta2df(path = "/file/to/fasta/species.fasta", ID_col = "Species")
fasta2df <- function(path, ID_col = 'fasta_header', ...) {
  
  # -- 1 -- Check input
  check_fasta(path = path)
  
  # -- 2 -- Check that fasta file is not empty
  fasta_li <- read.fasta(file = path, forceDNAtolower = FALSE, as.string = TRUE,
                         set.attributes = TRUE, seqonly = FALSE, ...,
                         strip.desc = FALSE, whole.header = FALSE ) |> unlist()
  
  # -- 3 -- Build a data.frame with given columns.
  fasta_df <- data.frame(ID_col = names(fasta_li), Sequence = fasta_li,
                         row.names = NULL, check.names = TRUE, 
                         fix.empty.names = TRUE, stringsAsFactors = FALSE) |>
    setNames( nm = c(ID_col, 'Sequence') ) 
  
  # -- 4 -- the end.
  return(fasta_df)
}

#' Align an input fasta file with Muscle
#'
#' @param input_path A path to a fasta file you want to import and align with MUSCLE.
#' @param output_path A file path to where you want to save the aligned fasta file. 
#' If omitted a `DNAStringSet` is returned.
#' @param ... other parameters passed to `msa::msaMuscle()`.
#'
#' @return A `DNAStringSet`. If `output_path` is specified it writes to file.
#' @importFrom Biostrings writeXStringSet DNAStringSetList readDNAStringSet
#' @importFrom msa msaMuscle
#' @export
#' 
#' @description This function is a wrapper to circumvent 2 annoying issues that 
#' the \href{http://dx.doi.org/10.1093/bioinformatics/btv494}{`msa`} package has 
#' when aligning a fasta sequence, namely not being able to sort the output 
#' fasta gapped file as the input sequences and not being able to write to fasta file.
#' 
#' @details This function only aligns DNA fasta files using an `R` implementation of the muscle 
#' multiple sequences aligner. \href{https://www.drive5.com/muscle/}{MUSCLE} 
#' is claimed to achieve both better average accuracy and 
#' better speed than ClustalW2 or T-Coffee, depending on the chosen options 
#' (as stated \href{https://www.ebi.ac.uk/Tools/msa/muscle/}{here}). 
#'
#' @examples
#' # Input file
#' dna_fasta_path <- file.path(seqs_dir, "DNA/All_Seq.fasta")
#' # Output
#' dna_alg_path <- file.path(seqs_dir, "Alignments/All_Seq_Ordered_as_Input.afa")
#' 
#' align_DNA_fasta(input_path = dna_fasta_path, output_path = dna_alg_path)
align_DNA_fasta <- function(input_path, output_path, ...) {
  
  check_fasta(path = input_path)
  
  dna_fasta <- readDNAStringSet(filepath = input_path, format = 'fasta', use.names = T)
  alg_Mscl <- msaMuscle(inputSeqs = dna_fasta, type = 'dna', verbose = F, ...)
  # sort Muscle aligned sequences as in input
  input_fasta_headers <- names(dna_fasta)
  aligned_fasta_headers <- names(alg_Mscl@unmasked)
  # get the index of the aligned sequences ordered as in the input
  indx_order <- match(input_fasta_headers, aligned_fasta_headers)
  # order the muscle aligned sequences DNAStringSet object as input
  ordered_alg_Mscl <- alg_Mscl@unmasked[indx_order]
  # unlist the DNAStringSet for exporting it to fasta file
  out_algn <- unlist(DNAStringSetList(ordered_alg_Mscl))
  
  if ( missing(output_path) ) {
    return(out_algn)
  } 
  
  if ( !file.exists(output_path) ) {
    writeXStringSet(x = out_algn, filepath = output_path, format = "fasta", width = 200)
  } else if (file.exists(output_path)){
    stop('Output file already exists!\n', 
         'To remove the existing file run:\n',
         'file.remove("', output_path, '")')
  } else {
    stop("Something is wrong with the outout_path... please provide a valid file path")
  }
}

#' Align a fasta file with MUSCLE and calculate the percentage of sequence identity of the aligned sequence
#'
#' @param input_path A path to a DNA fasta file for which you want the pairwise sequence identity
#' @param percentage_identity Logical, return the percentage of sequence identity. 
#' Default `TRUE`. If `FALSE` the squared root of the pairwise identity is returned. See details.
#' @param ... Other parameters passed to `msa::msaMuscle()`.
#'
#' @return A matrix with aligned DNA sequence identity
#' @importFrom Biostrings readDNAStringSet
#' @importFrom msa msaMuscle msaConvert
#' @importFrom seqinr dist.alignment
#' @export
#' 
#' @description This function first aligns a DNA fasta file with multiple sequences 
#' using MUSCLE using `msa::msaMuscle()`, where you can pass additional parameters thanks to `...`. 
#' Then the pairwise sequence indentity is calculated using `seqinr::dist.alignment()`. 
#' Gaps in the alignment will be counted in the identity measure. The output 
#' matrix is always ordered as the input fasta sequences.
#' 
#' @details By default this function returns the percentage of sequence identity 
#' from 0 to 100. By setting `percentage_identity = FALSE`, the sqrt(1 - identity) 
#' is returned. So if the identity between 2 sequences is 19% the squared root 
#' of (1.0 - 0.19) i.e. 0.9.
#'
#' @note If you want to inspect the multiple sequence alignment used to calculate 
#' the percentage of sequence identity you can write it to a fasta file 
#' with `align_DNA_fasta()` using the same MUSCLE parameters.
#'
#' @examples
#' pim <- fasta2pim(input_path = dna_fasta_path)
fasta2pim <- function(input_path, percentage_identity = TRUE, ...) {
  
  check_fasta(path = input_path)
  # input fasta
  dna_fasta <- readDNAStringSet(filepath = input_path, format = 'fasta', use.names = TRUE)
  # align with MUSCLE
  alg_Mscl <- msaMuscle(inputSeqs = dna_fasta, type = 'dna', verbose = FALSE, ...)
  # convert the msa MUSCLE to an object that seqinr can recognise.
  alg_seqinr <- msaConvert(alg_Mscl, type = "seqinr::alignment")
  
  # Calculate the percentage of identity matrix
  d <- dist.alignment(alg_seqinr, matrix = "identity", gap = TRUE)
  d <- as.matrix(d)
  
  # order distance matrix as as in input
  stopifnot(all(rownames(d) == colnames(d)))
  input_fasta_headers <- names(dna_fasta)
  mat_names <- rownames(d)
  # get the index of the aligned sequences ordered as in the input
  indx_order <- match(input_fasta_headers, mat_names)
  
  d <- d[indx_order, indx_order]
  
  # by default dist.alignment returns sqrt(1-dist).
  if (percentage_identity == TRUE) {
    # convert it back to the percentage of sequence identity
    pim <- (-(d ^ 2) + 1)*100
    stopifnot(all(dim(d) == dim(pim)))
    return(pim)
  } else if (percentage_identity == FALSE) {
    return(d)
  } else {
    stop("The parameter 'percentage_identity' must be either TRUE or FALSE!")
  }
}