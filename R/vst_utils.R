#' Grep the PSI values of an as event from a vast-tools output table. 
#'
#' @param inclusion_tbl path to vast-tools inclusion table that contains a vst_id event
#' @param vst_id vast-tools alternative splicing event to grep in the inclusion_tbl
#' @param tmp_dir path to a temporary location folder
#' @param clear_tmp Logical, remove the temporary files after the function is done?
#' @param fast_grep Logical, if TRUE use 'grep -m1' to stop the search after first match found. If FALSE all vst_expression_tbl is read.
#' @param split_temp Logical, if TRUE first create a temp header file then a temp AS event file and merge them. If FALSE create one temp file with both the header and the AS event PSI and quality scores.
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' #' grep_psi(path_to_vst_tbl, vst_id) |>
#'     tidy_vst_psi()
grep_psi <- function(inclusion_tbl, vst_id, tmp_dir = tempdir(), verbose = FALSE, 
                     clear_tmp = TRUE, fast_grep = TRUE, split_temp = TRUE) {
  
  # 0 -- Check params and set commands
  if (!file.exists(inclusion_tbl)) { 
    stop("Can't find an inclusion table at:\n", inclusion_tbl )
  }
  
  if (!dir.exists(tmp_dir)) { 
    if ( verbose ) { message("No temporary directory found, creating one." ) }
    dir.create(tmp_dir)
  }
  
  # the -m1 instructs grep to stop reading a file after 1 matching line. 
  # It is done to speed up the process
  # The number of searches is calibrated to the num of events to grep
  num_events <- length(vst_id)
  if (verbose) { message("Looking for ", num_events, " events to grep") }
  
  if ( fast_grep ) {
    grep_command <- paste0("grep -m", num_events, " ")
  } else if ( fast_grep == FALSE ){
    grep_command <- "grep "
  } else {
    stop("fast_grep must be logical, either TRUE or FALSE")
  }
  
  # Check vst_id 
  if ( missing(vst_id) ) { stop("You didn't specified a vst_id!") } 
  
  # If there's more than one vst_id in the query
  if ( length(vst_id) >= 2 ) {
    
    # Turn the vst_id into a "greppable" multiple vst_id
    vst_id <- paste(vst_id, collapse = "|")
    
    # grep searches an extended regular expression
    grep_command <- paste0(grep_command, "-E ")
  }
  
  # Either first get the file header and then the AS event line. split_temp TRUE
  # or get header and AS line with one command only with process substitution < ( ). split_temp FALSE.
  if ( split_temp == TRUE ) {
    
    # 1 -- Get table header
    get_header_line = paste0("head -n1 ", inclusion_tbl, " > ", 
                             tmp_dir, "/tmp_inclusion_tbl_hearder.tab")
    
    system(command = get_header_line)
    
    
    # 2 -- Get AS event PSI
    get_exon_line <- paste0(grep_command, "'", vst_id, "' ", 
                            inclusion_tbl, " > ",
                            tmp_dir, "/tmp_vst_id_PSI.tab")
    system(command = get_exon_line)
    
    
  } else if ( split_temp == FALSE ) {
    # 1, 2 -- Get table header and AS event PSI
    get_header_n_exon_line <- paste0( "cat <(head -n1 ", inclusion_tbl,
                                      ") <(", grep_command, "'", vst_id, "' ", 
                                      inclusion_tbl, ") > ", tmp_dir,
                                      "/tmp_vst_id_PSI.tab")
    
    # the R system() command calls the bourne shell (sh) instead of the bourne again shell (bash)
    # to overcome this print the command in the bourne shell and pipe the output into bash
    bash_header_AS_line <- paste('echo "', get_header_n_exon_line, '" | bash')
    
    system( command = bash_header_AS_line )
    
  } else {
    stop("split_tmp must be logical, either TRUE or FALSE")
  }
  
  # 3 -- Check that temp file is not empty
  exon_PSI_tmp_info <- file.info(file.path(tmp_dir, "tmp_vst_id_PSI.tab"))
  
  if ( exon_PSI_tmp_info$size == 0) {
    stop("Couldn't find event ", vst_id, " in table ",
         inclusion_tbl, " with the command: ", grep_command)
  }
  
  if ( split_temp == TRUE ) {
    # 4 -- Read in table header
    tmp_header <- read.delim(file = file.path(tmp_dir,
                                              "tmp_inclusion_tbl_hearder.tab"),
                             header = T, check.names = F, stringsAsFactors = F)
    
    # 5 -- Read in AS event vast table 
    tmp_exon <- read.delim(file = file.path(tmp_dir, "tmp_vst_id_PSI.tab"),
                           header = F, check.names = F, stringsAsFactors = F)
    
    # 6 -- Merge header and exon processed vast table
    colnames(tmp_exon) <- colnames(tmp_header)
    vst_psi_tbl <- rbind(tmp_header, tmp_exon)
    
  } else if ( split_temp == FALSE ) {
    # 4, 5, 6 -- Read in AS event vast table 
    vst_psi_tbl <- read.delim(file = file.path(tmp_dir, "tmp_vst_id_PSI.tab"),
                              header = T, check.names = F, stringsAsFactors = F)
    
  } else {
    stop("Something went wrong when reading tmp_vst_id_PSI.tab in\n", tmp_dir)
  }
  
  # 7 -- Info on what was found
  if ( verbose ) {
    
    if ( num_events == 1) {
      message("Found ", vst_id ) 
    }
    
    if ( num_events >= 2 ) {
      # Inform only on events that were found
      message("Found ", paste(vst_psi_tbl$EVENT, collapse = ", ") ) 
      vst_id <- unlist(strsplit(x = vst_id, split = "|", fixed = T) )
      missing_ids <- vst_id[ which( !vst_id %in% vst_psi_tbl$EVENT ) ]
      
      if ( identical(missing_ids, character(0) ) ) { 
        message("Found all queried events!")
      } else {
        warning("Missing event(s) not foud: ", paste(missing_ids, " ") )
      }
    }
  }
  
  # 8 -- Clean up
  if ( clear_tmp ) {
    
    if ( split_temp == TRUE ) {
      file.remove( file.path(tmp_dir, "tmp_inclusion_tbl_hearder.tab") )
    }
    
    file.remove( file.path(tmp_dir, "tmp_vst_id_PSI.tab") )
    if ( verbose ) { message("Temporary files have been removed.") } 
  }
  # 9 -- Return a tbl
  return(vst_psi_tbl)
}

#' Grep the expression values of a gene from a vast-tools output table. 
#'
#' @param vst_expression_tbl Path to a vast-tools expression table (either cRPKM or TPM)
#' @param ensembl_gene_id an ensembl gene id that is used to grep a line
#' @param tmp_dir a temporary directory where to same the intermediate files
#' @param verbose Logical, print some info?
#' @param clear_tmp Logical, remove the temporary files after the function is done?
#' @param fast_grep Logical, if TRUE use 'grep -m1' to stop the search after first match found. If FALSE all vst_expression_tbl is read.
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' grep_gene_expression(path_to_vst_Expression_tbl, ensembl_gene_id) |>
#'     tidy_vst_expr(ID_cols = "ensembl_geneid", expression_unit = "TMP")
grep_gene_expression <-  function(vst_expression_tbl, ensembl_gene_id, 
                                  tmp_dir = tempdir(), verbose = FALSE, 
                                  clear_tmp = TRUE, fast_grep = TRUE ) {
  
  if ( verbose ) { message("Looking for the expression levels of: ", 
                           ensembl_gene_id)
  }
  
  # 1 -- Get gene expression table header
  get_header_line <- paste0("head -n1 ", vst_expression_tbl, " > ", 
                            tmp_dir, "/tmp_expression_tbl_hearder.tab")
  
  system(command = get_header_line)
  
  # 2 -- Read in table header
  tmp_header <- read.delim(file = file.path(tmp_dir,
                                            "tmp_expression_tbl_hearder.tab"),
                           header = T, check.names = T, stringsAsFactors = F)
  
  # 3 -- Get gene EXPRESSION values for the gene
  # grep is always done using the ensembl gene id
  # the -m1 instructs grep to stop reading a file after 1 matching line. 
  if (fast_grep == TRUE) {
    grep_command <- "grep -m1 "
  } else if( fast_grep == FALSE ){
    grep_command <- "grep "
  } else {
    stop("fast_grep must be logical, either TRUE or FALSE")
  }
  
  get_expression_line <- paste0(grep_command, "'", ensembl_gene_id, "' ",
                                vst_expression_tbl, " > ",
                                tmp_dir, "/tmp_gene_expression.tab")
  system(command = get_expression_line)
  
  # 4 -- Read in GENE expression line
  gene_expr_tmp_info <- file.info(file.path(tmp_dir, "tmp_gene_expression.tab"))
  
  if ( gene_expr_tmp_info$size == 0 ) {
    stop("Couldn't find gene ", ensembl_gene_id, " in table ",
         vst_expression_tbl)
  } else { if ( verbose ) { message("Found ", ensembl_gene_id ) } }
  
  tmp_expression <- read.delim(file = file.path(tmp_dir, "tmp_gene_expression.tab"),
                               header = F, check.names = F,
                               stringsAsFactors = F)
  
  # If file has more than one row, take just the last one.
  # This is necessary in case the ensembl_gene_id provided is in one of the headers
  # but with the fast_grep = TRUE parameter this would probably not work.
  if ( nrow( tmp_expression ) > 1 ) {
    warning("When grep gene expression values more than one line was found!\n",
            "Returning only the last line, which is probably the what you want?")
    tmp_expression <- tail(tmp_expression, 1)
  }
  
  # 5 -- Merge header and exon processed vast table
  colnames(tmp_expression) <- colnames(tmp_header)
  vst_expression_tbl <- rbind(tmp_header, tmp_expression)
  
  # 6 -- Clean up
  if ( clear_tmp ) {
    file.remove( file.path(tmp_dir, "tmp_expression_tbl_hearder.tab") )
    file.remove( file.path(tmp_dir, "tmp_gene_expression.tab") )
    if ( verbose ) { message("Temporary files have been removed.") } 
  }
  
  # 7 -- Fix first column name 
  # If the value in the first column start as an ensembl gene id
  if ( grepl(pattern = "^ENS", x = as.character(vst_expression_tbl[, 1])) ) {
    if ( verbose ) { message("First column is an ENSEMBL gene id") }
    colnames(vst_expression_tbl)[1] <- "ensembl_gene_id"
  }
  
  # 7 -- Return a tbl
  return( vst_expression_tbl )
}

#' Reshape a wide vast-tools inclusion table into a tidyverse-friendly long table. 
#'
#' @param vst_psi_tbl A dataframe generated with `grep_psi()` from a vast-tools inclusion table.
#' @param num_id_cols How many first num_id_cols to consider as info/metadata/IDs of the AS event ID in the table.
#' @param num_of_Score_Types How many quality scores are present in the columns with headers ending with `quality_col_suffix`.
#' @param quality_col_suffix Suffix identifying the quality control columns. Default is "-Q".
#' @param return_quality_scores Logical. Do you want the individual scores to be returned in the output data.frame? Defaul `TRUE`.
#' @param return_S1_only Logical. Return only the first quality score, default `TRUE`, use `FALSE` for returning ALL `num_of_Score_Types` quality scores for each event.
#' @param verbose Print out information
#' @param add_ID_col Logical. Do you want to add an extra `col_ID_name` to the output data.frame? Default `FALSE`.
#' @param col_ID_name Extra column that could be used to add a new identifier to the data. Default is "banana".
#'
#' @description This function works with inclusion tables generated both by vast-tools combine or compare. In the latter case, if there's a 'dPSI' column it will be automatically detect and included in the reshaped table.
#'
#' @details This function works well in conjunction with `grep_psi()` function.
#'
#' @return A reshaped `data.frame` in long format as `tibble`.
#' @import dplyr
#' @import stringr
#' @import tidyr
#' 
#' @export
#'
#' @examples
#' grep_psi(path_to_vst_PSI_tbl, vst_id = "HsaEX0000001") |>
#'     tidy_vst_psi() -> psi_tbl
tidy_vst_psi <- function(vst_psi_tbl, num_id_cols = 6, num_of_Score_Types = 5,
                         quality_col_suffix = "-Q", return_quality_scores = TRUE,
                         return_S1_only = TRUE, verbose = FALSE, add_ID_col = FALSE,
                         col_ID_name = 'banana') {

  if (return_S1_only == TRUE & return_quality_scores == F) {
    return_S1_only <- FALSE
    # stop("If return_S1_only is TRUE return_quality_scores must be true as well")
  }
  
  # ID columns
  ID_cols <- 1:num_id_cols
  rgx_Q <- paste0(".*(?<!", quality_col_suffix, ")$")
  # All columns with PSI info: regex with negative lookhead
  PSI_cols <- grep(pattern = rgx_Q, x = colnames(vst_psi_tbl),
                   perl = T)
  PSI_cols <- PSI_cols[-ID_cols]
  
  # If it's a vast-tools compare table with dPSI remove it:
  if ( any(colnames(vst_psi_tbl[PSI_cols]) == "dPSI") ) {
    # check that dPSI is the last column
    if (colnames(vst_psi_tbl[ncol(vst_psi_tbl)]) == "dPSI") {
      if (verbose) { 
        message("I see that the last column is called dPSI. ",
                "I'm going to process this table as a vast-tools compare output table.")
      }
      # remove the dPSI from the PSI_cols
      PSI_cols <- PSI_cols[-which(colnames(vst_psi_tbl[PSI_cols]) == 'dPSI' )]
      ncol_dPSI <- which(colnames(vst_psi_tbl) == 'dPSI')
      vst_compare_tbl <- TRUE
    } else {
      # There's a column called dPSI but is not the last one so I'm not gonna consider this
      # as a vast-tools compare table
      # remove this non-last dPSI column from the PSI_cols
      PSI_cols <- PSI_cols[-which(colnames(vst_psi_tbl[PSI_cols]) == 'dPSI' )]
      vst_compare_tbl <- FALSE
    }
  } else {
    # This is not a vast-tools compare table
    vst_compare_tbl <- FALSE
  }
  
  # All columns with Q information
  rgx_Q <- paste0(".*", quality_col_suffix, "$")
  Qual_cols <- grep(pattern = rgx_Q, x = colnames(vst_psi_tbl),
                    perl = T)
  
  # Check that there are the quality score columns
  if ( identical(Qual_cols, integer(0)) ) {
    stop("The input inclusion table does NOT have quality score tables! ",
         "Can't process further if there are no columns with ending with: ",
         quality_col_suffix, "\n",
         "Maybe in the future I'll implement a way to handle this case. ",
         "The tables that have these columns are the output of vast-tools combine and compare.",
         "The tables without this columns are the output of vast-tools tidy.")
  }
  
  check_Q_col_suffix <- grepl(pattern = quality_col_suffix, 
                              x = colnames(vst_psi_tbl)[Qual_cols], 
                              perl = T)
  
  if ( all(check_Q_col_suffix) ) {
    
    if (verbose) { message("Correctly identified quality columns") }
    
  } else if ( !all(check_Q_col_suffix) ) {
    stop("Can't correctly identify quality columns. The suffix", 
         quality_col_suffix, " Try chanding the quality_col_suffix parameter.")
  } else {
    stop("There's a weird error. Debug this and/or ask Nicco")
  }
  
  
  veryverbose <- FALSE
  # Initialise dataframe as vector to populate with results from different pages
  li_vast_Q_cols <- list()
  # Loop through the Q columns and turn them into a dataframe
  for (i in 1:length(Qual_cols)) {
    
    # Use this if statement only for debugging
    if (veryverbose) {
      message("Column ", i, " ", colnames(vst_psi_tbl)[Qual_cols][i] )
    }
    
    # Split list of 5 elements by comma in quality score columns 
    # Merge all dataframes by populating a list
    vst_psi_tbl |>
      # Extract the quality score column i
      select(colnames(vst_psi_tbl)[Qual_cols][i]) |>
      # turn values to character
      pull() |>
      str_split_fixed(pattern = ',', n = num_of_Score_Types) |> 
      as.data.frame() |>
      setNames(paste(colnames(vst_psi_tbl)[Qual_cols][i], 
                        c("S1", "S2", "S3", "S4", "S5"), 
                        sep = "_")) -> li_vast_Q_cols[[i]]
  }
  
  # Turn into dataframe
  df_vast_Q <- as.data.frame(li_vast_Q_cols)
  
  # Re-define the column names to have the code more robust in case for example
  # the column names start with a number. In this way I have the col names well set
  # And the columns can be processed well with regexes later on. 
  colnames(df_vast_Q) <- unlist(lapply(li_vast_Q_cols, function(x) names(x)))
  
  # Coerc factors col to character col
  df_vast_Q <- mutate_if(df_vast_Q, is.factor, as.character)
  
  if ( nrow(vst_psi_tbl) == nrow(df_vast_Q) ) {
    # Merge data
    vst_psi_tbl_Q <- cbind(vst_psi_tbl, df_vast_Q) 
    if ( ncol(vst_psi_tbl_Q) == ncol(vst_psi_tbl) + ncol(df_vast_Q) ) {
      if (verbose) { message("Tidying up passed sanity check point 1") }  
    } else {
      stop("Something went wrong in reshaping the data around check point 1B")
    }
    
  } else {
      stop("Something went wrong in reshaping the data around check point 1A")
  }
    
  Qual_cols_Score <- grep(pattern = "-Q_S[1-5]$", x = colnames(vst_psi_tbl_Q),
                          perl = T)
  
  # Reshape to long format only the PSI quality score info
  pivot_longer(data = vst_psi_tbl_Q[, c(ID_cols, Qual_cols_Score )],
               cols = colnames(vst_psi_tbl_Q)[Qual_cols_Score],
               names_to = "Quality_Score",
               values_to = "Quality_Score_Value") |>
    # Change cols to make it compatible with joining
    mutate(Quality_Score_Type = gsub(".*Q_", "", Quality_Score, perl = T)) |>
    mutate(Sample = str_remove(Quality_Score, "\\-Q_S[1-5]$")) |> 
    select(-Quality_Score) -> lng_vst_psi_quality_tbl
  
  # Reshape to long format the PSI columns without the Quality info
  pivot_longer(data = vst_psi_tbl[, c(ID_cols, PSI_cols)],
               cols = colnames(vst_psi_tbl)[PSI_cols],
               names_to = "Sample",
               values_to = "PSI") -> lng_vst_psi_tbl
  
  # Sanity check
  check2 <- nrow(lng_vst_psi_tbl) * num_of_Score_Types == nrow(lng_vst_psi_quality_tbl)
  
  if( check2 ) {
    if (verbose) { message("Tidying up passed sanity check point 2")} 
  } else {
    stop("Something went wrong in reshaping the data around check point 2")
  }
  
  # Coerce factors col to character col
  lng_vst_psi_tbl <- mutate_if(lng_vst_psi_tbl, is.factor, as.character)
  
  # Merge datasets
  tidy_vst_psi_tbl <- left_join(lng_vst_psi_tbl, lng_vst_psi_quality_tbl,
                                       by = c("GENE", "EVENT", "COORD", "LENGTH",
                                              "FullCO", "COMPLEX","Sample") )
  
  # Add back the dPSI column if present
  if(vst_compare_tbl == TRUE) {
    tidy_vst_psi_tbl <- left_join(tidy_vst_psi_tbl, 
                                  vst_psi_tbl[, c(ID_cols, ncol_dPSI)],
                                  by = c("GENE", "EVENT", "COORD", "LENGTH",
                                         "FullCO", "COMPLEX"))
  }

  if (is.numeric(tidy_vst_psi_tbl$PSI)) {
    if (verbose) { message("Tidying up passed sanity check point 3") }
  } else{
    stop("PSI is not numeric")
  }
  
  if ( any( duplicated(tidy_vst_psi_tbl) ) ) {
    stop("Something went wrong in reshaping the data around, there are duplicates in the data.frame!")
  } else if ( !any( duplicated(tidy_vst_psi_tbl) ) ) {
    if (verbose) { message("Tidying up passed sanity check point 4") }
  } else {
    stop("Cannot check is there are duplicates is the data")
  }
  
  if ( return_quality_scores == TRUE ) {
    
    # Do nothing output already has scores
    if ( return_S1_only ) {
      tidy_vst_psi_tbl <- subset(tidy_vst_psi_tbl, Quality_Score_Type == "S1")
      
      # Coerc quality scores of Score 1 to factors
      Quality_Score_Values <- c("N", "VLOW", "LOW", "OK", "SOK")
      Quality_Score_Values <- factor(Quality_Score_Values,
                                     levels = Quality_Score_Values)
      tidy_vst_psi_tbl$Quality_Score_Value <- factor(
        tidy_vst_psi_tbl$Quality_Score_Value, 
        Quality_Score_Values
      )
      
    }
    
  } else if ( return_quality_scores == FALSE ) {
    
    cols_to_return <- c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX", "Sample", "PSI")
      
    if (vst_compare_tbl) {
      cols_to_return <- c(cols_to_return, "dPSI")
    }
    select(tidy_vst_psi_tbl, cols_to_return) |>
      unique() -> tidy_vst_psi_tbl
  } else {
    stop("return_quality_scores must me a logical (TRUE or FALSE)!")
  }

  # Add some sort of grouping/experiment name to the reshaped dataset
  if ( add_ID_col == TRUE ) {
    tidy_vst_psi_tbl$Group_Name <- col_ID_name
  }
  
  # Coerce the AS event type 
  tidy_vst_psi_tbl <- mutate(.data = tidy_vst_psi_tbl,
                             COMPLEX = factor(COMPLEX,
                                              levels = c('S', 'C1', 'C2', 'C3',
                                                         'MIC', 'ANN', 'Alt3',
                                                         'Alt5', 'IR') ) )
  return(tidy_vst_psi_tbl)
}

#' Reshape a wide vast-tools expression table to a long format
#'
#' @param data A dataframe generated with `grep_gene_expression` from a vast-tools gene expression table.
#' @param expression_unit Character: either `TPM` or `cRPKM` based on the gene expression unit in `data`.
#' @param ID_cols The columns names to use as IDs in `data`
#'
#' @return A reshaped `data.frame` in long format as `tibble`.
#' @import tidyr
#' @export
#'
#' @examples
#' grep_gene_expression(path_to_vst_Expression_tbl, ensembl_gene_id) |>
#'     tidy_vst_expr(ID_cols = "ensembl_geneid", expression_unit = "TMP")
tidy_vst_expr <- function(data, expression_unit = c("TPM", "cRPKM"),
                          ID_cols = c("ID", "NAME", "Names", "ensembl_gene_id")) {
  
  if (! any(expression_unit %in% c("TPM", "cRPKM")) ) {
    stop("expression_unit must be equal to TPM or cRPKM")
  }
  
  # If there are no ID_Cols in the table look for anther column IDs
  if ( !any( colnames(data) %in% ID_cols) ) {
    # This happens cause the header of the vast-tools expression table NORM
    # has only ensembl gene id in a column called "Names"
    # While the not normalised expression tables as 2 columns called ID and NAME
    ID_cols <- "Names"
  }
  
  # Select the expression columns in the data table
  sample_cols <- which(!colnames(data) %in% ID_cols)
  # Reshape into a longer format
  pivot_longer(data = data,
               cols = all_of(colnames(data)[sample_cols] ),
               names_to = "Sample",
               values_to = "Gene_Expr") -> lng_vst_expr_tbl
  
  # Coerc to numeric if not already
  if ( !is.numeric(lng_vst_expr_tbl$Gene_Expr) ) {
    lng_vst_expr_tbl$Gene_Expr <- as.numeric(lng_vst_expr_tbl$Gene_Expr)
  }
  
  return(lng_vst_expr_tbl)
}

#' Given a vast-tools AS event ID, return the species name.
#'
#' @param vst_id An alternative splicing event name from vast-tools.
#' @param latin_name Logical, whether or not you want the scientific name of the species.
#'
#' @return A character
#'
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' guess_species(vst_id = "HsaEX000001")
#' # Human
guess_species <- function(vst_id, latin_name = F) {
  
  guessed_species <- str_extract(string = vst_id, pattern = "^[A-Z][a-z][a-z]")
  
  if ( latin_name == FALSE) {
    if (guessed_species == "Hsa") { guessed_species <- "Human" } 
    if (guessed_species == "Mma") { guessed_species <- "Macaque" } 
    if (guessed_species == "Mmu" ) { guessed_species <- "Mouse"}
    if (guessed_species == "Rno" ) { guessed_species <- "Rat"}
    if (guessed_species == "Bta" ) { guessed_species <- "Cow"}
    if (guessed_species == "Cmi" ) { guessed_species <- "Elephant_Shark"}
    if (guessed_species == "Ath" ) { guessed_species <- "Arabidopsis"}
    
  } else if ( latin_name == TRUE ) {
    if (guessed_species == "Hsa") { guessed_species <- "hsapiens" } 
    if (guessed_species == "Mma") { guessed_species <- "mmulatta" } 
    if (guessed_species == "Mmu" ) { guessed_species <- "mmusculus"}
    if (guessed_species == "Rno" ) { guessed_species <- "rnovergicus"}
    if (guessed_species == "Bta" ) { guessed_species <- "btaurus"}
    if (guessed_species == "Cmi" ) { guessed_species <- "cmilli"}
    if (guessed_species == "Ath" ) { guessed_species <- "athaliana"}
    
  } else {
    
    stop("latin_name must be a logical, either TRUE or FALSE!\n")
  }
  
  if (nchar(guessed_species) == 0) { 
    stop("Couldn't guess the species for the vastID: ", vst_id)
  }
  
  return(guessed_species)
}

#' Wrapper for `readr::read_delim` for a tab delimited table. Fastest way to read this kind of data into `R`.
#'
#' @param path A character vector providing the 'path/to/the/table/to/read'.
#' @param verbose Logical. Use `TRUE` to know the table file size being read into `R`.
#' @param ... extra info to pass to `read_delim`.
#'
#' @return A tibble
#' @importFrom readr read_delim locale
#' @export
#'
#' @examples
#' read_vst_tbl(path = "location/to/tab/delimited/inclusion/table/file.tab") |>
#'    tidy_vst_psi()
#'  
read_vst_tbl <- function(path, verbose = FALSE, ...) {
  check_file(path)
  
  if (verbose) {
    message('Reading a vast-tools table of: ',
            round(file.info(path)$size / 10e3, 2), ' Kbytes')
  }
  return(read_delim(file = path, delim = '\t', col_names = T,
                    locale = locale(decimal_mark = "."), na = "NA", ...) )
}


#' Wrapper function to import pre-processed vast-tool analysis tables of Mouse Development ENCODE data.
#'
#' @param inclusion_tbl Path to vast-tools inclusion table that contains a `vst_id` event. Use default (`NULL`) to use pre-defined PSI table.
#' @param vst_expression_tbl Path to a vast-tools expression table (either `cRPKM` or `TPM`) that contains a `ensembl_gene_id` gene. Use default (`NULL`) to use pre-defined normalised TPMs table.
#' @param metadata_path A path to a metadata table. Use default (`NULL`) to use pre-defined metadata table.
#' @param ensembl_gene_id A valid Mouse ENSEMBL gene ID.
#' @param vst_id A valid Mouse Vast-tools ID.
#' @param filter_tbl Whether you want the full table to the filtered one with all the data.
#'
#' @return A tibble
#' @importFrom readr read_delim
#' @importFrom dplyr mutate select relocate left_join group_by
#' @export
#' 
#' @seealso [plot_mouse_tissue_devel()]
#'
#' @examples
#' # Get Pax6 gene expression and Pax6 exon 6 PSI across mouse tissues development.
#' get_mouse_tissue_devel_tbl(ensembl_gene_id = "ENSMUSG00000027168", vst_id = "MmuEX0033804") 
get_mouse_tissue_devel_tbl <- function(inclusion_tbl = NULL, 
                                       vst_expression_tbl = NULL,
                                       metadata_path = NULL,
                                       ensembl_gene_id, vst_id, 
                                       filter_tbl = TRUE) {
  # Check params 
  if ( missing(ensembl_gene_id) ) { stop("You didn't specified an ENSEMBL gene ID!") } 
  if ( missing(vst_id) ) { stop("You didn't specified a vst_id!") } 
  if ( !is.logical(filter_tbl) ) { stop("filter_tbl must be TRUE or FALSE")}
  
  # check if CRG cluster is mounted 
  if ( dir.exists('~/mnt/narecco/projects/07_Suz12AS/data' ) )  {
    mounted <- TRUE
  } else if ( dir.exists('/users/mirimia/narecco/projects/07_Suz12AS/data') ) {
    mounted <- FALSE
  } else {
    stop("Can't figure out if the CRG cluster is mounted")
  }
  
  if ( mounted == TRUE) {
    data_dir <- file.path('~/mnt/narecco/projects/07_Suz12AS/data')
  } else if ( mounted == FALSE ) {
    data_dir <- file.path('/users/mirimia/narecco/projects/07_Suz12AS/data')
  } else {
    stop("Can't specify the data dir")
  }

  vast_out <- file.path(data_dir, 'INCLUSION_tbl/Mouse_Development/vast_tools/vast_out')
    
    if ( is.null(inclusion_tbl) ) {
        inclusion_tbl <- file.path(vast_out, "INCLUSION_LEVELS_FULL-mm10-469-v251.tab")
    }
    
    if ( is.null(vst_expression_tbl) ) {
        vst_expression_tbl <- file.path(vast_out, "TPM-mm10-469-NORM.tab")
    }
    
    if( is.null(metadata_path) ) {
        metadata_dir <- file.path(data_dir, 'SRA_tbls/Mouse_Development')
        metadata_path <- file.path(metadata_dir, 'All_ENCODE_Mouse_Dev_Tissues.tab')
    }
    
    mtdt <- read_delim(file = metadata_path, delim = "\t", show_col_types = FALSE,
                       col_names = c("Run", "Sample", "Group", "Tissue", "Stage", "Sex"))
    
    # Import PSI and normalised gene expression TPMs
    grep_psi(inclusion_tbl = inclusion_tbl, vst_id = vst_id, verbose = F) |>
        tidy_vst_psi() -> psi_tbl
    
    grep_gene_expression(vst_expression_tbl = vst_expression_tbl,
                         ensembl_gene_id = ensembl_gene_id, verbose = F) |>
        tidy_vst_expr() -> expr_tbl
    
    # Split individual replicates and merged samples
    left_join(psi_tbl, expr_tbl, by = "Sample") |>
        mutate(Type =  case_when(Sample %in% mtdt$Sample ~ "Replicate",
                                 !Sample %in% mtdt$Sample ~ "Merge")) -> parsed_tbl
    
    subset(parsed_tbl, Type == "Replicate") |>
        left_join(mtdt, by = "Sample") -> individual_reps_tbl
    
    # Merge
    merged_reps_tbl <- subset(parsed_tbl, Type == "Merge") 
    colnames(merged_reps_tbl)[colnames(merged_reps_tbl) == "Sample"] <- "Group"
    merged_reps_tbl <- left_join(merged_reps_tbl, mtdt, by = "Group")
    
    parsed_tbl <- rbind(individual_reps_tbl, merged_reps_tbl)
    
    # Calculate mean PSI between replicates.
    parsed_tbl |>
        group_by(Tissue, Stage) |>
        mutate(mean_PSI = mean(PSI, na.rm = T),
               mean_Gene_Expr = mean(Gene_Expr)) |>
        relocate(mean_PSI, .after = PSI)  -> parsed_tbl
    
    # Clean up names and factorise tissues
    parsed_tbl |>
        mutate(Tissue = ifelse(Tissue == "CentralNervousSystem", yes = "CNS", no = Tissue),
               Tissue = ifelse(Tissue == "CraniofacialProminence", yes = "Craniofacial\nProminence", no = Tissue) ) |>
        mutate(Tissue = factor(Tissue, 
                               levels = c("Craniofacial\nProminence", "CNS", "Brain", 
                                          "Forebrain", "Midbrain", "Hindbrain", 
                                          "NeuralTube", "Heart", "SkeletalMuscle", "Limb",
                                          "Liver", "Stomach", "Intestine", "Kidney",
                                          "Lung", "Bladder",
                                          "Spleen", "AdrenalGland", "Thymus") ) ) |>
        relocate(Tissue, .after = "Sample") -> parsed_tbl
    
    if ( filter_tbl == TRUE ) {
        
        subset(parsed_tbl, !Tissue %in% c("CNS", "Brain", "SkeletalMuscle", "Bladder",
                                          "Spleen", "AdrenalGland") ) |>
            subset(! Stage %in% c("E14", "E18") ) |>
            select( c("GENE", "EVENT", "Tissue", "Stage", "mean_PSI", "mean_Gene_Expr")) |>
            unique() -> tidy_tbl
    } else if ( filter_tbl == FALSE) {
        tidy_tbl <- parsed_tbl
    } else {
        stop("filter_tbl must be logical!")
    }
    return(tidy_tbl)
}

#' Extract the PSI from an inclusion table and return it as a matrix
#'
#' @param inclusion_tbl path to vast-tools inclusion table that contains a vst_id event.
#' @param vst_id vast-tools alternative splicing event to grep in the `inclusion_tbl`.
#' @param quality_thrshld vast-tools event quantification quality score threshold. Must be one of "N", "VLOW", "LOW", "OK", "SOK". For more info read the official documentation [here](https://github.com/vastgroup/vast-tools#combine-output-format) under "Column 8, score 1".
#' @param verbose Print out extra info.
#' 
#' @return A matrix
#' @import dplyr 
#' @import tibble
#' @export
#' @description The matrix contains the samples as row names and the vastID as column name
#'
#' @examples
#' gimme_psi_mat(inclusion_tbl = /path/to/inclusion/table/file.tab,
#'               vst_id = "HsaEX0000001")
gimme_psi_mat <- function(inclusion_tbl, vst_id, quality_thrshld = "N", 
                          verbose = FALSE) {
  # 1 -- Check input parameters
  if ( missing(vst_id) ) { stop("You didn't specified a vst_id!") } 
  if ( !any(quality_thrshld == c("N", "VLOW", "LOW", "OK", "SOK")) ) {
    stop("Can't understand the filtering option:\t", quality_thrshld,
         "\nThe parameter quality_thrshld must be one of N, VLOW, LOW, OK, or, SOK")
  }
  
  
  # 2 -- GET AS EVENT PSI AND INFO INTO A TIDY TABLE
  grep_psi(inclusion_tbl = inclusion_tbl, vst_id = vst_id, verbose = verbose ) |>
    tidy_vst_psi(verbose = verbose) -> psi_tbl
  
  # 3 -- KEEP SAMPLES WITH AS QUALITY >= quality_thrshld 
  # Set the Score 1 quality values as factors 
  Quality_Score_Values <- c("N", "VLOW", "LOW", "OK", "SOK")
  psi_tbl$Quality_Score_Value <- factor(psi_tbl$Quality_Score_Value,
                                        Quality_Score_Values)
  
  num_quality_thrshld <- as.numeric(factor(quality_thrshld, 
                                           levels = Quality_Score_Values))
  psi_tbl <- subset(psi_tbl, 
                    as.numeric(Quality_Score_Value) >= num_quality_thrshld)
  # 4 -- Turn it into a matrix
  psi_tbl |>
    select(Sample, PSI) |>
    column_to_rownames("Sample") |>
    setNames(vst_id) |>
    as.matrix() -> psi_mat
  
  return(psi_mat)
}

#' Extract the gene expression counts from an vast-tools expression table and return them as a matrix
#'
#' @param vst_expression_tbl Path to a vast-tools expression table (either cRPKM or TPM).
#' @param min_mean_count Filter out low expressed genes in the table read from `inclusion_tbl`. Defines the minimum row mean expression value across all samples that a gene must have to be selected. 
#' @param verbose Print out extra info.
#'
#' @return A matrix
#' @import dplyr 
#' @import tibble
#' @export
#' @description The matrix contains the samples as row names and the gene ID as column name. The gene ID used are the one present in the first column of the table (usually ENSEMBL gene ID).
#'
#' @examples
#' gimme_expr_mat(vst_expression_tbl = /path/to/expr/table/file.tab)
gimme_expr_mat <- function(vst_expression_tbl, min_mean_count = 5, 
                           verbose = FALSE) {
  
  # 1 -- Check input parameters
  if ( missing(vst_expression_tbl) ) { 
    stop("You didn't specified a path to the expression table!") 
    } 
  
  # 2-- DUMP ALL GENE EXPRESSIONS IN A TABLE
  dump <- read_vst_tbl(path = vst_expression_tbl, show_col_types = FALSE)
  
  # 3 -- CHECK IF FIRST COLUMN IS AN ENSEMBL GENE ID
  # If the value in the first column and first row start as an ENSEMBL gene ID
  if ( grepl(pattern = "^ENS", x = as.character(dump[1, 1]) ) ) {
    if ( verbose ) { message("First column is an ENSEMBL gene ID") }
    colnames(dump)[1] <- "ensembl_gene_id"
  }
  
  # 4 -- Parse the gene expression table into a matrix
  ## Fist check that first column names are all unique names.
  num_rows <- nrow(dump)
  num_uniq_names <- dump |> pull(var = 1) |> unique() |> length()
  
  if (num_rows != num_uniq_names) {
    stop("The vast-tools gene expression table ", 
         "(specified with vst_expression_tbl) does NOT have unique ",
         "gene IDs! There are ", num_rows, " rows in the table, but only ",
         num_uniq_names, " unique gene name IDs. Do something about!")
  }
  # Turn the first column into row names
  dump <- column_to_rownames(dump, var = colnames(dump)[1])
  # Keep only columns that are numeric
  dump <- select(dump, where(is.numeric) )
  
  # 5 -- Filter out low expressed genes and transpose
  gene_expr_mat <- t(dump[which(rowMeans(dump, na.rm = T) >= min_mean_count), ])
  
  return(gene_expr_mat)
}


#' Calculate correlation of one AS event PSI vs many genes expression levels (1 vs many).
#' 
#' @param inclusion_tbl path to vast-tools inclusion table that contains a vst_id event.
#' @param vst_id vast-tools alternative splicing event to grep in the `inclusion_tbl`.
#' @param quality_thrshld vast-tools event quantification quality score threshold. Must be one of "N", "VLOW", "LOW", "OK", "SOK". For more info read the official documentation [here](https://github.com/vastgroup/vast-tools#combine-output-format) under "Column 8, score 1".
#' @param vst_expression_tbl Path to a vast-tools expression table (either cRPKM or TPM).
#' @param min_mean_count Filter out low expressed genes in the table read from `inclusion_tbl`. Defines the minimum row mean expression value across all samples that a gene must have to be selected. 
#' @param corr_method Either `spearman`, `pearson`, or `kendall` passed to the function `cor()`. 
#' @param num_genes Return only the top and bottom number of genes. Integer number greater or equal than 1. Default is `NULL` returning all genes.
#' @param map_ID_2_names Logical. Whether or not to map the ENSEMBL gene IDs to gene names. Can be used only if `num_genes` is specified and the table contains ENSEMBL gene ID (check automatically).
#' @param species Species character to use to map the ENSEMBL gene ID. Used by `gimme_mart()` to built a bioMaRt object. Default is guessed from `vst_id`.
#' @param verbose Print out information
#' @param ... Extra parameters passed to `cor()` like `use = "complete.obs"`.
#'
#' @return A tibble
#' @import dplyr 
#' @import tibble
#' @export
#' @description This function calculates the correlation between an alternative splicing event Percentage of Sequence Inclusion (PSI) with the expression levels (i.e. counts) of all genes present in the vast-tools expression table, after some basic filtering. The type of correlation can be either Spearman, Pearson, or Kendall. If using this type of correlation extra parameters can be passed to the `cor()` R function used to calculate the correlations. The mapping of ENSEMBL gene IDs to gene names is available only if selecting best correlating genes. This is done in order to avoid to query ENSEMBL bioMart servers with a massive gene list.
#' 
#'
#' @examples
#' # Return one PSI to all gene correlation
#' gimme_PSI_expr_corr(inclusion_tbl = psi_path, vst_id = "HsaEX0000001", 
#'                     vst_expression_tbl = expr_path, corr_method = "spearman", 
#'                     use = "complete.obs", verbose = TRUE ) -> corr_df
#'                     
#' # Return one PSI to all genes correlation filtered for the top and bottom correlating genes and map the IDs to names.
#' gimme_PSI_expr_corr(inclusion_tbl = psi_path, vst_id = "HsaEX0000001",
#'                     quality_thrshld = "VLOW", 
#'                     vst_expression_tbl = expr_path, min_mean_count = 100,
#'                     corr_method = "spearman", use = "complete.obs"
#'                     num_genes = 10, map_ID_2_names = T, species = "hsapiens",
#'                     verbose = TRUE ) -> best_corr_df                
gimme_PSI_expr_corr <- function(inclusion_tbl, vst_id, quality_thrshld = "N",
                                vst_expression_tbl, min_mean_count = 5,
                                corr_method = c("spearman", "pearson", "kendall"),
                                num_genes = NULL, 
                                map_ID_2_names = FALSE, species,
                                verbose = FALSE, ...) {
  # 1 ---- CHECK INPUT PARAMETERS ----
  if ( missing(vst_id) ) { stop("You didn't specified a vst_id!") } 
  if ( missing(corr_method) ) { 
    stop("You didn't specified an correlation method!",
         "Use '?cor' read more about what method to use.") 
  } 
  if ( ! any( corr_method %in%  c("spearman", "pearson", "kendall") ) ) {
    message("corr_method must be either:",
            "'spearman', 'pearson' or, 'kendall'.",
            "Use '?cor' read more about what method to use.")
  }  
  if (! is.numeric(min_mean_count) ) { stop("min_mean_count must be a number!")}
  if (map_ID_2_names == TRUE) {
    # Guess the species
    if ( missing(species) ) { 
      species <- guess_species(vst_id = vst_id, latin_name = TRUE)
    }
    
    if ( is.null(species) ) {
      stop("You need to specify the species to use for mapping the ENSEMBL",
           "gene IDs! Use '?gimme_mart()' to check which species are supported")
    }
  }
 
  # 2 -- GET AS EVENT PSI AND INFO INTO A MATRIX
  gimme_psi_mat(inclusion_tbl = inclusion_tbl, vst_id = vst_id, 
                quality_thrshld = quality_thrshld, verbose = verbose) -> psi_mat
  
  # 3-- GET GENE EXPRESSION COUNTS INTO A MATRIX
  gimme_expr_mat(vst_expression_tbl, min_mean_count = min_mean_count, 
                 verbose = verbose) -> gene_expr_mat
  
  # 4 -- CHECK IF COLUMN NAMES OF THE EXPRESSION MATRIX IS AN ENSEMBL GENE ID
  # If the value in the first column and first row start as an ENSEMBL gene ID
  if ( grepl(pattern = "^ENS", x = as.character(colnames(gene_expr_mat)[1]) )) {
    if ( verbose ) { message("First column is an ENSEMBL gene ID") }
    # Specify that these gene ID are mappale
    mappable <- TRUE
  }
  
  # 7 -- FILTER OUT SAMPLES FOR WHICH THE PSI WAS DISCARDED
  gene_expr_mat <- gene_expr_mat[rownames(gene_expr_mat) %in%  rownames(psi_mat), ]
  
  # 8 -- PRINT MATRIXES DIMENSIONS
  if (verbose) {
    message("Splicing inclusion table: ", nrow(psi_mat), " samples x ",
            ncol(psi_mat), " alternatively spliced event")
    message("Gene expression table: ", nrow(gene_expr_mat), " samples x ",
            ncol(gene_expr_mat), " genes")
  }
  
  # quick check that all samples are there
  stopifnot(all(rownames(psi_mat) %in% rownames(gene_expr_mat)))
  
  # 8 -- COMPUTE CORRELATION OF 1 AS EVENT PSI VS ALL GENES
  cor(x = psi_mat[, 1], y = gene_expr_mat, ..., method = corr_method) |>
    t() |>
    as.data.frame() |>
    setNames("Correlation") |>
    arrange(desc(Correlation)) |>
    rownames_to_column("gene_name") |>
    as_tibble() -> genes_corr_df
  
  # 9 -- RETURN ONLY TOP AND BOTTOM CORRELATING GENES 
  if ( all(is.numeric(num_genes), num_genes >= 1) ) {
    genes_corr_df |> (\(x) {
      rbind( head(x, num_genes), tail(x, num_genes))
    })() -> genes_corr_df
    
    # 10 -- TRY TO MAP ENSEMBL GENE IDs TO GENE NAMES
    if ( map_ID_2_names ) {
      if (mappable == TRUE) {
        
        # Use ENSEMBL's BioMart to map ensembl gene IDs to gene names
        ensembl <- gimme_mart(verbose = verbose, species = species)
        
        # Query for external gene names
        lapply(genes_corr_df$gene_name, function(x) {
          ensembl_id_2_gene_name(ensembl_gene_id = x,
                                 only_gene_name = T,
                                 mRt_objct = ensembl,
                                 verbose = verbose)
        }) |> unlist() -> mapped_gene_names
        
        # Add gene names
        genes_corr_df$external_gene_name <- mapped_gene_names
        genes_corr_df |>
          # If gene name not available use ENSEMBL gene ID.
          mutate(external_gene_name = ifelse(test = is.na(mapped_gene_names),
                                             yes = gene_name, 
                                             no = mapped_gene_names)) |>
          # Order external gene name as factors. Useful for plotting
          mutate(external_gene_name = fct_reorder(.f = external_gene_name, 
                                                  .x = rev(Correlation))
                 ) |>
          relocate(external_gene_name, .before = Correlation)-> genes_corr_df

      } else if ( mappable == FALSE){
        warning("These genes ID do not seem to be ENSEMBL gene IDs",
                "Skipping mapping")
      } else {
        stop("Can't figure out if the gene IDs are mappable...")
      }
    }
  }
  return(genes_corr_df)
}

#' Create a bed file from events in a vast-tools inclusion table 
#'
#' @param path Path to vast-tools inclusion file or a tab delimited file with the same first ID columns as vast-tools inclusion tables.
#' @param header Does the vast-tools inclusion file have a header? If `FALSE` consider first row as an event.
#' @param remove_chr If `TRUE` remove 'chr' from the first column of the bed files.
#' @param out_bed_name Name of output bed file, without the `.bed` extension. If missing use input file basename.
#' @param out_path Path to directory where to save bed file. If folder doesn't exist create it.
#' @param verbose If `TRUE` print now many lines are in the bed file.
#'
#' @return a sorted bed file without a header with EVENT ID in the 4th column and gene name in the 5th column.
#' @importFrom readr read_delim write_delim
#' @importFrom dplyr across arrange select
#' @importFrom stringr str_extract
#' 
#' @description
#' The input file (specified with path), doesn't need to be specifically a vast-tools output file. 
#' As long as the first five columns of the input file match the structure of a vast-tools inclusion table, this function will work.
#' 
#' @details
#' This function can also report the strand for exons, in the future I might implement this also for introns.
#'
#' @examples
#' inclusion_tbl2bed(path = vst_tbl_path, header = T, remove_chr = T, out_path = bed_dir)
inclusion_tbl2bed <- function(path, header = TRUE, remove_chr = FALSE, 
                              out_bed_name, out_path, verbose = TRUE) {
  
  if ( !file.exists(path) ) { stop('Input file does not exist!') }
  input_basename <- basename(path)
  
  if ( missing(out_bed_name) ) {
    out_bed_name <- str_extract(string = input_basename, pattern = '^*.*(?=\\.tab$)')
  }
  
  if ( missing(out_path) ) { stop('Specify an output directory!') }
  
  if ( !dir.exists(out_path) ) { dir.create(path = out_path, recursive = T) }
  
  input <- read_delim(file = path, col_select = c(3,2,1,5), progress = F,
                      delim = "\t", escape_double = FALSE, show_col_types = FALSE,
                      col_names = header, na = "empty", trim_ws = TRUE)  |>
    setNames(c('COORD', 'EVENT', 'GENE', 'FullCO')) 
  
  # if gene name is not there use the AS event ID
  input |>
    mutate(GENE = ifelse(test = GENE == "", yes = EVENT, no = GENE ) ) |>
    mutate(chr = str_extract(string = COORD, pattern = '^chr[A0-Z9]+(?=:)')) |>
    mutate(start = str_extract(string = COORD, pattern = '(?<=:)[0-9]+(?=-)')) |>
    mutate(end = str_extract(string = COORD, pattern = '(?<=-)[0-9]+$')) |>
    mutate(across(c(start, end), as.integer)) |>
    select(chr, start, end, EVENT, GENE, FullCO) -> input
  
  # extract strand for exons the others are just non-stranded
  indx_exons <- which(grepl(pattern = 'EX[0-9]+', x = input$EVENT))
  indx_others <- which(!grepl(pattern = 'EX[0-9]+', x = input$EVENT))
  # sanity check
  stopifnot( length(indx_exons) + length(indx_others) == nrow(input) )
  
  # TO DO: implement strand detection also for introns
  # For INT: chromosome:C1exon=C2exon:strand.
  ex_in <- input[indx_exons, ]
  
  # From FullCO column extract the C1donor and C2acceptor splice sites coord
  ex_in |>
    mutate(C1donor = str_extract(pattern = '(?<=:)[0-9]+', string = FullCO), .before = EVENT) |>
    mutate(C2acceptor = str_extract(pattern = '[0-9]+$', string = FullCO), .before = EVENT) |>
    mutate(across(c(C1donor, C2acceptor), as.integer)) |>
    # If the downstream acceptor is bigger than the upstream donor then it's on the positive strand.
    mutate(strand = ifelse(test = C1donor <= C2acceptor, yes = '+', no = '-'), .before = FullCO ) |>
    select(chr, start, end, EVENT, GENE, strand) -> ex_in
  
  # Add non-strand ('.') to others
  other_in <- input[indx_others, ]
  other_in$strand <- '.'
  other_in <- other_in |> select(chr, start, end, EVENT, GENE, strand)
  
  # sort
  rbind(ex_in, other_in) |> arrange(chr, start) -> input
  
  if ( remove_chr ) {
    input$chr <- gsub(pattern = '^chr', replacement = '', x = input$chr)
  }
  
  write_delim(x = input, delim = '\t', append = F, col_names = F, quote = 'none', 
              progress = F, escape = 'none', 
              file = file.path(out_path, paste0(out_bed_name, '.bed')) )
  
  if (verbose == TRUE) { 
    message('Saved a bed file with ', nrow(input), ' lines') 
  } else if (verbose == FALSE) { 
    # do nothing 
  } else { 
    stop('Parameter verbose must be a either TRUE or FALSE!') 
  }
}
