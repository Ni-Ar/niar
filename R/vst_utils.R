#' Reshape a wide vast-tools inclusion table into a tidyverse-friendly long table. This function works with inclusion tables generated both by vast-tools combine or compare. In the latter case, if there's a 'dPSI' column it will be automatically detect and included in the reshaped table.
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
#' @return A reshaped data.frame in long format as tibble.
#' @import dplyr
#' @importFrom magrittr extract extract2 %>%
#' @import stringr
#' @import tidyr
#' 
#' @export
#'
#' @examples
#' grep_psi(path_to_vst_PSI_tbl, vst_id) %>%
#'     tidy_vst_psi()
tidy_vst_psi <- function(vst_psi_tbl, num_id_cols = 6, num_of_Score_Types = 5,
                         quality_col_suffix = "-Q", return_quality_scores = TRUE,
                         return_S1_only = TRUE, verbose = TRUE, add_ID_col = FALSE,
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
      # This is not a vast-tools compare table with a dPSI column at the end
      vst_compare_tbl <- FALSE
    }
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
    vst_psi_tbl %>% 
      extract2(colnames(vst_psi_tbl)[Qual_cols][i]) %>% 
      str_split_fixed(string = ., pattern = ',', n = num_of_Score_Types) %>% 
      as.data.frame() %>% 
      setNames(., paste(colnames(vst_psi_tbl)[Qual_cols][i], 
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
               values_to = "Quality_Score_Value") %>%
    # Change cols to make it compatible with joining
    dplyr::mutate(Quality_Score_Type = gsub(".*Q_", "", Quality_Score, perl = T)) %>%
    dplyr::mutate(Sample = str_remove(Quality_Score, "\\-Q_S[1-5]$")) %>% 
    dplyr::select(-Quality_Score) -> lng_vst_psi_quality_tbl
  
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
  
  # Coerc factors col to character col
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
    select(tidy_vst_psi_tbl, cols_to_return) %>%
      unique() ->  tidy_vst_psi_tbl
  } else {
    stop("return_quality_scores must me a logical (TRUE or FALSE)!")
  }

  # Add some sort of grouping/experiment name to the reshaped dataset
  if ( add_ID_col == TRUE ) {
    tidy_vst_psi_tbl$Group_Name <- col_ID_name
  }
  return(tidy_vst_psi_tbl)
}

#' Reshape a wide vast-tools expression table to a long format
#'
#' @param data A dataframe generated with `grep_gene_expression` from a vast-tools gene expression table.
#' @param expression_unit Character: either `TPM` or `cRPKM` based on the gene expression unit in `data`.
#' @param ID_cols The columns names to use as IDs in `data`
#'
#' @return
#' @import tidyr
#' @export
#'
#' @examples
#' grep_gene_expression(path_to_vst_Expression_tbl, ensembl_gene_id) %>%
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
  
  # suppressMessages( require('tidyr') )
  
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
#' @param vastid An alternative splicing event name from vast-tools.
#' @param latin_name Logical, whether or not you want the scientific name of the species.
#'
#' @return A character
#'
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' guess_species(vastid = "HsaEX000001")
#' # Human
guess_species <- function(vastid, latin_name = F) {
  
  # suppressMessages( require('stringr') )
  guessed_species <- stringr::str_extract(string = vastid, 
                                          pattern = "^[A-Z][a-z][a-z]")
  
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
  
  return(guessed_species)
}

#' Function to grep the PSI values of an as event from a vast-tools output table. 
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
#' #' grep_psi(path_to_vst_tbl, vst_id) %>%
#'     tidy_vst_psi()
grep_psi <- function(inclusion_tbl, vst_id, tmp_dir = tempdir(), verbose = TRUE, 
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
  if ( missing(vst_id) ) { stop("You didn't specified a vastID!") } 
  
  # If there's more than one vastID in the query
  if ( length(vst_id) >= 2 ) {
    
    # Turn the vastID into a "greppable" multiple vastID
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
                            tmp_dir, "/tmp_vastID_PSI.tab")
    system(command = get_exon_line)
    
    
  } else if ( split_temp == FALSE ) {
    # 1, 2 -- Get table header and AS event PSI
    get_header_n_exon_line <- paste0( "cat <(head -n1 ", inclusion_tbl,
                                      ") <(", grep_command, "'", vst_id, "' ", 
                                      inclusion_tbl, ") > ", tmp_dir,
                                      "/tmp_vastID_PSI.tab")
    
    # the R system() command calls the bourne shell (sh) instead of the bourne again shell (bash)
    # to overcome this print the command in the bourne shell and pipe the output into bash
    bash_header_AS_line <- paste('echo "', get_header_n_exon_line, '" | bash')
    
    system( command = bash_header_AS_line )
    
  } else {
    stop("split_tmp must be logical, either TRUE or FALSE")
  }
  
  # 3 -- Check that temp file is not empty
  exon_PSI_tmp_info <- file.info(file.path(tmp_dir, "tmp_vastID_PSI.tab"))
  
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
    tmp_exon <- read.delim(file = file.path(tmp_dir, "tmp_vastID_PSI.tab"),
                           header = F, check.names = F, stringsAsFactors = F)
    
    # 6 -- Merge header and exon processed vast table
    colnames(tmp_exon) <- colnames(tmp_header)
    vst_psi_tbl <- rbind(tmp_header, tmp_exon)
    
  } else if ( split_temp == FALSE ) {
    # 4, 5, 6 -- Read in AS event vast table 
    vst_psi_tbl <- read.delim(file = file.path(tmp_dir, "tmp_vastID_PSI.tab"),
                              header = T, check.names = F, stringsAsFactors = F)
    
  } else {
    stop("Something went wrong when reading tmp_vastID_PSI.tab in\n", tmp_dir)
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
    
    file.remove( file.path(tmp_dir, "tmp_vastID_PSI.tab") )
    if ( verbose ) { message("Temporary files have been removed.") } 
  }
  # 9 -- Return a tbl
  return(vst_psi_tbl)
}

#' Function to grep the expression values of a gene from a vast-tools output table. 
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
#' grep_gene_expression(path_to_vst_Expression_tbl, ensembl_gene_id) %>%
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
    if ( verbose ) { message("First column is an ensebl gene id") }
    colnames(vst_expression_tbl)[1] <- "ensembl_gene_id"
  }
  
  # 7 -- Return a tbl
  return( vst_expression_tbl )
}


#' Wrapper for readr::read_delim for a tab delimited table. Fastest way to read this kind of data into `R`.
#'
#' @param path A character vector providing the 'path/to/the/table/to/read'.
#' @param verbose Logical. Use `TRUE` to know the table file size being read into `R`.
#' @param ... extra info to pass to `read_delim`.
#'
#' @return A tibble
#' @importFrom readr read_delim
#' @export
#'
#' @examples
#' read_vst_tbl(path = "location/of/tab/delimited/file/inclusion/table.tab") %>%
#'    tidy_vst_psi()
#'  
#'    
read_vst_tbl <- function(path, verbose = FALSE, ...) {
  if (!file.exists(path)) { stop('Cannot find vast-tools table!\n') }
  
  if (verbose) {
    message('Reading a vast-tools table of: ',
            round(file.info(cmpr_vst_tbl)$size / 10e3, 2), ' Kbytes')
  }
  return(read_delim(file = path, delim = '\t', col_names = T,
                    locale = locale(decimal_mark = "."), na = "NA", ...) )
}

