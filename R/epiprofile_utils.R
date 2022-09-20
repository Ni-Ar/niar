#' Auxiliary function to modify the EpiProfile columns names on the fly while importing data.
#'
#' @param col the columns to rename
#'
Rename_EpiProfile_colnames <- function(col) {
  # Fix sample names columns
  fixed <- gsub(pattern = "^[0-9]+,", "", col) 
  # Add a name to first column
  fixed[1] <- "First_Col"
  # find empty col names
  fixed[ which( grepl(x = fixed, pattern = "^$", perl = T )) ] <- "Empty"
  # Select cols that are not First_Col nor empty
  indx_Num_Col <- which( !grepl(x = fixed, pattern = "First_Col|Empty", perl = T ))
  # Identify the number of samples in the imported table
  num_samples <- length(indx_Num_Col) / 3
  fixed[indx_Num_Col] <- paste0(c(rep("Ratio_", num_samples), 
                                  rep("Area_", num_samples),
                                  rep("RT_", num_samples)),
                                fixed[indx_Num_Col])
  return(fixed)
}

#' Import EpiProfile version2.1 output file `histone_ratio.xls` into R
#'
#' @param histone_ratio_path Path to where the xls file is
#'
#' @return A tidy tibble
#' @import dplyr
#' @importFrom readr read_delim locale
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom rlang as_function
#' @importFrom vctrs vec_fill_missing
#' @export
#'
#' @examples
#' read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
#'     tidy_hPTMs()   
read_EpiProfile_histone_ratio <- function(histone_ratio_path) {
  
  # This could come from a txt file maybe.
  c("TKQTAR(H3_3_8)", "KSTGGKAPR(H3_9_17)", "KQLATKAAR(H3_18_26)", 
    "KSAPATGGVKKPHR(H3_27_40)", 'KSAPSTGGVKKPHR(H33_27_40)', 
    "YRPGTVALR(H3_41_49)", 'YQKSTELLIR(H3_54_63)', 'EIAQDFKTDLR(H3_73_83)', 
    'VTIMPKDIQLAR(H3_117_128)',
    # Histone H4 peptides
    'GKGGKGLGKGGAKR(H4_4_17)', 'KVLR(H4_20_23)', 
    "DNIQGITKPAIR(H4_24_35)", 'RGGVKR(H4_40_45)', "KTVTAMDVVYALKR(H4_79_92)",
    # Histone H1
    'KAAGGAKR(H14_25_32)', 'KASGPPVSELITKAVAASKER(H12_33_53)', 
    'KATGPPVSELITKAVSASKER(H15_33_53)', 'unmod(H1_1_35)', 'unmod(H1_54_81)', 
    # H2A
    'KGNYAER(H2A1_36_42)', 'KGNYSER(H2A3_36_42)', 'KGHYAER(H2AX_36_42)', 
    'GKQGGKAR(H2A1_4_11)', 'GKQGGKVR(H2AJ_4_11)', 'GKTGGKAR(H2AX_4_11)', 
    'SGRGKQGGKAR(H2A1_1_11)', 'AGGKAGKDSGKAKAKAVSR(H2AV_1_19)', 
    'AGGKAGKDSGKAKTKAVSR(H2AZ_1_19)', 'AKAKTR(H2A1_12_17)', 'AKAKSR(H2A3_12_17)',
    'DNKKTR(H2A1_72_77)', 'unmod(H2A_1_88)',  "HLQLAIR(H2A_82_88)",
    # H2B
    'unmod(H2B_1_29)', 'LAHYNKR(H2B1B_80_86)'
  ) -> list_of_peptides
  
  
  # Import data and fix col names during importing
  data <- read_delim(file = histone_ratio_path, quote = "", delim = "\t", 
                     escape_double = FALSE, locale = locale(decimal_mark = "."),
                     trim_ws = TRUE, show_col_types = FALSE, 
                     name_repair = as_function( ~ Rename_EpiProfile_colnames(col = .x) ),
                     col_select = c(!matches("Empty")) ) |>
    # remove first line that contains info included in the header
    filter(row_number() != 1) |> 
    # Create a new column with the peptide type for that line
    mutate(Peptide_Type = case_when(First_Col %in% list_of_peptides ~ First_Col)) |> 
    relocate(Peptide_Type, .after = First_Col ) |>
    mutate(Peptide_Type = vec_fill_missing(Peptide_Type, direction = c("down"), max_fill = NULL ) ) |> 
    subset(First_Col != Peptide_Type) |> 
    # Turn data into long format and use positive lookbehind regex for splitting column names
    pivot_longer(cols = c(starts_with("Ratio"), starts_with("Area"), starts_with("RT") ), 
                 names_to = c("Feature", "Sample"),
                 names_sep = '(?<=Ratio|Area|RT)_',
                 values_to = "Value") |>
    # Coerc Value column to numeric
    mutate(Value = as.numeric(Value)) |>
    # Turn data into wide format to have Ratio Area and RT in different cols
    pivot_wider(id_cols = c("First_Col", "Peptide_Type", "Sample"), 
                names_from = Feature, values_from = "Value")
  return(data)
}

#' Parse the information about histone PTMs from an EpiProfile output import
#'
#' @param data A tibble imported with `read_EpiProfile_histone_ratio()`.
#' @param split_double_PTMs Logical. Split the modification on peptides with 2 lysine co-occurring PTMs in `Sub_PTM1` and `Sub_PTM2` new columns. 
#' @param max_nchar_PTM Integer. Longest character length of the modification text in column `PTM`.
#'
#' @return A tidy tibble
#' @import dplyr
#' @import stringr
#' @export
#'
#' @examples
#' read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
#'     tidy_hPTMs()   
tidy_hPTMs <- function(data, split_double_PTMs = FALSE, max_nchar_PTM = 20) {
  
  data |>
    # Parse the Histone PTM info
    mutate(Histone = gsub("_[0-9]+_[0-9]+\\s", " ", First_Col)) |>
    mutate(Histone = str_split_fixed(string = Histone, pattern = " ", n = 2)[,1]) |>
    mutate(Modification = str_split_fixed(string = First_Col, pattern = " ", n = 2)[,2]) |>
    # Parse the Histone Peptide info
    mutate(Peptide_Boundaries = str_extract(string = Peptide_Type, pattern = "(?<=\\()[A-Z].*(?=\\))") ) |>
    mutate(Peptide_Boundaries = str_extract(string = Peptide_Boundaries, pattern = "_[0-9]+_[0-9]+(?=$)") ) |>
    mutate(Peptide_Boundaries = str_replace_all(string = Peptide_Boundaries, pattern = "_", " ") ) |>
    mutate(Peptide_Start = as.integer(str_extract(string = Peptide_Boundaries, pattern = "(?<=^\\s)[0-9]+") ) ) |>
    mutate(Peptide_End = as.integer(str_extract(string = Peptide_Boundaries, pattern = "(?<=\\s)[0-9]+(?=$)") ) ) |>
    mutate(Peptide_Sequence = str_remove(string = Peptide_Type, pattern = "\\([A-Z].*\\)") ) |>
    relocate(Peptide_Sequence, Peptide_Start, Peptide_End, .after = Peptide_Type) |>
    relocate(Histone, Modification, .after = First_Col) |>
    select(!c(Peptide_Boundaries)) |>
    # Label better Histone variants like H3.3
    mutate(Histone = gsub("^H33", "H3.3", Histone ) ) %>%
    mutate(PTM = case_when(Modification == "unmod" ~ paste0(Histone, "un ", Peptide_Start, "-", Peptide_End),
                           Modification != "unmod" ~ paste0(Histone, Modification) ) )  %>%
    mutate(PTM = str_trunc(string = PTM, width = max_nchar_PTM, side = "right" )) %>%
    relocate(PTM, .after = Modification) -> data
  
  # Split peptides with 2 modified lysines into 2 new columns called:
  #  - Sub_PTM1 with the first modification 
  #  - Sub_PTM2 with the second modification
  if(split_double_PTMs) {
    ## Regex break-down:
    # ^K[1-9]+[m|a][e|c]([1-3]?)
    # '^K' starts with K
    # '[1-9]+' lysine residue after the K is one number or more 2 numbers.
    # '[m|a][e|c]' me or ac modification
    # '([1-3]?)' me1-2-3 or ac (with no number). ? = zero or 1 occurrance of a number
    # repeated 2 times to get only the modifications with 2 K co-modifications and ends with $.
    rgx_1st_K <- "^K[1-9]+[m|a][e|c]([1-3]?)"
    rgx_2nd_K <- "K[1-9]+[m|a][e|c]([1-3]?)$"
    
    data |>
      pull(Modification) |>
      grep(pattern = paste0(rgx_1st_K, rgx_2nd_K), value = T, perl = T) |> 
      unique() -> dbl_cmb_PTM
    
    data |>
      mutate(Sub_PTM1 = case_when(
        Modification %in% dbl_cmb_PTM ~ str_extract(pattern = rgx_1st_K, string = Modification)
      )) |>
      mutate(Sub_PTM2 = case_when(
        Modification %in% dbl_cmb_PTM ~ str_extract(pattern = rgx_2nd_K, string = Modification)
      )) |>
      relocate(Sub_PTM1, Sub_PTM2, .after = PTM) -> data
  } 
  return(data)
}

#' Calculate summary stats across all samples for each hPTM and calculate hPTM are log2 Fold Change relative to control sample
#'
#' @param data Tibble imported with `read_EpiProfile_histone_ratio()`.
#' @param Pseudocount Number. By default adds 1 pseudo count to the `Area` values to avoid diving by zero.
#' @param Reference_Sample_Name name of the reference sample to be used for area log2 FC (e.g. WT/untreated condition)
#'
#' @return A tibble
#' @import dplyr
#' @export
#' 
#' @examples
#'read_EpiProfile_histone_ratio("path/to/histone_ratio.xls") |>
#'     tidy_hPTMs() |>
#'     get_hPTM_FC(Reference_Sample_Name = 'WT_1_DIA')
get_hPTM_FC <- function(data, Pseudocount = 1, Reference_Sample_Name) {
  data |>
    # Add Area pseudocount
    mutate(Area = Area + Pseudocount) |>
    # Calculate summary stats for Area and RT across all samples
    group_by(First_Col, Peptide_Type) |>
    mutate( across(.cols = c("Area","RT"), ~ mean(.x, na.rm = TRUE), .names = '{.col}_mean' )  ) |>
    mutate( across(.cols = c("Area", "RT"), ~ sd(.x, na.rm = TRUE), .names = '{.col}_sd' )  ) |> 
    # Calculate log2 FC of Area to reference sample
    mutate(Area_log2FC = log2( (Area ) / (Area[Sample == Reference_Sample_Name] ) ) ) -> data
  return(data)
}

