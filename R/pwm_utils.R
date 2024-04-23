#' Helper function to use as a tidyverse verb
#'
#' @param column 
#' @param alphabet 
#' @param verbose Logical. Whether or not to print warning messages about case-sensitivity. Default F
#'
#' @return a column frequency
#' @export
#'
col_freq <- function(column, alphabet, verbose = FALSE) {
  
  # check case of the alphabet
  if ( all(grepl(pattern = "[A-Z]", x = alphabet)) ) {
    alphabet_case <- 'UPPER'
  } else if ( all(grepl(pattern = "[a-z]", x = alphabet)) ){
    alphabet_case <- 'lower'
  } else if ( all( c( any(grepl(pattern = "[a-z]", x = alphabet)), any(grepl(pattern = "[A-Z]", x = alphabet) ) ) ) ) {
    alphabet_case <- 'MiXeD'
  } else if ( all(grepl(pattern = "[A-Z]|\\-", x = alphabet)) ){
    alphabet_case <- 'UPPER_gapped'
  } else if ( all(grepl(pattern = "[a-z]|\\-", x = alphabet)) ){
    alphabet_case <- 'lower_gapped'
  } else {
    warning("I'm not sure I understand the case of the alphabet.")
  }
  
  # check case of the input data
  input_column_chars <- unlist(strsplit(column, split = ""))
  if ( all(grepl(pattern = "[A-Z]", x = input_column_chars)) ) {
    column_case <- 'UPPER'
  } else if ( all(grepl(pattern = "[a-z]", x = input_column_chars)) ){
    column_case <- 'lower'
  } else if ( all( c( any(grepl(pattern = "[a-z]", x = input_column_chars) ), any( grepl(pattern = "[A-Z]", x = input_column_chars) ) ) ) ) {
    column_case <- 'MiXeD'
  } else if ( all(grepl(pattern = "[A-Z]|\\-", x = input_column_chars)) ){
    column_case <- 'UPPER_gapped'
  } else if ( all(grepl(pattern = "[a-z]|\\-", x = input_column_chars)) ){
    column_case <- 'lower_gapped'
  } else {
    warning("I'm not sure I understand the case of the input data letters")
  }
  
  if ( alphabet_case != column_case) {
    if ( verbose == TRUE) {
      warning("Case-sensitive warning! Alphabet is in ", alphabet_case, 
              "case, the input data is in ", column_case, "case.")
    }
  }
  # Get frequencies for each letter
  letter_freq <- table(column)
  # get index position of letter frequencies as they appear in the alphabet
  indx <- match(alphabet, names(letter_freq) )
  # sort frequencies according to the alphabet order
  letter_freq <- letter_freq[indx]
  # Set NA missing frequencies to zero
  letter_freq[is.na(letter_freq)] <- 0
  # coerc to integer numbers
  letter_freq <- as.integer(letter_freq)
  # set the names to the letters in the alphabet (probably redundant)
  names(letter_freq) <- alphabet
  return(letter_freq)
}

#' Turn a dataframe into a matrix
#'
#' @param data 
#' @param ID_col 
#'
#' @return a matrix
#' @export
df2mat <- function(data, ID_col) {
  
  # -- 1 -- check that all required columns are there.
  data_required_cols <- c(ID_col, "Sequence")
  if ( !all(data_required_cols %in% colnames(data) ) ) {
    
    missing_col <- which( ! data_required_cols %in% colnames(data) )
    stop("Input dataframe is missing the required column(s): ", 
         paste0(data_required_cols[missing_col], collapse = " ") )
  }
  
  # -- 2 -- check that the ID columns are unique
  if( length(data$ID_col) != length(unique(data$ID_col)) ) {
    warning('The names in the "ID_col" are NOT unique!')
  }
  
  # -- 3 -- Check that Sequence length is the same for all rows
  seq_len_uniq <- unique(nchar(data$Sequence))
  
  if ( length(seq_len_uniq) >= 2) { 
    stop('The sequences in ID_col are NOT of the same length! ',
         'The sequences have a length of: ', 
         paste0(seq_len_uniq, collapse = ', '), '. ',
         'Fix your input and try again.')
  }
  
  # -- 4 -- Turn the dataframe into a matrix, splitting every letter into a column
  data |>
    select(Sequence) |>
    # split every letter in the Sequence as an individual column of a matrix
    mutate(str_split_fixed(string = Sequence, pattern = '', n = seq_len_uniq)) |> 
    # this dataframe as 2 columns: the original Sequence and N which is the matrix itself
    setNames(nm = c('Sequence', 'N' ) ) |>
    # extract the matrix
    pull(N) -> mat
  
  # -- 5 -- Annotate the matrix
  colnames(mat) <- paste0('N', 1:seq_len_uniq)
  rownames(mat) <- data |> pull(ID_col)
  return(mat)
}

#' Handy wrapper to grep regexes on multiple vectors
#'
#' @param chrctr A character vector
#' @param regex A character with PERL compatible regexes, case sensitive.
#'
#' @return A named logical vector
#' @export
test_regexes <- function(chrctr, regex) {
  
  stopifnot(is.character(chrctr))
  stopifnot(is.character(regex))
  
  logic_res <- grepl(x = chrctr, pattern = regex, perl = TRUE, ignore.case = FALSE)
  names(logic_res) <- names(chrctr)
  return(logic_res)
}

#' Guess an alphabet type in a dataframe sequences.
#'
#' @param data A data.frame with a minimum of 2 columns. One named `Sequence`, the other named as you prefer that will be specified with `ID_col`.
#' @param ID_col The name of the column in `data` to be used as the identifier of the `Sequence` column.
#' @param include_gaps Logical, include also the gap character '-' (dash) in the guess alphabet regex. Default `TRUE`.
#' @param name Logical, return only the guessed alphabet name instead of the letters in that alphabet. Default `FALSE`.
#' @param verbose Logical, print some info messages. Default `FALSE`.
#'
#' @return A character vector with the letters of that alphabet
#' @export
#' 
#' @description It can can recognise sequences that are DNA/RNA; upper, lower, or both UPPER&lower case; gapped or non-gapped.
#' 
guess_alphabet <- function(data, ID_col, include_gaps = TRUE, name = FALSE, 
                           verbose = FALSE) {
  
  # -- 1 -- Define the regexes to search. Add proteins in the future
  regex_DNA_UPPER <- '[ACGT]+'
  regex_RNA_UPPER <-  gsub(pattern = 'T', replacement = 'U', x = regex_DNA_UPPER)
  
  regex_DNA_lower <- tolower(regex_DNA_UPPER)
  regex_RNA_lower <- tolower(regex_RNA_UPPER)
  
  # UPPER and lower case regex
  # The look ahead ensures that there is at least one UPPER and lowercase letter in any order in the DNA
  # regex_DNA_UPlowER <- '(?=[acgt]+[ACGT]+|[ACGT]+[acgt]+)[aAcCgGtT]+'
  # regex_RNA_UPlowER <- '(?=[acgu]+[ACGU]+|[ACGU]+[acgu]+)[aAcCgGuU]+'
  regex_DNA_UPlowER <- '(?=.*[acgt])(?=.*[ACGT])[aAcCgGtT]+'
  regex_RNA_UPlowER <- '(?=.*[acgu])(?=.*[ACGU])[aAcCgGuU]+'
  
  if (include_gaps == TRUE) {
    # add the gap character (dash '-') to the letters to search.
    regex_DNA_UPPER_gapped <- '[ACGT\\-]+'
    regex_RNA_UPPER_gapped <-  gsub(pattern = 'T', replacement = 'U', x = regex_DNA_UPPER_gapped)
    
    regex_DNA_lower_gapped <- tolower(regex_DNA_UPPER_gapped)
    regex_RNA_lower_gapped <- tolower(regex_RNA_UPPER_gapped)
    
    regex_DNA_UPlowER_gapped <- '(?=.*[actg])(?=.*[ACTG])(?=.*\\-)[aAcCgGtT\\-]+'
    regex_RNA_UPlowER_gapped <- '(?=.*[acug])(?=.*[ACUG])(?=.*\\-)[aAcCgGuU\\-]+'
  }

  # -- 2 -- Make a vector with the content of all variables that start with 'regex_'
  regexes <- ls(pattern = "^regex_")
  alphabets_to_search <- sapply(regexes, function(x) { get(x) }, simplify = TRUE )

  names_alphabets_to_search <- gsub(x = names(alphabets_to_search), 
                                     pattern = "^regex_", replacement = "")
  
  # change name of gapped regexes
  names_alphabets_to_search <- gsub("_gapped", "\\-gapped", names_alphabets_to_search)
  
  # make regex search the whole line
  alphabets_to_search <- paste0('^', alphabets_to_search, '$')
  names(alphabets_to_search) <- names_alphabets_to_search
  
  # -- 3 -- Turn the dataframe into a letter matrix
  df2mat(data, ID_col) |>
    # concatenate together the columns of the letter matrix
    apply(2, paste, collapse = "") -> col_letter
  
  # collapse together all letters (by colum) into one single string.
  all_letters_by_col <- paste0(col_letter, collapse = '')
  
  # -- 4 -- Guess the alphabet
  # apply each regex in the list 'alphabets_to_search' to the concatenated string
  # same as: sapply(alphabets_to_search, function(x) { test_regexes(column_characters, regex = x)} )
  
  guess <- sapply(X = alphabets_to_search, test_regexes, chrctr = all_letters_by_col) |> 
    # summaries the logical results column-wise with all.
    # apply(2, all) |> 
    # select only the TRUE regex 
    which() |> names()
  
  # -- 5 -- Make sure only one guess is found
  if ( length(guess) != 1 ) {
    if ( verbose ) {
      message("I found more than one of the default alphabets could match the ", 
              "input sequences! Namely: \n", 
              paste0(alphabets_to_search[guess], collapse = ', '), 
              "\nMaybe specify the alphabet yourself?")
    }
    
    alphabet_subtype <- unique(gsub(pattern = "^[D|R]NA_", replacement = "", guess))
    seq_type <- unique(gsub(pattern = "_UP*.*$", replacement = "", guess))
    
    if ( seq_type == 'DNA' ) { 
      if(verbose) { message("It's a DNA seq type alphabet")  }
    } else if (seq_type == 'RNA' ) {
      if(verbose) { message("It's a RNA seq type alphabet") }
    } else {
      if(verbose) { message("I can't even understand if it's a DNA or RNA seq type alphabet!") }
    }
    
    if ( any(grepl(pattern = "\\-gapped", x = alphabet_subtype)) ) {
      if(verbose) { message("It is likely a ", seq_type, " gapped alphabet") }
      case_subtype <- unique(gsub(pattern = "\\-gapped", replacement = "", alphabet_subtype))
      if( length(case_subtype) > 1) {
        if(verbose) { message('The alphabet case is either ', paste0(case_subtype, collapse = " or ") ) }
      }
    }
    # -- for simplicity's sake select the last guessed alphabet type
    if(verbose) { message('Trying ', guess[length(guess)], ' as a wild guess...') }
    guess <- guess[length(guess)]
  } else {
    if(verbose) { message('The found alphabet is: ', guess) }
  }
  
  # -- 6 -- Set the letters based on the guessed alphabet
         if ( guess == 'DNA_UPPER' ) { alphabet <- c('A', 'C', 'G', 'T') 
  } else if ( guess == 'DNA_lower' ) { alphabet <- c('a', 'c', 'g', 't') 
  } else if ( guess == 'DNA_UPlowER' ) { alphabet <- c('A', 'C', 'G', 'T', 'a', 'c', 'g', 't') 
  } else if ( guess == 'RNA_UPPER') { alphabet <- c('A', 'C', 'G', 'U') 
  } else if ( guess == 'RNA_lower' ) { alphabet <- c('a', 'c', 'g', 'u') 
  } else if ( guess == 'RNA_UPlowER' ) { alphabet <- c('A', 'C', 'G', 'U', 'a', 'c', 'g', 'u')
  } else if ( guess == 'DNA_UPPER-gapped' ) { alphabet <- c('A', 'C', 'G', 'T', '-') 
  } else if ( guess == 'DNA_lower-gapped' ) { alphabet <- c('a', 'c', 'g', 't', '-') 
  } else if ( guess == 'DNA_UPlowER-gapped' ) { alphabet <- c('A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '-') 
  } else if ( guess == 'RNA_UPPER-gapped') { alphabet <- c('A', 'C', 'G', 'U', '-') 
  } else if ( guess == 'RNA_lower-gapped' ) { alphabet <- c('a', 'c', 'g', 'u', '-') 
  } else if ( guess == 'RNA_UPlowER-gapped' ) { alphabet <- c('A', 'C', 'G', 'U', 'a', 'c', 'g', 'u', '-')
  } else {
    stop("I wasn't able to guess the type of alphabet of the input sequences. ",
         "Please specify the alphabet! ")
  }
  
  # -- 7 -- Decide what to return
  if (name == TRUE) { return(guess) } else { return(alphabet) }
}

#' Turn a data.frame with sequences into a positional frequency matrix. 
#'
#' @param data A data.frame with a minimum of 2 columns. One named `Sequence`, the other named as you prefer that will be specified with `ID_col`.
#' @param ID_col The name of the column in `data` to be used as the identifier of the `Sequence` column.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default.
#' 
#' @return A data frame with positional frequencies.
#' @importFrom dplyr summarise across everything
#' @export
#' 
#' @description Given a dataframe with a column called 'Sequence' it transforms the sequences into a matrix where the frequency of each letter in the alphabet is counted at every position
#' @details For a great song click \href{https://youtu.be/tqvvaY2LvuI?t=110}{here}.
#' 
df2pfm <- function(data, ID_col, alphabet) {
  
  if( missing(alphabet) ) { alphabet <- guess_alphabet(data, ID_col) }
  
  pfm <- data |>
    df2mat(ID_col) |>
    as.data.frame() |>
    # calculate for each column the frequency of each letter in the alphabet
    summarise(across(.cols = everything(), .fns = ~col_freq(column = .x, alphabet) ) ) 
  
  # set the alphabet as row names of the position frequency matrix
  rownames(pfm) <- alphabet
  return(pfm)
}

#' Turn a data.frame with sequences into a positional probability matrix. 
#'
#' @param data A data.frame with a minimum of 2 columns. One named `Sequence`, the other named as you prefer that will be specified with `ID_col`.
#' @param ID_col The name of the column in `data` to be used as the identifier of the `Sequence` column.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default.
#' @param c an integer used to round up the numbers 
#'
#' @return A data frame with positional probabilities.
#' @export
#'
#' @description Generates a positional probability matrix (PPM).
#'
#' @examples
#' df2ppm(data, ID_col = 'Species', alphabet = c("A", "C", "G", "T")) 
df2ppm <- function(data, ID_col, alphabet, c = 5) {
  
  if( missing(alphabet) ) { alphabet <- guess_alphabet(data, ID_col) }
  
  df2pfm(data, ID_col, alphabet) |>
    # calculate the relative frequency, aka the probability
    summarise(across(.cols = everything(), .fns = ~.x/sum(.x) ) ) -> ppm
  # set the alphabet as row names of the position probability matrix
  rownames(ppm) <- alphabet
  # round to c numbers after the comma
  ppm <- signif(ppm, c)  
  return(ppm)
}

#' Create a PWM from a dataframe with sequences
#'
#' @param data A data.frame with a minimum of 2 columns. One named `Sequence`, the other named as you prefer that will be specified with `ID_col`.
#' @param ID_col The name of the column in `data` to be used as the identifier of the `Sequence` column.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default.
#' @param bg_prob_letter Explain well in details.
#' @param pseudocount_letter Explain well in details.
#' @param long_format Logical. If `TRUE` reshape the PWM into a tidy long data.frame format. Default `FALSE`.
#'
#' @return A data.frame or a tidy long format data.frame
#' @importFrom dplyr mutate arrange relocate
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export
#' @description Generate a positional weight matrix (PWM) from a set of sequences. This function comes with batteries included.
#' @details Here I'll explain more things, ideally with formulas if LaTeX syntax was supported.
#'
#' @examples
#' df2pwm(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), long_format = T) 
df2pwm <- function(data, ID_col, alphabet, bg_prob_letter, 
                   pseudocount_letter, long_format = FALSE) {
  
  # -- 1 -- check input params
  if ( !is.data.frame(data) ) { stop("'data' must be a data.frame!") }
  
  if ( !is.logical(long_format)) { 
    stop("'long_format' must be a logical, either TRUE or FALSE")
  }
  
  if( missing(alphabet) ) { alphabet <- guess_alphabet(data, ID_col) }
  
  # -- 2 -- check background probabilities
  # If background probabilities for each letter in the alphabet are not defined
  # Set them to equal background probability.
  if( missing(bg_prob_letter) ) {
    bg_prob_letter <- rep(1/length(alphabet), length(alphabet))
    names(bg_prob_letter) <- alphabet
  }
  
  if( length(bg_prob_letter) == length(alphabet) ) {
    # If background probabilities are instead defined 
    indx_bg <- match(alphabet, names(bg_prob_letter) )
    # order them according to the alphabet order
    bg_prob_letter <- bg_prob_letter[indx_bg]
  } else if( length(bg_prob_letter) == 1 ) {
    # If only one background probability is set
    # distribute it uniformly across the alphabet
    bg_prob_letter <- rep(bg_prob_letter, length(alphabet))
    names(bg_prob_letter) <- alphabet
  } else {
    stop("Something is wrong with the 'bg_prob_letter'...?" )
  }
  
  # Check that background probabilities sum up to 1.
  if ( round( sum(bg_prob_letter), 4) != 1 ) {
    stop('The sum of the background probabilities for all letters in the alphabet is not 1!')
  }
  
  # -- 3 -- Check pseudo-counts
  num_seq_input <- nrow(df2mat(data, ID_col = ID_col))
  
  # If pseudo counts for each letter in the alphabet are not defined
  if ( missing(pseudocount_letter) ) {
    pseudocount_letter <- ( sqrt(num_seq_input) * 1/length(alphabet) ) 
  } 
  
  if( length(pseudocount_letter) == length(alphabet) ) {
    # If pseudo counts are instead defined 
    indx_pc <- match(alphabet, names(pseudocount_letter) )
    # order them according to the alphabet order
    pseudocount_letter <- pseudocount_letter[indx_pc]
  } else if ( length(pseudocount_letter) == 1 ) {
    # Distribute the pseudocounts uniformly across the alphabet
    pseudocount_letter <- rep(pseudocount_letter, length(alphabet))
    names(pseudocount_letter) <- alphabet
  } else {
    stop("Something is wrong with the 'bg_prob_letter'...?" )
  }
  
  #  -- 4 -- Calculate the summation of the pseudocounts
  sum_pseudocounts <- sum( rep(pseudocount_letter, length(alphabet)) )
  
  # -- 5 -- Calculate c(b,i) + s(b) / N + ∑s(b′)
  pfm <- df2pfm(data, ID_col, alphabet)
  corrected_probability <- ( pfm + pseudocount_letter ) / ( num_seq_input + sum_pseudocounts )
  # -- 6 -- Calculate PWM
  pwm <- log2( corrected_probability / bg_prob_letter )
  
  # -- 7 -- Reshape data into long format if needed.
  if( long_format == TRUE) {
    pwm |>
      rownames_to_column('Alphabet') |>
      pivot_longer(cols = starts_with("N"), names_to = "Position", 
                   values_to = "Weights") |>
      mutate(Position = as.integer( gsub("N", "", Position ) ) ) |>
      relocate(Position, .before = Alphabet) |> arrange(Position) |>
      mutate(Power_2_Weights = 2^Weights) -> tidy_pwm
    return(tidy_pwm)
  } else {
    return(pwm)
  }
}

#' Visualise a positional weight matrix as a heatmap
#'
#' @param long_df A data.frame obtained with `df2pwd(long_format = TRUE)`
#' @param cell_txt_size Text size of the heatmap cells that corresponds to the letter weight in that position. Default 3.
#' @param axis_txt_size Text size of axis labels. Default 10.
#'
#' @return A ggplot plot
#' @import ggplot2
#' @export
#'
#' @examples
#' df2pwm(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), long_format = T) |>
#'   plot_pwm_as_hm() 
plot_pwm_as_hm <- function(long_df, cell_txt_size = 3, axis_txt_size = 10) {
  
  num_pos <- max(long_df$Position)
  len_alphabet <- length(unique(long_df$Alphabet))
  
  ggplot(long_df)+
    aes(x = Position, y = Alphabet, fill = Weights ) +
    geom_tile( colour = "black", linewidth = 0.15) +
    geom_text(aes(label = round(Weights, 2) ), size = cell_txt_size, family = 'Arial' ) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick3", midpoint = 0) +
    scale_x_continuous(n.breaks = num_pos,
                       expand = expansion(mult = 0, add = 0) ) +
    coord_fixed(ratio = ((len_alphabet * 4) -1) / num_pos ) +
    guides(fill = guide_colorbar(barheight = unit(2, 'mm'), barwidth = num_pos * 0.75) ) +
    labs(x = "Sequence position", y = "Alphabet") +
    theme_minimal(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(panel.grid = element_blank(), 
          plot.background = element_blank(),
          legend.position = 'bottom', 
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.margin = margin(t = -1, unit = "mm"),
          legend.title = element_text(vjust = 1),
          axis.line = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(vjust = 1, margin = margin(t = -0.5, unit = "mm")),
          axis.text.y = element_text(hjust = 1, margin = margin(r = -0.25, unit = "mm")),
          axis.ticks =  element_blank()
    ) -> p_heatmap_pwm
  return(p_heatmap_pwm)
}

#' Visualise a PWM as a stacked barplot
#'
#' @param long_df A data.frame obtained with `df2pwd(long_format = TRUE)`
#' @param axis_txt_size axis_txt_size Text size of axis labels. Default 10.
#' @param alphabet_palette A colour palette vector with names matching the PWM alphabet
#'
#' @return A ggplot2 plot
#' @import ggplot2
#' @export
#' @description This function is useful to visually explore the different contribution of pseudocounts and weighted background probabilities to a PWM. See exaples.
#'
#' @examples
#' df2pwm(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), long_format = T) |>
#'   plot_pwm_as_bars() + ggtitle('uniform bg prob & default pseudocounts')
#'   
#' df2pwm(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), long_format = T, pseudocount_letter = 0.05) |>
#'   plot_pwm_as_bars() + ggtitle('Uniform bg prob & small pseudocounts')
#'   
#' mm10_bg_prob <- c('a' = 0.2917, 't' = 0.2917, 'c' = 0.2083, 'g' = 0.2083 )
#' df2pwm(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), long_format = T, 
#'        pseudocount_letter = 0.05, bg_prob_letter = mm10_bg_prob ) |>
#'        plot_pwm_as_bars() + ggtitle('mm10 bg prob & small pseudocounts')
#' 
#' # Get the alphabet and build a palette        
#' alphabet <- guess_alphabet(data = data, ID_col = 'Species')
#' plt <- palette.colors(n = length(alphabet))
#' names(plt) <- alphabet
#' 
#' df2pwm(data, ID_col = 'Species', long_format = T) |>
#'    plot_pwm_as_bars(alphabet_palette = plt)
plot_pwm_as_bars <- function(long_df, axis_txt_size = 10,
                             alphabet_palette = c('t' = "forestgreen", 
                                                  'c' = "dodgerblue",
                                                  'g' = "orange",
                                                  'a' = "firebrick3" ) ) {
  
  if ( ! all(names(alphabet_palette) %in% long_df$Alphabet) ){
    stop("The letters in the alphabet palette do not match the letters in the pwm!")
  }
  
  num_pos <- max(long_df$Position)
  len_alphabet <- length(unique(long_df$Alphabet))
  
  ggplot(long_df) +
    aes(x = Position, y = Power_2_Weights, fill = Alphabet, group = Position ) +
    geom_col(colour = "black", linewidth = 0.15, width = 0.85, position = position_stack() ) +
    scale_fill_manual(values = alphabet_palette ) +
    scale_x_continuous(n.breaks = num_pos, expand = expansion(mult = 0, add = 0) ) +
    scale_y_continuous(n.breaks = 8, expand = expansion(mult = 0, add = c(0, 0.05) ) ) +
    labs(x = "Sequence position", y = expression(paste(2^{'Weights'})) ) +
    theme_minimal(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(panel.grid = element_blank(),
          plot.background = element_blank(),
          legend.position = 'bottom', 
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key.size = unit(5, 'points'),
          legend.key = element_blank(),
          legend.margin = margin(t = -1, unit = "mm"),
          legend.title = element_text(vjust = 1),
          axis.line = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.y = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.length.y = unit(1, "mm"),
          axis.text = element_text(colour = "black"))
}

#' Calculate the information content expressed in bits of sequences stored in a data.frame
#'
#' @param data A data.frame with a minimum of 2 columns. One named `Sequence`, the other named as you prefer that will be specified with `ID_col`.
#' @param ID_col The name of the column in `data` to be used as the identifier of the `Sequence` column.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default.
#' @param small_n_correction Apply a small correction to the Shannon Entropy. See details. Default `FALSE`.
#' @param long_format Logical. If `TRUE` reshape the bits into a tidy long data.frame format. Default `FALSE`.
#' @param ignore_case Logical. If `TRUE` the length of the alphabet is calculated ignoring the case of the alphabet.
#' Meaning that the maximum bits height will calculated on the case-insensitive length of the alphabet. See notes for more explanation. Default `FALSE`.
#'
#' @return A data.frame or a tidy long format data.frame
#' @export
#' @description This function calculates the information content expressed in bits using the Shannon entropy. 
#' Check the details for full explanation and formulas. 
#' However, currently there's no support for LaTeX syntax for subscript text and fractions. 
#' To display them properly once could copy-paste the details section in \href{https://www.overleaf.com/}{Overleaf}.
#' 
#' @details Given an alphabet of letters of length \emph{W} where every letter 
#' defined as \emph{l} for which \emph{l} belongs to \emph{W}, 
#' we can represent the DNA alphabet as \eqn{l'} belongs to \emph{A,C,G,T} where \emph{W = 4}. 
#' With a multiple sequence alignment of \eqn{N} sequences of length \emph{I} we denote 
#' the information content expressed in bits of the letter \emph{l} at position \emph{i} 
#' with \eqn{bits_l_,_i} we define the following 
#' formula \deqn{bits_l_,_i = R(l,i) \times ( log_2(W) - (H_i + \epsilon) )} 
#' where \eqn{H_i} is the Shannon entropy representing the uncertainty of 
#' position \emph{i} is defined as: 
#' \deqn{ -\sum_{i = 1}^{W} { p_l_i \times log_2 p_l_,_i }} where \eqn{p_l_i} is the
#' relative frequency (a.k.a. probability) of letter \eqn{l} at position \eqn{i}; 
#' \eqn{\epsilon} is the approximation for small-sample corrections, 
#' i.e. a correction for an alignment of \emph{N} sequences in the alignment defined as 
#' \deqn{\epsilon = \frac{1}{log_e{2}} \times \frac{W-1}{2N}} and \eqn{R(l,i)} sequences 
#' position probability matrix containing the \eqn{p_l_i} for \emph{N} sequences.
#' 
#' @note When having an upper and lower case DNA sequence, with an `alphabet` that 
#' as both 'ATGC' and 'atgc' one case force the maximum information 
#' content to `log2(4)` instead of `log2(8)` by doing `ignore_case = TRUE`.
#' 
#' @examples
#' df2bits(data, ID_col = 'Species', 
#'         alphabet = c('a', 'c', 'g', 't'), 
#'         small_n_correction = F, 
#'         long_format = T)
df2bits <- function(data, ID_col, alphabet, small_n_correction = FALSE, 
                    long_format = FALSE, ignore_case = FALSE) {
  
  # -- 1 -- Check input params
  if( missing(alphabet) ) { alphabet <- guess_alphabet(data, ID_col) }
  
  # -- 2 -- Calculate the Shannon entropy for each column.
  # If in a column there's only one letter of the alphabet the entropy will be zero.
  H_i <- df2ppm(data, ID_col, alphabet) |>
    reframe(across(.cols = everything(), .fns = ~ -sum( (.x * log2(.x)), na.rm = T) ) )
  
  # Get also the ppm that will be multiplied by R_i
  ppm <- df2ppm(data, ID_col, alphabet)
  
  # -- 3 -- Ignore case, if needed. 
  # Trick to have DNA upper and lower case letters in the same plot with values 
  # as if the alphabet was case insentitive (i.e. 4 letters instead of 8)
  if (ignore_case == TRUE) {
    alphabet <- unique(toupper(alphabet))
  }
  
  # -- 4 -- Perform small number of sequences correction, if needed
  if (small_n_correction) {
    num_seq <- nrow(df2mat(data, ID_col))
    e_n = ( 1 / log(2) ) * ( length(alphabet) - 1 ) / ( 2 * num_seq ) 
  } else {
    e_n = 0
  }
  
  # -- 5 -- Calculate information content of each position
  R_i = log2( length(alphabet) ) - (H_i  + e_n)
  R_i <- as.numeric(R_i)
  
  # -- 6 -- Calculate the product of vector (R_i) with a matrix (ppm) by row
  # First a quick check
  stopifnot(ncol(ppm) == length(R_i))
  
  # sweep base function is used to multiply a vector with a matrix
  # x is the matrix, MARGIN = 2 means "by row" 
  # STAT is the vector, and FUN is the function (i.e. multiplication)
  bits_heights <- sweep(x = ppm, MARGIN = 2, STATS = R_i, FUN = "*")

  # -- last -- Reshape data into long format if needed.
  if( long_format == TRUE) {
    bits_heights |>
      rownames_to_column('Alphabet') |>
      pivot_longer(cols = starts_with("N"), names_to = "Position", 
                   values_to = "bits") |>
      mutate(Position = as.integer( gsub("N", "", Position ) ) ) |>
      relocate(Position, .before = Alphabet) |> arrange(Position) -> tidy_bits
    return(tidy_bits)
  } else {
    return(bits_heights)
  }
}

#' Quick visualisation of a logo. Not a fancy plot, just useful for quick checks.
#'
#' @param long_df A data.frame obtained with `df2bits(long_format = TRUE)`.
#' @param cell_txt_size Text size of the heatmap cells that corresponds to the letter weight in that position. Default 3.
#' @param axis_txt_size Text size of axis labels. Default 10.
#'
#' @return A ggplot plot
#' @import ggplot2
#' @export
#' @description If you're looking a good way to make logos check out the package \href{https://omarwagih.github.io/ggseqlogo/}{ggseqlogo}
#'
#' @examples
#' df2bits(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), small_n_correction = FALSE, long_format = T) |>
#' plot_rudimentary_logo()
#' 
#' df2bits(data, ID_col = 'Species', alphabet = c('a', 'c', 'g', 't'), small_n_correction = TRUE, long_format = T ) |>
#' plot_pwm_as_bars()
plot_rudimentary_logo <- function(long_df, axis_txt_size = 10,
                                  alphabet_palette = c('t' = "forestgreen", 
                                                       'c' = "dodgerblue",
                                                       'g' = "orange",
                                                       'a' = "firebrick3") ) {
  
  num_pos <- max(long_df$Position)
  len_alphabet <- length(unique(long_df$Alphabet))
  
  ggplot(long_df) +
    aes(x = Position, y = bits, fill = Alphabet) +
    geom_col(colour = 'black', linewidth = 0.15, width = 0.85, position = position_stack() ) +
    scale_fill_manual(values = alphabet_palette)+
    scale_x_continuous(n.breaks = num_pos,
                       expand = expansion(mult = 0, add = 0) )  +
    scale_y_continuous(n.breaks = 5, limits = c(0, log2(len_alphabet)),
                       expand = expansion(mult = 0, add = c(0, 0.05) ) ) +
    coord_cartesian(clip = 'off') +
    labs(x = "Sequence Position", y = "Information Content (bits)") +
    theme_minimal(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = 'grey84', linewidth = 0.15),
          plot.background = element_blank(),
          legend.position = 'bottom', 
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key.size = unit(5, 'points'),
          legend.key = element_blank(),
          legend.margin = margin(t = -1, unit = "mm"),
          legend.title = element_text(vjust = 1),
          axis.line = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.y = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.length.y = unit(1, "mm"),
          axis.text = element_text(colour = "black"))
}

#' Given a sequence of upper and lower cases creates the best breaks and labels for a ggplot sequence logo
#'
#' @param seq A character vector with UPPER and lower cases
#' @param lowercase_spacer An integer indicating the distance between the label numbers in the lowercase sequence.
#' @param uppercase_spacer An integer indicating the distance between the label numbers in the uppercase sequence.
#' @param uppercase_adjustment An integer defining the actual length of the region in the coord2seq.sh script. 
#' If the upper-case region has been trimmed upstream (e.g. -u -39) this parameter 
#' corrects the X-axis numbering to match the real length before trimming. Can be omitted.
#'
#' @return A list the best breaks and labels
#' @export
#' 
#' @note If proving a non-integer number for `*case_spacer` parameters it will be rounded up to the closer integer.
#' 
#' @description This is a helper function used in `plot_bits_logo()` and `plot_JSD_logo()`. 
#'
renumber_logo_seq_breaks <- function(seq, lowercase_spacer = 5, 
                                     uppercase_spacer = 10, 
                                     uppercase_adjustment) {
  # Extract UPPER- & lower-case transition boundaries in the sequences
  # I often use it as the 3'ss and 5'ss positions respectively
  
  # -- 1 -- Check input params
  if (lowercase_spacer < 2) { stop('lowercase_spacer must be >= 2') }
  if (uppercase_spacer < 2) { stop('uppercase_spacer must be >= 2') }
  # If input is not an integer number round it up, to avoid errors.
  lowercase_spacer <- round(lowercase_spacer, 0)
  uppercase_spacer <- round(uppercase_spacer, 0)
  
  # Transition positions in sequence in seq
  lwrcs_l_up <- nchar(str_extract(string = seq, pattern = "^[a-z]*(?=[A-Z])"))
  lwrcs_l_do <- nchar(str_extract(string = seq, pattern = "(?<=[A-Z])[a-z]*$"))
  last_UPPERCASE_l <- nchar(str_extract(string = seq, pattern = "[A-Z]+"))
  seq_len <- nchar(seq)
  
  ## Set X-axis ggplot based on the UPPER- lower-case transitions
  # breaks for upstream lower-case numbers
  if ( lwrcs_l_up > 0 ) { 
    first_lowercase_breaks <- seq(from = 1, to = lwrcs_l_up, by = lowercase_spacer)
  } else {
    first_lowercase_breaks <- NULL
  }
  # breaks for downstream lower-case breaks
  if ( lwrcs_l_do > 0 ) {
    last_lowercase_breaks <- seq(from = (seq_len - lwrcs_l_do + lowercase_spacer), to = seq_len, by = lowercase_spacer) 
  } else {
    last_lowercase_breaks <- NULL
  }
  
  c( first_lowercase_breaks,
     lwrcs_l_up + 1, # breaks first exon number,
     seq(from = lwrcs_l_up + uppercase_spacer, to = (seq_len - lwrcs_l_do), by = uppercase_spacer), # upper-case breaks
     as.integer(seq_len - lwrcs_l_do), # last Upper-case letter position
     last_lowercase_breaks
  ) -> best_breaks
  best_breaks <- unique(best_breaks)
  
  # ggplot X-axis labels
  if ( lwrcs_l_up > 0 ) { 
    first_lowercase_labels <- seq(from = -lwrcs_l_up, to = -1, by = lowercase_spacer)
    middle_uppercase_labels <- seq(from = 1, to = last_UPPERCASE_l / uppercase_spacer, by = 1 )*uppercase_spacer
    middle_uppercase_labels <- c(1, middle_uppercase_labels, last_UPPERCASE_l)
  } else {
    # If the sequence does not start with lowercase sequences but just uppercase letters
    # this should be length of the original coord2seq region before slopping - last_UPPERCASE_l
    # For example and exon of length 50nt: 50 - last_UPPERCASE_l
    
    if( missing(uppercase_adjustment) ) {
      first_lowercase_labels <- 1
      middle_uppercase_labels <- seq(from = 1, to = last_UPPERCASE_l / uppercase_spacer, by = 1 )*uppercase_spacer
      middle_uppercase_labels <- c(middle_uppercase_labels, last_UPPERCASE_l)
      
    } else if ( is.numeric(uppercase_adjustment) ) {
      first_lowercase_labels <- uppercase_adjustment - last_UPPERCASE_l
      last_uppercase_label <- uppercase_adjustment
      middle_uppercase_labels <- seq(from = first_lowercase_labels, to = last_uppercase_label, by = uppercase_spacer )
      
    } else {
      first_lowercase_labels <- 1
      warning("There's something wrong with `uppercase_adjustment`.", 
              "It must be a number. Setting to 1")
    }
  }
  
  if ( lwrcs_l_do > 0 ) {
    last_lowercase_labels <- paste0("+", seq(from = lowercase_spacer, to = lwrcs_l_do, by = lowercase_spacer) )
  } else {
    last_lowercase_labels <- NULL
  }
  
  c( first_lowercase_labels, 
     middle_uppercase_labels, 
     last_lowercase_labels
  ) -> best_labels
  
  best_labels <- unique(best_labels)
  
  # sanity check
  if( length(best_breaks) != length(best_labels) ){
    stop("Something went wrong when calculating the ideal breaks & labels. ",
         "Try setting `uppercase_spacer` or `lowercase_spacer` to another value.")
  }
  
  # prepare a list for output
  output <- list(best_breaks, best_labels)
  names(output) <- c('Breaks', 'Labels')
  
  return(output)
}

#' Given an alphabet name set a palette for the plots.
#'
#' @param alphabet A character vector with the alphabet letters.
#' @param alphabet_name The alphabet name returned by `guess_alphabet()`.
#'
#' @return An palette for `ggseqlogo`
#' @importFrom ggseqlogo make_col_scheme
#' @export
get_alphabet_palette <- function(alphabet, alphabet_name) {
  
  supported_alphabets <- c('DNA_UPlowER', 'DNA_UPPER', 'RNA_UPlowER', 
                           'DNA_UPPER-gapped', 'RNA_UPPER-gapped', 
                           'user_specified_alphabet')
  
  if ( !any(alphabet_name %in% supported_alphabets) ){
    stop('The alphabet name ', alphabet_name, ' is not a supported alphabet!')
  }

  if (alphabet_name == 'DNA_UPlowER') {
    
    palette <- make_col_scheme(
      chars = c(unique(toupper(sort(alphabet))), unique(tolower(sort(alphabet)))), 
      groups = rep( unique( toupper( sort(alphabet) ) ), 2), 
      cols = rep( c('A' = 'firebrick3', 'C' = 'dodgerblue',
                    'G' = 'darkgoldenrod1', 'T' = 'darkgreen'), 2),
      name = 'UPPER/lower-case DNA'
    )
    # seq_type <- 'dna' # Define seq type for ggseqlogo, probably not required
    
  } else if (alphabet_name == 'DNA_UPPER') {
    palette <- make_col_scheme(
      chars = sort(alphabet),
      groups = sort(alphabet),
      cols = c('A' = 'firebrick3', 'C' = 'dodgerblue',
               'G' = 'darkgoldenrod1', 'U' = 'darkgreen'),
      name = 'UPPER-case DNA'
    )
    # seq_type <- 'rna' # Define seq type for ggseqlogo, probably not required
  } else if (alphabet_name == 'RNA_UPlowER') {
    palette <- make_col_scheme(
      chars = c(unique(toupper(sort(alphabet))), unique(tolower(sort(alphabet)))), 
      groups = rep( unique( toupper( sort(alphabet) ) ), 2), 
      cols = rep( c('A' = 'firebrick3', 'C' = 'dodgerblue',
                    'G' = 'darkgoldenrod1', 'U' = 'darkgreen'), 2),
      name = 'UPPER/lower-case RNA'
    )
    # seq_type <- 'rna' # Define seq type for ggseqlogo, probably not required
    
  } else if (alphabet_name ==  'DNA_UPPER-gapped') {
    
    palette <- make_col_scheme(
      chars = unique( sort(alphabet) ),
      groups = unique( sort(alphabet) ), 
      cols = c('-' = 'gray16', 'A' = 'firebrick3', 'C' = 'dodgerblue',
               'G' = 'darkgoldenrod1', 'T' = 'darkgreen'),
      name = 'UPPER-gapped DNA'
    )
  } else if (alphabet_name ==  'RNA_UPPER-gapped') {
    palette <- make_col_scheme(
      chars = unique( sort(alphabet) ),
      groups = unique( sort(alphabet) ), 
      cols = c('-' = 'gray16', 'A' = 'firebrick3', 'C' = 'dodgerblue',
               'G' = 'darkgoldenrod1', 'U' = 'darkgreen'),
      name = 'UPPER-gapped RNA'
    )
  } else if ( alphabet_name == 'user_specified_alphabet' ) {
    # set a default palette
    palette <- make_col_scheme(
      chars = alphabet,
      cols = rev( palette.colors(n = length(alphabet), palette = 'Alphabet') ),
      name = 'default_palette'
    )
    # seq_type <- 'other' # Define seq type for ggseqlogo, probably not required
    
  } else {
    stop('Problem setting up the colour palette for the letters. ',
         'Alphabet type: ', alphabet_name, ' not yet supported.',
         'The Alphabets currently supported are: ', 
         paste0(supported_alphabets, collapse = ", ") )
  }
  return(palette)
}

#' Plot sequence information content as bits.
#'
#' @param df A data.frame with a column called `Sequence` and anther defined in `ID_col`. Use the wrapper `fasta2df()` to import a fasta file to a dataframe that is compatible with this function.
#' @param ID_col A character specifying the column name in `df` to be used as the sequence ID.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default from `df`.
#' @param y_lims A numeric vector of length 2 specifying the Y-axis min and max value. Default `c(0, NA)`
#' @param InfoContent_thrshld The information content (bits) threshold to consider a letter position 
#' to be highlighted. Position whose letter information sum is lower than `InfoContent_thrshld` are highlithed by a vertical bar in the plot.
#' @param anno_width A small number that defines how wide the vertical bar should be. Default `0.75`.
#' @param highlight_colour A colour name to fill the letter highlighting rectangle. Default `grey74` with `alpha = 0.5`.
#' @param axis_txt_size A number specifying the size of the axis text in the plot. Default 10.
#' @param ttl_txt Some text in quotes specifying the plot title.
#' @param small_n_correction Logical for applying a small correction to the Shannon Entropy for low number of input sequences. Parameter for \link[niar]{df2bits}, type `?df2bits` for details. Default `FALSE`.
#' @param ... Advanced parameters passed to \link[niar]{renumber_logo_seq_breaks}.
#'
#' @return A ggplot sequence logo
#' @import ggplot2
#' @importFrom ggseqlogo make_col_scheme geom_logo
#' @importFrom dplyr between
#' @export
#' 
#' @description This function makes a publication quality logo. 
#' If the sequences contains an alphabet with letters in upper AND lower case, 
#' the Shannon entropy is calculated in a case-insensitive way. 
#' Meaning that the maximum bit value will be `log2` the number of in the alphabet. 
#' Instead in the plot the letters preserve their case.
#' 
#' @details This function uses \link[ggseqlogo]{geom_logo} to plot the logo.
#'
#' @examples
#' plot_bits_logo(df = df_w_seqs, ID_col = 'Species')
#' 
#' plot_bits_logo(df = df_w_seqs, ID_col = 'Species', 
#'                InfoContent_thrshld = 0.5, anno_width = 0.5, 
#'                highlight_colour = 'lightblue', 
#'                uppercase_spacer = 5, lowercase_spacer = 6)
#'                
#' # compare the plot when a small number epsilon is added to the Shannon's Entropy formula 
#' # This is the default               
#' plot_bits_logo(df = suz12_ex4_eutheria, ID_col = 'Species', y_lims = c(0, 2),
#'                small_n_correction = F, ttl_txt = 'Without small correction')    
#'
#' # This is optional, but recommended when having few input sequences
#' plot_bits_logo(df = suz12_ex4_eutheria, ID_col = 'Species', y_lims = c(0, 2),
#'                small_n_correction = F, ttl_txt = 'Without small correction')                               
plot_bits_logo <- function(df, ID_col, alphabet, small_n_correction = FALSE,
                           y_lims = c(0, NA), anno_width = 0.75, 
                           InfoContent_thrshld, highlight_colour = 'grey74',
                           axis_txt_size = 10, ttl_txt = NULL, ...) {
  
  # -- 1 -- Calculate the information content bits using the Shannon Entropy formula
  InfoContent <- df2bits(data = df, ID_col = ID_col, 
                         small_n_correction = small_n_correction, 
                         long_format = FALSE, ignore_case = TRUE )
  
  # Coerce to matrix
  if(is.data.frame(InfoContent)){ InfoContent <- as.matrix(InfoContent) }
  
  # -- 2 -- Set colour palette
  if ( missing(alphabet) ) {
    alphabet <- guess_alphabet(data = df, ID_col = ID_col, verbose = TRUE)
    alphabet_name <- guess_alphabet(data = df, ID_col = ID_col, name = TRUE)
    # message('Alphabet type: ', alphabet_name, '\nLetters: ', paste0(alphabet, collapse = ', ') )
  } else {
    alphabet_name <- 'user_specified_alphabet'
  }
  
  palette <- get_alphabet_palette(alphabet = alphabet,
                                  alphabet_name = alphabet_name)
  
  # -- 3 -- Check and set Y-axis ranges
  # Get the theoretical and observed maximum information content
  num_letter_in_alphabet <- length( unique( toupper( alphabet ) ) )
  theoretical_max <- log2(num_letter_in_alphabet)
  observed_max <- max(InfoContent)
  
  # By default the Y-axis scales to the maximum bits value 
  if ( is.na(y_lims[2]) ) {
    y_lims[2] <- observed_max
  }
  
  # check Y-axis limits
  stopifnot(between(x = y_lims[1], left = 0, right = theoretical_max))
  stopifnot(between(x = y_lims[2], left = 0, right = theoretical_max))
  
  if( y_lims[1] > y_lims[2] ) {
    stop("The Shannon's information content bits (Y-axis) must be between 0 ", 
         "and ", theoretical_max, " for an input sequence with ", 
         num_letter_in_alphabet, " letters.\nThe first 'y_lims' is the lower ", 
         " bound limit and must be less than the second 'y_lims' value.")
  }
  
  # Check that the maximal JSD is set between 0 and 1
  if( !missing(InfoContent_thrshld) ) {
    if( !between(x = InfoContent_thrshld, left = 0, right = theoretical_max) ) {
      stop("The minimum Shannon's information content bits (Y-axis) must be ",
           "between 0 and ", theoretical_max, ", not ", 
           InfoContent_thrshld, "!")
    }
  }
  
  ## -- 4 -- Extract UPPER- & lower-case transition boundaries in the sequences
  # (I often use it for exons 3'ss and 5'ss positions respectively)
  
  # pick a random representative sequence from each dataset
  seq <- sample(df$Sequence, size = 1)
  Xaxis_numbers <- renumber_logo_seq_breaks(seq, ...)
  
  # if (alphabet_name %in% c('DNA_UPlowER', 'RNA_UPlowER') ) {
  #   # set X-axis breaks & labels based on the UPPER- lower-case transitions
  #   message('Trying to find the best breaks based on letters case' )
  #   Xaxis_numbers <- renumber_logo_seq_breaks(seq, ...)
  # } else {
  #   message('Trying nothing with the x-axis breaks' )
  #   # if the sequence doesn't have UPPER or lower cases just start from 1
  #   Xaxis_numbers <- list()
  #   x_breaks_spacing <- round(nchar(seq) * 0.05, 0) + 1
  #   Xaxis_numbers$Breaks <- seq(from = 1, to = nchar(seq), by = x_breaks_spacing)
  #   Xaxis_numbers$Labels <- seq(from = 1, to = nchar(seq), by = x_breaks_spacing)
  # }

  breaks <- Xaxis_numbers$Breaks
  labels <- Xaxis_numbers$Labels
  
  # -- 5 -- Get positions with low information content
  if ( missing(InfoContent_thrshld) ) {
    # If letters highlighting is not required make an empty one df
    annotate_df <- data.frame(Position = max(breaks),
                              ymin = y_lims[1], ymax = y_lims[2])
    highlight_colour <- NA
  } else {
    # If highlighting is required 
    lowInfo_pos <- which(signif(apply(InfoContent, 2, function(x) { sum(x) }), 2) <= InfoContent_thrshld)
    
    if( length(lowInfo_pos) < 1 ) {
      message("There are NO positions in the logo with information ", 
              "content (bits) lower or equal to ", InfoContent_thrshld, 
              ". Ignoring the letter highlighting request.")
      
      annotate_df <- data.frame(Position = max(breaks),
                                ymin = y_lims[1], ymax = y_lims[2])
      highlight_colour <- NA
      
    } else {
      # -- Annotate letters with low information content
      annotate_df <- data.frame(Position = lowInfo_pos, 
                                ymin = y_lims[1], ymax = y_lims[2])
    }
  } 
  
  ## -- 6 -- Plot logo of Shannon entropy bits
  ggplot() + 
    geom_rect(data = annotate_df,
              aes(xmin = Position - anno_width, xmax = Position + anno_width, ymin = ymin, ymax = ymax),
              colour =  NA, linewidth = 0, alpha = 0.5, fill = highlight_colour) +
    geom_logo(data = InfoContent, method = 'custom', col_scheme = palette,
              stack_width = 0.98, font = "roboto_regular") +
    scale_x_continuous(expand = expansion(mult = 0, add = 0), 
                       breaks = breaks, labels = labels ) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0, 0.01) ),
                       n.breaks = 3, limits = y_lims) +
    labs(y = 'bits', title = ttl_txt) +
    theme_classic(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_line(linewidth = 0.15, colour = 'grey84'), 
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = 'none', 
          # axis.line = element_line(linewidth = 0.15, colour = 'black'),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 0, unit = "mm"), angle = 90),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(vjust = 1, margin = margin(t = -0.5, unit = "mm")),
          axis.text.y = element_text(hjust = 1, margin = margin(r = 0, unit = "mm")),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.ticks.y = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.length.y = unit(1, 'mm'),
          plot.title = element_text(size = axis_txt_size - 1.1, hjust = 0.05, vjust = 0.75, margin = margin(b = 0)) 
    ) -> p_bits_logo
  return(p_bits_logo)
}

#' Helper function to calculate the Jensen-Shannon divergence 
#'
#' @param p A data.frame of positional probabilities
#' @param q Another data.frame of positional probabilities
#'
#' @return A numeric list 
#' @description The 2 data.frames `p` and `q` must have the same dimensions. The calculated data.frame `m` is the mean between `p` and `q`
#' @export
#' 
#' @examples
#' H_i <- mapply(FUN = col_JS_divergence, p = ppm_p, q = ppm_q) 
col_JS_divergence <- function(p, q) {
  m <- (p + q) / 2
  H_i = 0.5*sum(p*log2(p/m), na.rm = TRUE) + 0.5*sum(q*log2(q/m), na.rm = TRUE)
  return(H_i)
}

#' Helper function to calculate the sum of the absolute difference between 2 position probability matrices
#'
#' @param p 
#' @param q 
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' sum_abs_delta_prob_l <- mapply(FUN = col_sub_abs, p = ppm_p, q = ppm_q) 
col_sub_abs <- function(p, q) {
  sum_abs_delta <- sum(abs(p - q), na.rm = T)
  return(sum_abs_delta)
}

#' Calculate Jensen-Shannon divergence (JSD) between 2 set of sequences. 
#'
#' @param df1 
#' @param df2 
#' @param ID_col 
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default.
#'
#' @return A data.frame or a tidy long formatted data.frame
#'
#' @export
#' @description Letters in the alphabet that are less abundant in `df2` sequences relative to `df1` sequences have positive values.
#' More abundant letters in `df2` sequences relative to `df1` have nevative values. 
#'
#' @examples
#' df2JSD(df1 = data_p, df2 = data_q, ID_col, long_format = F)
df2JSD <- function(df1, df2, ID_col, alphabet, long_format = FALSE) {
  
  # -- 0 -- Get alphabet if missing
  if( missing(alphabet) ) { 
    alphabet1 <- guess_alphabet(df1, ID_col) 
    alphabet2 <- guess_alphabet(df2, ID_col) 
    
    if (identical(alphabet1, alphabet2) ) { 
      alphabet <- alphabet1 
    } else{
        stop("I wasn't able to understand if the 2 input dataframes have ", 
             "the same alphabet. Please specify it with alphabet = c(...) ")
    }
  }
  
  # -- 1 -- Extract ppms from sequence dataframes
  ppm_p <- df2ppm(data = df1, ID_col, alphabet)
  ppm_q <- df2ppm(data = df2, ID_col, alphabet)
  
  # -- 2 -- Apply the JS function col-wise on 2 dataframes 
  H_i <- mapply(col_JS_divergence, p = ppm_p, q = ppm_q) 
  
  # -- 3 -- Calculate normalised probability difference at position i for letter l
  if ( !all(dim(ppm_p) == dim(ppm_q) ) ) {
    stop('The 2 input sequences do NOT have the same length!',
         '\n', 'df1 sequences are ', unique(nchar(df1$Sequence)), 
         ' letters long\n', 'df2 sequences are ', unique(nchar(df2$Sequence)),
         ' letters long\n')
  }
  delta_prob_l_i <- ppm_p - ppm_q
  sum_abs_delta_prob_i <- mapply(FUN = col_sub_abs, p = ppm_p, q = ppm_q) 
  
  # -- 4 -- Normalise the probability difference using the sum of the absolute difference at position i
  probability_difference_l_i <- sweep(x = delta_prob_l_i, MARGIN = 2, STATS = sum_abs_delta_prob_i, FUN = "/")
  
  # since sum_abs_delta_prob_i is zero if p == q, I set all NaN to zeros
  probability_difference_l_i[is.na(probability_difference_l_i)] <- 0
  
  # -- 6 -- Calculate Jensen-Shannon divergence
  JSdivergence <- sweep(x = probability_difference_l_i, MARGIN = 2, STATS =  H_i, FUN = "*") 
  
  # -- 7 -- Reshape data into long format if needed.
  if( long_format == TRUE) {
    JSdivergence |>
      rownames_to_column('Alphabet') |>
      pivot_longer(cols = starts_with("N"), names_to = "Position", 
                   values_to = "JS_divergence") |>
      mutate(Position = as.integer( gsub("N", "", Position ) ) ) |>
      relocate(Position, .before = Alphabet) |> arrange(Position) -> tidy_JSdiv
    return(tidy_JSdiv)
  } else {
    return(JSdivergence)
  }
}

#' Quick visualisation of the Jensen-Shannon divergence 
#'
#' @param long_df A data.frame created with `df2JSD(long_format = T)`.
#' @param axis_txt_size 
#' @param alphabet_palette 
#'
#' @return A ggplot plot
#' @import ggplot2
#' @export
#'
#' @examples
#' df2JSD(df1 = data_p, df2 = data_q, ID_col, alphabet, long_format = T) |>
#'   plot_rudimentary_JSdiv_logo()
#' # alphabet is case sensitive, in this example the lowercase letters will be in grey.
#' df2JSD(df1 = df_clade1, df2 = df_clade2, ID_col = 'Species', 
#'                 alphabet = c('A', 'C', 'T', 'G', 'a', 'c', 't', 'g'), long_format = T) |>
#'    plot_rudimentary_JSD_logo(alphabet_palette = c('A' ="firebrick3", 'C' = "dodgerblue", 'T' = "forestgreen", 'G' ="orange") )
plot_rudimentary_JSD_logo <- function(long_df, axis_txt_size = 10,
                                      alphabet_palette = c('t' = "forestgreen", 
                                                           'c' = "dodgerblue",
                                                           'g' = "orange", 
                                                           'a' = "firebrick3")
                                      ) {
  num_pos <- max(long_df$Position)
  len_alphabet <- length(unique(long_df$Alphabet))
  
  ggplot(long_df) +
    aes(x = Position, y = JS_divergence, fill = Alphabet) +
    geom_hline(yintercept = 0, linewidth = 0.75, colour = 'black', linetype = 'solid') +
    geom_col(colour = 'black', linewidth = 0.15, width = 0.85, position = position_stack() ) +
    scale_fill_manual(values = alphabet_palette)+
    scale_x_continuous(n.breaks = num_pos,
                       expand = expansion(mult = 0, add = 0) )  +
    scale_y_continuous(n.breaks = 5,
                       expand = expansion(mult = 0.05 ) ) +
    coord_cartesian(clip = 'off') +
    labs(x = "Sequence Position", y = "Jensen-Shannon divergence") +
    theme_minimal(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(panel.grid = element_blank(),
          panel.grid.major = element_line(colour = 'grey84', linewidth = 0.15),
          plot.background = element_blank(),
          legend.position = 'bottom', 
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key.size = unit(5, 'points'),
          legend.key = element_blank(),
          legend.margin = margin(t = -1, unit = "mm"),
          legend.title = element_text(vjust = 1),
          axis.line = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.y = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.length.y = unit(1, "mm"),
          axis.text = element_text(colour = "black")) -> p_JS_div_logo
  return(p_JS_div_logo)
} 

#' Calculate and plot the Jensen-Shannon divergence between 2 sets of sequences contained in 2 data.frames
#'
#' @param df1 A data.frame with a column called `Sequence` and anther defined in `ID_col`. Use the wrapper `fasta2df()` to import a fasta file to a dataframe that is compatible with this function.
#' @param df2 A different data.frame with the same properties as `df1`.
#' @param ID_col A character specifying the column name in `df1` and `df2` to be used as the sequence ID.
#' @param alphabet A character vector containing the alphabet letters present in `Sequence`. Guessed by default from `df1`.
#' @param y_lims A numeric vector of length 2 specifying the Y-axis min and max value. Default `c(-0.5, 0.5)`
#' @param max_JS_div_thrshld The maximal absolute JS divergence threshold to consider a letter position to be highlighted by a vertical bar in the plot.
#' @param anno_width A small number that defines how wide the vertical bar should be. Default `0.75`.
#' @param highlight_colour A colour name to fill the letter highlighting rectangle. Default `grey74` with `alpha = 0.5`.
#' @param axis_txt_size A number specifying the size of the axis text in the plot. Default 10.
#' @param ttl_txt Some text in quotes specifying the plot title.
#' @param ... Advanced parameters passed to \link[niar]{renumber_logo_seq_breaks}.
#'
#' @return A ggplot sequence logo with the Jensen-Shannon divergence on the Y-axis.
#' @import ggplot2
#' @importFrom ggseqlogo make_col_scheme geom_logo
#' @importFrom stringr str_extract
#' @importFrom dplyr between
#' @export
#'
#' @examples
#' plot_JSD_logo(df1 = clade1, df2 = clade2f, ID_col = 'Species',
#'               max_JS_div_thrshld = 0.2, anno_width = 0.5)
plot_JSD_logo <- function(df1, df2, ID_col, alphabet, 
                          y_lims = c(-0.5, 0.5), max_JS_div_thrshld, 
                          anno_width = 0.75, highlight_colour = 'grey74',
                          axis_txt_size = 10, ttl_txt = NULL, 
                          ...) {
  
  # -- 1 -- Check input parameters
  # Check that the maximal JSD is set between 0 and 1
  if( !missing(max_JS_div_thrshld) ) {
    if( !between(x = max_JS_div_thrshld, left = 0, right = 1) ) {
      stop("The maximal absolute Jensen-Shannon divergence must be between 0 and 1, ",
           "not ", max_JS_div_thrshld, "!")
    }
  }
  
  # check Y-axis limits
  stopifnot(between(x = y_lims[1], left = -0.5, right = 0.5))
  stopifnot(between(x = y_lims[2], left = -0.5, right = 0.5))
  
  if( y_lims[1] > y_lims[2] ) {
    stop("The Jensen-Shannon divergence bits (Y-axis) must be between -0.5 ", 
         "and 0.5, the first 'y_lims' is the lower bound limit and must be ", 
         "less than the second 'y_lims' value")
  }
  
  # -- 2 -- Set colour palette
  if ( missing(alphabet) ) {
    alphabet <- guess_alphabet(data = df1, ID_col = ID_col)
    alphabet_name <- guess_alphabet(data = df1, ID_col = ID_col, name = TRUE)
    # message('Alphabet type: ', alphabet_name, '\nLetters: ', paste0(alphabet, collapse = ', ') )
  } else {
    alphabet_name <- 'user_specified_alphabet'
  }

  palette <- get_alphabet_palette(alphabet = alphabet,
                                  alphabet_name = alphabet_name)
  
  # -- 3 -- Calculate Jensen-Shannon divergence
  JSdiv <- df2JSD(df1, df2, ID_col, alphabet, long_format = FALSE)
  
  # Coerc to matrix
  if(is.data.frame(JSdiv)){ JSdiv <- as.matrix(JSdiv) }
  
  ## -- 4 -- Extract UPPER- & lower-case transition boundaries in the sequences
  # (I often use it for exons 3'ss and 5'ss positions respectively)
  
  # pick a random representative sequence from each dataset
  seq_p <- sample(df1$Sequence, size = 1)
  seq_q <- sample(df2$Sequence, size = 1)
  
  # set X-axis breaks & labels based on the UPPER- lower-case transitions
  Xaxis_numbers_p <- renumber_logo_seq_breaks(seq_p, ...)
  Xaxis_numbers_q <- renumber_logo_seq_breaks(seq_q, ...)
  
  # if (alphabet_name %in% c('DNA_UPlowER', 'RNA_UPlowER') ) {
  #   Xaxis_numbers_p <- renumber_logo_seq_breaks(seq_p, ...)
  #   Xaxis_numbers_q <- renumber_logo_seq_breaks(seq_q, ...)
  # } else {
  #   # if the sequences don't have UPPER or lower cases just start from 1
  #   Xaxis_numbers_p <- list()
  #   Xaxis_numbers_p$Breaks <- 1:length(seq_p)
  #   Xaxis_numbers_p$Labels <- 1:length(seq_p)
  #   
  #   Xaxis_numbers_q <- list()
  #   Xaxis_numbers_q$Breaks <- 1:length(seq_q)
  #   Xaxis_numbers_q$Labels <- 1:length(seq_q)
  # }
  
  breaks_p <- Xaxis_numbers_p$Breaks
  labels_p <- Xaxis_numbers_p$Labels
  
  breaks_q <- Xaxis_numbers_q$Breaks
  labels_q <- Xaxis_numbers_q$Labels
  
  # Check that the 2 representative sequences have case-transitions in the same spot
  if ( !all( c(breaks_p == breaks_q, labels_p == labels_q) ) ) {
    warning("The positions of the UPPER to lower case transitions are NOT ", 
            "the same between the 2 sets of inputs!")
  }
  
  # -- 5 -- Get maximal Jensen-Shannon divergence
  if ( missing(max_JS_div_thrshld) ) {
    # If letters highlighting is not required make an empty one df
    annotate_df <- data.frame(Position = max(breaks_p),
                              ymin = y_lims[1], ymax = y_lims[2])
    highlight_colour <- NA
  } else {
    # If highlighting is required 
    divergent_pos <- which(signif(apply(JSdiv, 2, function(x) { sum(abs(x)) }), 2) >= max_JS_div_thrshld)
    
    if( length(divergent_pos) < 1 ) {
      message("There are NO positions in the sequences with maximal JS ",
              "divergence greater or equal to ", max_JS_div_thrshld, 
              ". Ignoring the highlighting request.")
      
      annotate_df <- data.frame(Position = max(breaks_p),
                                ymin = y_lims[1], ymax = y_lims[2])
      highlight_colour <- NA
      
    } else {
      # -- Annotate highly divergent letters
      annotate_df <- data.frame(Position = divergent_pos, 
                                ymin = y_lims[1], ymax = y_lims[2])
    }
  } 
  
  ## -- 6 -- Plot logo of JSD
  ggplot() + 
    geom_rect(data = annotate_df,
              aes(xmin = Position - anno_width, xmax = Position + anno_width, ymin = ymin, ymax = ymax),
              colour =  NA, linewidth = 0, alpha = 0.5, fill = highlight_colour) +
    geom_logo(data = JSdiv, method = 'custom', col_scheme = palette,
              stack_width = 0.98, font = "roboto_regular") + # , seq_type = seq_type
    scale_x_continuous(expand = expansion(mult = 0, add = 0), 
                       breaks = breaks_p, labels = labels_p) +
    scale_y_continuous(expand = expansion(mult = 0, add = c(0, 0.01)),
                       n.breaks = 3, limits = y_lims) +
    labs(y = 'JSD (bits)', title = ttl_txt) +
    theme_classic(base_size = axis_txt_size, base_family = 'Arial') %+replace% 
    theme(#panel.grid.major.x = element_line(linewidth = 0.15, colour = 'grey84'), 
          panel.grid.major.x = element_blank(), 
          panel.grid.major.y = element_line(linewidth = 0.15, colour = 'grey84'), 
          panel.background = element_blank(),
          plot.background = element_blank(),
          legend.position = 'none', 
          # axis.line = element_line(linewidth = 0.15, colour = 'black'),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 0, unit = "mm"), angle = 90),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(vjust = 1, margin = margin(t = -0.5, unit = "mm")),
          axis.text.y = element_text(hjust = 1, margin = margin(r = 0, unit = "mm")),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          # axis.ticks.y = element_line(linewidth = 0.15, colour = 'black'),
          axis.ticks.length.y = unit(1, 'mm'),
          plot.title = element_text(size = axis_txt_size - 1.1, hjust = 0.05, vjust = 0.75, margin = margin(b = 0)) 
          ) -> p_JSD_logo
  return(p_JSD_logo)
}
