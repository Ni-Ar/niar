#' Handy function to create a biomaRt R object based on the ENSEMBL database to be used when querying biomart
#'
#' @param species a character specifying the species to work with. Currently the supported species are:
#' \itemize{
#' \item{\code{hsapiens} - Human (Homo Sapiens) Default}
#' \item{\code{mmusculus} - Mouse (Mus Musculus)}
#' \item{\code{rnorvegicus} - Rat (Rattus Norvegicus)}
#' \item{\code{ggallus} - Chicken (Gallus Gallus)}
#' } 
#' @param latest Logical. Do you want to use the latest biomart version? Default `TRUE`. If `FALSE` penultimate version is used.
#' @param use_mirror Logical. Do you want to use a mirror for creating the biomart object? Default `FALSE`.
#' @param which_mirror Specify which biomart mirror to use. Chose from:
#' \itemize{
#' \item{\code{www} - UK (Sanger Institute). in Cambridge, default}
#' \item{\code{uswest} - US West (Amazon AWS) on West Coast of US of America}
#' \item{\code{useast} - US East (Amazon AWS) on East Coast of US of America}
#' \item{\code{asia} - Asia (Amazon AWS) in Singapore}
#' } 
#' @param out_dir Either `FALSE` (default) or a character specifying a path where to save the biomaRt object. If directory does not exist it will be created.
#' @param verbose Logical. Print info messages? Default `FALSE`.
#'
#' @return
#' @import biomaRt
#' @export
#'
#' @examples
#' # Human ENSEMBL biomaRT object
#' ensembl <- gimme_mart()
gimme_mart <- function( species = "hsapiens", latest = TRUE, 
                        use_mirror = FALSE, which_mirror, out_dir = FALSE, 
                        verbose = FALSE) {
  supported_speciess <- c("hsapiens", "mmusculus", "rnorvegicus", "ggallus")
  if (! species %in% supported_speciess ) {
    stop("The species you selected ", species, " is not yet supported. ",
         "Use first letter and second full name of the species scientific ",
         " name like 'hsapiens' for Homo Sapiens or 'mmusculus' for ",
         "Mus Musculus.")
  }
  
  if (is.null(out_dir)) {
    stop("You didn't specify the out_dir parameter! It must be FALSE to ",
         "keep the biomart object in the environment, or a character (path) ",
         " to specify where to save the biomart object as a R object.")
    
  }
  
  if (out_dir == TRUE) {
    stop("out_dir parameter is either FALSE or a character specifying the path.")
  }
  
  if (use_mirror) {
    if (is.null(which_mirror)) {
      stop("You want to use a mirror, but need to specify one from: ",
           "'www', 'uswest', 'useast', 'asia'.")
    }
  }
  
  # ---- Make Biomart object
  ensembl <- useMart(biomart = "ensembl")
  datasets <- listDatasets(ensembl)
  DB_versions <- listEnsemblArchives()
  
  # If retrieve Ensembl Archive fails because of Mart connection
  if ( !exists("DB_versions") ) {
    # Create a fake database version for human
    if (species == "hsapiens") {
      DB_versions <- data.frame(name = c("Ensembl 104", "Ensembl 103"), 
                                date = c("May 2021", "Feb 2021"),
                                url = c("https://may2021.archive.ensembl.org", 
                                        "https://feb2021.archive.ensembl.org"),
                                version = c(104, 103),
                                current_release = c("*", "")
      )
    }
  }
  
  # ---- Get species dataset name
  species_ensembl <- datasets[grep(species, datasets$dataset), "dataset" ]
  
  if ( length(species_ensembl) == 0 ) {
    stop("The species you selected ", species, " is wrong or ", 
         "not supported. Use first letter and second full name of the ", 
         "species scientific name like 'hsapiens' for Homo Sapiens or ", 
         " 'mmusculus' for Mus Musculus.")
  }
  
  if (!use_mirror) {
    # -- If using a mirror you cannot specify the version
    
    # -- Get db version
    if ( latest == TRUE ) {
      
      DB_v <- DB_versions[grep("\\*" , DB_versions$current_release), "version"]
      DB_v <- as.integer(DB_v)
      
      # -- Get db version release date
      DB_v_date <- DB_versions[grep("\\*" , DB_versions$current_release), "date"]
      
      if (verbose) {
        message("Creating Mart object with latest version ", DB_v, 
                " released on ", DB_v_date )
        
      }
      
      DB_v_date <- sub(" ", "_", DB_v_date)
      
    } else if (latest == FALSE ) {
      
      DB_v <- DB_versions[grep("\\*" , DB_versions$current_release) + 1, "version"]
      DB_v <- as.integer(DB_v)
      
      # -- Get db version release date
      DB_v_date <- DB_versions[grep("\\*" , DB_versions$current_release) + 1, "date"]
      if (verbose) {
        message("Creating Mart object with penultimate version ", DB_v, 
                " released on ", DB_v_date )
      }
      DB_v_date <- sub(" ", "_", DB_v_date)
      
    } else { stop("latest argument must be logical either T or F")
    }
    
  }
  
  if (use_mirror) {
    
    if (which_mirror %in% c('www', 'uswest', 'useast', 'asia') ) {
      # ---- Ensembl Dataset on a mirro
      ensembl <- useEnsembl(biomart = "ensembl", 
                            dataset = species_ensembl,
                            mirror  = which_mirror)
      
    } else {
      stop("The which_mirror argument must be one of: ",
           "'www', 'uswest', 'useast', 'asia'.",
           which_mirror, " is not supported")
    }
    
  } else if (!use_mirror) {
    # ---- Ensembl Dataset on main host
    ensembl <- useEnsembl(biomart = "ensembl", 
                          dataset = species_ensembl, 
                          version = DB_v, 
                          host = "https://www.ensembl.org")
  }
  
  # -- Write output or keep in the environment
  if ( is.character(out_dir) ) {
    
    # ---- Output Directory
    out_biomart_dir <- file.path(out_dir, species)
    # ---- Check if output dir exists, if not create it
    if (!dir.exists(out_biomart_dir)) { dir.create(out_biomart_dir) }
    
    if (!use_mirror) {
      output_ensembl_name <- paste0(out_biomart_dir, "/", 
                                    format(Sys.Date(), "%Y_%m_%d"), "_",
                                    species_ensembl, "_v", DB_v, "_",
                                    DB_v_date, ".rds")
      
    } else if (use_mirror) {
      output_ensembl_name <- paste0(out_biomart_dir, "/", 
                                    format(Sys.Date(), "%Y_%m_%d"), "_",
                                    species_ensembl, "_mirror_", which_mirror,
                                    ".rds")
    } else {
      stop("Something is wrong when creating a name for the ENSEMBL Object.")
    }
    
    
    saveRDS(object = ensembl, file = output_ensembl_name)
    if ( verbose ) { message("ENSEMBL Mart Object saved in ", out_biomart_dir) }
    
  } else if ( out_dir == FALSE ) {
    return(ensembl)
  } 
}

#' Convenient converter from ensembl ID to gene name
#'
#' @param ensembl_gene_id A character specifying the ENSEMBL gene ID
#' @param only_gene_name Logical, whether to return only the gene name or a `data.frame` with more info. Default `TRUE`
#' @param verbose Logical. Do you want me to be chatty?
#' @param mRt_objct A `biomaRt` object created with the `gimme_mart()` function.
#'
#' @return
#' @import biomaRt
#' @export
#'
#' @examples
#' ensembl <- gimme_mart()
#' ensembl_id_2_gene_name("ENSG00000000003")
#' "TSPAN6"
#' # Example with mouse ensembl
#' ensembl <- gimme_mart(species = "mmusculus", latest = T, use_mirror = F, 
#'                       out_dir = F, verbose = T)
#' ensembl_id_2_gene_name("ENSMUSG00000067377", only_gene_name = F)
#'      ensembl_gene_id external_gene_name chromosome_name start_position   gene_biotype
#' 1 ENSMUSG00000067377             Tspan6               X      132791817 protein_coding
#'
ensembl_id_2_gene_name <- function(ensembl_gene_id, only_gene_name = TRUE,
                                    verbose = FALSE, mRt_objct = ensembl) {
  
  # -- 0 -- Check input parameters ----
  stopifnot(is.character(ensembl_gene_id))
  ensembl_gene_id_regex <- "ENS[A-Z]+[0-9]{11}"
  valid_ensembl_id <- grepl(pattern = ensembl_gene_id_regex, 
                            x = ensembl_gene_id, ignore.case = F)
  
  if ( !all(valid_ensembl_id) )  {
    stop( paste(ensembl_gene_id[which(!valid_ensembl_id)], collapse = ", "), 
          " not valid ensembl gene id!\n")
  }
  
  # -- 1 -- Check ENSEMBL BiomaRt R Object ----
  if( verbose ) {
    message("Dataset: ", mRt_objct@dataset, "\n",
            "ENSEMBL host: ", mRt_objct@host) 
  }

  # -- 2 -- Retrieve some basic info for this gene ----
  things_to_retrieve <- c("ensembl_gene_id", "external_gene_name",
                          "chromosome_name", "start_position", "gene_biotype")
  
  goi_info <- getBM(attributes = things_to_retrieve,
                    filters = 'ensembl_gene_id',
                    values = ensembl_gene_id,
                    mart = mRt_objct,
                    verbose = verbose,
                    uniqueRows = T,
                    quote = "\'", 
                    useCache = TRUE)
  
  if ( nrow( goi_info ) == 0 ) {
    warning("Could not retrieve info for gene '", ensembl_gene_id, "'.\n",
            "Maybe it is not a valid gene name?")
    return("Error")
  } 
  
  # -- 3 -- Return some or all info ----
  if ( only_gene_name == TRUE ) {
    gene_name <- goi_info$external_gene_name
    return(gene_name)
  } else if ( only_gene_name == FALSE){
    return(goi_info)
  } else {
    stop("only_gene_name must be a logical, i.e. either TRUE or FALSE")
  }
}

#' Convenient converter from gene name to ensembl ID
#'
#' @param gene_name A character specifying the gene name
#' @param only_ensembl_id Logical, whether to return only the gene name or a `data.frame` with more info. Default `TRUE`
#' @param verbose Logical. Do you want me to be chatty?
#' @param mRt_objct A `biomaRt` object created with the `gimme_mart()` function.
#'
#' @return
#' @import biomaRt
#' @export
#'
#' @examples
#' ensembl <- gimme_mart()
#' gene_name_2_ensembl_id("TSPAN6")
#' "ENSG00000000003"
#' 
#' # Example with mouse ensembl
#' ensembl <- gimme_mart(species = "mmusculus", latest = T, use_mirror = F, 
#'                       out_dir = F, verbose = T)
#' gene_name_2_ensembl_id(gene_name = "TSPAN6")
#' "ENSMUSG00000067377"
gene_name_2_ensembl_id <- function(gene_name, only_ensembl_id = TRUE,
                                    verbose = FALSE, mRt_objct = ensembl) {
  
  stopifnot(is.character(gene_name))
  
  # -- 0 -- Check ENSEMBL BiomaRt R Object ----
  if( verbose ) {
    message("Dataset: ", mRt_objct@dataset, "\n",
            "ENSEMBL host: ", mRt_objct@host) 
  }
  
  # -- 1 -- Retrieve some basic info for this gene ----
  things_to_retrieve <- c("ensembl_gene_id", "external_gene_name",
                          "chromosome_name", "start_position", "gene_biotype")

  goi_info <- getBM(attributes = things_to_retrieve,
                    filters = 'external_gene_name',
                    values = gene_name,
                    mart = mRt_objct,
                    verbose = verbose,
                    uniqueRows = T,
                    quote = "\'", 
                    useCache = TRUE)
  
  # }
  
  
  if ( nrow( goi_info ) == 0 ) {
    warning("Could not retrieve info for gene '", ensembl_gene_id, "'.\n",
            "Maybe it is not a valid gene name?")
    return("Error")
  } 
  
  # -- 2 -- Return some or all info ----
  if ( only_ensembl_id == TRUE ) {
    
    ensembl_id <- goi_info$ensembl_gene_id
    ensembl_gene_id_regex <- "ENS[A-Z]+[0-9]{11}"
    valid_ensembl_id <- grepl(pattern = ensembl_gene_id_regex, x = ensembl_id,
                              ignore.case = F)
    
    if (!valid_ensembl_id) {
      warning(ensembl_gene_id, " is not a valid ensembl gene id!\n")
    }
    
    return(ensembl_id)
    
  } else if ( only_ensembl_id == FALSE){
    return(goi_info)
  } else {
    stop("only_ensembl_id must be a logical, i.e. either TRUE or FALSE")
  }
}


#' Convenient converter from ensembl ID to entrez ID
#'
#' @param ensembl_gene_id A character specifying the ENSEMBL gene ID
#' @param only_entrez Logical, whether to return only the gene name or a `data.frame` with more info. Default `TRUE`
#' @param verbose Logical. Do you want me to be chatty?
#' @param mRt_objct A `biomaRt` object created with the `gimme_mart()` function.
#'
#' @return
#' @import biomaRt
#' @export
#'
#' @examples
#' # Example with human ensembl
#' ensembl <- gimme_mart()
#' ensembl_id_2_gene_name("ENSG00000000003")
#' "7105"
#' # Example with mouse ensembl
#' ensembl <- gimme_mart(species = "mmusculus", latest = T, use_mirror = F, 
#'                       out_dir = F, verbose = T)
#' ensembl_id_2_entrez("ENSMUSG00000067377")
#' "56496"
ensembl_id_2_entrez <- function(ensembl_gene_id, only_entrez = TRUE,
                                 verbose = FALSE, mRt_objct = ensembl) {
  # -- 0 -- Check input parameters ----
  stopifnot(is.character(ensembl_gene_id))
  ensembl_gene_id_regex <- "ENS[A-Z]+[0-9]{11}"
  valid_ensembl_id <- grepl(pattern = ensembl_gene_id_regex, 
                            x = ensembl_gene_id, ignore.case = F)
  
  if ( !all(valid_ensembl_id) )  {
    stop( paste(ensembl_gene_id[which(!valid_ensembl_id)], collapse = ", "), 
          " not valid ensembl gene id!\n")
  }
  
  # -- 1 -- Check ENSEMBL BiomaRt R Object ----
  if( verbose ) {
    message("Dataset: ", mRt_objct@dataset, "\n",
            "ENSEMBL host: ", mRt_objct@host) 
  }
  
  
  # -- 2 -- Retrieve some basic info for this gene ----
  things_to_retrieve <- c("ensembl_gene_id", "external_gene_name",
                          "entrezgene_id", "chromosome_name", "gene_biotype" )
  
  goi_info <- getBM(attributes = things_to_retrieve,
                    filters = 'ensembl_gene_id',
                    values = ensembl_gene_id,
                    mart = mRt_objct,
                    verbose = verbose,
                    uniqueRows = T,
                    quote = "\'", 
                    useCache = TRUE)
  
  if ( nrow( goi_info ) == 0 ) {
    warning("Could not retrieve info for gene '", ensembl_gene_id, "'.\n",
            "Maybe it is not a valid gene name?")
    return("Error")
  } 
  
  # coerc entrez gene id  (now NCBI gene id)
  goi_info$entrezgene_id <- as.character(goi_info$entrezgene_id)
  
  # -- 3 -- Return some or all info ----
  if ( only_entrez == TRUE ) {
    entrezgene_id <- goi_info$entrezgene_id
    return(entrezgene_id)
  } else if ( only_entrez == FALSE){
    return(goi_info)
  } else {
    stop("only_entrez must be a logical, i.e. either TRUE or FALSE")
  }
}

