#' Handy function to create a biomaRt R object based on the ENSEMBL database to be used when querying biomart
#'
#' @description A function to get a biomaRt object during an R analysis that can later on be used for retrieving gene IDs, coordinates and more.
#' 
#' @param species a character specifying the species to work with. Currently the supported species are:
#' \itemize{
#' \item{\code{hsapiens} - Human (Homo Sapiens) Default}
#' \item{\code{mmusculus} - Mouse (Mus Musculus)}
#' \item{\code{rnorvegicus} - Rat (Rattus Norvegicus)}
#' \item{\code{ggallus} - Chicken (Gallus Gallus)}
#' } 
#' @param version Which database version of biomart do you want? Use:
#' \itemize{
#' \item{\code{latest} - for the latest version of the database. }
#' \item{\code{penultimate} - for the one before the latest version of the database. }
#' \item{\code{number} - a specific version of the biomart database. If the version you query is not available an error is returned.}
#' }  
#' @param use_mirror Logical. Do you want to use a mirror for creating the biomart object? Default `FALSE`. If `TRUE` the parameter version will be converted to `TRUE`.
#' @param which_mirror Specify which biomart mirror to use. Chose from:
#' \itemize{
#' \item{\code{www} - UK (Sanger Institute). in Cambridge, default}
#' \item{\code{uswest} - US West (Amazon AWS) on West Coast of US of America}
#' \item{\code{useast} - US East (Amazon AWS) on East Coast of US of America}
#' \item{\code{asia} - Asia (Amazon AWS) in Singapore}
#' } 
#' @param out_dir Either `FALSE` (default) or a character specifying a path where to save the biomaRt object. If directory does not exist it will be created.
#' @param verbose Logical. Print the version and genome assembly for the species queried. Default `TRUE`.
#'
#' @return A biomaRt ensembl object
#' @importFrom biomaRt useMart listDatasets listEnsemblArchives useEnsembl
#' @export
#'
#' @examples
#' # Human ENSEMBL biomaRt object
#' human_ensembl <- gimme_mart()
gimme_mart <- function(species = "hsapiens", version = 'latest', verbose = TRUE,
                       out_dir = FALSE, use_mirror = FALSE, which_mirror) {
  # ---- Check input parameters ----
  supported_speciess <- c("hsapiens", "mmusculus", "rnorvegicus", "ggallus")
  if (! species %in% supported_speciess ) {
    stop("The species you selected ", species, " is not yet supported. ",
         "Use first letter and second full name of the species scientific ",
         " name like 'hsapiens' for Homo Sapiens or 'mmusculus' for ",
         "Mus Musculus.")
  }
  
  if ( is.null(out_dir) ) {
    stop("You didn't specify the out_dir parameter! It must be FALSE to ",
         "keep the biomart object in the environment, or a character (path) ",
         " to specify where to save the biomart object as a R object.")
  }
  
  if (out_dir == TRUE) {
    stop("out_dir parameter is either FALSE or a character specifying the path.")
  }
  
  # -- If using a mirror server the version can only be the latest
  if ( use_mirror == TRUE ) {
    if ( is.null(which_mirror) ) {
      stop("You want to use a mirror, but need to specify one from: ",
           "'www', 'uswest', 'useast', 'asia'.")
    } else if ( which_mirror %in% c('www', 'uswest', 'useast', 'asia') ) {
      if (version != 'latest') {
        warning("When specifying a mirror the database version is set to 'latest'!")
        version <- 'latest'
      }
    } else {
      stop("use_mirror must be either TRUE or FALSE" )
    }
  }
  
  # ---- Get available archived versions of Ensembl ----
  DB_versions <- listEnsemblArchives()
  
  # If retrieve Ensembl Archive fails because of Mart connection
  if ( !exists("DB_versions") ) {
    stop("Can't retrieve ENSEMBL Database Archives info! ",
         "This is usually caused by a connection error or biomart downtime")
  }
  
  # ---- Get db version ----
  if( is.numeric(version) == TRUE | is.integer(version) == TRUE ) {
    
    DB_v_available <- as.integer(DB_versions$version[-1])
    
    if( version < min(DB_v_available) ) {
      stop("The version you required (v", version, ") is no longer supported! ",
           "Increase to at least version", min(DB_v_available))
    }
    
    if( version > max(DB_v_available) ) {
      stop("The version you required (v", version, ") is not yet released! ",
           "Decrease to at least version ", max(DB_v_available) )
    }
    
    if ( version %in% DB_v_available ) {
      DB_v <- as.integer(version)
      # -- Get db version release date
      DB_v_date <- DB_versions[which(DB_versions$version == DB_v), "date"]
      
    } else{
      stop("The specified version is not currently available! ",
           "Pick one of these available versions: ",
           paste0(DB_v_available, sep = ", "))
    }
    
  } else if ( version == 'latest' ) {
    # -- latest version
    DB_v <- DB_versions[grep("\\*" , DB_versions$current_release), "version"]
    DB_v <- as.integer(DB_v)
    # -- Get db version release date
    DB_v_date <- DB_versions[grep("\\*" , DB_versions$current_release), "date"]
    
  } else if ( version == 'penultimate' ) {
    # -- penultimate version
    DB_v <- DB_versions[grep("\\*" , DB_versions$current_release) + 1, "version"]
    DB_v <- as.integer(DB_v)
    # -- Get db version release date
    DB_v_date <- DB_versions[grep("\\*" , DB_versions$current_release) + 1, "date"]
  } else { 
    stop("The 'version' argument must be either 'latest' or 'penultimate' or an integer", 
         "number between ", max(DB_v_available), " and ", min(DB_v_available),
         " And 'use_mirror' can be used only with version = 'latest' ")
  }
  
  if ( !exists("DB_v") ) {
    stop("I couldn't get the database version you wanted! ",
         "The 'version' argument must be either 'latest' or FALSE or an integer", 
         "number between ", max(DB_v_available), " and ", min(DB_v_available),
         " And 'use_mirror' can be used only with version = 'latest' ")
    
  }
  
  DB_v_host_url <- DB_versions[which(DB_versions$version == DB_v), "url"]
  
  ensembl <- useMart(biomart = "ensembl", host = DB_v_host_url)
  datasets <- listDatasets(ensembl)
  # ---- Get species dataset name
  species_ensembl <- datasets[grep(species, datasets$dataset), "dataset" ]
  
  if ( length(species_ensembl) == 0 ) {
    stop("The species you selected ", species, " is wrong or ", 
         "not supported!")
  }
  # ---- Get genome assembly name version
  genome_assembly_version <- datasets[grep(species, datasets$dataset), "version" ]
  
  # ---- Print info ----
  if (verbose == TRUE) {
    message("Creating Mart object:\n",
            "version: ", DB_v, "\n", 
            "released in: ", DB_v_date, "\n",
            "dataset: ", species_ensembl, "\n",
            "genome assembly: ", genome_assembly_version )
  }
  
  # ---- Create ENSEMBL BioMart Object ----
  if (use_mirror == TRUE) {
    # -- Ensembl Dataset from a mirrored server
    ensembl <- useEnsembl(biomart = "ensembl", 
                          dataset = species_ensembl,
                          mirror  = which_mirror)
    
  } else if ( use_mirror == FALSE) {
    # -- Ensembl Dataset on main host
    ensembl <- useEnsembl(biomart = "ensembl", 
                          dataset = species_ensembl, 
                          version = DB_v, 
                          host = DB_v_host_url)
  }
  
  # ---- Write output or keep in the environment ---
  if ( is.character(out_dir) ) {
    
    # ---- Output Directory
    out_biomart_dir <- file.path(out_dir, species)
    # ---- Check if output dir exists, if not create it
    if (!dir.exists(out_biomart_dir)) { dir.create(out_biomart_dir) }
    
    if (use_mirror == FALSE) {
      DB_v_date <- sub(" ", "_", DB_v_date)
      output_ensembl_name <- paste0(out_biomart_dir, "/", 
                                    format(Sys.Date(), "%Y_%m_%d"), "_",
                                    species_ensembl, "_v", DB_v, "_",
                                    DB_v_date, ".rds")
      
    } else if ( use_mirror == TRUE) {
      output_ensembl_name <- paste0(out_biomart_dir, "/", 
                                    format(Sys.Date(), "%Y_%m_%d"), "_",
                                    species_ensembl, "_mirror_", which_mirror,
                                    ".rds")
    } else {
      stop("Error occurred when preparing to write to file the biomart object.")
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
#' @return A character corresponding to an ensembl gene name or a `data.frame`.
#' @importFrom biomaRt getBM 
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
                    verbose = FALSE,
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
#' @return A character corresponding to an ensembl gene ID or a `data.frame`.
#' @importFrom biomaRt getBM 
#' @export
#' @description The ENSEMBL gene IDs are returned in the same order as the input gene names.
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
                    verbose = FALSE,
                    uniqueRows = T,
                    quote = "\'", 
                    useCache = TRUE)
  
  if ( nrow( goi_info ) == 0 ) {
    warning("Could not retrieve info for gene '", ensembl_gene_id, "'.\n",
            "Maybe it is not a valid gene name?")
    return("Error")
  } 
  
  # -- 2 -- Sort the genes using the input order ----
  # use factors to define the input order
  goi_info$external_gene_name <- factor(goi_info$external_gene_name, gene_name) 
  # order by external gene name
  goi_info <- goi_info[order(goi_info$external_gene_name), ]
  # turn external gene name back to character
  goi_info$external_gene_name <- as.character(goi_info$external_gene_name)
  
  # -- 3 -- Return some or all info ----
  if ( only_ensembl_id == TRUE ) {
    
    ensembl_id <- goi_info$ensembl_gene_id
    ensembl_gene_id_regex <- "ENS[A-Z]+[0-9]{11}"
    valid_ensembl_id <- grepl(pattern = ensembl_gene_id_regex, x = ensembl_id,
                              ignore.case = F)
    
    if (!all(valid_ensembl_id)) {
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
#' @return A character corresponding to an entrez ID or a `data.frame`.
#' @importFrom biomaRt getBM 
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
                    verbose = FALSE,
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

