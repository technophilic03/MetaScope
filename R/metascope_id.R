
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("source_name", "species", "genus", "family", "order", 
                           "class", "phylum", "domain", "count", ".", "v", "n", ":=", 
                           "lowest_unique_taxonomy_level", "final_genome_name", "Genome", 
                           "read_count", "TaxonomyID", "Proportion", "readsEM", "EMProportion"))
}


get_max_index_matrix <- function(mat) {
  stopifnot(inherits(mat, "Matrix"))  # Ensure sparse Matrix
  # Convert to triplet format for efficient row-wise max
  trip <- Matrix::summary(mat)
  row_max <- tapply(trip$x, trip$i, max)
  # Keep only entries that match the row max
  keep <- trip$x == row_max[as.character(trip$i)]
  # Build a new sparse logical matrix marking max positions
  is_max <- Matrix::sparseMatrix(i = trip$i[keep],
               j = trip$j[keep],
               x = TRUE,
               dims = dim(mat))
  return(is_max)
}

obtain_reads <- function(input_file, input_type, aligner, quiet) {
  to_pull <- c("qname", "rname", "cigar", "qwidth", "pos")
  if (identical(input_type, "bam")) {
    if (!quiet) message("Reading .bam file: ", input_file)
    if (identical(aligner, "bowtie2")) {
      params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("AS"))
    } else if (identical(aligner, "subread")) {
      params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("NM"))
    } else if (identical(aligner, "other")) {
      params <- Rsamtools::ScanBamParam(what = to_pull)
    }
    reads <- Rsamtools::scanBam(input_file, param = params)
  } else if (identical(input_type, "csv.gz")) {
    if (!quiet) message("Reading .csv.gz file: ", input_file)
    reads <- data.table::fread(input_file, sep = ",", header = FALSE) %>%
      magrittr::set_colnames(c("qname", "rname", "cigar", "qwidth", "pos", "tag")) %>% as.list() %>% list()
    if (identical(aligner, "bowtie2")) {
      reads[[1]]$tag <- list("AS" = reads[[1]]$tag)
    } else if (identical(aligner, "subread")) {
      reads[[1]]$tag <- list("NM" = reads[[1]]$tag)
    }
  }
  return(reads)
}

identify_rnames <- function(reads, unmapped = NULL) {
  reads_in <- reads[[1]]$rname
  if(!is.null(unmapped)) reads_in <- reads[[1]]$rname[!unmapped]
  # FOR 2015 - if >50% have 8 pipe symbols
  if (mean(stringr::str_count(reads_in, "\\|") == 8) > 0.5) {
    mapped_rname <- stringr::str_split(reads_in, "\\|") |>
      plyr::laply(function(x) magrittr::extract(x, 8))
    return(mapped_rname)
  }
  # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole
  prefixes <- c("AC", "NC", "NG", "NT", "NW", "NZ", "NR") %>%
    paste0("_") %>%
    paste(collapse = "|")
  this_read <- stringr::str_split_i(reads_in, prefixes, -1) |>
    stringr::str_split_i("[| ]", 1)
  read_prefix <- stringr::str_extract(reads_in, prefixes)
  mapped_rname <- paste0(read_prefix, this_read)
  return(mapped_rname)
}


find_accessions <- function(accessions, accession_path, quiet) {
  # Convert accessions to taxids and get genome names
  if (!quiet) message("Obtaining taxonomy and genome names")
  taxids <- taxonomizr::accessionToTaxa(accessions, sqlFile = accession_path)
  genome_names <- apply(taxonomizr::getTaxonomy(taxids, sqlFile = accession_path,
                                                desiredTaxa = c("superkingdom", "phylum", "class",
                                                                "order", "family", "genus", "species", "strain")),
                        1,function(x) paste0(x, collapse = ";"))
  return(genome_names)
}

get_alignscore <- function(aligner, cigar_strings, count_matches, scores,
                           qwidths) {
  #Subread alignment scores: CIGAR string matches - edit score
  if (identical(aligner, "subread")) {
    num_match <- unlist(vapply(cigar_strings, count_matches,
                               USE.NAMES = FALSE, double(1)))
    alignment_scores <- num_match - scores
    scaling_factor <- 1.0 / max(alignment_scores)
    relative_alignment_scores <- alignment_scores - min(alignment_scores)
    exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
  } else if (identical(aligner, "bowtie2")) {
    # Bowtie2 alignment scores: AS value + read length (qwidths)
    alignment_scores <- scores + qwidths
    scaling_factor <- 1.0 / max(alignment_scores)
    relative_alignment_scores <- alignment_scores - min(alignment_scores)
    exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
  } else if (identical(aligner, "other")) {
    # Other alignment scores: No assumptions
    exp_alignment_scores <- 1
  }
  return(exp_alignment_scores)
}


get_assignments <- function(combined, convEM, maxitsEM, unique_taxids,
                            unique_genome_names, update_bam = TRUE, input_file, priors_df, quiet) {
  combined$index <- seq.int(1, nrow(combined))
  input_distinct <- dplyr::distinct(combined, .data$qname, .data$rname,
                                    .keep_all = TRUE)
  qname_inds_2 <- input_distinct$qname
  rname_tax_inds_2 <- input_distinct$rname
  scores_2 <- input_distinct$scores
  non_unique_read_ind <- unique(combined[[1]][(
    duplicated(input_distinct[, 1]) | duplicated(input_distinct[, 1],
                                                 fromLast = TRUE))])
  # 1 if read is multimapping, 0 if read is unique
  y_ind_2 <- as.numeric(unique(input_distinct[[1]]) %in% non_unique_read_ind)
  gammas <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = scores_2)
  pi_old <- rep(1 / nrow(gammas), ncol(gammas))
  pi_new <- theta_new <- Matrix::colMeans(gammas)
  if (!is.null(priors_df)) {
    if (colnames(priors_df)[1] != "species") {
      stop("priors_df must have two columns. \n",
           "The first column is named 'species' and contains the species names. \n",
           "The second column is named 'prior_weights' and contains the weights.")
    }
    colnames(priors_df)[2] <- "prior_weights"
    weights <- tidyr::as_tibble(unique_genome_names) |>
      dplyr::rename("species" = "value") |>
      dplyr::left_join(priors_df, by = "species") |>
      tidyr::replace_na(replace = list(prior_weights = 0))
    unique_reads <- weights$prior_weights * max(combined$qname)
    exp_weights <- exp(weights$prior_weights)
    #weighted_pi_new <- pi_new + weights / (1 + weights)
    posterior <- pi_new * exp_weights
    pi_new <- posterior
  } else {
    unique_reads = 0
  }
  epsilon = 1e-8
  conv <- max(abs(pi_new - pi_old) / (pi_old + epsilon))
  it <- 0
  if (!quiet) message("Starting EM iterations")
  while (conv > convEM && it < maxitsEM) {
    # Expectation Step: Estimate expected value for each read to ea genome
    pi_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                   x = pi_new[rname_tax_inds_2])
    theta_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                      x = theta_new[rname_tax_inds_2])
    weighted_gamma <- gammas * pi_mat * theta_mat
    weighted_gamma_sums <- Matrix::rowSums(weighted_gamma)
    gammas_new <- weighted_gamma / weighted_gamma_sums
    # Maximization step: proportion of reads to each genome
    pi_new <- (Matrix::colSums(gammas_new) + unique_reads) / (sum(Matrix::colSums(gammas_new)) + sum(unique_reads))
    theta_new_num <- (Matrix::colSums(y_ind_2 * gammas_new) + 1)
    theta_new <- theta_new_num / (nrow(gammas_new) + 1)
    # Check convergence
    it <- it + 1
    conv <- max(abs(pi_new - pi_old) / (pi_old + epsilon), na.rm = TRUE)
    pi_old <- pi_new
    if (!quiet) message(c(it, " ", conv))
  }
  if (!quiet) message("DONE! Converged in ", it, " iterations.")
  hit_which <- get_max_index_matrix(gammas_new)
  combined <- combined |>
    dplyr::mutate("hit" = hit_which[cbind(!!rlang::sym("qname"), !!rlang::sym("rname"))])
  combined_single <- combined |>
    dplyr::group_by(!!rlang::sym("qname"), !!rlang::sym("rname")) |>
    dplyr::mutate("best_hit" = !!rlang::sym("hit") &
                    !!rlang::sym("scores") == max(!!rlang::sym("scores"))) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!rlang::sym("qname")) |>
    dplyr::mutate("single_hit" = dplyr::row_number() == min(dplyr::row_number()[!!rlang::sym("best_hit")]))
  combined_distinct <- dplyr::distinct(combined, !!rlang::sym("qname"),
                                       !!rlang::sym("rname"), .keep_all = TRUE) |>
    dplyr::filter(!!rlang::sym("hit") == TRUE)
  best_hit <- Matrix::colSums(hit_which)
  names(best_hit) <- seq_along(best_hit)
  best_hit <- best_hit[best_hit != 0]
  hits_ind <- as.numeric(names(best_hit))
  final_taxids <- unique_taxids[hits_ind]
  final_genomes <- unique_genome_names[hits_ind]
  proportion <- best_hit / sum(best_hit)
  gammasums <- Matrix::colSums(gammas_new)
  readsEM <- round(gammasums[hits_ind], 1)
  propEM <- gammasums[hits_ind] / sum(gammas_new)
  results_tibble <- dplyr::tibble(TaxonomyID = final_taxids, Genome = final_genomes,
                                  read_count = best_hit, Proportion = proportion,
                                  readsEM = readsEM, EMProportion = propEM,
                                  hits_ind = hits_ind) %>%
    dplyr::arrange(dplyr::desc(.data$read_count))
  if (!quiet) message("Found reads for ", nrow(results_tibble), " genomes")

  return(list(results_tibble, combined_distinct, combined_single))
}

#' Count the number of base lengths in a CIGAR string for a given operation
#'
#' The 'CIGAR' (Compact Idiosyncratic Gapped Alignment Report) string is how the
#' SAM/BAM format represents spliced alignments. This function will accept a
#' CIGAR string for a single read and a single character indicating the
#' operation to be parsed in the string. An operation is a type of column that
#' appears in the alignment, e.g. a match or gap. The integer following the
#' operator specifies a number of consecutive operations. The
#' \code{count_matches()} function will identify all occurrences of the operator
#' in the string input, add them, and return an integer number representing the
#' total number of operations for the read that was summarized by the input
#' CIGAR string.
#'
#' This function is best used on a vector of CIGAR strings using an apply
#' function (see examples).
#'
#' @param x Character. A CIGAR string for a read to be parsed. Examples of
#' possible operators include "M", "D", "I", "S", "H", "=", "P", and "X".
#' @param char A single letter representing the operation to total for the
#' given string.
#'
#' @return an integer number representing the total number of alignment
#' operations for the read that was summarized by the input CIGAR string.
#'
#' @export
#'
#' @examples
#' # A single cigar string: 3M + 3M + 5M
#' cigar1 <- "3M1I3M1D5M"
#' count_matches(cigar1, char = "M")
#'
#' # Parse with operator "P": 2P
#' cigar2 <- "4M1I2P9M"
#' count_matches(cigar2, char = "P")
#'
#' # Apply to multiple strings: 1I + 1I + 5I
#' cigar3 <- c("3M1I3M1D5M", "4M1I1P9M", "76M13M5I")
#' vapply(cigar3, count_matches, char = "I",
#'        FUN.VALUE = numeric(1))
#'

count_matches <- function(x, char = "M") {
  if (length(char) != 1) {
    stop("Please provide a single character ",
         "operator with which to parse.")
  } else if (length(x) != 1) {
    stop("Please provide a single CIGAR string to be parsed.")
  }
  pattern <- paste0("\\d+", char)
  ind <- gregexpr(pattern, x)[[1]]
  start <- as.numeric(ind)
  end <- start + attr(ind, "match.length") - 2
  out <- cbind(start, end) %>% apply(
    1, function(y) substr(x, start = y[1], stop = y[2])) %>%
    as.numeric() %>% sum()
  return(data.table::fifelse(is.na(out[1]), yes = 0, no = out[1]))
}

#' Helper Function for MetaScope ID
#'
#' Used to create plots of genome coverage for any number of accession numbers
#'
#' @param which_taxid Which taxid to plot
#' @param which_genome Which genome to plot
#' @param accessions List of accessions from \code{metascope_id()}
#' @param taxids List of accessions from \code{metascope_id()}
#' @param reads List of reads from input file
#' @param out_base The basename of the input file
#' @param out_dir The path to the input file
#'
#' @return A plot of the read coverage for a given genome

locations <- function(which_taxid, which_genome,
                      accessions, taxids, reads, out_base, out_dir) {
  plots_save <- file.path(out_dir, paste(out_base, "cov_plots",
                                         sep = "_"))
  if (!dir.exists(plots_save)) dir.create(plots_save)
  # map back to accessions
  choose_acc <- paste(accessions[which(as.numeric(taxids) %in% which_taxid)])
  # map back to BAM
  map2bam_acc <- which(reads[[1]]$rname %in% choose_acc)
  # Split genome name to make digestible
  this_genome <- strsplit(which_genome, " ")[[1]][c(1, 2)]
  use_name <- paste(this_genome, collapse = " ") %>% stringr::str_replace(",", "")
  coverage <- round(mean(seq_len(338099) %in% unique(
    reads[[1]]$pos[map2bam_acc])), 3)
  # Plotting
  dfplot <- dplyr::tibble(x = reads[[1]]$pos[map2bam_acc])
  ggplot2::ggplot(dfplot, ggplot2::aes(.data$x, fill = ggplot2::after_stat(count))) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::scale_fill_gradient(low = 'red', high = 'yellow') +
    ggplot2::theme_classic() +
    ggplot2::labs(title = bquote("Positions of reads mapped to"~italic(.(use_name))),
                  xlab = "Aligned position across genome (leftmost read position)",
                  ylab = "Read Count",
                  caption = paste0("Accession Number: ", choose_acc))

  ggplot2::ggsave(paste0(plots_save, "/",
                         stringr::str_replace(use_name, " ", "_"), ".png"),
                  device = "png")
}

#' Identify which genomes are represented in a processed sample
#'
#' This function will read in a .bam or .csv.gz file, annotate the taxonomy and
#' genome names, reduce the mapping ambiguity using a mixture model, and output
#' a .csv file with the results. Currently, it assumes that the genome
#' library/.bam files use NCBI accession names for reference names (rnames in
#' .bam file).
#'
#' @param input_file The .bam or .csv.gz file of sample reads to be identified.
#' @param input_type Extension of file input. Should be either "bam" or
#'   "csv.gz". Default is "csv.gz".
#' @param aligner The aligner which was used to create the bam file. Default is
#'   "bowtie2" but can also be set to "subread" or "other".
#' @param db Currently accepts one of \code{c("ncbi", "silva", "other")} Default
#'   is \code{"ncbi"}, appropriate for samples aligned against indices compiled
#'   from NCBI whole genome databases. Alternatively, usage of an alternate
#'   database (like Greengenes2) should be specified with \code{"other"}.
#' @param db_feature_table If \code{db = "other"}, a data.frame must be supplied
#'   with two columns, "Feature ID" matching the names of the alignment indices,
#'   and a second \code{character} column supplying the taxon identifying
#'   information.
#' @param accession_path (character) Filepath to NCBI accessions SQL database.
#'   See \code{taxonomzr::prepareDatabase()}.
#' @param out_dir The directory to which the .csv output file will be output.
#'   Defaults to \code{dirname(input_file)}.
#' @param convEM The convergence parameter of the EM algorithm. Default set at
#'   \code{1/10000}.
#' @param maxitsEM The maximum number of EM iterations, regardless of whether
#'   the convEM is below the threshhold. Default set at \code{50}. If set at
#'   \code{0}, the algorithm skips the EM step and summarizes the .bam file 'as
#'   is'.
#' @param num_species_plot The number of genome coverage plots to be saved.
#'   Default is \code{NULL}, which saves coverage plots for the ten most highly
#'   abundant species.
#' @param update_bam Whether to update BAM file with new read assignments.
#'   Default is \code{FALSE}. If \code{TRUE}, requires \code{input_type = "bam"}
#'   such that a BAM file is the input to the function.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#' @param tmp_dir Path to a directory to which bam and updated bam files can be
#'   saved. Required.
#' @param priors_df data.frame containing priors data. The data.frame consists
#'   of two columns, 'species' containing species name, and 'prior_weights'
#'   containing the prior weights (as a percent; integer).
#' @param group_by_taxa Character. Taxonomy level at which accessions should be
#'   grouped. Defaults to \code{"species"}
#' @param force_calls Boolean. Determines whether or not to force down to 
#'   group_by_taxa calls or to make ambiguous no calls when uncertain. Recommended
#'   FALSE for 16S data and TRUE for metagenomic data. 
#'
#' @return This function exports a .csv file with annotated read counts to
#'   genomes with mapped reads to the location returned by the function.
#'   Depending on the parameters specified, can also output an updated BAM
#'   file, and fasta files for additional analysis downstream.
#'
#' @importFrom rlang :=
#' 
#' @export
#'
#' @examples
#' #### Align reads to reference library and then apply metascope_id()
#' ## Assuming filtered bam files already exist
#'
#' ## Create temporary directory
#' file_temp <- tempfile()
#' dir.create(file_temp)
#'
#' ## Get temporary accessions database
#' tmp_accession <- system.file("extdata", "example_accessions.sql", package = "MetaScope")
#'
#' #### Subread aligned bam file
#'
#' ## Create object with path to filtered subread csv.gz file
#' filt_file <- "subread_target.filtered.csv.gz"
#' bamPath <- system.file("extdata", filt_file, package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' ## Run metascope id with the aligner option set to subread
#' metascope_id(input_file = file.path(file_temp, filt_file),
#'              aligner = "subread", num_species_plot = 0,
#'              input_type = "csv.gz", accession_path = tmp_accession)
#'
#' #### Bowtie 2 aligned .csv.gz file
#'
#' ## Create object with path to filtered bowtie2 bam file
#' bowtie_file <- "bowtie_target.filtered.csv.gz"
#' bamPath <- system.file("extdata", bowtie_file, package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' ## Run metascope id with the aligner option set to bowtie2
#' metascope_id(file.path(file_temp, bowtie_file), aligner = "bowtie2",
#'              num_species_plot = 0, input_type = "csv.gz",
#'              accession_path = tmp_accession)
#'
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'

metascope_id <- function(input_file, input_type = "csv.gz",
                         aligner = "bowtie2",
                         db = "ncbi",
                         db_feature_table = NULL,
                         accession_path = NULL,
                         priors_df = NULL,
                         tmp_dir = dirname(input_file),
                         out_dir = dirname(input_file),
                         convEM = 1 / 10000,
                         maxitsEM = 25,
                         update_bam = FALSE,
                         num_species_plot = NULL,
                         group_by_taxa = "species",
                         force_calls = TRUE, 
                         quiet = TRUE)  {

  out_base <- input_file %>% base::basename() %>% strsplit(split = "\\.") %>%
    magrittr::extract2(1) %>% magrittr::extract(1)
  out_file <- file.path(out_dir, paste0(out_base, ".metascope_id.csv"))
  # Check to make sure valid aligner is specified
  if (!(aligner %in% c("bowtie2", "subread", "other"))) {
    stop("Please make sure aligner is set to either 'bowtie2', 'subread',",
         " or 'other'")
  }
  if (db == "other" && is.null(db_feature_table)) {
    stop("Please supply a data.frame for db_feature_table if 'db = other'")
  }
  if (input_type == "csv.gz") {
    if (!quiet) message("Note, cannot generate updated_bam from csv.gz file")
    update_bam <- FALSE
    if (update_bam) stop("`update_bam` is set to TRUE but input is a `csv.gz`")
  }
  # Check that tmp_dir is specified
  if (is.null(tmp_dir)) stop("Please supply a directory for 'tmp_dir' argument.")
  reads <- obtain_reads(input_file, input_type, aligner, quiet)
  unmapped <- is.na(reads[[1]]$rname)
  if (db == "ncbi") reads[[1]]$rname <- identify_rnames(reads)
  mapped_rname <- as.character(reads[[1]]$rname[!unmapped])
  mapped_qname <- reads[[1]]$qname[!unmapped]
  mapped_cigar <- reads[[1]]$cigar[!unmapped]
  mapped_qwidth <- reads[[1]]$qwidth[!unmapped]
  if (aligner == "bowtie2") {
    # mapped alignments used
    map_edit_or_align <- reads[[1]][["tag"]][["AS"]][!unmapped]
  } else if (aligner == "subread") {
    map_edit_or_align <-
      reads[[1]][["tag"]][["NM"]][!unmapped] # mapped edits used
  }
  read_names <- unique(mapped_qname)
  accessions <- as.character(unique(mapped_rname))

  # Let user decide which taxa level to group by
  if (is.null(group_by_taxa)) {
    group_by_taxa <- "species"
  }
  if (db == "ncbi") {
    if (is.null(accession_path)) stop("Please provide a valid accession_path argument")
    taxonomy_indices <- tidyr::tibble(accessions) |>
      dplyr::mutate("taxids" = taxonomizr::accessionToTaxa(accessions, sqlFile = accession_path)) |>
      dplyr::mutate("taxa_names" = taxonomizr::getTaxonomy(taxids, sqlFile = accession_path,
                                                         desiredTaxa = group_by_taxa)[,]) |>
      dplyr::mutate("taxa_names" = ifelse(is.na(!!dplyr::sym("taxa_names")), paste0("unknown genome; accession ID is", accessions), !!dplyr::sym("taxa_names"))) |>
      dplyr::mutate("taxa_index" = match(!!dplyr::sym("taxa_names"), unique(!!dplyr::sym("taxa_names"))))

  } else if (db == "silva") {
    tax_id_all <- stringr::str_split(accessions, ";", n =2)
    taxids <- sapply(tax_id_all, `[[`, 1)
    genome_names <- sapply(tax_id_all, `[[`, 2)
    # Fix names
    mapped_rname <- stringr::str_split(mapped_rname, ";", n = 2) %>%
      sapply(`[[`, 1)
    accessions <- as.character(unique(mapped_rname))
  } else if (db == "other") {
    tax_id_all <- dplyr::tibble(`Feature ID` = accessions) %>%
      dplyr::left_join(db_feature_table, by = "Feature ID")
    taxids <- tax_id_all %>% dplyr::select(1) %>% unlist() %>% unname()
    genome_names <- tax_id_all %>% dplyr::select(2) %>% unlist() %>%
      unname()
  }


  # Make an alignment matrix (rows: reads, cols: unique taxids)
  if (!quiet) message("Setting up the EM algorithm")
  qname_inds <- match(mapped_qname, read_names)
  if (db == "ncbi") {
    rname_tax_inds <- taxonomy_indices$taxa_index[match(mapped_rname, taxonomy_indices$accessions)]
    if (!quiet) message("\tFound ", length(taxonomy_indices$taxa_index),
                        " unique taxa")
    unique_tax_data <- taxonomy_indices |>
      dplyr::mutate("taxid_rank" = suppressWarnings(as.numeric(taxids))) |>
      dplyr::group_by(!!dplyr::sym("taxa_index")) |>
      dplyr::slice(which.min(ifelse(is.na(!!dplyr::sym("taxid_rank")),
                                    Inf, !!dplyr::sym("taxid_rank")))) |>
      dplyr::ungroup()

    unique_taxids <- unique_tax_data$taxids
    unique_genome_names <- unique_tax_data$taxa_names
  }
  # Order based on read names
  rname_tax_inds <- rname_tax_inds[order(qname_inds)]
  cigar_strings <- mapped_cigar[order(qname_inds)]
  qwidths <- mapped_qwidth[order(qname_inds)]

  if (aligner == "bowtie2") {
    # mapped alignments used
    scores <- map_edit_or_align[order(qname_inds)]
  } else if (aligner == "subread") {
    # mapped edits used
    scores <- map_edit_or_align[order(qname_inds)]
  } else if (aligner == "other") scores <- 1
  qname_inds <- sort(qname_inds)
  exp_alignment_scores <- get_alignscore(aligner, cigar_strings,
                                         count_matches, scores, qwidths)
  combined <- dplyr::bind_cols("qname" = qname_inds,
                               "rname" = rname_tax_inds,
                               "scores" = exp_alignment_scores)
  
  results <- get_assignments(combined, convEM, maxitsEM, unique_taxids,
                             unique_genome_names, quiet = quiet, priors_df = priors_df)
  
  ## If you don't want to force species calls, then compute which no calls:
  if (!force_calls) {
    combined_dt <- data.table::as.data.table(combined)
    combined_uniques <- combined_dt[, .SD[which.max(scores)], by = .(qname, rname)]

    rnames_tax_table <- taxonomizr::getTaxonomy(unique_taxids, accession_path) |>
      as.data.frame()
    rnames_tax_table$rname <- seq.int(nrow(rnames_tax_table))
    
    combined_tax <- dplyr::left_join(combined_uniques, rnames_tax_table, by = "rname") 
    
    
    lowest_unique_level <- function(dt, q_col = "qname",
                                    levels = c("species", "genus", "family", "order", "class", "phylum", "domain")) {
      data.table::setDT(dt)
      res <- rep(NA_character_, nrow(dt))
      for (lev in levels) {
        map <- unique(dt[, .(q = get(q_col), v = get(lev))])[!is.na(v)]
        if (nrow(map) == 0) next
        cnt <- map[, .(n = data.table::uniqueN(v)), by = q]
        uq_q <- cnt[n == 1, q]
        uniq_vals <- map[q %in% uq_q, unique(v)]
        idx <- is.na(res) & dt[[lev]] %in% uniq_vals
        res[idx] <- lev
      }
      dt[, lowest_unique_taxonomy_level := res]
      
      
      return(dt)
    }
    
    combined_tax_levels <- lowest_unique_level(combined_tax) |> 
      dplyr::select(-!!dplyr::sym("qname"), -!!dplyr::sym("rname"), -!!dplyr::sym("scores")) |>
      unique() |>
      dplyr::mutate(final_genome_name = dplyr::case_when(
        lowest_unique_taxonomy_level == "species" ~ species,
        lowest_unique_taxonomy_level == "genus"   ~ genus,
        lowest_unique_taxonomy_level == "family"  ~ family,
        lowest_unique_taxonomy_level == "order"   ~ order,
        lowest_unique_taxonomy_level == "class"   ~ class,
        lowest_unique_taxonomy_level == "phylum"  ~ phylum,
        lowest_unique_taxonomy_level == "domain"  ~ domain,
        TRUE ~ NA_character_
      ))
    
    results_tibble <- results[[1]]
    data.table::setDT(results_tibble)
    data.table::setDT(combined_tax_levels)
    
    # Construct a single 'source_name' from dt (species preferred, falling back)
    combined_tax_levels[, source_name := data.table::fcoalesce(species, genus, family, order, class, phylum, domain)]
    
    name_map <- unique(combined_tax_levels[!is.na(source_name), .(source_name, final_genome_name)])
    
    # Join by text name and replace
    results_tibble[name_map, Genome := data.table::fifelse(!is.na(final_genome_name), final_genome_name, Genome),
                   on = .(Genome = source_name)]
    
    
    metascope_id_file <- results_tibble[
        order(-read_count), 
        .(
          TaxonomyID   = TaxonomyID[1],  # one with max read_count
          read_count   = sum(read_count,   na.rm = TRUE),
          Proportion   = sum(Proportion,   na.rm = TRUE),
          readsEM      = sum(readsEM,      na.rm = TRUE),
          EMProportion = sum(EMProportion, na.rm = TRUE)
        ),
        by = .(Genome)
      ][order(-read_count)]
    
  }
  else {
    metascope_id_file <- results[[1]] %>% dplyr::select("TaxonomyID", "Genome",
                                                        "read_count", "Proportion",
                                                        "readsEM", "EMProportion")
  }
  
  
  utils::write.csv(metascope_id_file, file = out_file, row.names = FALSE)
  if (!quiet) message("Results written to ", out_file)

  if (update_bam) {
    # Convert accessions back into taxids. Some taxids may be from different accessions
    # This step may produce false alignments to accessions of the same taxid
    if (db == "ncbi") {
      accessions_taxids <- taxonomy_indices |>
        dplyr::rename(rname = "taxa_index",
                      rname_names = accessions)
    } else {
      accessions_taxids <- tidyr::tibble(rname_names = accessions, taxids, rname = match(taxids, unique(taxids)))
    }
    combined_distinct <- results[[2]]
    combined_distinct <- combined_distinct |>
      dplyr::right_join(accessions_taxids, combined_distinct, by = "rname", relationship = "many-to-many") |>
      dplyr::mutate(qname_names = read_names[.data$qname])

    if (db == "ncbi") {
      bam_index_df <- data.frame("original_bam_index" = seq_along(reads[[1]]$qname),
                                 "qname_names" = reads[[1]]$qname,
                                 "rname_names" = as.character(reads[[1]]$rname))
    }
    if (db == "silva"){
      bam_index_df <- data.frame("original_bam_index" = seq_along(reads[[1]]$qname),
                                 "qname_names" = reads[[1]]$qname,
                                 "rname_names" = sub(';.*$','', as.character(reads[[1]]$rname)))
    }

    combined_bam_index <- dplyr::right_join(bam_index_df, combined_distinct,
                                            by = (c("qname_names", "rname_names"))) |>
      dplyr::mutate("qname_rname" = paste(.data$qname, .data$rname, sep = "_"),
                    "first_qname_rname" = !duplicated(.data$qname_rname)) |>
      dplyr::filter(.data$first_qname_rname == TRUE)

    filter_which <- rep(FALSE, nrow(bam_index_df))
    filter_which[combined_bam_index$original_bam_index] <- TRUE
    # Filter based on batches and row 
    i <- 0L
    keep_rule <- S4Vectors::FilterRules(list(
      keep = function(rec) {
        n <- NROW(rec$qname)
        idx <- seq.int(i + 1L, length.out = n)
        i <<- i + n
        filter_which[idx]
      }
    ))
    param <- Rsamtools::ScanBamParam(what = "qname")  # just to know chunk size
    
    Rsamtools::indexBam(files = input_file)
    input_bam <- Rsamtools::BamFile(input_file, yieldSize = 1e6)
    
    bam_out <- file.path(tmp_dir, paste0(out_base, ".updated.bam"))
    
    Rsamtools::filterBam(
      input_bam,
      destination = bam_out,
      param = param,
      filter = keep_rule,
      indexDestination = TRUE
    )
    
  }
  # Plotting genome locations
  num_plot <- num_species_plot
  if (is.null(num_plot)) num_plot <- min(nrow(metascope_id_file), 10)
  if (num_plot > 0) {
    if (!quiet) message("Creating coverage plots at ",
                        out_base, "_cov_plots")
    lapply(seq_along(metascope_id_file$TaxonomyID)[seq_len(num_plot)],
           function(x) {
             locations(as.numeric(metascope_id_file$TaxonomyID)[x],
                       which_genome = metascope_id_file$Genome[x],
                       accessions, taxids, reads, out_base, out_dir)})
  } else if (!quiet) message("No coverage plots created")
  return(out_file)
}
