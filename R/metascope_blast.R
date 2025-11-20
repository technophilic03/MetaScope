#' Gets Taxonomy IDS from accessions names

#' Returns a dataframe with the original accession, the associated taxonomy ID
#' and associated taxonomy species

#' @param all_fastas Biostrings loaded object of all fastas
#' @param accession_path Path to taxonomizr accessions database
#' @keywords internal
#' @return Data.frame of taxonomy Ids and Species name
#'

find_taxids <- function(all_fastas, accession_path) {
  if (grepl("ti|",names(all_fastas)[1], fixed = TRUE)){
    accessions_taxids <- tibble::tibble(
      "accessions" = unique(names(all_fastas))) |>
      dplyr::mutate(
        "TaxonomyID" = as.numeric(stringr::str_match(!!dplyr::sym("accessions"),
                                                     "ti\\|(\\d+)\\|")[,2]),
        "species" = taxonomizr::getTaxonomy(!!dplyr::sym("TaxonomyID"),
                                            sqlFile = accession_path,
                                            desiredTaxa = "species")[,]
      )
  } else {
    accessions_taxids <- tibble::tibble("accessions" = unique(names(all_fastas))) |>
      dplyr::mutate("TaxonomyID" = taxonomizr::accessionToTaxa(!!dplyr::sym("accessions"), accession_path),
                    "species" = taxonomizr::getTaxonomy(!!dplyr::sym("TaxonomyID"),
                                                        sqlFile = accession_path,
                                                        desiredTaxa = "species")[,])
  }
  return(accessions_taxids)
}

#' Adds in taxa if silva database

#' Returns MetaScope Table with silva taxa in separate columns

#' @param metascope_id_in MetaScope ID file with silva taxa
#' @param caching Boolean for if all_silva_headers.rds is already downloaded
#' @param path_to_write Path to save all_silva_headers.rds
#' @keywords internal
#' @return Data.frame of taxonomy information

add_in_taxa <- function(metascope_id_in, caching, path_to_write) {
  location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
  filename <- "all_silva_headers.rds"
  if (!caching) {
    if (!dir.exists(path_to_write)) dir.create(path_to_write)
    destination <- paste(path_to_write, filename, sep = "/")
    utils::download.file(location, destination)
  } else if (caching) {
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, filename, "rname")$rid
    if (!length(rid)) {
      rid <- names(BiocFileCache::bfcadd(bfc, filename, location))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
      BiocFileCache::bfcdownload(bfc, rid)
      BiocFileCache::bfcrpath(bfc, rids = rid)
    } else {
      message("Caching is set to TRUE, ",
              "and it appears that this file is already downloaded ",
              "in the cache. It will not be downloaded again.")
    }
    destination <- BiocFileCache::bfcrpath(bfc, rids = rid)
  }
  all_silva_headers <- readRDS(destination)
  tax_table_pre <- metascope_id_in |>
    dplyr::mutate("TaxonomyID" = as.character(!!dplyr::sym("TaxonomyID"))) |>
    dplyr::distinct(!!dplyr::sym("TaxonomyID"), .keep_all = TRUE) |>
    dplyr::left_join(all_silva_headers, by = c("TaxonomyID")) |>
    dplyr::relocate("genus", "species", .after = "family") |>
    dplyr::relocate("read_count", "Proportion", "readsEM", "EMProportion", .after = "species") |>
    dplyr::select(-"Genome")
  return(tax_table_pre)
}

#' Adds in taxa if input used NCBI database
#'
#' Returns MetaScope Table with NCBI taxa in separate columns
#'
#' @param metascope_id_in MetaScope ID file with NCBI taxa qnames
#' @param BPPARAM An optional BiocParallelParam instance determining the
#'   parallel back-end to be used during evaluation.
#' @inheritParams metascope_blast
#' @keywords internal
#' @return data.frame or tibble of taxonomy information

add_in_taxa_ncbi <- function(metascope_id_in, accession, BPPARAM) {
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class",
                   "order", "family", "genus", "species")
  all_ncbi <- plyr::llply(metascope_id_in$TaxonomyID, .fun = class_taxon,
                          accession = accession, .progress = "text")
  # fix unknowns
  unk_tab <- data.frame(name = "unknown", rank = taxon_ranks, id = 0)
  all_ncbi <- lapply(all_ncbi, function(x) if (all(is.na(x))) unk_tab else x)
  # Create table
  taxonomy_table <- BiocParallel::bplapply(all_ncbi, mk_table,
                                           taxon_ranks = taxon_ranks,
                                           BPPARAM = BPPARAM) |>
    dplyr::bind_rows() |> as.data.frame() |>
    magrittr::set_colnames(taxon_ranks) |>
    dplyr::mutate("TaxonomyID" = metascope_id_in$TaxonomyID) |>
    dplyr::relocate("TaxonomyID") |>
    dplyr::distinct(!!dplyr::sym("TaxonomyID"), .keep_all = TRUE) |>
    dplyr::left_join(metascope_id_in, by = c("TaxonomyID")) |>
    dplyr::relocate("genus", "species", .after = "family") |>
    dplyr::relocate("read_count", "Proportion", "readsEM", "EMProportion", .after = "species") |>
    dplyr::select(-"Genome")
  return(taxonomy_table)
}

#' Gets sequences from bam file
#'
#' Returns fasta sequences from a bam file with a given taxonomy ID
#'
#' @param id Taxonomy ID of genome to get sequences from
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param n Number of sequences to retrieve
#' @param bam_seqs A list of the sequence IDs from the bam file
#'
#' @return Biostrings format sequences

get_seqs <- function(id, bam_file, n = 10, bam_seqs) {
  rlang::is_installed("GenomicRanges")
  rlang::is_installed("IRanges")
  allGenomes <- stringr::str_detect(id, bam_seqs)
  # Get sequence info (Genome Name) from bam file
  seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
    tibble::rownames_to_column("seqnames")
  # Sample one of the Genomes that match
  Genome <- sample(seq_info_df$seqnames[allGenomes], 1)
  # Scan Bam file for all sequences that match genome
  param <- Rsamtools::ScanBamParam(what = c("rname", "seq"),
                                   which = GenomicRanges::GRanges(
                                     seq_info_df$seqnames[allGenomes],
                                     IRanges::IRanges(1, 1e+07)))
  allseqs <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  n <- min(n, length(allseqs$seq))
  seqs <- sample(allseqs$seq, n)
  return(seqs)
}


#' Gets multiple sequences from different accessions in a bam file
#'
#' Returns fasta sequences from a bam file with given taxonomy IDs
#'
#' @param ids_n List of vectors with Taxonomy IDs and the number of sequences to get from each
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param seq_info_df Dataframe of sequence information from metascope_blast()
#' @param metascope_id_tax Data.frame of taxonomy information
#' @param sorted_bam_file Filepath to sorted bam file
#'
#' @return Biostrings format sequences
#' @keywords internal
get_multi_seqs <- function(ids_n, bam_file, seq_info_df,
                           metascope_id_tax,
                           sorted_bam_file) {
  ids <- ids_n[1]
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  allGenomes <- seq_info_df |>
    dplyr::filter(!!dplyr::sym("seqnames") %in% ids) |>
    dplyr::pull("original_seqids")
  # Define BAM file with smaller yield size
  YS <- 1000
  bf <- Rsamtools::BamFile(sorted_bam_file, yieldSize = YS)
  open(bf)
  param <- Rsamtools::ScanBamParam(what = c("rname", "seq"),
                                   which = GenomicRanges::GRanges(
                                     allGenomes,
                                     IRanges::IRanges(1, 1e+07)))
  max_n <- Rsamtools::countBam(bf, param = param)$records
  close(bf)
  n_seqs <- min(as.integer(ids_n[2]), sum(max_n))

  YIELD <- function(bf) {
    value <- Rsamtools::scanBam(bf, param = param)
    flattened_value <- list()
    flattened_value$rname <- lapply(value, function(x) x$rname)
    flattened_value$seq <- lapply(value, function(x) x$seq |> as.character())
    return(flattened_value)
  }
  MAP <- function(value) {
    if(length(value$rname) == 0) return(dplyr::tibble(NULL))
    y <- dplyr::tibble("UID" = value$rname |> unlist(),
                       "Sequence" = value$seq |> unlist())
    return(y)
  }
  REDUCE <- function(x, y) {
    all_join <- dplyr::bind_rows(x, y)
    # Print all current reads to file
    return(all_join)
  }

  reduceByYield_multiseqs <- function(X, YIELD, MAP, REDUCE, n_end) {
    if (!Rsamtools::isOpen(X)) {
      Rsamtools::open.BamFile(X)
      on.exit(Rsamtools::close.BamFile(X))
    }
    result <- dplyr::tibble(NULL)

    while(nrow(result) < n_end) {
      data <- YIELD(X)
      result <- REDUCE(result, MAP(data))
      message(result)
    }
    return(result)
  }

  seqs <- reduceByYield_multiseqs(bf, YIELD, MAP, REDUCE, n_end = n_seqs)
  #seqs <- seqs[sample(nrow(seqs), n_seqs), ]

  final_seqs <- list(rname = seqs$UID, seq = Biostrings::DNAStringSet(seqs$Sequence))
  return(final_seqs)
}

#' Converts NCBI taxonomy ID to scientific name
#'
#' @param taxids List of NCBI taxids to convert to scientific name
#' @importFrom rlang .data
#' @inheritParams metascope_blast
#' @return Returns a dataframe of blast results for a metascope result
#'

taxid_to_name <- function(taxids, accession_path) {
  taxids <- stringr::str_replace(taxids, ";(.*)$", "") |> as.integer()
  out <- taxonomizr::getTaxonomy(taxids, accession_path)
  out_df <- as.data.frame(out) |> tibble::rownames_to_column("staxid") |>
    dplyr::select("staxid", "genus", "species") |>
    dplyr::mutate("species" = stringr::str_replace(!!dplyr::sym("species"), paste0(!!dplyr::sym("genus"), " "), "")) |>
    dplyr::mutate("staxid" = gsub("\\.", "", !!dplyr::sym("staxid"))) |>
    dplyr::mutate("staxid" = gsub("X", "", !!dplyr::sym("staxid"))) |>
    dplyr::mutate("staxid" = as.integer(!!dplyr::sym("staxid")))
  return(out_df)
}

#' Checks if blastn is installed
#'
#' @title Check if blastn exists on the system
#' @description This is an internal function that is not meant to be used
#'   outside of the package. It checks whether blastn exists on the system.
#' @return Returns TRUE if blastn exists on the system, else FALSE.
#'

check_blastn_exists <- function() {
  if (file.exists(Sys.which("blastn"))) return(TRUE)
  return(FALSE)
}

# BLASTn sequences


blastn_seqs <- function(db_path, fasta_path, res_path, hit_list, num_threads) {
  if (!check_blastn_exists()) {
    stop("BLAST executable not found. Please install before running.")
  }
  sys::exec_wait(
    "blastn", c("-db", db_path, "-query", fasta_path, "-out", res_path,
                "-outfmt", paste("10 qseqid sseqid pident length mismatch gapopen",
                                 "qstart qend sstart send evalue bitscore staxid"),
                "-max_target_seqs", hit_list, "-num_threads", num_threads,
                "-task", "megablast"),
    std_out = FALSE, 
    std_err = FALSE)
}



#' blastn_single_result
#'
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param which_result Index in results_table for which result to Blast search
#' @param num_reads Number of reads to blast per result
#' @param hit_list Number of how many blast results to fetch per read
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#' @param bam_seqs A list of the sequence IDs from the bam file
#' @param out_path Path to output results.
#' @param fasta_dir Path to where fasta files are stored.
#' @inheritParams metascope_blast
#'
#' @return Returns a dataframe of blast results for a metascope result
#'

blastn_single_result <- function(results_table, bam_file, which_result,
                                 num_reads = 100, hit_list = 10,
                                 num_threads, db_path, quiet,
                                 accession_path, bam_seqs, out_path,
                                 sample_name, fasta_dir = NULL) {
  res <- tryCatch({ #If any errors, should just skip the organism
    genome_name <- results_table$Genome[which_result]
    if (!quiet) message("Current id: ", genome_name)
    tax_id <- results_table$TaxonomyID[which_result]
    if (!quiet) message("Current ti: ", tax_id)
    # Generate sequences to blast
    if (!is.null(fasta_dir)) {
      fasta_path <- list.files(path = fasta_dir, full.names = TRUE)[which_result]
    } else {
      fasta_path <- get_seqs(id = tax_id, bam_file = bam_file, n = num_reads,
                             bam_seqs = bam_seqs)
    }

    res_path <- file.path(out_path, paste0(sprintf("%05d", which_result), "_",
                                          sample_name, ".csv"))
    blastn_seqs(db_path, fasta_path, res_path = res_path, hit_list, num_threads)
    blast_res <- utils::read.csv(res_path, header = FALSE) |>
      magrittr::set_colnames(c("qseqid", "sseqid", "pident", "length", "mismatch",
                               "gapopen","qstart", "qend", "sstart", "send",
                               "evalue", "bitscore", "staxid")) |>
      dplyr::mutate("MetaScope_Taxid" = tax_id,
                    "MetaScope_Genome" = as.character(genome_name))
      taxize_genome_df <- taxid_to_name(unique(blast_res$staxid),
                                        accession_path = accession_path)
    blast_res <- dplyr::left_join(blast_res, taxize_genome_df, by = "staxid")
    blast_res
  },
  error = function(e) {
    message(conditionMessage(e), "\n")
    this_format <- paste("qseqid sseqid pident length mismatch gapopen",
                         "qstart qend sstart send evalue bitscore staxid")
    all_colnames <- stringr::str_split(this_format, " ")[[1]]
    blast_res <- matrix(NA, ncol = length(all_colnames),
                        dimnames = list(c(), all_colnames)) |> as.data.frame() |>
      dplyr::mutate("MetaScope_Taxid" = results_table[which_result, 1],
                    "MetaScope_Genome" =
                      as.character(results_table[which_result, 2]),
                    "name" = NA)
    blast_res
  }
  )
  return(res)
}

#' Reformat BLASTn results
#'
#' @param results_table data.frame containing the MetaScope results.
#' @param bam_file \code{Rsamtools::bamFile} instance for the given sample.
#' @param num_results Integer; maximum number of Metascope results to BLAST.
#'   Default is 10.
#' @param num_reads_per_result Integer; number of reads to BLAST per result.
#'   Default is 100.
#' @param hit_list Integer; how many BLAST results to fetch for each read.
#'   Default is 10.
#' @param num_threads Integer; how many threads to use if multithreading.
#'   Default is 1.
#' @param db_path Character string; filepath for the location of the
#'   pre-installed BLAST database.
#' @param out_path Character string; Output directory to save CSV output files,
#'   including base name of files. For example, given a sample "X78256",
#'   filepath would be \code{file.path(directory_here, "X78256")} with extension
#'   omitted.
#' @param fasta_dir Character string; Directory where fastas from
#'   \code{metascope_id} are stored.
#' @param BPPARAM An optional BiocParallelParam instance determining the
#'   parallel back-end to be used during evaluation.
#' @inheritParams metascope_blast
#'
#' @return Creates and exports num_results number of csv files with blast
#'   results from local blast

blastn_results <- function(results_table, bam_file, num_results = 10,
                           num_reads_per_result = 100, hit_list = 10,
                           num_threads = 1, db_path, out_path, db = NULL,
                           sample_name = NULL, quiet = quiet,
                           accession_path, fasta_dir = NULL,
                           BPPARAM) {
  # Grab all identifiers
  if (is.null(fasta_dir)) {
    seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
      tibble::rownames_to_column("seqnames")
    bam_seqs <- find_accessions(seq_info_df$seqnames,
                                quiet = quiet) |>
      plyr::aaply(1, function(x) x[1])
  } else {
    bam_seqs <- NULL
  }
  # Grab results
  num_results2 <- min(num_results, nrow(results_table))
  run_res <- function(i) {
    df <- blastn_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result,
                               hit_list = hit_list, num_threads = num_threads,
                               db_path = db_path, quiet = quiet, bam_seqs = bam_seqs,
                               out_path = out_path, sample_name = sample_name,
                               fasta_dir = fasta_dir, accession_path = accession_path)
    utils::write.csv(df, file.path(out_path,
                                   paste0(sprintf("%05d", i), "_", sample_name, ".csv")),
                     row.names = FALSE)
  }
  BiocParallel::bplapply(seq_len(num_results2), run_res, BPPARAM = BPPARAM)
}

#' Calculates result metrics from a blast results table
#'
#' This function calculates the best hit (genome with most blast read hits),
#' uniqueness score (total number of genomes hit), species percentage hit
#' (percentage of reads where MetaScope species also matched the blast hit
#' species), genus percentage hit (percentage of reads where blast genus matched
#' MetaScope aligned genus) and species contaminant score (percentage of reads
#' that blasted to other species genomes) and genus contaminant score
#' (percentage of reads that blasted to other genus genomes)
#'
#' @param blast_results_table_path path for blast results csv file
#' @inheritParams metascope_blast
#'
#' @return a vector with best_hit, uniqueness_score, species_percentage_hit
#'   genus_percentage_hit, species_contaminant_score, and
#'   genus_contaminant_score

blast_result_metrics <- function(blast_results_table_path, accession_path, db = NULL){
  tryCatch({
    # Load in blast results table
    blast_results_table <- utils::read.csv(blast_results_table_path, header = TRUE)

    # Clean species results
    blast_results_table <- blast_results_table |>
      dplyr::filter(!grepl("sp.", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!grepl("uncultured", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!is.na(!!dplyr::sym("genus")))

    # Remove any empty tables
    if (nrow(blast_results_table) < 2) {
      error_vector <- c(rep(0, 5), rep(NA, 3)) |>
        magrittr::set_names(c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                            "species_contaminant_score", "genus_contaminant_score",
                            "best_hit_genus", "best_hit_species", "best_hit_strain"))
    return(error_vector)
    }
    # Clean species results
    blast_results_table <- blast_results_table |>
      dplyr::filter(!grepl("sp.", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!grepl("uncultured", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!is.na(!!dplyr::sym("genus")))
    # Adding MetaScope Species and Genus columns
    if (db == "silva") {
      # Cleaning up silva genome names
      #silva_genome <- gsub(";([a-z0-9])", " \\1", blast_results_table$MetaScope_Genome[1])
      silva_genome <- gsub(" uncultured", ";uncultured", blast_results_table$MetaScope_Genome[1])
      # Creating silva genus and species names
      silva_genus <- magrittr::extract(strsplit(silva_genome, split = "_"), 1)[[1]]
      silva_species <- magrittr::extract(strsplit(silva_genome, split = "_"), 2)[[1]]
      blast_results_table_2 <- blast_results_table |>
        dplyr::mutate("MetaScope_genus" = rep(silva_genus, nrow(blast_results_table)),
                      "MetaScope_species" = rep(silva_species, nrow(blast_results_table))) |>
        dplyr::rename("query_genus" = "genus",
                      "query_species" = "species") |>
        # Remove rows with NA
        tidyr::drop_na("query_genus", "query_species") |>
        # Getting best hit per read
        dplyr::group_by(!!dplyr::sym("qseqid")) |>
        dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
        # Removing duplicate query num and query species
        dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("query_species"), .keep_all = TRUE)

    } else {
      meta_tax <- taxid_to_name(unique(blast_results_table$MetaScope_Taxid),
                                accession_path = accession_path) |>
        dplyr::select(-"staxid")
      blast_results_table_2 <- blast_results_table |>
        dplyr::mutate("MetaScope_genus" = meta_tax$genus[1],
                      "MetaScope_species" = meta_tax$species[1]) |>
        dplyr::rename("query_genus" = "genus",
                      "query_species" = "species") |>
        # Remove rows with NA
        tidyr::drop_na() |>
        # Getting best hit per read
        dplyr::group_by(!!dplyr::sym("qseqid")) |>
        dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
        # Removing duplicate query num and query species
        dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("query_species"), .keep_all = TRUE)
    }
    # Calculating Metrics
    best_hit_genus <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_genus")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = FALSE)

    best_hit_species <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_species")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = FALSE)

    #best_hit_strain <- blast_results_table_2 |>
    #  tidyr::separate("sseqid", c(NA, "gi", NA, NA, NA), sep = "\\|")  |>
    #  dplyr::group_by(!!dplyr::sym("gi")) |>
    #  dplyr::summarise("num_reads" = dplyr::n()) |>
    #  dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = TRUE)
    #if (nrow(best_hit_strain) > 1 | db == "silva") {
    #  best_strain <- NA
    #} else if (nrow(best_hit_strain) == 1) {
    #  res <- taxonomizr::getTaxonomy(best_hit_strain$gi, sqlFile = accession_path)
    #  best_strain <- attr(res[[1]], "name")
    #}

    uniqueness_score <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_species")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      nrow()

    species_percentage_hit <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_species") == !!dplyr::sym("query_species")) |>
      nrow() / length(unique(blast_results_table_2$qseqid))

    genus_percentage_hit <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_genus") == !!dplyr::sym("query_genus")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = FALSE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))

    species_contaminant_score <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_species") != !!dplyr::sym("query_species")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = TRUE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))

    genus_contaminant_score <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_genus") != !!dplyr::sym("query_genus")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = TRUE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))


    all_results <- c(uniqueness_score, species_percentage_hit, genus_percentage_hit,
                     species_contaminant_score, genus_contaminant_score,
                     best_hit_genus$query_genus, best_hit_species$query_species)
    names(all_results) <- c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                            "species_contaminant_score", "genus_contaminant_score",
                            "best_hit_genus", "best_hit_species")
    return(all_results)
  },
  error = function(e)
  {
    message("Error", conditionMessage(e), "/n")
    error_vector <- c(0,0,0,0,0, NA, NA)
    names(error_vector) <- c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                             "species_contaminant_score", "genus_contaminant_score",
                             "best_hit_genus", "best_hit_species")

    return(error_vector)
  }
  )
}


#' Blast reads from MetaScope aligned files
#'
#' This function allows the user to check a subset of identified reads against
#' NCBI BLAST and the nucleotide database to confirm or contradict results
#' provided by MetaScope. It aligns the top `metascope_id()` results to NCBI
#' BLAST database. It REQUIRES that command-line BLAST and a separate nucleotide
#' database have already been installed on the host machine. It returns a csv
#' file updated with BLAST result metrics.
#'
#' This function assumes that you used the NCBI nucleotide database to process
#' samples, with a download date of 2021 or later. This is necessary for
#' compatibility with the bam file headers.
#'
#' This is highly computationally intensive and should be ran with multiple
#' cores, submitted as a multi-threaded computing job if possible.
#'
#' Note, if best_hit_strain is FALSE, then no strain was observed more often
#' among the BLAST results.
#'
#' @param metascope_id_path Character string; path to a csv file output by
#'   `metascope_id()`.
#' @param bam_file_path Character string; full path to bam file for the same
#'   sample processed by `metascope_id`. Note that the `filter_bam` function
#'   must have output a bam file, and not a .csv.gz file. See
#'   `?filter_bam_bowtie` for more details. Defaults to
#'   \code{list.files(file_temp, ".updated.bam$")[1]}.
#' @param tmp_dir Character string, a temporary directory in which to host
#'   files.
#' @param out_dir Character string, path to output directory.
#' @param sample_name Character string, sample name for output files.
#' @param num_results Integer, the maximum number of taxa from the metascope_id
#'   output to check reads. Takes the top n taxa, i.e. those with largest
#'   abundance. Default is 10.
#' @param num_reads Integer, the maximum number of reads to blast per microbe.
#'   If the true number of reads assigned to a given taxon is smaller, then the
#'   smaller number will be chosen. Default is 100. Too many reads will involve
#'   more processing time.
#' @param hit_list Integer, number of blast hit results to keep. Default is 10.
#' @param num_threads Integer, number of threads if running in parallel
#'   (recommended). Default is 1.
#' @param db_path Character string. The database file to be searched (including
#'   basename, but without file extension). For example, if the nt database is
#'   in the nt folder, this would be /filepath/nt/nt assuming that the database
#'   files have the nt basename. Check this path if you get an error message
#'   stating "No alias or index file found".
#' @param quiet Logical, whether to print out more informative messages. Default
#'   is FALSE.
#' @param db Currently accepts one of \code{c("ncbi", "silva", "other")} Default
#'   is \code{"ncbi"}, appropriate for samples aligned against indices compiled
#'   from NCBI whole genome databases. Alternatively, usage of an alternate
#'   database (like Greengenes2) should be specified with \code{"other"}.
#' @param fasta_dir Directory where fasta files for blast will be stored.
#' @param accession_path (character) Filepath to NCBI accessions SQL
#'   database. See \code{taxonomzr::prepareDatabase()}.
#'
#' @returns This function writes an updated csv file with metrics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ### Create temporary directory
#' file_temp <- tempfile()
#' dir.create(file_temp)
#'
#' bamPath <- system.file("extdata", "bowtie_target.bam",
#'                        package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' metascope_id(file.path(file_temp, "bowtie_target.bam"), aligner = "bowtie2",
#'              input_type = "bam", out_dir = file_temp, num_species_plot = 0,
#'              update_bam = TRUE)
#'
#' ## Run metascope blast
#' ### Get export name and metascope id results
#' out_base <- bamPath |> basename() |> tools::file_path_sans_ext() |>
#'   tools::file_path_sans_ext()
#' metascope_id_path <- file.path(file_temp, paste0(out_base, ".metascope_id.csv"))
#'
#' # NOTE: change db_path to the location where your BLAST database is stored!
#' db <- "/restricted/projectnb/pathoscope/data/blastdb/nt/nt"
#'
#' tmp_accession <- system.file("extdata", "example_accessions.sql", package = "MetaScope"
#'
#' metascope_blast(metascope_id_path,
#'                 bam_file_path = file.path(file_temp, "bowtie_target.bam"),
#'                 tmp_dir = file_temp,
#'                 out_dir = file_temp,
#'                 sample_name = out_base,
#'                 db_path = db,
#'                 num_results = 10,
#'                 num_reads = 5,
#'                 hit_list = 10,
#'                 num_threads = 3,
#'                 db = "ncbi",
#'                 quiet = FALSE,
#'                 fasta_dir = NULL,
#'                 accession_path = tmp_accession)
#'
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'}
#'

metascope_blast <- function(metascope_id_path,
                            bam_file_path = list.files(tmp_dir, ".updated.bam$",
                                                       full.names = TRUE)[1],
                            tmp_dir, out_dir, sample_name, fasta_dir = NULL,
                            num_results = 10, num_reads = 100, hit_list = 10,
                            num_threads = 1, db_path, quiet = FALSE,
                            db = NULL, accession_path = NULL) {
  # Parallelization settings
  if (!is.numeric(num_threads)) num_threads <- 1
  BPPARAM <- BiocParallel::MulticoreParam(workers = num_threads)

  # Sort and index bam file
  sorted_bam_file_path <- file.path(tmp_dir, paste0(sample_name, "_sorted"))
  message("Sorting bam file to ", sorted_bam_file_path, ".bam")
  Rsamtools::sortBam(bam_file_path, destination = sorted_bam_file_path,
                     overwrite = TRUE)
  sorted_bam_file <- paste0(sorted_bam_file_path, ".bam")
  Rsamtools::indexBam(sorted_bam_file)
  bam_file <- Rsamtools::BamFile(sorted_bam_file, index = sorted_bam_file)

  # Generate fasta file from sorted bam file
  command_string <- paste0("samtools view ", file.path(tmp_dir, sample_name),
                           "_sorted.bam | awk \'{print \">\"$3 \"\\n\" $10}\' > ",
                           file.path(tmp_dir, sample_name), "_aligned.fasta")
  system(command_string)
  all_fastas <- Biostrings::readDNAStringSet(paste0(file.path(tmp_dir, sample_name), "_aligned.fasta"))
  accessions_taxids <- find_taxids(all_fastas, accession_path)


  # Load in metascope id file and clean unknown genomes
  metascope_id_in <- utils::read.csv(metascope_id_path, header = TRUE)

  # Group metascope id by species and create metascope species id
  if (db == "silva") {
    metascope_id_species <- add_in_taxa(metascope_id_in, caching = FALSE,           # TODO Add initialization steps to install required databases
                                    path_to_write = tmp_dir)
  } else if (db == "ncbi") {
    tax_df <- taxonomizr::getTaxonomy(metascope_id_in$TaxonomyID, sqlFile = accession_path)
    rownames(tax_df) <- trimws(rownames(tax_df))
    tax_df <- as.data.frame(tax_df) |> 
      tibble::rownames_to_column("TaxonomyID") |>
      dplyr::mutate("TaxonomyID" = ifelse(grepl("NA", !!dplyr::sym("TaxonomyID"), fixed=TRUE), 
                                        NA, 
                                        gsub("^X[.]?", "", !!dplyr::sym("TaxonomyID")))) |>
      dplyr::filter(!is.na(!!dplyr::sym("TaxonomyID"))) |>
      dplyr::mutate("TaxonomyID" = as.integer(!!dplyr::sym("TaxonomyID")))
      
    metascope_id_species <- dplyr::left_join(metascope_id_in, tax_df,
                                             by = "TaxonomyID")
  }

  # Create fasta directory in tmp directory to save fasta sequences
  fastas_tmp_dir <- file.path(tmp_dir, "fastas")
  if(!dir.exists(fastas_tmp_dir)) dir.create(fastas_tmp_dir,
                                             recursive = TRUE)
  unlink(paste0(fastas_tmp_dir, "/*"), recursive = TRUE)
  message("Generating fasta sequences from bam file")
  # How many taxa
  num_taxa_loop <- min(nrow(metascope_id_species), num_results)

  # Sampling num_reads number of fastas per genome
  sample_fastas <- function(i) {
    current_species <- metascope_id_species$species[i]
    if (!is.na(current_species)) {
      current_accessions <- accessions_taxids |>
        dplyr::filter(!!dplyr::sym("species") == current_species) |>
        dplyr::pull("accessions")
    } else {
      accession_id <- gsub("unknown genome; accession ID is", "", 
                           metascope_id_species$Genome[i])
      current_accessions <- accessions_taxids$accessions[
        which(grepl(accession_id, accessions_taxids$accessions))[1]
      ]
    }
    seqs <- all_fastas[names(all_fastas) %in% current_accessions]
    num_seqs <- min(length(seqs), num_reads) 
    num_digits <- nchar(as.character(num_reads))
    seqs <- sample(seqs, size = num_seqs)
    if (length(seqs) > 0) {
      names(seqs) <- paste0(
        "ms_index", i, "|",
        names(seqs), 
        "|read_", 
        sprintf(paste0("%0", num_digits, "d"), seq_along(seqs)))
    }
    return(seqs)
  }

  # Parallelize sampling of fastas
  all_seqs <- BiocParallel::bplapply(seq_len(num_taxa_loop), sample_fastas, BPPARAM = BPPARAM)
  fasta_path <- file.path(fastas_tmp_dir, "blast_input_fastas.fa")
  
  Biostrings::writeXStringSet(
    do.call(c,all_seqs), # flattens list  
    filepath = fasta_path)
  
  # Create blast directory in tmp directory to save blast results in
  blast_tmp_dir <- file.path(tmp_dir, "blast")
  if(!dir.exists(blast_tmp_dir)) dir.create(blast_tmp_dir, recursive = TRUE)
  unlink(paste0(blast_tmp_dir, "/*"), recursive = TRUE)

  # Run rBlast on all metascope microbes
  if (!check_blastn_exists()) {
    stop("BLAST executable not found. Please install before running.")
  } else {
    message("BLASTN executable found on system.")
  }
  res_path <- file.path(blast_tmp_dir, "blast_output.csv")
  message("Running BLASTN on all sequences")
  sys::exec_wait(
    "blastn", c("-db", db_path, "-query", fasta_path, "-out", res_path,
                "-outfmt", paste("10 qseqid sseqid pident length mismatch gapopen",
                                 "qstart qend sstart send evalue bitscore staxid"),
                "-max_target_seqs", hit_list, "-max_hsps", 1, 
                "-num_threads", num_threads,
                "-task", "megablast"))

  # Calculate Blast metrics
  message("Running BLAST metrics on all blast results")
  blast_result_metrics_df <- utils::read.csv(res_path, header = FALSE)
  colnames(blast_result_metrics_df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                         "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid")
  ## Extract MetaScope Taxid
  blast_result_metrics_df <- blast_result_metrics_df |>
    dplyr::mutate(
      ms_index = stringr::str_extract(!!dplyr::sym("qseqid"), "(?<=ms_index)\\d+"), 
      metascope_taxid = purrr::map_chr(stringr::str_split(!!dplyr::sym("qseqid"), pattern = "\\|"), ~ .x[3]),
      read_id = purrr::map_chr(stringr::str_split(!!dplyr::sym("qseqid"), pattern = "\\|"), ~ .x[8]))
  ## Calculate MetaScope species and Blast Species
  ms_tax <- taxonomizr::getTaxonomy(unique(blast_result_metrics_df$metascope_taxid), sqlFile = accession_path) |>
    as.data.frame() |>
    tibble::rownames_to_column("metascope_taxid") |>
    dplyr::mutate(metascope_taxid = stringr::str_trim(!!dplyr::sym("metascope_taxid")))
  colnames(ms_tax) <- c("metascope_taxid", "ms_domain", "ms_phylum", "ms_class", "ms_order", "ms_family", "ms_genus", "ms_species")
  b_tax <- taxonomizr::getTaxonomy(unique(blast_result_metrics_df$staxid), sqlFile = accession_path) |>
    as.data.frame() |>
    tibble::rownames_to_column("staxid") |>
    dplyr::mutate(staxid = as.integer(stringr::str_trim(!!dplyr::sym("staxid"))))
  colnames(b_tax) <- c("staxid", "b_domain", "b_phylum", "b_class", "b_order", "b_family", "b_genus", "b_species")
  blast_result_metrics_df <- dplyr::left_join(blast_result_metrics_df, ms_tax, by = "metascope_taxid")
  blast_result_metrics_df <- dplyr::left_join(blast_result_metrics_df, b_tax, by = "staxid")
  # Calculate Blast Metrics grouped by qseqid
  ms_blast_results <- data.frame(
    ms_index = as.character(seq(1, num_results))
  )
  ## Best Hit Genus
  best_hit_genus <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_genus"))) |> 
    dplyr::group_by(!!dplyr::sym("qseqid"), 
                    !!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id"),
                    !!dplyr::sym("b_genus")) |>
    dplyr::summarize(genus_counts = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id")) |>
    dplyr::slice_max(!!dplyr::sym("genus_counts")) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::count(!!dplyr::sym("ms_index"), !!dplyr::sym("b_genus"), name = "genus_counts") |> 
    dplyr::slice_max(order_by = !!dplyr::sym("genus_counts"), n = 1, with_ties = FALSE) |>
    dplyr::rename("best_hit_genus" = "b_genus")
  ms_blast_results <- dplyr::left_join(ms_blast_results, best_hit_genus, by = "ms_index")
  ## Best Hit Species 
  best_hit_species <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_species"))) |> 
    dplyr::group_by(!!dplyr::sym("qseqid"), 
                    !!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id"),
                    !!dplyr::sym("b_species")) |>
    dplyr::summarize(species_counts = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id")) |>
    dplyr::slice_max(!!dplyr::sym("species_counts")) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::count(!!dplyr::sym("ms_index"), !!dplyr::sym("b_species"), name = "species_counts") |> 
    dplyr::slice_max(order_by = !!dplyr::sym("species_counts"), n = 1, with_ties = FALSE) |>
    dplyr::rename("best_hit_species" = "b_species")
  ms_blast_results <- dplyr::left_join(ms_blast_results, best_hit_species, by = "ms_index")
  ## Uniqueness Score
  uniqueness_score <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_species"))) |> 
    dplyr::group_by(!!dplyr::sym("qseqid"), 
                    !!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id"),
                    !!dplyr::sym("b_species")) |>
    dplyr::summarize(species_counts = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index"),
                    !!dplyr::sym("read_id")) |>
    dplyr::slice_max(!!dplyr::sym("species_counts")) |> # Up to here is best hit species
    dplyr::ungroup() |>
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::summarise(uniqueness_score = dplyr::n_distinct(!!dplyr::sym("b_species")))
  ms_blast_results <- dplyr::left_join(ms_blast_results, uniqueness_score, by = "ms_index")
  ## Species Percentage Hit
  species_percentage_hit <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_species"))) |> 
    dplyr::mutate(species_match = !!dplyr::sym("b_species") == !!dplyr::sym("ms_species")) |> 
    dplyr::group_by(!!dplyr::sym("ms_index"), !!dplyr::sym("read_id")) |> 
    dplyr::summarize(has_match = any(!!dplyr::sym("species_match")), .groups = "drop") |> # Check if read has at least one match
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::summarize(species_percentage_hit = mean(!!dplyr::sym("has_match")), .groups = "drop")
  ms_blast_results <- dplyr::left_join(ms_blast_results, species_percentage_hit, by = "ms_index")
  ## Genus Percentage Hit
  genus_percentage_hit <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_genus"))) |> 
    dplyr::mutate(genus_match = !!dplyr::sym("b_genus") == !!dplyr::sym("ms_genus")) |> 
    dplyr::group_by(!!dplyr::sym("ms_index"), !!dplyr::sym("read_id")) |> 
    dplyr::summarize(has_match = any(!!dplyr::sym("genus_match")), .groups = "drop") |> # Check if read has at least one match
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::summarize(genus_percentage_hit = mean(!!dplyr::sym("has_match")), .groups = "drop")
  ms_blast_results <- dplyr::left_join(ms_blast_results, genus_percentage_hit, by = "ms_index")
  ## Species Contaminant Score
  species_contaminant_score <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_species"))) |>
    dplyr::mutate(species_match = !!dplyr::sym("b_species") == !!dplyr::sym("ms_species")) |>
    dplyr::group_by(!!dplyr::sym("ms_index"), !!dplyr::sym("read_id")) |>
    dplyr::summarize(has_mismatch = any(!(!!dplyr::sym("species_match"))), .groups = "drop") |>
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::summarize(species_contaminant_score = mean(!!dplyr::sym("has_mismatch")), .groups = "drop")
  ms_blast_results <- dplyr::left_join(ms_blast_results, species_contaminant_score, by = "ms_index")
  ## Genus Contaminant Score
  genus_contaminant_score <- blast_result_metrics_df |>
    dplyr::filter(!is.na(!!dplyr::sym("b_genus"))) |>
    dplyr::mutate(genus_match = !!dplyr::sym("b_genus") == !!dplyr::sym("ms_genus")) |>
    dplyr::group_by(!!dplyr::sym("ms_index"), !!dplyr::sym("read_id")) |>
    dplyr::summarize(has_mismatch = any(!(!!dplyr::sym("genus_match"))), .groups = "drop") |>
    dplyr::group_by(!!dplyr::sym("ms_index")) |>
    dplyr::summarize(genus_contaminant_score = mean(!!dplyr::sym("has_mismatch")), .groups = "drop")
  ms_blast_results <- dplyr::left_join(ms_blast_results, genus_contaminant_score, by = "ms_index")
  
  # Append Blast Metrics to MetaScope results
  metascope_id_species <- metascope_id_species |>
    dplyr::mutate(ms_index = as.character(dplyr::row_number()))
  metascope_blast_df <- dplyr::left_join(metascope_id_species, ms_blast_results, by = "ms_index") |>
    dplyr::select(-!!dplyr::sym("ms_index"))
  metascope_blast_df[is.na(metascope_blast_df)] <- 0
  
  # Export File.
  print_file <- file.path(out_dir, paste0(sample_name, ".metascope_blast.csv"))
  utils::write.csv(metascope_blast_df, print_file, row.names = FALSE)
  message("Results written to ", print_file)
  
  utils::write.csv(blast_result_metrics_df, file.path(blast_tmp_dir, "blast_output_taxonomy.csv"), row.names = FALSE)
  return(metascope_blast_df)
}

#' Reassign reads from MetaScope BLASTn alignment
#'
#' Using the output from \code{metascope_blast()}, the
#' \code{blast_reassignment()} function takes the results and alters the
#' original \code{metascope_id()} output to reassign reads that were invalidated
#' by the BLAST findings. Currently, the implementation of this function only
#' reassigns reads to a taxon that was already found in the sample at a higher
#' abundance.
#'
#' @param metascope_blast_path Character string. The filepath to a
#'   \code{metascope_blast} CSV output file.
#' @param species_threshold Numeric. A number between 0 and 1 indicating
#'   the minimum proportion of reads needed for a taxon to be considered
#'   validated from the BLAST results. Default is 0.2, or 20\%.
#' @param num_hits Integer. The number of hits for which to assess validation.
#'   Default is 10, i.e., only the top 10 taxa will be assessed.
#' @param blast_tmp_dir Character string. Filepath of the directory where
#'   BLAST results were output from the \code{metascope_blast} function.
#'   Referencing the arguments from \code{metascope_blast}, this would be
#'   \code{file.path(tmp_dir, "blast")}
#' @param reassign_validated Logical. Should reads from validated accessions be
#'   reassigned to other validated accessions. Defaults to \code{FALSE}
#' @param reassign_to_validated Logical. Should reads only be re-assigned to
#'   validated accessions or to any accession with more reads than the current
#'   accession. Defaults to \code{TRUE}
#' @inheritParams metascope_blast
#'
#' @returns Returns a \code{data.frame} with the reassigned taxa and read counts.
#'
#' @export

blast_reassignment <- function(metascope_blast_path, 
                               species_threshold, num_hits,
                               blast_tmp_dir, out_dir, sample_name,
                               reassign_validated = FALSE,
                               reassign_to_validated = TRUE) {
  # Create Validated Column based on species thresholds
  metascope_blast_df <- data.table::fread(metascope_blast_path) |>
    tibble::as_tibble() |>
    dplyr::mutate(blast_validated = (!!dplyr::sym("species_percentage_hit") > species_threshold))
  # Create Index for later reassignment
  metascope_blast_df$ms_index <- rownames(metascope_blast_df) |> as.numeric()
  # Tidy columns read_counts to read_count
  if ("read_counts" %in% colnames(metascope_blast_df)) {
    ind <- colnames(metascope_blast_df) == "read_counts"
    colnames(metascope_blast_df)[ind] <- "read_count"
  }

  reassigned_metascope_blast <- metascope_blast_df

  blast_file <- list.files(blast_tmp_dir, pattern = "blast_output_taxonomy.csv", full.names = TRUE)
  blast_result_metrics_df <- data.table::fread(blast_file)
  # Create vector of indices that have been reassigned
  drop_indices <- c()

  # test all indices if reassign_validated, otherwise only test unvalidated indices
  if (reassign_validated) {
    test_indices <- c(2:min(num_hits, nrow(metascope_blast_df)))
  } else {
    test_indices <- which(!metascope_blast_df$blast_validated)
    test_indices <- test_indices[test_indices <= num_hits]
  }

  # define helper function get_blast_summary
  get_blast_summary <- function(i) {
    blast_summary <-
      tryCatch(
        {
          blast_summary <- data.table::fread(blast_file) |>
            tibble::as_tibble() |>
            dplyr::filter(!!dplyr::sym("ms_index") == i) |>
            dplyr::group_by(!!dplyr::sym("qseqid")) |>
            # Select lowest evalues
            dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
            # Removing duplicate query num and query species
            dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("b_species"), .keep_all = TRUE) |>
            dplyr::ungroup() |>
            dplyr::group_by(!!dplyr::sym("b_species")) |>
            dplyr::summarise("num_reads" = dplyr::n(), .groups="keep") |>
            dplyr::slice_max(order_by = !!dplyr::sym("num_reads"), with_ties = TRUE) |>
            dplyr::left_join(metascope_blast_df[seq_len(num_hits), ],
                             by = dplyr::join_by("b_species" == "species")) |>
            # Do not reassign below metascope's current reassignment
            dplyr::filter(!is.na(!!dplyr::sym("ms_index"))) |>
            dplyr::filter(!!dplyr::sym("ms_index") < i)
          # If reassign_to_validated, then only reassign reads to blast_validated accesions
          if (reassign_to_validated) {
            blast_summary <- blast_summary |>
              dplyr::filter(!!dplyr::sym("blast_validated") == TRUE)
          }
          blast_summary <- blast_summary |>
            dplyr::ungroup() |>
            dplyr::mutate(reassignment_proportion = !!dplyr::sym("num_reads") / sum(!!dplyr::sym("num_reads")),
                          reassigned_read_count = metascope_blast_df$read_count[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_Proportion = metascope_blast_df$Proportion[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_readsEM = metascope_blast_df$readsEM[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_EMProportion = metascope_blast_df$EMProportion[i] * !!dplyr::sym("reassignment_proportion"))
          blast_summary
        },
        error = function(cond) {
          message("Problem with blast summary for current file: ")
          message(blast_file[i])
          message(cond)
          # Choose a return value in case of error
          blast_summary <- dplyr::tibble(
            "genus" = numeric(),
            "species" = numeric(),
            "num_reads" = numeric()
          )
          blast_summary
        }
      )
    return(blast_summary)
  }
  for (i in test_indices){
    if (!metascope_blast_df$blast_validated[i]){
      blast_summary <- get_blast_summary(i)
      if (nrow(blast_summary) > 0) {
        # Updated reassigned_metascope_blast df for every row in blast_summary
        for (n in 1:nrow(blast_summary)) {
          metascope_index <- blast_summary$ms_index[n]
          #print(reassigned_metascope_blast$read_count[metascope_index])
          #print(blast_summary$reassigned_read_count[n])
          reassigned_metascope_blast$read_count[metascope_index] <-
            reassigned_metascope_blast$read_count[metascope_index] + blast_summary$reassigned_read_count[n]
          reassigned_metascope_blast$Proportion[metascope_index] <-
            reassigned_metascope_blast$Proportion[metascope_index] + blast_summary$reassigned_Proportion[n]
          reassigned_metascope_blast$readsEM[metascope_index] <-
            reassigned_metascope_blast$readsEM[metascope_index] + blast_summary$reassigned_readsEM[n]
          reassigned_metascope_blast$EMProportion[metascope_index] <-
            reassigned_metascope_blast$EMProportion[metascope_index] + blast_summary$reassigned_EMProportion[n]
        }
        drop_indices <- append(drop_indices, i)
      }
    }
  }
  
  if (length(drop_indices) > 0) {
    reassigned_metascope_blast <- reassigned_metascope_blast[-drop_indices, ]
  }
  reassigned_metascope_blast <- reassigned_metascope_blast |>
    dplyr::select(-"ms_index")

  print_file <- file.path(out_dir, paste0(sample_name, ".metascope_blast_reassigned.csv"))
  data.table::fwrite(reassigned_metascope_blast, print_file, row.names = FALSE)
  message("Results written to ", print_file)
  return(reassigned_metascope_blast)
}
