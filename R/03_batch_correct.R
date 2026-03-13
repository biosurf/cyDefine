#' Perform batch correction via cyCombine
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @param num.threads Number of threads to parallelize with
#' @inheritParams cyCombine::batch_correct
#' @param verbose Verbosity
#' @param ... Additional optional parameters for batch correction.
#'  See `?cyCombine::batch_correct` for options.
#'
#' @return Named list of tibbles of batch corrected reference and query data
#' @export
#'
batch_correct <- function(
    reference,
    query,
    markers,
    covar = NULL,
    label = NULL,
    norm_method = "scale",
    ref.batch = NULL,
    seed = 332,
    xdim = 6,
    ydim = 6,
    rlen = 10,
    mc.cores = NULL,
    num.threads = 1,
    verbose = TRUE,
    ...) {
  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c(markers))

  # Use mc.cores if given (for compatibility with cyCombine)
  if (!is.null(mc.cores)) num.threads <- mc.cores
  # add batch and sample columns if not provided
  if ("sample" %!in% colnames(reference)) {
    if (verbose) {
      message("Sample information is not provided for reference - assuming one sample.")
    }
    reference$sample <- "ref_sample"
  }

  if ("sample" %!in% colnames(query)) {
    if (verbose) {
      message("Sample information is not provided for query - assuming one sample.")
    }
    query$sample <- "dat_sample"
  }

  if ("batch" %!in% colnames(reference)) {
    if (verbose) {
      message("Batch information is not provided for reference - assuming one batch.")
    }
    reference$batch <- "ref_batch"
  }

  if ("batch" %!in% colnames(query)) {
    if (verbose) {
      message("Batch information is not provided for query - assuming one batch.")
    }
    query$batch <- "dat_batch"
  }

  # check uniqueness of samples
  common_samples <- intersect(
    unique(reference$sample),
    unique(query$sample))
  if (length(common_samples) > 0) {
    warning(
      "Overlapping sample ID(s) found between reference and query. ",
      "Assuming that these represent different samples. Adding '_ref' ",
      "and '_query', respectively, to the end of overlapping sample ID(s)"
    )

    for (s in common_samples) {
      reference <- reference |>
        dplyr::mutate(sample = replace(
          sample,
          sample == s,
          paste0(s, "_ref")
        ))
      query <- query |>
        dplyr::mutate(sample = replace(
          sample,
          sample == s,
          paste0(s, "_query")
        ))
    }
  }

  if (length(intersect(
    unique(reference$batch),
    unique(query$batch)
  )) > 0 &
    verbose) {
    message(
      "Overlapping batch IDs found between reference and query. ",
      "If originating from different batches, make the IDs unique and ",
      "rerun batch correction."
    )
  }


  # combine query and reference
  uncorrected <- dplyr::full_join(reference,
    query,
    by = intersect(
      colnames(reference),
      colnames(query)
    )
  )
  uncorrected <- check_id(uncorrected)

  # run batch correction
  corrected <- cyCombine::batch_correct(
    uncorrected,
    markers = markers,
    label = label,
    xdim = xdim,
    ydim = ydim,
    ref.batch = ref.batch,
    norm_method = norm_method,
    covar = covar,
    rlen = rlen,
    seed = seed,
    mc.cores = num.threads,
    ...
  )

  rm(uncorrected)

  corrected_ref <- corrected |>
    dplyr::filter(sample %in% unique(reference$sample)) |>
    dplyr::select(colnames(reference))
  corrected_query <- corrected |>
    dplyr::filter(sample %in% unique(query$sample)) |>
    dplyr::select(colnames(query))

  return(list(
    "reference" = corrected_ref,
    "query" = corrected_query
  ))
}
