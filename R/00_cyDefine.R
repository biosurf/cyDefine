
#' Run the full cyDefine pipeline
#' @inheritParams adapt_reference
#' @inheritParams batch_correct
#' @inheritParams classify_cells
#' @inheritParams identify_unassigned
#' @param reference Tibble of reference data (cells in rows, markers in
#' columns) if you want to apply a custom reference. If NULL, the Seurat
#' PBMC atlas will be adapted and used as reference.
#' @param batch_correct Boolean indicating whether or not you want to perform
#' batch correction via cyCombine
#' @param norm_method Normalization method for cyCombine batch correction.
#' Should be either 'rank', 'scale', or 'qnorm'. Default is 'scale'.
#' @param identify_unassigned Boolean indicating whether or not you want to
#' identify unassigned cells after classifying the cells.
#' @param using_pbmc Boolean indicating whether the Seurat PBMC atlas is
#' being supplied as reference. This info is needed to enable automated
#' screening of common marker names in `map_marker_names()` and tree-based
#' merging of cell type labels.
#' @param exclude_redundant Exclude celltypes in the reference that are not
#' present in the query.
#'
#' @param adapt_reference Boolean indicating whether the provided reference
#' should be adapted to match the query marker panel. Default: FALSE if
#' `using_pbmc` is FALSE, otherwise TRUE.
#'
#' @param save_adapted_reference Default: NULL. Save the adapted reference in a
#'  specified object name or as "adapted_reference.rds" in a specified path.
#'  Usage: FAlSE/NULL - not stored.
#'  TRUE/path - stored as "adapted_reference.rds" at "./" or path.
#'  path/to/filename.rds - stored as "filename.rds" at path/to/.
#' @return Tibble of query data with added columns: "model_prediction"
#' indicating the canonical cell type predicted by the model, and
#' "predicted_celltype" indicating the final predicted cell type after
#' identifying unassigned cells, i.e. will either be "unassigned" or the same
#' as "model_prediction".
#' @export
#'
cyDefine <- function(
    reference,
    query,
    markers,
    using_pbmc = FALSE,
    adapt_reference = ifelse(
      using_pbmc,
      TRUE,
      FALSE),
    min_f1 = 0.7,
    batch_correct = TRUE,
    ref.batch = NULL,
    xdim = 6, ydim = 6,
    identify_unassigned = TRUE,
    MAD_factor = 2.5,
    norm_method = "scale",
    covar = NULL,
    load_model = NULL,
    save_model = NULL,
    mtry = floor(length(markers)/3),
    splitrule = "gini",
    probability = TRUE,
    use.weights = TRUE,
    min.node.size = 1,
    num.trees = 600,
    num.threads = 1,
    mc.cores = NULL,
    save_adapted_reference = NULL,
    exclude_redundant = FALSE,
    exclude_celltypes = c(
      "Doublet",
      "Platelet",
      "Eryth",
      "HSPC"
    ),
    unassigned_name = "unassigned",
    identify_type = "probability",
    probability_threshold = 0.5,
    train_on_unassigned = FALSE,
    seed = 332,
    verbose = TRUE,
    ...) {

  if (inherits(reference, "character")) {
    if (reference == "pbmc" | reference == "pbmc_reference") {
      reference <- get_reference("pbmc", verbose = verbose)
    }
  }
  # Use mc.cores if given (for compatibility with cyCombine)
  if (!is.null(mc.cores) & num.threads == 1) num.threads <- mc.cores
  if (train_on_unassigned) identify_type <- "trained"


  if (adapt_reference) {

    if (verbose) message("Adapting reference to query marker panel")

    t <- system.time({
    reference <- adapt_reference(
      reference = reference,
      markers = markers,
      num.threads = num.threads,
      mtry = mtry,
      min_f1 = min_f1,
      using_pbmc = using_pbmc,
      exclude_celltypes = exclude_celltypes,
      verbose = verbose
      )
    })
    if (verbose) message(
      "Reference adaptation took ", round(t[[3]], 2), " seconds")
    # Store adapted reference
    if (!is(save_adapted_reference, "NULL")) {
      if (is(save_adapted_reference, "logical")){
        if (save_adapted_reference) saveRDS(reference, "adapted_reference.rds")
      } else {
        if (grepl("rds", tolower(save_adapted_reference))) {
          saveRDS(reference, save_adapted_reference)
        } else {
          saveRDS(reference, file.path(save_adapted_reference, "adapted_reference.rds"))
        }
      }
    }
  }

  if (batch_correct) {
    t <- system.time({
    # Batch correction via cyCombine
    corrected <- batch_correct(
      reference = reference,
      query = query,
      markers = markers,
      xdim = xdim,
      ydim = ydim,
      norm_method = norm_method,
      ref.batch = ref.batch,
      seed = seed,
      verbose = verbose,
      num.threads = ifelse(is.null(mc.cores), num.threads, mc.cores),
      ...
      )
    })
    if (verbose) message(
      "Batch correction took ", round(t[[3]], 2), " seconds")
    reference <- corrected$reference
    query <- corrected$query

  }

  # Optionally exclude redundant cells in the reference
  if (exclude_redundant) {
    reference <- excl_redundant_populations(
      reference = reference,
      query = query,
      markers = markers,
      min_cells = 50,
      min_pct = 0.05,
      num.threads = num.threads,
      mtry = mtry,
      seed = seed,
      verbose = verbose
    )

  }

  # Canonical cell type assignment
  t <- system.time({
  classified <- classify_cells(
    reference = reference,
    query = query,
    markers = markers,
    unassigned_name = unassigned_name,
    load_model = load_model,
    save_model = save_model,
    mtry = mtry,
    splitrule = splitrule,
    min.node.size = min.node.size,
    num.trees = num.trees,
    use.weights = use.weights,
    probability = probability,
    num.threads = num.threads,
    seed = seed,
    verbose = verbose
    )
  })
  if (verbose) message(
    "Classification took ", round(t[[3]], 2), " seconds")
  rm(query)


  # Identification of unassigned cells
  if (identify_unassigned) {
    t <- system.time({
    classified$query <- identify_unassigned(
      query = classified$query,
      reference = reference,
      markers = markers,
      model = classified$model,
      mtry = mtry,
      num.threads = num.threads,
      unassigned_name = unassigned_name,
      identify_type = identify_type,
      probability_threshold = probability_threshold,
      train_on_unassigned = train_on_unassigned,
      MAD_factor = MAD_factor,
      seed = seed,
      verbose = verbose
      )
    })
    if (verbose) message(
      "Outlier detection took ", round(t[[3]], 2), " seconds")
  }

  return (list(
    "query" = classified$query,
    "reference" = reference,
    "model" = classified$model
    ))
}


