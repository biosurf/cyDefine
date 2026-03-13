#' Modify query marker names to match reference marker names
#'
#' @param using_pbmc Boolean indicating whether or not the Universal PBMC atlas is
#' used as reference
#' @param ref_markers Character vector of reference markers
#' @param map_specific_from Optional: used together with `map_specific_to`.
#' Character vector of specific marker(s) of the query to manually map to
#' specific marker(s) of the reference (given by `map_specific_to`).
#' @param map_specific_to Optional: used together with `map_specific_from`.
#' Markers of the reference that the markers given by `map_specific_from` should
#' be mapped to.
#' @param verbose Verbosity
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param query_markers Character vector of query markers
#' @param return_unrecognized Boolean specifying whether to return names of
#' unrecognized markers.
#'
#' @importFrom dplyr rename_at vars all_of select
#' @importFrom stringr str_replace_all
#'
#' @return Tibble of query data modified to make marker names match reference.
#' If return_unrecognized is TRUE, returns a list with this tibble and a
#' character vector with the unrecognized markers.
#' @family adapt
#'
#' @examples
#'
#' map_marker_names(
#'   using_pbmc = TRUE, query = example_query,
#'   query_markers = example_markers, map_specific_from = c("CXCR5", "TCRgd"),
#'   map_specific_to = c("CD185", "gdTCR")
#' )
#'
#' @export
map_marker_names <- function(
    query,
    query_markers,
    ref_markers = cyDefine::pbmc_markers,
    using_pbmc = TRUE,
    map_specific_from = NULL,
    map_specific_to = NULL,
    return_unrecognized = FALSE,
    verbose = TRUE) {
  # check required columns
  check_colnames(
    avail_colnames = colnames(query),
    req_colnames = c(query_markers)
  )

  if (!is.null(map_specific_from)) {
    if (verbose) {
      message("Mapping specified markers to reference markers:\n",
        paste(map_specific_from,
          "to",
          map_specific_to,
          collapse = "\n"
        ),
        sep = "\n"
      )
    }

    query <- query |>
      dplyr::rename_at(
        dplyr::vars(dplyr::all_of(map_specific_from)),
        ~map_specific_to
      )

    # update data markers
    query_markers <- c(
      query_markers[query_markers %!in% map_specific_from],
      map_specific_to
    )
  }

  avail_markers_in_ref <- query_markers[query_markers %in% ref_markers]


  # further investigation if not all markers are found in reference
  if (length(avail_markers_in_ref) != length(query_markers)) {
    not_in_ref <- setdiff(query_markers, ref_markers)

    if (using_pbmc) {
      # common marker names
      common_names <- list(
        "CXCR5" = "CD185", "CCR4" = "CD194", "CCR6" = "CD196",
        "TIM3" = "CD366", "TIM-3" = "CD366", "PDL1" = "CD274",
        "PD-L1" = "CD274", "PD1" = "CD279", "PD-1" = "CD279",
        "CTLA4" = "CD152", "TCRgd" = "gdTCR", "TCRab" = "abTCR",
        "CCR2" = "CD192", "HER2" = "CD340", "MDR1" = "CD243",
        "IL-7R" = "CD127", # "B220"="CD45", "CD45R"="CD45",
        "IL-7Ra" = "CD127", "NKP46" = "CD335", "GYPA" = "CD235A",
        "GPA" = "CD235A", "Ter119" = "CD235A", "TRAIL" = "CD253",
        "TNFSF10" = "CD253", "CXCR4" = "CD184", "OX40" = "CD134", #"CCR7" = "CD197",
        "TNFRSF4" = "CD134", "ACT35" = "CD134", "CCR5" = "CD195", #"CXCR3" = "CD183",
        "CXCR6" = "CD186", "TCRVa7.2" = "TCR-V-7.2", "Va7.2" = "TCR-V-7.2",
        "TCRV7.2" = "TCR-V-7.2", "TRAV1-2" = "TCR-V-7.2", "Vg9" = "TCR-V-9", "TCRVg9" = "TCR-V-9",
        "TCRVa24-JaQ" = "TCR-V-24-J-18", "TCRVa24-Ja18" = "TCR-V-24-J-18",
        "CCR9" = "CD199", "CDw199" = "CD199", "ICOS" = "CD278",
        "LAG3" = "CD223", "CD40L" = "CD154", "CD40LG" = "CD154", "CD41A" = "CD41",
        "PDL2" = "CD273", "CD3" = "CD3E"#, "NKG2" = "CD159", "NKG2A" = "CD159a"
      )

      old_names <- intersect(names(common_names), not_in_ref)
      new_names <- as.character(common_names[old_names])

      # update not_in_ref
      not_in_ref <- not_in_ref[not_in_ref %!in% old_names]
    } else {
      old_names <- c()
      new_names <- c()
    }

    # if still un-mapped markers, check minor differences in names
    for (name in not_in_ref) {
      # get index of reference marker matching name
      # ignore non-alphabetic and non-numeric signs + case
      ref_marker_idx <- which(toupper(stringr::str_replace_all(
        string = ref_markers,
        pattern = "[^[:alnum:]]",
        replacement = ""
      )) ==
        toupper(stringr::str_replace_all(
          string = name,
          pattern = "[^[:alnum:]]",
          replacement = ""
        )))
      if (length(ref_marker_idx) == 1) {
        old_names <- c(old_names, name)
        new_names <- c(new_names, ref_markers[ref_marker_idx])
      }
    }

    if (length(old_names) > 0) {
      # rename markers
      if (verbose) {
        message(paste("Renaming",
          old_names,
          "to",
          new_names,
          sep = " ",
          collapse = "\n"
        ))
      }

      query <- query |>
        dplyr::rename_at(dplyr::all_of(old_names), ~new_names)

      not_in_ref <- not_in_ref[not_in_ref %!in% old_names]
    }

    # remove markers not found in reference
    if (length(not_in_ref) > 0) {
      warning(
        "The following markers were not detected to be in the reference and were excluded from the data:\n  ",
        paste(not_in_ref, collapse = ", "), "\n  ",
        "Markers of the Universal PBMC reference can be found in 'pbmc_markers'. ",
        "If needed, markers can be manually mapped to reference markers using the 'map_specfic_from' and 'map_specfic_to' arguments"
      )

      query <- query |>
        dplyr::select(-dplyr::all_of(not_in_ref))
    }
  } else {
    not_in_ref <- c()
  }
  if (return_unrecognized) {
    return(list(query, not_in_ref))
  } else {
    return(query)
  }
}




#' Merge groups of similar cell populations
#'
#' @inheritParams cyDefine
#' @param populations_to_merge A data frame of clustered populations
#' @param reference Data frame of reference data (cells in rows, markers in columns)
#' @note param name_unassigned_if_merged Cell types of the labeling tree, where if merged, the merged population should rather be named 'unassigned'
#' @family adapt
merge_populations <- function(
    populations_to_merge,
    reference,
    # name_unassigned_if_merged = "unassigned",
    using_pbmc = TRUE
    ) {
  # get labels of clusters of merged populations
  populations_to_merge <- populations_to_merge |>
    dplyr::group_by(cluster) |>
    dplyr::reframe(
      popu = popu,
      merged_label = get_merged_label(popu, using_pbmc = using_pbmc)
    )
  # |>
  #   dplyr::mutate(
  #     merged_label = dplyr::case_when(
  #       stringr::str_detect(merged_label, name_unassigned_if_merged) ~ name_unassigned_if_merged,
  #       TRUE ~ merged_label
  #     )
  #   )

  # modify cell type labels in reference to correspond to clusters
  reference <- dplyr::left_join(reference,
    populations_to_merge,
    by = c("celltype" = "popu")
  ) |>
    dplyr::mutate(
      celltype = dplyr::coalesce(
        merged_label,
        celltype
    )) |>
    dplyr::select(-merged_label, -cluster)

  return(reference)
}



#' Get label of a merged population
#'
#' @param populations Character vector of individual cell types in merged population
#' @inheritParams merge_populations
#' @importFrom stats na.omit
#' @return Label of merged population
#' @family adapt
get_merged_label <- function(populations, using_pbmc = TRUE) {
  if (using_pbmc) {
    check_package("reshape2")

    levels <- c(stringr::str_c(
      "L",
      seq_len(ncol(labeling_tree) - 1)
    ), "leaf")
    i <- 1

    while (i < length(levels) - 1) {
      parent_popus <- unique(na.omit(labeling_tree[[levels[i]]]))

      for (parent in parent_popus) {
        leaf_popus <- dplyr::filter(
          labeling_tree,
          get(levels[i]) == parent
        ) |>
          dplyr::pull(leaf)

        # if all popus of parent node are to be merged
        if (length(setdiff(leaf_popus, populations)) == 0) {
          populations <- c(populations[(populations %!in% leaf_popus)], parent)
        }
      }
      i <- i + 1
    }
  }

  # concatenate remaining populations
  label <- paste(populations,
    collapse = " / "
  )

  return(label)
}




#' Identify indistinguishable cell type populations in the reference, based on
#' query marker panel
#'
#' @inheritParams adapt_reference
#' @inheritParams ranger::ranger
#' @note param population_pairs Data frame of two character columns with one cell
#' type in each, and rows representing all pairs of cell types to be tested.
#' @param reference desc
#' @param min_f1 desc
#' @param seed desc
#' @param verbose desc
#' @importFrom dplyr slice_sample group_by group_modify ungroup n
#' @importFrom igraph components graph_from_data_frame
#'
#' @return A tibble of similar populations
#'
identify_similar_populations <- function(
    reference,
    markers,
    num.threads = 4,
    mtry = 4,
    min_f1 = 0.7,
    seed = 332,
    verbose = TRUE
    ) {
  check_colnames(colnames(reference), c(markers, "celltype"))

  # check number of cells of each cell type
  celltype_counts <- table(reference$celltype)
  small_populations <- names(celltype_counts[celltype_counts < 5])
  if (length(small_populations) > 0) {

    if (verbose) {
      message(
        "OBS: Less than 5 cells of type(s):\n",
        paste(small_populations,
          collapse = "\n"
        ),
        "\nare present in the reference. These will not be ",
        "considered for population merging."
      )
    }

    reference <- reference[reference$celltype %!in% small_populations, ]
  }


  if (verbose) {
    message("Running classification to identify similar populations")
  }

  # split data 50% by celltype
  reference <- reference |>
    dplyr::group_by(celltype) |>
    dplyr::mutate(train_idx = sample(c(TRUE, FALSE), size = dplyr::n(), replace = TRUE, prob = c(0.5, 0.5))) |>
    dplyr::ungroup()

  tmp_ref <- reference[reference$train_idx, ]
  tmp_query <- reference[!reference$train_idx, ]


  y_pred <- classify_cells(
    reference = tmp_ref,
    query = tmp_query,
    markers = markers,
    num.threads = num.threads,
    mtry = mtry,
    return_pred = TRUE,
    seed = seed,
    verbose = FALSE
  )

  # ensure same levels for pred and obs
  y_test <- factor(tmp_query$celltype)
  y_pred <- factor(y_pred,
    levels = levels(y_test)
  )


  # get celltypes with low one-vs-rest F1 score
  f1 <- compute_f1(y_pred, y_test)
  merging_candidates <- names(f1[f1 < min_f1])

  if (length(merging_candidates) > 0) {
    # get confusion matrix of merging candidates
    similar_population_pairs <- create_confusion_matrix(y_pred, y_test) |>
      dplyr::filter(
        observed %in% merging_candidates,
        observed != predicted
      ) |>
      dplyr::group_by(observed) |>
      dplyr::reframe(predicted = predicted[n == max(n)])

    population_clusters <- igraph::components(
      igraph::graph_from_data_frame(similar_population_pairs)
    )$membership

    similar_populations <- tibble::tibble(
      popu = names(population_clusters),
      cluster = population_clusters
    )
  } else {
    return(tibble::tibble())
  }


  return(similar_populations)
}

#' Exclude Redundant Cell Populations from Reference
#'
#' This function filters out redundant cell types from a reference dataset based on their counts and proportions in the query dataset. It uses an initial projection to classify the cell types and then excludes those with counts below a specified minimum or a proportion below a specified percentage.
#'
#' @param min_cells An integer specifying the minimum number of cells for a cell type to be retained. Default is 50.
#' @param min_pct A numeric value specifying the minimum percentage of cells for a cell type to be retained. Default is 0.001.
#' @inheritParams adapt_reference
#' @inheritParams classify_cells
#' @return A data frame containing the filtered reference dataset with redundant cell types excluded.
#' @export
#'
excl_redundant_populations <- function(
    reference,
    query,
    markers,
    num.threads = 4,
    mtry = 22,
    min_cells = 50,
    min_pct = 0.001,
    seed = 332,
    verbose = TRUE
    ) {

  if (verbose) {
    message("Making initial projection to filter out redundant cell types of the reference")
  }

  y_pred <- classify_cells(
    reference = reference,
    query = query,
    markers = markers,
    num.threads = num.threads,
    mtry = mtry,
    return_pred = TRUE,
    seed = seed,
    verbose = FALSE
  )


  celltype_counts <- summary(factor(y_pred,
    levels = unique(reference$celltype)
  ))

  excl_celltypes <- tibble::tibble(
    celltype = names(celltype_counts),
    n = celltype_counts
  ) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::filter(n < min_cells | pct < min_pct) |>
    dplyr::pull(celltype)

  if (verbose) {
    if (length(excl_celltypes) == 0) {
      message("No cell types found to be redundant, continuing with all populations")
    } else {
      message(
        "Excluding the following redundant celltypes from the reference: \n",
        paste(excl_celltypes, collapse = ", ")
      )
    }
  }

  reference <- reference[reference$celltype %!in% excl_celltypes,]

  return(reference)
}




#' Adapt PBMC reference to query marker panel
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param markers Character vector of available markers present in both query
#' and reference.
#' @param initial_project Boolean indicating whether cyDefine should perform an
#' initial projection to remove redundant cell types of the reference. Set to
#' FALSE, if you want to keep all cell types of the reference. Defaults to TRUE.
#' @param min_f1 Minimum F1 score required for two cell types to be separated
#' @param exclude_celltypes Only relevant if `using_pbmc = TRUE`. Character
#' vector of cell types of the Universal PBMC reference to NOT consider during cell
#' type classification. Defaults to non-PBMCs.
#' @param celltype_col Column with celltype information.
#' @inheritParams map_marker_names
#' @inheritParams merge_populations
#' @inheritParams identify_similar_populations
#' @inheritParams classify_cells
#' @inheritParams ranger::ranger
#' @importFrom stats aggregate
#'
#' @return Tibble of adapted reference
#' @export
#'
adapt_reference <- function(
    reference,
    markers,
    using_pbmc = TRUE,
    initial_project = FALSE,
    exclude_celltypes = c(
      "Doublet",
      "Platelet",
      "Eryth",
      "HSPC"
      ),
    celltype_col = ifelse(
      using_pbmc,
      "celltype.l2",
      "celltype"
      ),
    min_f1 = 0.7,
    seed = 332,
    num.threads = 4,
    mtry = 22,
    verbose = TRUE
    ) {
  # remove excluded cell types
  names(reference)[names(reference) == celltype_col] <- "celltype"
  reference <- reference[reference$celltype %!in% exclude_celltypes, ]
  reference$celltype_original <- reference$celltype

  reference <- check_id(reference)

  # check for similar cell types for as long as new merges are performed
  set.seed(seed)
  new_queries <- TRUE
  i <- 1
  while (new_queries) {
    if (verbose) {
      message("# ------- Population merging - Round ", i, " ------- #")
    }

    # find similar populations for all vs all remaining populations
    similar_populations <- identify_similar_populations(
      reference = reference,
      markers = markers,
      mtry = mtry,
      num.threads = num.threads,
      min_f1 = min_f1
    )

    # merge groups of similar populations
    if (nrow(similar_populations) > 0) {
      if (verbose) message("Merging groups of similar populations")

      reference <- merge_populations(
        populations_to_merge = similar_populations,
        reference = reference,
        using_pbmc = using_pbmc
      )


      if (verbose) {
        message(
          "Merging:\n\t",
          stats::aggregate(
            popu ~ cluster,
            data = similar_populations,
            FUN = function(x) paste(x, collapse = ", ")
            )[, 2] |>
            paste(collapse = "\n\t")
          , "\n"
        )
      }
    } else {
      new_queries <- FALSE
    }
    # reference$cluster <- NULL

    i <- i + 1
  }

  if (verbose) {
    message("\nReference adapted!")
  }

  return(reference)
}
