#' Compute Mahalanobis distance from each cell to a given population
#'
#' @param cells Tibble of cells (cells in rows, markers in columns) from each of which the Mahalanobis distance to the population should be computed. ???
#' @param population Tibble of a population of cells (cells in rows, markers in columns) to which the Mahalanobis distance should be computed. ???
#' @param weights Numeric vector of marker weights for computing the weighted Mahalanobis distance. ???
#'
#' @return Numeric vector of Mahalanobis distance to population for each cell.
#'
mahalanobis_distances <- function(cells,
                                  population,
                                  weights = NULL) {
  check_package("matrixcalc")

  if (nrow(population) < 3) {
    warning(
      "Less than 3 observations are available for population ",
      population
    )
    return(rep(NA, nrow(cells)))
  }

  # covariance and mean of population
  covar <- stats::cov(population)
  mu <- apply(population, 2, mean)

  if (!(matrixcalc::is.positive.definite(covar))) {
    warning(
      "Decomposition impossible for population ", population,
      " as the covariance matrix is not positive definite"
    )
    return(rep(NA, nrow(cells)))
  }

  # Mahalanobis distance from each cell to population
  if (is.null(weights)) {
    check_package("Rfast")
    distances <- apply(cells,
      1,
      Rfast::mahala,
      mu = mu,
      sigma = covar
    )
  } else {
    distances <- apply(as.matrix(cells),
      1,
      wmahalanobis,
      center = mu,
      cov = covar,
      weight = weights
    )
  }
  return(distances)
}



#' Use median absolute deviation to find the upper boundary for outlier detection
#'
#' @param distances Numeric vector of Mahalanobis distances
#' @param MAD_factor Coefficient used to define the boundary. Generally, 3 is
#' considered very conservative, 2.5 is moderately conservative, and 2 is poorly
#' conservative.
#'
#' @return Upper limit
#'
MAD_max_distance <- function(distances, MAD_factor = 3) {
  # remove outliers and find maximum distance
  MAD_score <- stats::mad(distances, na.rm = TRUE)
  upper_thres <- stats::median(distances, na.rm = TRUE) + MAD_factor * MAD_score

  return(upper_thres)
}





#' Identify unassigned cells
#' @inheritParams classify_cells
#' @param query Tibble of classified query data (cells in rows, markers in columns)
#' @param train_on_unassigned Boolean indicating whether unassigned cells should be included in model training.
#' Recommended when reference and query are samples stemming from the same experiment and unassigned
#' cells of the query are assumed to be representative of those found in the reference.
#' @param pct_expl_var The percentage threshold of explained variance for PC selection
#' @param MAD_factor Median average distance threshold of the outlier to its assigned population
#' @param model Random forest model for `identify_type = "probability`.
#' @param probability_threshold Threshold for random forest probability threshold. Use in conjunction with `unassigned_type = "probability"`.
#' @param identify_type (Default: probability) Type of outlier detection.
#'  Options are (multiple allowed in conjunction):
#'  "probability" - Cut cells with a lower probability than the `probability_threshold`.
#'  "trained" - Train on unassigned in the reference.
#'  "mad" - MAD threshold on mahalanobis distance in PC space per cell type.
#'  "maha" - Same as "mad", but in marker space.
#'  "sd" - Unassign cells with more than 20% of markers that lay more than two standard deviations from reference expression
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats quantile sd
#' @importFrom dplyr group_map
#' @return A tibble of unassigned cells
#' @export
#'
identify_unassigned <- function(reference,
                                query = NULL,
                                markers,
                                model = NULL,
                                num.threads = 1,
                                mtry = ceiling(length(markers)/3),
                                identify_type = "probability",
                                probability_threshold = 0.8,
                                train_on_unassigned = FALSE,
                                unassigned_name = "unassigned",
                                pct_expl_var = 0.95,
                                MAD_factor = 2.5,
                                seed = 332,
                                verbose = TRUE) {

  if (is(reference, "list")) {
    if (!is.null(reference$model)) model <- reference$model
    query <- reference$query
    reference <- reference$reference
  }

  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c("model_prediction", markers))

  query$predicted_celltype <- NULL

  # keep track of ids
  query <- check_id(query)

  # if (verbose) {message("Identifying unassigned cells using ", n_components,
  #                       " principal components and MAD-factor ", MAD_factor)}

  if ("trained" %in% identify_type) {
    unassigned_name <- reference$celltype[grepl(unassigned_name, reference$celltype)][1]
    # check dimensions
    n_unassigned <- reference |>
      dplyr::filter(celltype == !!unassigned_name) |>
      nrow()

    if (n_unassigned < 20) {
      warning(
        "Too few cells labeled '", unassigned_name, "' are present in the reference (",
        n_unassigned, ") to train on these cells. ",
        "Please modify the 'unassigned_name' or use the unsupervised analysis type.",
        "Skipping identification. Rerun with `identify_unassigned()`."
      )
      return(query)
    }

    if (verbose) {
      message("Identifying unassigned cells per predicted cell type")
    }
    preds <- pbmcapply::pbmclapply(mc.cores = 1,
      unique(query$model_prediction),
      function(popu) {
        # if (verbose) message("Filtering unassigned cells from ", popu)

        # cell type specific query and reference
        celltype_query <- dplyr::filter(
          query,
          model_prediction == popu
        )

        celltype_ref <- dplyr::filter(
          reference,
          celltype == popu | celltype == !!unassigned_name
        )

        # check dimension
        n_popu <- celltype_ref |>
          dplyr::filter(celltype == !!popu) |>
          nrow()

        if (n_popu < 30) {
          warning(
            "Too few cells are available for celltype '", popu,
            "' for modelling, thus no cells predicted to belong to ", popu,
            " will be identified as 'unassigned'"
          )

          return(dplyr::tibble(
            "id" = celltype_query$id,
            "predicted_celltype" = popu
          ))
        }

        y_pred <- classify_cells(
          reference = celltype_ref,
          query = celltype_query,
          markers = markers,
          num.threads = num.threads,
          mtry = mtry,
          unassigned_name = FALSE,
          return_pred = TRUE,
          # n_trees = 50,
          verbose = FALSE
        )


        return(dplyr::tibble(
          "id" = celltype_query$id,
          "predicted_celltype" = y_pred
        ))
      }
    )
    preds <- do.call(rbind, preds)

    query <- dplyr::left_join(query,
      preds,
      by = "id"
    )
  }
  if ("mad" %in% identify_type){
    if (verbose) {
      message("Identifying unassigned cells per predicted cell type")
    }
    distances <- pbmcapply::pbmclapply(
      mc.cores = num.threads,
      unique(query$model_prediction),
      function(popu) {
        # cell type specific query and reference
        celltype_query <- dplyr::filter(
          query,
          model_prediction == popu
        )

        celltype_ref <- dplyr::filter(
          reference,
          celltype == popu
        )

        # check dimensions
        if (nrow(celltype_ref) <= length(markers)) {
          warning(
            "Fewer observations than number of markers are available for celltype '", popu, "'\n",
            "Mahalanobis distance cannot be computed, thus no cells predicted to belong to ", popu, " will be identified as 'unassigned'"
          )

          return(dplyr::tibble(
            "id" = celltype_query$id,
            "max_distance" = rep(Inf, nrow(celltype_query)),
            "distance" = rep(-Inf, nrow(celltype_query))
          ))
        }

        celltype_all <- dplyr::bind_rows(
          celltype_query,
          celltype_ref
        )

        # compute cell type specific PCA
        pca_embed <- stats::prcomp(
          x = celltype_all[, markers],
          retx = TRUE,
          center = TRUE,
          scale. = TRUE
        )

        # find no of PCs
        cum_expl_var <- summary(pca_embed)$importance["Cumulative Proportion", ]
        n_components <- which(cum_expl_var >= pct_expl_var)[1]
        # n_components <- ncol(pca_embed$x)

        pca_all <- dplyr::as_tibble(pca_embed$x[, 1:n_components])

        # split back into data and reference
        query_pcs <- pca_all[1:nrow(celltype_query), ]
        ref_pcs <- pca_all[(nrow(celltype_query) + 1):nrow(pca_all), ]


        # compute weighted Mahalanobis distance from each reference cell to its own population
        weights <- summary(pca_embed)$importance["Proportion of Variance", ][1:n_components]

        ref_distances <- mahalanobis_distances(
          cells = ref_pcs,
          population = ref_pcs,
          weights = diag(weights, length(weights))
        )

        # get max distance, using MAD to sort out outliers
        max_distance <- MAD_max_distance(ref_distances,
          MAD_factor = MAD_factor
        )

        # get Mahalanobis distances to predicted populations for all cells
        query_distances <- mahalanobis_distances(
          cells = query_pcs,
          population = ref_pcs,
          weights = diag(weights, length(weights))
        )

        return(dplyr::tibble(
          "id" = celltype_query$id,
          "max_distance" = max_distance,
          "distance" = query_distances
        ))
      }
    )
    distances <- do.call(rbind, distances)

    query <- query |>
      dplyr::select(-dplyr::any_of(c("distance", "max_distance"))) |>
      dplyr::left_join(distances, by = "id")

    # change predicted population to "unassigned" for cells with distance > max_distance to their predicted celltype
    if (!"predicted_celltype" %in% colnames(query)) query$predicted_celltype <- as.character(query$model_prediction)
    query$predicted_celltype[query$distance > query$max_distance] <- unassigned_name
    # query <- query |>
    #   dplyr::mutate(predicted_celltype = case_when(
    #     distance > max_distance ~ unassigned_name,
    #     TRUE ~ predicted_celltype
    #   ))
  }
  if ("probability" %in% identify_type) {
    stopifnot("Please provide the RF model." = !is.null(model))
    stopifnot("Please run cyDefine with 'probability = TRUE'" = !is.null(dim(model$predictions)))

    # Identify cells with low certainty
    probs <- model$predictions
    max_prob <- apply(probs, 1, max)

    if (!"predicted_celltype" %in% colnames(query)) query$predicted_celltype <- as.character(query$model_prediction)
    query$predicted_celltype[max_prob < probability_threshold] <- unassigned_name
  }
  if ("cov" %in% identify_type) {
    if (verbose) {
      message("Identifying unassigned cells per predicted cell type")
    }
    distances <- pbmcapply::pbmclapply(
      mc.cores = num.threads,
      unique(query$model_prediction),
      function(popu) {
        # cell type specific query and reference
        celltype_query <- dplyr::filter(
          query,
          model_prediction == popu
        )

        celltype_ref <- dplyr::filter(
          reference,
          celltype == popu
        )

        # check dimensions
        if (nrow(celltype_ref) <= length(markers)) {
          warning(
            "Fewer observations than number of markers are available for celltype '", popu, "'\n",
            "Mahalanobis distance cannot be computed, thus no cells predicted to belong to ", popu, " will be identified as 'unassigned'"
          )

          return(dplyr::tibble(
            "id" = celltype_query$id,
            "max_distance" = rep(Inf, nrow(celltype_query)),
            "distance" = rep(-Inf, nrow(celltype_query))
          ))
        }

        celltype_all <- dplyr::bind_rows(
          celltype_query,
          celltype_ref
        )

        # compute cell type specific PCA
        pca_embed <- stats::prcomp(
          x = celltype_ref[, markers],
          retx = TRUE,
          center = TRUE,
          scale. = TRUE
        )


        # find no of PCs
        cum_expl_var <- summary(pca_embed)$importance["Cumulative Proportion", ]
        n_components <- which(cum_expl_var >= pct_expl_var)[1]
        # n_components <- ncol(pca_embed$x)

        pca_all <- dplyr::as_tibble(pca_embed$x[, 1:n_components])

        # split back into data and reference
        query_pcs <- predict(pca_embed, newdata = celltype_query)[, 1:n_components]
        # query_pcs <- pca_all[1:nrow(celltype_query), ]
        ref_pcs <- pca_all#pca_all[(nrow(celltype_query) + 1):nrow(pca_all), ]


        # compute weighted Mahalanobis distance from each reference cell to its own population
        weights <- summary(pca_embed)$importance["Proportion of Variance", ][1:n_components]

        ref_distances <- mahalanobis_distances(
          cells = ref_pcs,
          population = ref_pcs,
          weights = diag(weights, length(weights))
        )

        # get max distance, using MAD to sort out outliers
        max_distance <- MAD_max_distance(ref_distances,
                                         MAD_factor = MAD_factor
        )

        # get Mahalanobis distances to predicted populations for all cells
        query_distances <- mahalanobis_distances(
          cells = query_pcs,
          population = ref_pcs,
          weights = diag(weights, length(weights))
        )

        return(dplyr::tibble(
          "id" = celltype_query$id,
          "max_distance" = max_distance,
          "distance" = query_distances
        ))
      }
    )
    distances <- do.call(rbind, distances)

    query <- query |>
      dplyr::select(-dplyr::any_of(c("distance", "max_distance"))) |>
      dplyr::left_join(distances, by = "id")

    # change predicted population to "unassigned" for cells with distance > max_distance to their predicted celltype
    if (!"predicted_celltype" %in% colnames(query)) query$predicted_celltype <- as.character(query$model_prediction)
    query$predicted_celltype[query$distance > query$max_distance] <- unassigned_name



  }
  if ("sd" %in% identify_type) {
    check_package("purrr")
    check_package("tidyr")
    if (verbose) {
      message("Identifying unassigned cells per predicted cell type")
    }


    sd_threshold  <- 2    # how many SDs from mean = outlier marker
    marker_cutoff <- 0.2  # proportion of markers that must be outliers to flag cell

    ref_stats <- reference |>
      group_by(celltype) |>
      summarise(across(all_of(markers), list(mean = mean, sd = sd))) |>
      tidyr::pivot_longer(-celltype,
                   names_to = c("marker", ".value"),
                   names_sep = "_")

    # --- Step 2: Compute outlier scores on the REFERENCE itself ---
    # This tells us what a "normal" outlier score looks like per cell type
    ref_thresholds <- reference |>
      mutate(cell_id = row_number()) |>
      tidyr::pivot_longer(all_of(markers), names_to = "marker", values_to = "value") |>
      left_join(ref_stats, by = c("celltype", "marker")) |>
      group_by(cell_id, celltype) |>
      summarise(prop_outlier_markers = mean(abs(value - mean) > sd_threshold * sd),
                .groups = "drop") |>
      group_by(celltype) |>
      summarise(threshold = quantile(prop_outlier_markers, 0.95))  # 95th percentile of reference = threshold

    # --- Step 3: Score query cells and apply per-class thresholds ---
    final_labels <- query |>
      mutate(cell_id = row_number()) |>
      tidyr::pivot_longer(all_of(markers), names_to = "marker", values_to = "value") |>
      left_join(ref_stats, by = c("model_prediction" = "celltype", "marker")) |>
      group_by(cell_id, model_prediction) |>
      summarise(prop_outlier_markers = mean(abs(value - mean) > sd_threshold * sd),
                .groups = "drop") |>
      left_join(ref_thresholds, by = c("model_prediction" = "celltype")) |>
      mutate(label = if_else(prop_outlier_markers > threshold,
                             "unassigned", model_prediction))

    if (!"predicted_celltype" %in% colnames(query)) query$predicted_celltype <- as.character(query$model_prediction)
    query$predicted_celltype[final_labels$label == "unassigned"] <- unassigned_name
  }
  if ("maha" %in% identify_type) {
    check_package("purrr")
    check_package("tidyr")
    check_package("Rfast")

    # Compute per-celltype covariance matrices
    ref_params <- reference |>
      group_by(celltype) |>
      dplyr::group_map(~ {
        mat     <- .x |> select(all_of(markers)) |> as.matrix()
        cov_fit <- list(center = colMeans(mat), cov = cov(mat))
      }) |>
      rlang::set_names(group_keys(group_by(reference, celltype))$celltype)

    # Calibrate mahalanobis thresholds per celltype
    ref_thresholds <- purrr::imap_dfr(ref_params, ~ {
      mat   <- reference |>
        filter(celltype == .y) |>
        select(all_of(markers)) |>
        as.matrix()

      dists <- Rfast::mahala(mat, .x$center, .x$cov)
      tibble(model_prediction = .y, threshold = MAD_max_distance(dists, MAD_factor = MAD_factor))#quantile(dists, 0.95))#
    })

    # Score query cells

    # Preallocate distance vector
    mahal_dists <- numeric(nrow(query))

    for (cls in names(ref_params)) {
      idx <- which(query$model_prediction == cls)
      if (length(idx) == 0) next
      mahal_dists[idx] <- Rfast::mahala(as.matrix(query[idx, markers]),
                                        ref_params[[cls]]$center,
                                        ref_params[[cls]]$cov)
    }

    # Apply threshold on scores
    query_scores <- query |>
      mutate(
        cell_id    = row_number(),
        mahal_dist = mahal_dists
      ) |>
      left_join(ref_thresholds, by = "model_prediction") |>
      mutate(
        outlier_score = mahal_dist / threshold,
        label         = if_else(mahal_dist > threshold, "unassigned", model_prediction)
      )

    if (!"predicted_celltype" %in% colnames(query)) query$predicted_celltype <- as.character(query$model_prediction)
    query$predicted_celltype[query_scores$label == "unassigned"] <- unassigned_name
  }

  return(dplyr::arrange(query, id))
}
