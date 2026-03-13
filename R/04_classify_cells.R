
#' Compute macro F1
#' @noRd
compute_f1 <- function(y_pred, y_test) {
  y_pred <- as.character(y_pred)
  y_test <- as.character(y_test)
  # Get the unique labels
  labels <- unique(c(y_pred, y_test))

  # Initialize a named vector to store the F1 score for each label
  f1_scores <- numeric(length(labels))
  names(f1_scores) <- labels

  # Compute precision, recall, and F1 score for each label
  for (i in seq_along(labels)) {
    label <- labels[i]

    # True positives, false positives, false negatives
    tp <- sum(y_pred == label & y_test == label)
    fp <- sum(y_pred != label & y_test == label)
    fn <- sum(y_pred == label & y_test != label)

    # Calculate precision and recall
    precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
    recall <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)

    # Calculate F1 score
    f1_scores[i] <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  }

  return(f1_scores)
}

#' Create confusion matrix
#' @importFrom dplyr rename
#' @noRd
create_confusion_matrix <- function(y_pred, y_test) {
  table(y_pred, y_test) |>
    as.data.frame() |>
    dplyr::rename(predicted = y_pred, observed = y_test, n = Freq)
}



#' Summary function for using macro-average F1-score as metric in caret::trainControl
#'
#' @noRd
#' @return Macro-average F1 score
#'
macroF1_summary <- function(data) {
  # multiclass classification performance
  f1 <- compute_f1(data$pred, data$obs)
  f1[is.na(f1)] <- 0
  f1 <- mean(f1)
  names(f1) <- "macroF1"
  return(f1)
}



#' Classify canonical cell types
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @param subsample Number of cells to use for training from each cell type
#' @param unassigned_name Name used for unassigned cells
#' @param use.weights (Default: FALSE) Boolean whether to use celltype frequency as class weights
#' @param load_model Optional: Path to an rda file of a previously trained model
#' @param save_model Optional: Path to save an rda file of the trained model
#' @param return_pred Boolean indicating if only predictions should be returned as a character vector
#' @param seed Random seed
#' @param verbose Verbosity
#' @inheritParams ranger::ranger
#' @importFrom stats predict
#' @return Tibble of query data with an added column of the predicted cell type, 'model_prediction'
#' @export
#'
classify_cells <- function(
    reference,
    query,
    markers,
    subsample = 2000,
    mtry = floor(length(markers)/3),
    min.node.size = 1,
    splitrule = "gini",
    num.trees = 300,
    use.weights = TRUE,
    probability = TRUE,
    unassigned_name = "unassigned",
    load_model = NULL,
    save_model = NULL,
    num.threads = 1,
    return_pred = FALSE,
    seed = 332,
    verbose = TRUE) {
  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c(markers))

  stopifnot("Please select an mtry smaller than the number of features" =
              mtry <= ncol(reference))

  # remove unassigned cells from reference prior to classification
  reference <- reference |>
    dplyr::filter(celltype != !!unassigned_name)

  # keep track of ids
  if ("id" %!in% colnames(query)) {
    query$id <- 1:nrow(query)
  }
  if ("id" %!in% colnames(reference)) {
    reference$id <- 1:nrow(reference)
  }


  if (!is.null(load_model)) {
    if (inherits(load_model) == "ranger") {
      rf_model <- load_model
      rm(load_model)
    } else {
      if (verbose) message("Loading saved model: ", load_model)
      if (endsWith(toupper(load_model), "RDS")) {
        rf_model <- readRDS(load_model)
      } else if (endsWith(toupper(load_model), "RDA")) {
        load(load_model)
        model_name <- gsub("\\..*$", "", basename(load_model))
        rf_model <- eval(parse(text = model_name))
      }
    }


  } else {
    if (verbose) message(
      "Training random forest model using ", num.threads, " threads")

    set.seed(seed)
    subset <- reference |>
      dplyr::group_by(celltype) |>
      dplyr::slice_sample(n = subsample)

    if (use.weights) {
      model_weights <- reference |>
        dplyr::group_by(celltype) |>
        dplyr::reframe(weight = dplyr::n() / nrow(reference)) |>
        dplyr::pull(weight)
        # dplyr::mutate(weight = nrow(reference) / dplyr::n()) |>
        # dplyr::ungroup() |>
        # dplyr::mutate(weight = weight / sum(weight)) |>
        # dplyr::filter(id %in% subset$id) |>
        # dplyr::pull(weight)
    } else {model_weights <- NULL}


    t <- system.time({
    # Set up your ranger model
    rf_model <- ranger::ranger(
      y = as.factor(subset$celltype),          # Formula interface, assuming 'data' includes all predictors
      x = subset[, markers],                   # Data frame containing the variables
      num.trees = num.trees,                   # Number of trees
      mtry = mtry,                             # Number of variables to possibly split at in each node
      importance = "impurity",       # Type of importance metric
      write.forest = TRUE,                     # Save the forest model
      probability = probability,               # If you need probability output for classification
      min.node.size = min.node.size,           # Minimum node size
      replace = TRUE,                          # Sampling with replacement
      sample.fraction = 1,                     # Fraction of observations to sample (1 for replacement)
      # case.weights = model_weights,
      class.weights = model_weights,           # Celltype weights
      splitrule = splitrule,                   # Splitting rule
      num.threads = num.threads,               # Number of threads to use
      seed = seed,                             # Seed for reproducibility
      save.memory = FALSE,                     # Save memory mode
      verbose = verbose,                       # Show progress and output
      oob.error = FALSE                        # Compute out-of-bag error estimate
    )
    })
    if (verbose) message(
      "Model training took ", round(t[[3]], 2), " seconds")

    # save model
    if (!is.null(save_model)) {
      saveRDS(rf_model, file = save_model)
    }
  }

  if (verbose) message("Predicting..")
  rf_model$predictions <- stats::predict(
    object = rf_model,
    data = query[, markers]
  )$predictions

  model_prediction <- rf_model$predictions
  if (probability) {
    model_prediction <- colnames(model_prediction)[apply(model_prediction, 1, which.max)]
  }

  # ARI: aricode::ARI(pred, query$celltype)
  if (return_pred) {
    return(model_prediction)
  }

  # add model predictions to query
  query <- dplyr::bind_cols(query, "model_prediction" = model_prediction) |>
    dplyr::arrange(id)

  return(list("query" = query, "model" = rf_model))
}
