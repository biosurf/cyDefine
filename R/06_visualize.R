#' Get distinct colors for a range of cell types
#'
#' @param populations Character vector of populations to get colors for
#' @param add_unassigned Boolean indicating whether the color 'black' should be
#' added for unassigned cells
#' @return Named list with values HEX colors and populations as names
#' @export
#'
get_distinct_colors <- function(populations, add_unassigned = TRUE) {
  # unassigned cells will be assigned to black
  populations <- populations[populations != "unassigned"]

  colors <- c(
    "#DC050C", "#7BAFDE", "#FDB462", "#33A02C", "#FB8072",
    "#1965B0", "#882E72", "#B2DF8A", "#B17BA6", "#A6761D",
    "#E7298A", "#55A1B1", "#E6AB02", "#7570B3", "#ba5ce3",
    "#FAE174", "#B5651D", "#E78AC3", "#aeae5c", "#FF7F00",
    "#56ff0d", "#AA0A0A", "#BEAED4", "#1e90ff", "#aa8282",
    "#0AC8D4", "#808000", "#7800FA", "#00FAFA", "#641A00",
    "#8DD3C7", "#666666", "#999999", "#d4b7b7", "#8600bf",
    "#00bfff", "#ffff00", "#D4E1C8", "#D470C8", "#64C870",
    "#64C80C", "#0C00FA", "#FA00FA", "#707A00"
  )

  colors <- colors[1:length(populations)]
  names(colors) <- populations

  # add unassigned cells
  if (add_unassigned) {
    colors <- c(colors, "unassigned" = "#000000")
  }

  return(colors)
}

#' Expand colors to new/merged labels
#'
#' @noRd
expand_colors <- function(adapted_celltypes, original_celltypes = NULL, colors = NULL) {

  unique_celltypes <- sort(unique(adapted_celltypes))
  if (is.null(colors)) {
    return(get_distinct_colors(unique_celltypes))
  }

  if (!all(unique_celltypes %in% names(colors))) {
    new_celltypes <- unique_celltypes[!unique_celltypes %in% names(colors)]
    cell_counts <- table(original_celltypes)
    colors[new_celltypes] <- sapply(new_celltypes, function(celltype) {
      sub_cells <- unique(original_celltypes[adapted_celltypes == celltype])
      sub_cell_counts <- cell_counts[sub_cells]
      major_celltype <- names(sub_cell_counts[sub_cell_counts == max(sub_cell_counts)])
      colors[major_celltype]
    })
  }
  return(colors)
}


#' Visualize reference and query by UMAPs
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @param colors Optional: Named list of colors per cell type. Values should be
#' HEX colors and names should be the unique cell types
#' @param col Column to color by.
#' @param query_col Column to color query by. Default "predicted_celltype"
#' @param ref_col Column to color reference by. Default "celltype"
#' @param build_umap_on Character describing what data should be applied for
#' performing the UMAP embedding. Options are "both" (default), "reference",
#' and "query". "Both" implies that the UMAP embedding is done on both the
#' reference and query data.
#' @param shuffle Boolean indicating whether rows should be shuffled before plotting
#' @param down_sample Boolean indicating whether cells should be down-sampled
#' @param sample_n Number of cells to sample from each of reference and query,
#' if down_sample = TRUE
#' @param return_data Boolean indicating if the data used for plotting should be
#' returned rather than a ggplot2 object
#' @param seed Random seed
#' @param verbose Verbosity
#' @inheritParams plot_embedding
#'
#' @return ggplot2 plot
#' @export
#'
plot_umap <- function(reference,
                      query = NULL,
                      markers = NULL,
                      colors = NULL,
                      col = NULL,
                      query_col = "predicted_celltype",
                      ref_col = "celltype",
                      build_umap_on = c("both", "reference"),
                      metric = "euclidean",
                      shuffle = TRUE,
                      down_sample = TRUE,
                      sample_n = 10000,
                      return_data = FALSE,
                      seed = 332,
                      verbose = TRUE,
                      title = c("Reference", "Query - predicted cell types"),
                      legend_title = "Cell type",
                      add_centroids =  c(FALSE, TRUE, "text", "label"),
                      highlight_labels = FALSE) {

  build_umap_on <- match.arg(build_umap_on)
  add_centroids <- as.character(add_centroids) |> match.arg(add_centroids)
  if (!is.null(col)) ref_col <- query_col <- col
  # check required packages
  check_package("ggplot2")
  check_package("uwot")
  check_package("patchwork")
  if (is(add_centroids, "character")) {
    check_package("ggforce")
    check_package("ggrepel")
  }


  # Extract data based on input type
  if (is(reference, "list")) {
    reference <- reference$reference
    query <- reference$query
  }
  if ("cell_id" %in% colnames(reference)) id <- "cell_id" else id <- "id"

  if (is.null(markers)) markers <- cyCombine::get_markers(reference)

  check_colnames(colnames(reference), c(ref_col, markers))
  if ("UMAP1" %in% colnames(reference)) {
    reference[, c("UMAP1", "UMAP2")] <- NULL
  }
  if (!is.null(query)) {
  check_colnames(colnames(query), c(query_col, markers))
    if ("UMAP1" %in% colnames(query)) {
      query[, c("UMAP1", "UMAP2")] <- NULL
    }
  }

  if (verbose) message("Generating UMAP")

  set.seed(seed)
  # shuffle data
  if (shuffle) {
    reference <- reference[sample(nrow(reference)), ]
    if (!is.null(query)) {
      query <- query[sample(nrow(query)), ]
    }
  }


  if (down_sample) {
    reference <- reference |>
      dplyr::slice_sample(n = min(sample_n, nrow(reference)))
    if (!is.null(query)) {
      query <- query |>
        dplyr::slice_sample(n = min(sample_n, nrow(query)))
    }
  }

  umap_data <- reference[, markers]

  if (!is.null(query)) {
    umap_data <- rbind(umap_data, query[, markers])
  }

  if (build_umap_on == "both" | is.null(query)) {
    if (verbose) {
      message("Computing UMAP embedding of all cells of reference and query")
    }

    # UMAP embedding of both reference and query
    full_umap <- uwot::umap(
        umap_data,
        n_neighbors = 15,
        min_dist = 0.2,
        metric = "euclidean",
        ret_model = FALSE
      )

    colnames(full_umap) <- c("UMAP1", "UMAP2")

    ref_umap <- full_umap[1:nrow(reference), ]
    if (!is.null(query)) {
      query_umap <- full_umap[(nrow(reference) + 1):nrow(full_umap), ]
    }

    # ref_umap <- full_umap[1:nrow(reference), ]
    # query_umap <- full_umap[(nrow(reference) + 1):nrow(full_umap), ]
  } else if (build_umap_on == "reference") {
    if (verbose) {
      message(
        "Computing UMAP embedding only of reference cells. ",
        "Be aware that this can hide potential novel populations in query!")
    }

    # UMAP embedding of reference
    ref_umap_embed <- reference[, markers] |>
      uwot::umap(
        n_neighbors = 15,
        min_dist = 0.2,
        metric = "euclidean",
        ret_model = TRUE
      )

    ref_umap <- ref_umap_embed$embedding
    colnames(ref_umap) <- c("UMAP1", "UMAP2")

    if (verbose) {
      message("Projecting query cells onto reference UMAP embedding")
    }

    # projection of new data
    query_umap <- query[, markers] |>
      uwot::umap_transform(model = ref_umap_embed)
    colnames(query_umap) <- c("UMAP1", "UMAP2")
  }

  reference <- cbind(reference, ref_umap)
  if (!is.null(query)) query <- cbind(query, query_umap)


  ref_plot <- plot_embedding(
    reference,
    add_centroids = add_centroids,
    col = ref_col,
    colors = colors,
    title = title[1],
    highlight_labels = highlight_labels,
    legend_title = legend_title)

  if (is.null(query)) {
    if (return_data) {
      return(list("data" = reference, "plot" = ref_plot))
    } else {
      return(ref_plot)
    }
  }

  query_plot <- plot_embedding(
    query,
    add_centroids = add_centroids,
    col = query_col,
    colors = colors,
    title = title[2],
    highlight_labels = highlight_labels,
    legend_title = legend_title)

  p <- ref_plot + ggplot2::theme(legend.position = "none") +
    query_plot + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  if (return_data) {
    list("reference" = reference, "query" = query, "plot" = p)
  } else {
    return(p)
  }
}


#' Add centroids to embedding plot
#'
#' Helper function that adds centroid labels or points to an existing ggplot2
#' embedding visualization. Centroids are calculated as the mean UMAP coordinates
#' for each group defined by the specified column.
#'
#' @param df Data frame containing UMAP coordinates and grouping variables
#' @param embedding_plot Existing ggplot2 plot object to add centroids to
#' @param add_centroids Character or logical indicating centroid display type.
#' Options are "text", "label", or TRUE (equivalent to "text")
#' @param col Character string specifying the column name to group by for centroids
#' @param highlight_labels Logical indicating whether to highlight labels with
#' ellipses for cell types containing "/" (uses ggforce::geom_mark_ellipse)
adding_centroids <- function(df, embedding_plot, add_centroids, col, highlight_labels) {
  check_package("ggforce")
  check_package("ggrepel")
  centroids <- df |>
    dplyr::group_by(.data[[col]]) |>
    dplyr::summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))


    if (add_centroids == "text" | add_centroids == "TRUE") {

      if (highlight_labels) {
        merge_list <- get_merge_list(df)
        embedding_plot <- embedding_plot +
          ggforce::geom_mark_ellipse(data = centroids |> dplyr::filter(celltype %in% names(merge_list)),
            ggplot2::aes(label = .data[[col]]), label.fontsize = 7.5,
              label.buffer = unit(1.5, "mm"),
            con.size = 0.3, con.type = "straight", con.cap = unit(2 , "mm"), con.border = "one",
                expand = unit(0.01, "mm"), label.fill = NA)

      } else {
        embedding_plot <- embedding_plot +
          ggrepel::geom_text_repel(
            data = centroids,
            aes(label = .data[[col]]),
            color = "black", size = 3, vjust = 1,
            max.overlaps = Inf)
        }

    } else {
      embedding_plot <- embedding_plot +
        ggrepel::geom_label_repel(
          data = centroids,
          ggplot2::aes(label = .data[[col]]),
          show.legend = FALSE,
          size = 2.5,
          box.padding = 0.5,
          point.padding = 1, max.overlaps = Inf)
    }
    embedding_plot <- embedding_plot +
      ggplot2::theme(legend.position = "none")


  return(embedding_plot)
}

#' Plot already generated UMAP coordinates
#'
#' Plot an already generated UMAP. Handles both categorical and continuous
#' coloring variables with appropriate scales and legends.
#'
#' @param embedding Data frame containing embedding coordinates (UMAP1, UMAP2) and
#' coloring variables
#' @param col Character string specifying the column name to color points by
#' @param colors Optional named vector of colors for categorical variables. If NULL,
#' colors are automatically generated using cyDefine::get_distinct_colors()
#' @param title Character string for plot title
#' @param add_centroids Character or logical indicating whether to add centroids.
#' Options are "text"/TRUE, "label", or FALSE (default)
#' @param highlight_labels Logical indicating whether to highlight centroid labels
#' with ellipses (default FALSE)
#' @param legend_title Manually assigned legend title.
#' @param slot used if embedding is a list. Specifies the slot with UMAPs. Default: "data"
#'
#' @return ggplot2 object showing the embedding visualization
#' @export
plot_embedding <- function(embedding, col, colors = NULL, title = "", add_centroids = c(FALSE, TRUE, "text", "label"), highlight_labels = FALSE,
                           legend_title = ggplot2::waiver(), slot = "data") {

  add_centroids <- as.character(add_centroids) |> match.arg(add_centroids)

  is_factor <- class(embedding[[col]]) != "numeric"

  if (is(embedding, "list")) embedding <- embedding[[slot]]

  if (is_factor) {

    # avoid too long cell type labels
    embedding <- embedding |>
      as.data.frame() |>
      dplyr::mutate("{col}" := stringr::str_wrap(embedding[[col]], 20))

    if (is.null(colors)) {
      colors <- cyDefine::get_distinct_colors(sort(unique(embedding[[col]])))
    } else {
      original_celltypes <- switch("celltype_original" %in% colnames(embedding)+1,
                                   embedding[[col]],
                                   embedding$celltype_original)
      colors <- cyDefine:::expand_colors(embedding[[col]], original_celltypes, colors = colors)
    }


    # get number of columns in legend
    n_legend_cols <- dplyr::case_when(
      length(colors) > 45 ~ 4,
      length(colors) > 30 ~ 3,
      length(colors) > 15 ~ 2,
      TRUE ~ 1
    )


    color_scale <- function(col, values, legend_title) {
      scale <- ggplot2::scale_colour_manual(col, values = values)
      guides <- ggplot2::guides(color = ggplot2::guide_legend(
        title = legend_title,
      override.aes = list(size = 2, alpha = 1),
      ncol = n_legend_cols
      ))
      return(list(scale, guides))
      }
  } else {
    color_scale <- function(col, ...) {}
  }


  embedding_plot <- embedding |>
    ggplot2::ggplot(ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      color = .data[[col]]
    )) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    color_scale(col, colors, legend_title)

  if (add_centroids != "FALSE") embedding_plot <-
    adding_centroids(embedding, embedding_plot, add_centroids, col, highlight_labels)

  return(embedding_plot)
}

#' Visualize multiple markers on UMAP embedding
#'
#' Create a multi-panel visualization showing the expression or intensity of
#' multiple markers overlaid on a UMAP embedding. Each
#' marker is displayed as a separate subplot with continuous color scaling.
#'
#' @param embedding Data frame containing embedding coordinates (UMAP1, UMAP2) and
#' marker expression values
#' @param markers Character vector of marker names (column names in embedding) to visualize
#' @param title Character string for the overall plot title (default "CyTOF UMAP")
#' @param ncol Integer specifying number of columns in the panel layout (default 4)
#' @param show_legend Logical indicating whether to display color legends for each
#' marker subplot (default TRUE)
#'
#' @return patchwork object containing multiple ggplot2 subplots arranged in a grid
#' @export
#'
plot_markers <- function(embedding, markers, title = "CyTOF UMAP", ncol = 4, show_legend = TRUE) {

  check_package("patchwork")
  # Create a list to store all plots
  plot_list <- list()

  # Generate a plot for each marker
  for (marker in markers) {
    # Create individual plot with marker name in title
    p <- plot_embedding(embedding,
                        col = marker,
                        title = marker,
                        highlight_labels = FALSE)

    # Add to plot list
    plot_list[[marker]] <- p
  }

  # Combine all plots using patchwork
  n_markers <- length(markers)
  nrow <- ceiling(n_markers / ncol)

  combined_plot <- patchwork::wrap_plots(plot_list,
                              ncol = ncol,
                              nrow = nrow)

  # Add overall title
  combined_plot <- combined_plot +
    patchwork::plot_annotation(
      title = paste(title, " - All Markers"),
      theme = ggplot2::theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  if (!show_legend) combined_plot <- combined_plot & ggplot2::theme(legend.position = "none")

  return(combined_plot)
}





#' Make bar plot of predicted cell type abundances
#'
#' @inheritParams plot_umap
#' @param predicted_populations Character vector of (predicted) populations to plot abundance for
#'
#' @return A ggplot2 bar plot of abundances
#' @export
#'
plot_abundance <- function(
  predicted_populations,
  colors = NULL,
  return_data = FALSE,
  verbose = TRUE
) {

  if (verbose) {
    message("Visualizing abundance of predicted cell types")
  }

  if (is.null(colors)) {
    colors <- get_distinct_colors(unique(predicted_populations))
  }

  freqs <- table(predicted_populations)  |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      prop = n / sum(n),
      label = paste0(round(100 * prop, 1), "%"),
      col = colors[predicted_populations]
    )

  if (return_data) {
    return(freqs)
  }

  p <- ggplot2::ggplot(freqs, ggplot2::aes(
      x = .data$predicted_populations,
      y = .data$prop,
      fill = .data$predicted_populations
    )) +
    ggplot2::geom_bar(
      stat = "identity",
      width = 0.9
    ) +
    ggplot2::geom_text(
      data = freqs,
      ggplot2::aes(
        # x = .data$predicted_populations,
        # y = .data$prop,
        label = paste0(round(100 * .data$prop, 1), "%")
      ),
      position = ggplot2::position_dodge(width = 1),
      vjust = -0.5,
      hjust = 0.5,
      size = 3
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 10
      ),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::ggtitle("Cell type abundance") +
    ggplot2::xlab("Predicted cell types") +
    ggplot2::ylab("Proportion of cells")

  return(p)
}


#' Compare cell type abundances between reference and query datasets
#'
#' Creates a grouped bar plot comparing the proportions of cell types between
#' reference and query datasets. The plot shows percentage abundance on the y-axis
#' and cell types on the x-axis, with bars grouped by dataset.
#'
#' @param reference Tibble of reference data or a list containing reference and query data.
#' If a list, should contain elements named 'reference' and 'query'
#' @param query Tibble of query data (cells in rows). Ignored if reference is a list
#' @param ref_col Character string specifying the column name containing cell types
#' in the reference dataset (default "celltype")
#' @param query_col Character string specifying the column name containing cell types
#' in the query dataset (default "model_prediction")
#' @param ref_name Character string for reference dataset label in legend (default "Reference")
#' @param query_name Character string for query dataset label in legend (default "Query")
#' @param colors Named vector of colors for datasets. If NULL, uses default green
#' for reference and orange for query
#' @param return_data Logical indicating whether to return the processed data instead
#' of the plot (default FALSE)
#'
#' @return ggplot2 object showing grouped bar plot of cell type abundances, or
#' data frame if return_data = TRUE
#' @export
#'
plot_abundance_comparison <- function(
    reference,
    query = NULL,
    ref_col = "celltype",
    query_col = "model_prediction",
    ref_name = "Reference",
    query_name = "Query",
    colors = NULL,
    return_data = FALSE
    ) {

  # Extract cell type columns based on input type
  if (is(reference, "list")) {
    ref_celltypes <- as.character(reference$reference[[ref_col]])
    query_celltypes <- as.character(reference$query[[query_col]])
  } else {
    ref_celltypes <- as.character(reference[[ref_col]])
    query_celltypes <- as.character(query[[query_col]])
  }

  # Calculate proportions for reference
  ref_freq <- table(ref_celltypes) |>
    tibble::as_tibble() |>
    dplyr::rename(celltype = ref_celltypes) |>
    dplyr::mutate(
      prop = n / sum(n),
      dataset = ref_name
    )

  # Calculate proportions for query
  query_freq <- table(query_celltypes) |>
    tibble::as_tibble() |>
    dplyr::rename(celltype = query_celltypes) |>
    dplyr::mutate(
      prop = n / sum(n),
      dataset = query_name
    )

  # Combine data
  combined_freq <- rbind(ref_freq, query_freq)

  # Get all unique cell types
  all_celltypes <- unique(c(ref_celltypes, query_celltypes))

  # Ensure all cell types are present in both datasets (with 0 if missing)
  complete_data <- tidyr::expand_grid(
    celltype = all_celltypes,
    dataset = c(ref_name, query_name)
  ) |>
    dplyr::left_join(combined_freq, by = c("celltype", "dataset")) |>
    dplyr::mutate(
      n = tidyr::replace_na(n, 0),
      prop = tidyr::replace_na(prop, 0)
    )

  if (return_data) {
    return(complete_data)
  }

  # Set colors if not provided
  if (is.null(colors)) {
    colors <- c("#228B22", "#FFA500")  # Green for reference, orange for query
    names(colors) <- c(ref_name, query_name)
  }

  # Create grouped bar plot
  p <- ggplot2::ggplot(complete_data) +
    ggplot2::aes(x = celltype, y = prop * 100, fill = dataset) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      x = "Cell Type",
      y = "Percentage (%)",
      fill = "Dataset"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "bottom",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"))

  return(p)
}



#' Plot heatmap of marker expressions per population
#'
#' @inheritParams plot_umap
#' @param data Tibble with data (cells in rows, markers in columns)
#' @param population_col Name of the cell type column in 'data' that you want to
#' plot expressions for
#' @param markers_to_plot Character vector of the markers that you want to
#' include in the heatmap
#' @param title Title of the plot
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @return a pheatmap heatmap
#' @export
#'
plot_heatmap <- function(data,
                         population_col = "predicted_celltype",
                         markers_to_plot,
                         title = "Average marker expression of predicted cell types\n",
                         return_data = FALSE,
                         verbose = TRUE) {
  check_package("pheatmap")
  check_colnames(colnames(data), c(population_col, markers_to_plot))

  if (verbose) {
    message("Visualizing expression of predicted cell types")
  }

  data <- data |>
    dplyr::group_by_at(population_col) |>
    dplyr::summarise_at(dplyr::all_of(markers_to_plot), mean) |>
    dplyr::mutate_at(dplyr::all_of(markers_to_plot), scale)

  plot_dat <- data |>
    dplyr::select(dplyr::all_of(markers_to_plot)) |>
    as.matrix() |>
    t()

  colnames(plot_dat) <- data[[population_col]]

  if (return_data) {
    return(plot_dat)
  }


  p <- pheatmap::pheatmap(
    mat = plot_dat,
    scale = "none",
    main = title,
    cluster_rows = TRUE,
    cluster_cols = TRUE,

    # breaks = my.breaks,
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
      n = 7, name =
        "RdBu"
    )))(100),
    legend = TRUE,
    border_color = "grey",
    treeheight_row = 30,
    treeheight_col = 30,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 90
  )

  return(p)
}

#' Plot a Diagram from Merged Cell Types
#'
#' This function generates a diagram that visualizes the merging of cell types based on
#' a provided list or data frame. The diagram illustrates how original cell types (nodes) merge into
#' final cell types.
#'
#' @param input A `list` or `data.frame`. If a `data.frame`, it should have been processed using
#'   `adapt_reference()` and must contain columns `celltype` and `celltype_original`.
#'   If a `list`, it should be structured such that each key represents a merged cell type and
#'   each value is a vector of original cell types.
#' @param colors A named vector where the names correspond to cell types (both original and merged),
#'   and the values are color codes to be used for filling the nodes in the diagram.
#' @param fontcolor_nodes A named vector (optional) that specifies font colors for specific nodes.
#'   If a node's font color is not provided in this vector, `default_fontcolor` will be used.
#' @param default_fontcolor A `character` string specifying the default font color to be used
#'   for nodes that are not included in `fontcolor_nodes`. Default is `"black"`.
#'
#' @return A rendered diagram produced by `DiagrammeR::grViz()` that visually represents the merging
#'   of cell types. The diagram includes nodes for both original and merged cell types, with edges
#'   showing the relationship between them.
#'
#' @details The function first checks if the input is a data frame or a list. If it is a data frame,
#'   it calls `get_merge_list()` to convert it into the required list structure. The function then
#'   generates a DOT language script to define the diagram, including node styling, edge creation,
#'   and graph layout. Finally, it uses `DiagrammeR::grViz()` to render the flowchart.
#'
#' @examples
#' \dontrun{
#' # Assuming `adapted_reference` is a preprocessed tibble from `adapt_reference()`
#' fontcolor_nodes <- c("NK" = "white", "NK Proliferating / NK" = "white")
#'
#' plot_diagram(adapted_reference, fontcolor_nodes)
#' }
#'
#' @seealso `adapt_reference()`, `get_merge_list()`
#'
#' @export
plot_diagram <- function(input, colors = NULL, fontcolor_nodes = c("unassigned" = "white"), default_fontcolor = "black") {

  check_package("DiagrammeR")

  if (inherits(input, "data.frame")) {
    if (!"celltype_original" %in% colnames(input)) {
      message("No populations merged in the adaptation.")
      return()
    }
    merge_list <- get_merge_list(input)
  } else if (inherits(input) == "list") {
    merge_list <- input
  }


  if (all(c("celltype", "celltype_original") %in% names(input)) & !is(colors, "NULL")) {
    colors <- expand_colors(input$celltype, input$celltype_original, colors)
  } else if (is.null(colors)) {
    colors <- cyDefine::get_distinct_colors(unique(unlist(merge_list)))
  }
  if (!all(names(merge_list) %in% names(colors))) {
    colors[names(merge_list)] <- sapply(names(merge_list), function(merged) {
      as.character(colors[stringr::str_split(merged, " /")[[1]][1]])
    })
  }

  rm(input)
  # Initialize the diagram
  diagram_code <- "
    digraph flowchart {
      # Global node styles
      node [fontname = Helvetica, style = filled, color = black, shape = box];
  "

  # To track the incrementing value for node names
  node_counter <- 0

  # Function to get a unique node name by appending an incrementing value
  get_unique_node_name <- function(node_name) {
    node_counter <<- node_counter + 1
    return(paste0(node_name, "_", node_counter))
  }

  # Loop through the merge list to create nodes and edges
  for (merge_into in names(merge_list)) {
    merged_nodes <- merge_list[[merge_into]]

    # Set the font color based on the input vector or default
    if (merge_into %in% names(fontcolor_nodes)) {
      fontcolor <- fontcolor_nodes[merge_into]
    } else if (colors[merge_into] %in% colors[names(fontcolor_nodes)]) {
      inherit_font <- colors[names(fontcolor_nodes)][colors[merge_into] == colors[names(fontcolor_nodes)]]
      fontcolor <- fontcolor_nodes[names(inherit_font)]
    } else {fontcolor <- default_fontcolor}

    # Create a unique identifier for the merge_into node
    merge_into_unique <- get_unique_node_name(merge_into)

    # Define the merged node
    diagram_code <- paste0(diagram_code, "
      \"", merge_into_unique, "\" [shape=ellipse, fillcolor='", colors[merge_into], "', fontcolor='", fontcolor, "', label = \"", merge_into, "\"];
    ")

    # Define the individual nodes and edges
    for (node in merged_nodes) {
      fontcolor <- ifelse(node %in% names(fontcolor_nodes),
                          fontcolor_nodes[node], default_fontcolor)

      # Create a unique identifier for the node
      node_unique <- get_unique_node_name(node)

      diagram_code <- paste0(diagram_code, "
        \"", node_unique, "\" [shape=box, fillcolor='", colors[node], "', fontcolor='", fontcolor, "', label = \"", node, "\"];
        \"", node_unique, "\" -> \"", merge_into_unique, "\";
      ")
    }
  }

  # Close the diagram definition
  diagram_code <- paste0(diagram_code, "
      # Set global graph attributes
      graph [layout = dot, rankdir = TB];
    }
  ")

  # Render the diagram
  DiagrammeR::grViz(diagram_code)
}

#' Get Merge List for create_diagram
#'
#' This function generates a list mapping merged cell types to their original cell types
#' based on the clustering results from an adapted reference dataset. The function assumes
#' that the dataset has been preprocessed using `adapt_reference()`.
#'
#' @param adapted_reference A `data.frame` or `tibble` that contains the adapted reference data.
#'   The data should include the columns `cluster`, `celltype`, and `celltype_original`.
#'   These columns are expected to be present in the dataset after running `adapt_reference()`.
#'
#' @return A named list where each element represents a merged cell type. The names of the list
#'   correspond to the unique cell types in the `celltype` column, and the values are lists of
#'   the original cell types (`celltype_original`) that were merged into the final cell type.
#'
#' @details The function filters the provided reference data to include only unique combinations
#'   of `cluster`, `celltype`, and `celltype_original`. It then groups the data by `cluster` and
#'   creates a list where each unique `celltype` is associated with its corresponding original
#'   cell types. The list is returned in a format suitable for use in generating flowcharts or
#'   diagrams that visualize the merging of cell types.
#'
#' @examples
#' \dontrun{
#' # Assuming `adapted_reference` is a preprocessed tibble from `adapt_reference()`
#' merge_list <- get_merge_list(adapted_reference)
#' }
#'
get_merge_list <- function(adapted_reference) {
  stopifnot(
    "Please run adapt_reference() before creating a diagram." =
      all(c("celltype", "celltype_original") %in%
            colnames(adapted_reference)))

  merged_populations <- adapted_reference |>
    dplyr::distinct(celltype_original, .keep_all = TRUE) |>
    dplyr::count(celltype) |>
    dplyr::filter(n > 1)

  # merge_list <- adapted_reference |>
  #   dplyr::distinct(celltype, celltype_original) |>
  #   dplyr::filter(celltype != celltype_original) |>
  #   dplyr::group_by(celltype) |>
  #   dplyr::reframe(
  #     celltype = unique(celltype),
  #     celltype_original = list(celltype_original)) |>
  #   tibble::deframe()

  merge_list <- adapted_reference |>
    dplyr::distinct(celltype, celltype_original) |>
    dplyr::filter(celltype %in% merged_populations$celltype) |>
    dplyr::group_by(celltype) |>
    dplyr::reframe(
      celltype = unique(celltype),
      celltype_original = list(celltype_original)) |>
    tibble::deframe()

  return(merge_list)
}

#' Create an Alluvial Plot
#'
#' This function creates an alluvial plot to compare true and predicted clusterings.
#'
#' @param df A data frame containing the classified data.
#' @param true A character string indicating the column name for true labels.
#' @param predicted A character string indicating the column name for predicted labels.
#' @param color A named vector of colors for the true labels.
#' @param color_by A character string indicating which variable to use for filling the alluvia ("true" or "predicted").
#' @param n An integer indicating the sample size.
#' @param seed An integer indicating the seed value for random sampling.
#' @param title Plot title
#' @param names Names of compared variables as printed on the plot
#'
#' @return A ggplot object representing the alluvial plot.
#'
#' @importFrom methods is
#' @importFrom dplyr sample_n
#' @import ggplot2
#'
#' @export
plot_alluvium <- function(
    df,
    true = "celltype",
    predicted = "model_prediction",
    color = NULL,
    color_by = c("true", "predicted"),
    names = c("True", "Predicted"),
    title = "True vs Predicted",
    n = 500,
    seed = 1286
    ) {
  check_package("ggalluvial")
  color_by <- match.arg(color_by)
  color_by <- if (color_by == "true") true else predicted
  stopifnot("Please set 'true' as a column in the input data." = true %in% colnames(df))
  stopifnot("Please set 'predicted' as a column in the input data." = predicted %in% colnames(df))
  if (is(color, "NULL")) {
    color <- get_distinct_colors(unique(df[[color_by]]), add_unassigned = F)
  }

  # Set the seed for reproducibility
  set.seed(seed)

  if (n < nrow(df)) {
    df <- dplyr::sample_n(df, size = n)
  }

  # Create the alluvial plot
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(axis1 = true, axis2 = predicted)) +
    ggalluvial::geom_alluvium(ggplot2::aes_string(fill = color_by), alpha = .7, reverse = TRUE) +
    ggalluvial::geom_stratum(width = 1/4) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, ggplot2::aes(label = ggplot2::after_stat(stratum))) +
    ggplot2::scale_x_discrete(limits = c("True", "Predicted"), expand = c(0.15, 0.05)) +
    ggplot2::scale_fill_manual(values = sort(color)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::ggtitle(title)

  return(plot)
}




#' Plot correlation of marker expression between reference and query by cell type
#'
#' Creates faceted scatter plots showing the correlation between mean marker expression
#' in reference and query datasets for each cell type. Each plot includes a linear
#' regression line and displays the R-squared correlation coefficient.
#'
#' @param reference Tibble of reference data or a list containing reference and query data.
#' If a list, should contain elements named 'reference' and 'query'
#' @param query Tibble of query data (cells in rows). Ignored if reference is a list
#' @param markers Character vector of marker names (column names) to include in analysis
#' @param ref_col Character string specifying the column name containing cell types
#' in the reference dataset (default "celltype")
#' @param query_col Character string specifying the column name containing cell types
#' in the query dataset (default "model_prediction")
#' @param ref_name Character string for reference dataset label (default "Reference")
#' @param query_name Character string for query dataset label (default "Query")
#' @param celltypes_to_plot Character vector of specific cell types to include. If NULL,
#' uses intersection of cell types present in both datasets
#' @param ncol Integer specifying number of columns in facet layout (default 3)
#' @param point_size Numeric value for point size in scatter plots (default 3)
#' @param marker_colors Named vector of colors for markers. If NULL, generates colors
#' automatically using cyDefine::get_distinct_colors()
#'
#' @importFrom stats cor
#'
#' @return ggplot2 object showing faceted correlation plots with R-squared values
#' @export
#'
plot_expression_correlation <- function(
    reference,
    query = NULL,
    markers,
    ref_col = "celltype",
    query_col = "model_prediction",
    ref_name = "Reference",
    query_name = "Query",
    celltypes_to_plot = NULL,
    ncol = 3,
    point_size = 3,
    marker_colors = NULL
    ) {

  # Extract data based on input type
  if (is(reference, "list")) {
    ref_data <- reference$reference
    query_data <- reference$query
  } else {
    ref_data <- reference
    query_data <- query
  }

  # Get cell types
  ref_celltypes <- as.character(ref_data[[ref_col]])
  query_celltypes <- as.character(query_data[[query_col]])

  # If specific cell types to plot are not provided, use common cell types
  if (is.null(celltypes_to_plot)) {
    celltypes_to_plot <- intersect(unique(ref_celltypes), unique(query_celltypes))
  }

  # Set colors for markers if not provided
  if (is.null(marker_colors)) {
    marker_colors <- cyDefine::get_distinct_colors(markers)#scales::hue_pal()(length(markers))
    names(marker_colors) <- markers
  }

  # Add cell type to data frames and select only markers
  ref_data_long <- ref_data |>
    dplyr::mutate(celltype = ref_celltypes) |>
    dplyr::select(celltype, dplyr::all_of(markers)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(markers),
                 names_to = "marker",
                 values_to = "ref_expr") |>
    dplyr::filter(celltype %in% celltypes_to_plot)

  query_data_long <- query_data |>
    dplyr::mutate(celltype = query_celltypes) |>
    dplyr::select(celltype, dplyr::all_of(markers)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(markers),
                 names_to = "marker",
                 values_to = "query_expr") |>
    dplyr::filter(celltype %in% celltypes_to_plot)

  # Calculate mean expression for each marker per cell type
  ref_means <- ref_data_long |>
    dplyr::group_by(celltype, marker) |>
    dplyr::summarise(ref_mean = mean(ref_expr, na.rm = TRUE), .groups = "drop")

  query_means <- query_data_long |>
    dplyr::group_by(celltype, marker) |>
    dplyr::summarise(query_mean = mean(query_expr, na.rm = TRUE), .groups = "drop")

  # Join the means
  combined_data <- ref_means |>
    dplyr::inner_join(query_means, by = c("celltype", "marker")) |>
    dplyr::rename(ref_expr = ref_mean, query_expr = query_mean)

  # Calculate correlations per cell type
  correlation_data <- combined_data |>
    dplyr::group_by(celltype) |>
    dplyr::summarise(
      r_squared = round(stats::cor(ref_expr, query_expr,
                                   method = "pearson",
                                   use = "complete.obs")^2, 2),
      .groups = "drop"
    )

  # Add correlation data back to combined data
  combined_data <- combined_data |>
    dplyr::left_join(correlation_data, by = "celltype")

  # Get ranges for annotation positioning
  annotation_data <- combined_data |>
    dplyr::group_by(celltype) |>
    dplyr::summarise(
      x_pos = max(ref_expr, na.rm = TRUE) * 0.95,
      y_pos = min(query_expr, na.rm = TRUE) +
        (max(query_expr, na.rm = TRUE) - min(query_expr, na.rm = TRUE)) * 0.05,
      .groups = "drop"
    ) |>
    dplyr::left_join(correlation_data, by = "celltype") |>
    dplyr::mutate(label = paste0("r² = ", r_squared))

  # Create faceted plot
  p <- ggplot2::ggplot(combined_data) +
    ggplot2::aes(y = ref_expr, x = query_expr) +
    ggplot2::geom_point(aes(color = marker),
               size = point_size,
               alpha = 0.8) +
    ggplot2::scale_color_manual(values = marker_colors) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black",
                linetype = "dashed", size = 0.5) +
    ggplot2::facet_wrap(~ celltype, ncol = ncol, scales = "free") +
    ggplot2::labs(
      y = paste(ref_name, "expression"),
      x = paste(query_name, "expression")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = ifelse(length(markers) > 10, "right", "none"),
      strip.text = ggplot2::element_text(size = 11, face = "bold"),
      axis.title = ggplot2::element_text(size = 9),
      axis.text = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, size = 0.5)
    )

  # Add correlation annotations
  p <- p +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      hjust = 1, vjust = 0,
      size = 3.5,
      fontface = "bold",
      inherit.aes = FALSE)

  # Add marker labels if there are few markers
  if (length(markers) <= 10) {
    check_package("ggrepel")
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = marker),
      size = 2.5,
      max.overlaps = 20,
      segment.size = 0.2,
      segment.alpha = 0.5
    )
  }

  return(p)
}



#' Extract and plot per-class marker importance
#'
#' Calculates marker importance for each cell type by analyzing how each marker
#' contributes to distinguishing that cell type from others using one-vs-all approach.
#'
#' @param reference Training data used to build the model
#' @param col Name of the cell type column
#' @param markers Vector of marker names
#' @param top_n Number of top variables to show per class (default 10)
#'
#' @return ggplot2 heatmap showing variable importance by cell type
#' @export
#'
plot_marker_importance_cell <- function(reference, col = "celltype", markers = get_markers(reference), top_n = length(markers)) {

  cell_types <- unique(reference[[col]])
  importance_matrix <- matrix(0, nrow = length(markers), ncol = length(cell_types))
  rownames(importance_matrix) <- markers
  colnames(importance_matrix) <- cell_types

  # Calculate importance for each cell type using one-vs-all
  for (ct in cell_types) {
    binary_labels <- ifelse(reference[[col]] == ct, ct, "Other")

    if (length(unique(binary_labels)) > 1) {  # Only if we have both classes
      binary_model <- ranger::ranger(
        y = as.factor(binary_labels),
        x = reference[, markers],
        importance = "impurity",
        num.trees = 100
      )
      importance_matrix[, ct] <- binary_model$variable.importance
    }
  }

  # Convert to long format for plotting
  importance_df <- importance_matrix |>
    as.data.frame() |>
    tibble::rownames_to_column("marker") |>
    tidyr::pivot_longer(cols = -marker, names_to = "celltype", values_to = "importance")

  # Get top markers overall
  top_markers <- importance_df |>
    dplyr::group_by(marker) |>
    dplyr::summarise(total_importance = sum(importance), .groups = "drop") |>
    dplyr::arrange(desc(total_importance)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(marker)

  # Filter to top markers
  plot_data <- importance_df |>
    dplyr::filter(marker %in% top_markers) |>
    dplyr::mutate(marker = factor(marker, levels = rev(top_markers)))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = celltype, y = marker, fill = importance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                                  name = "Importance") +
    ggplot2::labs(
      title = "Variable Importance by Cell Type",
      x = "Cell Type",
      y = "Markers"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )

  return(p)
}


#' Plot minimal marker sets for cell type discrimination
#'
#' Identifies the smallest set of markers needed to distinguish each cell type
#' from all others using recursive feature elimination approach.
#'
#' @param reference Training data
#' @param col Name of cell type column
#' @param markers Vector of marker names
#' @param n_markers Number of markers to consider per cell type (default 5)
#'
#' @return ggplot2 object showing minimal marker sets
#' @export
#'
plot_top_marker_importance <- function(
    reference,
    col = "celltype",
    markers = NULL,
    n_markers = 5) {

  if (is.null(markers)) markers <- get_markers(reference)

  cell_types <- unique(reference[[col]])
  minimal_sets <- list()

  for (ct in cell_types) {
    # Create binary classification problem
    binary_labels <- ifelse(reference[[col]] == ct, "Target", "Other")

    if (length(unique(binary_labels)) > 1) {
      # Build model and get importance
      binary_model <- ranger::ranger(
        y = as.factor(binary_labels),
        x = reference[, markers],
        importance = "impurity",
        num.trees = 100
      )

      # Get top markers for this cell type
      importance_scores <- binary_model$variable.importance
      top_markers <- names(sort(importance_scores, decreasing = TRUE))[1:n_markers]

      # Normalize importance scores for this cell type (0-1 scale)
      # normalized_importance <- importance_scores / max(importance_scores)
      # normalized_importance <- scale(importance_scores)[,1]
      normalized_importance <- (importance_scores - min(importance_scores)) /
        (max(importance_scores) - min(importance_scores))

      minimal_sets[[ct]] <- data.frame(
        celltype = ct,
        marker = top_markers,
        importance = importance_scores[top_markers],
        normalized_importance = normalized_importance[top_markers],
        rank = 1:n_markers
      )
    }
  }

  # Combine results
  minimal_df <- do.call(rbind, minimal_sets)

  p <- ggplot2::ggplot(minimal_df, ggplot2::aes(x = rank, y = celltype, fill = normalized_importance)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = marker), size = 2.5, color = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue",
                                 name = "Normalized\nImportance",
                                 labels = scales::percent_format()) +
    ggplot2::scale_x_continuous(breaks = 1:n_markers, labels = paste("Marker", 1:n_markers)) +
    ggplot2::labs(
      title = "Minimal Marker Sets for Cell Type Discrimination",
      subtitle = "Top markers needed to distinguish each cell type from others",
      x = "Marker Priority",
      y = "Cell Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )

  return(p)
}

#' Plot marker importance from Random Forest model
#'
#' Creates a horizontal bar plot showing the importance scores of markers used
#' in the Random Forest model. Markers are ordered by importance score.
#'
#' @param model Random Forest model object (from ranger package)
#' @param top_n Integer specifying number of top markers to display. If NULL,
#' shows all variables (default NULL)
#' @param title Character string for plot title (default "Marker Importance")
#' @param color Single color for bars or named vector of colors by marker
#'
#' @return ggplot2 object showing marker importance scores
#' @export
#'
plot_marker_importance <- function(model, top_n = NULL, title = "Marker Importance", color = "steelblue") {

  # Extract variable importance
  importance_data <- tibble::tibble(
    variable = names(model$variable.importance),
    importance = model$variable.importance
  ) |>
    dplyr::arrange(dplyr::desc(importance))

  # Filter to top_n if specified
  if (!is.null(top_n)) {
    importance_data <- importance_data |>
      dplyr::slice_head(n = top_n)
  }

  # Reorder factor levels for plotting
  importance_data$variable <- factor(importance_data$variable,
                                     levels = importance_data$variable)

  p <- ggplot2::ggplot(importance_data, ggplot2::aes(x = importance, y = variable)) +
    ggplot2::geom_col(fill = color, alpha = 0.8) +
    ggplot2::labs(
      title = title,
      x = "Importance Score",
      y = "Variables"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}



#' Analyze marker usage frequency across trees
#'
#' Shows how often each marker is used for splitting across all trees,
#' revealing which variables are consistently important vs occasionally useful.
#'
#' @param model Random Forest model object
#' @param top_n Number of top variables to display (default 15)
#'
#' @return ggplot2 object showing variable usage frequency
#' @export
#'
plot_marker_usage <- function(model, top_n = 15) {

  forest <- model$forest
  marker_names <- forest$independent.variable.names

  # Count variable usage across all trees
  all_splits <- unlist(forest$split.varIDs)
  variable_counts <- table(all_splits)

  # Convert to marker names
  usage_data <- data.frame(
    marker_index = as.numeric(names(variable_counts)),
    usage_count = as.numeric(variable_counts)
  ) |>
    dplyr::mutate(
      marker_name = marker_names[marker_index + 1],  # R is 1-indexed, split.varIDs are 0-indexed
      usage_percentage = usage_count / sum(usage_count) * 100
    ) |>
    dplyr::arrange(dplyr::desc(usage_count)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(marker_name = factor(marker_name, levels = rev(marker_name)))

  p <- ggplot2::ggplot(usage_data, ggplot2::aes(x = usage_percentage, y = marker_name)) +
    ggplot2::geom_col(fill = "darkgreen", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(usage_percentage, 1), "%")),
                       hjust = -0.1, size = 3) +
    ggplot2::labs(
      title = "Marker Usage",
      subtitle = "How often each marker is used for splitting decisions",
      x = "Percentage of Total Splits",
      y = "Markers"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    ) +
    ggplot2::xlim(0, max(usage_data$usage_percentage) * 1.1)

  return(p)
}

