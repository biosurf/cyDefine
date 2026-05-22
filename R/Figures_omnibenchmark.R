library(jsonlite)
library(tibble)
library(dplyr)
library(tidyr)
library(pbmcapply)

performance_files <- list.files(
  path = "../ob-pipeline-cytof/out/",
  pattern = "data_import.flow_metrics.json.gz",
  recursive = TRUE,
  full.names = TRUE
)

results_list <- pbmcapply::pbmclapply(performance_files, function(file) {
  con <- gzfile(file, "rb")
  on.exit(close(con))

  # 1. Read raw content as text
  raw_text <- paste(readLines(con, warn = FALSE), collapse = "\n")

  # 2. Replace invalid 'NaN' with valid JSON 'null'
  # We use regex to ensure we only replace the literal NaN, not parts of other words
  clean_text <- gsub(":\\s*NaN", ": null", raw_text)

  # 3. Parse the cleaned text
  json_data <- fromJSON(clean_text, simplifyVector = TRUE)

  # json_data <- fromJSON(con, simplifyVector = TRUE)

  rows_list <- list()
  row_idx <- 1

  # Check if results exist
  if (is.null(json_data$results)) return(NULL)

  # 1. Dynamically iterate over all runs (e.g., run0, run1, etc.)
  for (run_id in names(json_data$results)) {
    run_data <- json_data$results[[run_id]]

    if (is.null(run_data$per_sample)) next

    # 2. Iterate over each sample file within the run
    for (sample_name in names(run_data$per_sample)) {
      sample_data <- run_data$per_sample[[sample_name]]

      if (is.null(sample_data$per_population)) next

      # 3. Iterate over each population within the sample
      for (pop_id in names(sample_data$per_population)) {
        pop_data <- sample_data$per_population[[pop_id]]

        # Create a single row tibble for this population
        row_df <- tibble(
          file_source       = basename(file),
          method            = stringr::str_extract(file, "/analysis/(.*)/metrics/", group = 1) |> stringr::str_remove("/.*"),
          dataset           = stringr::str_extract(file, "dataset_name-(.*)_name-data_import.data_raw", group = 1),
          params            = stringr::str_extract(file, "(raw_potential-batches.*)/preprocessing/", group = 1),
          run_id            = run_id,          # Dynamic run ID
          sample_file       = sample_name,     # e.g., "label-10.csv"
          population_id     = as.integer(pop_id),
          population_name   = pop_data$population_name,
          accuracy          = pop_data$accuracy,
          precision         = pop_data$precision,
          recall            = pop_data$recall,
          f1                = pop_data$f1,
          tp                = pop_data$tp,
          fp                = pop_data$fp,
          fn                = pop_data$fn,
          tn                = pop_data$tn,
          scaling_rate      = pop_data$scaling_rate,
          support           = pop_data$support,
          drop_ungated_test = stringr::str_extract(file, pattern = "drop-ungated-test-(true|false)", group = 1) |> as.logical(),
          drop_ungated_training = stringr::str_extract(file, pattern = "drop-ungated-training-(true|false)", group = 1) |> as.logical(),
          n_cells           = pop_data$n_cells
        )
        row_df$method <- stringr::str_replace_all(row_df$method, "cydefine", "cyDefine")

        row_df$runtime <- sample_data$runtime_seconds
        rows_list[[row_idx]] <- row_df
        row_idx <- row_idx + 1
      }
    }
  }

  bind_rows(rows_list)
})


# Combine all results
final_metrics <- bind_rows(Filter(Negate(is.null), results_list))

final_metrics |>
  group_by(dataset, drop_ungated_test, drop_ungated_training) |>
  summarise(median(f1, na.rm = T))


# F1 ----
## Figure 3ab - Heatmap ----
# Calculate average F1 scores by dataset and method
df_heatmap <- final_metrics |>
  # mutate(method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))) |>
  group_by(dataset, method, drop_ungated_test, drop_ungated_training) |>
  summarise(
    avg_f1 = median(f1, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(dataset, drop_ungated_test, drop_ungated_training) |>
  mutate(
    is_max = avg_f1 == max(avg_f1, na.rm = TRUE),
    # Round to 3 decimal places for display
    avg_f1 = round(avg_f1, 3),
    # Create text labels
    label = as.character(avg_f1),
    # Add asterisk for missing data (if needed)
    label = ifelse(is.na(avg_f1), "*", label),
    is_max = ifelse(is.na(is_max), FALSE, is_max),
  )

# Create heatmaps for both conditions
create_heatmap <- function(data, title) {
  data |>
    ggplot(aes(x = method, y = dataset, fill = avg_f1)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = ifelse(.data$is_max, .data$label, ""), fontface = "bold"), size = 3.5) +
    geom_text(aes(label = ifelse(!.data$is_max, .data$label, "")), size = 3.5) +
    # geom_text(aes(label = label),
    #           color = "black",
    #           size = 3.5,
    #           fontface = "bold") +
    # scale_fill_gradient(low = "#e3f5e1", high = "#30a719", na.value = "white",
    scale_fill_gradient2(
      low = "#F5F5F5",
      mid = "lightgreen",
      high = "darkgreen",
      na.value = "white",
      midpoint = 0.4,
      limits = c(0.2, 1.0),
      name = "Average\nF1 score",
      breaks = seq(0.4, 1.0, 0.1)
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Dataset"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm")
    ) +
    coord_fixed()
}

# Create heatmap for including unassigned cells
fig3a <- create_heatmap(
  df_heatmap |> filter(!drop_ungated_test, !drop_ungated_training),
  "Average F1 score\n(including unassigned cells)"
)

# Create heatmap for excluding unassigned cells
fig3b <- create_heatmap(
  df_heatmap |> filter(!drop_ungated_test, drop_ungated_training),
  "Average F1 score\n(excluding unassigned cells)"
) +
  ylab("")

# Create heatmap for excluding unassigned cells
fig3z <- create_heatmap(
  df_heatmap |> filter(drop_ungated_test, drop_ungated_training),
  "Average F1 score\n(disregarding unassigned cells)"
) +
  ylab("")

# Combine heatmaps side by side
library(patchwork)

fig3ab <- fig3a + fig3b + fig3z +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a")

fig3ab

# Save plots
ggsave("figs/Figure3ab.png", fig3ab, width = 10, height = 5, dpi = 300)


## Figure 3c Boxplot ----
fig3c <- final_metrics |>
  # mutate(method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))) |>
  dplyr::filter(!drop_ungated_test, !drop_ungated_training) |>
  ggplot(aes(x = dataset, y = f1, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.size = 0,
               outlier.alpha = 0.7) +
  # facet_wrap(~unassigned) +
  geom_point(shape = 21, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "F1 score per cell type\n(including unassigned cells)",
    x = "Dataset",
    y = "F1 score per cell type",
    fill = ""
  ) +
  ggpubr::theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0, size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.2)
  )

# print(fig3c)
ggsave(plot = fig3c, filename = "figs/Figure3c.png", width = 10, height = 5, dpi = 300)



# performance <- readr::read_tsv("../ob-pipeline-cytof/out/performances.tsv")
# performance |> filter(module == "cydefine")

# Runtimes ----

# df_runtimes <- readRDS("results/runtimes.rds")

# Runtime bar plot by dataset and method
fig_runtime <- final_metrics |>
  group_by(method, dataset, drop_ungated_test, drop_ungated_training) |>
  summarise(runtime = sum(runtime), .groups = "keep") |>
  mutate(Case = case_when(
    !drop_ungated_test & !drop_ungated_training ~ "Include unassigned",
    !drop_ungated_test & drop_ungated_training ~ "Exclude unassigned",
    drop_ungated_test & drop_ungated_training ~ "Disregard unassigned"
  )) |>
  ggplot(aes(x = dataset, y = runtime, fill = method)) +
  geom_col(position = position_dodge(width = 1), alpha = 0.8) +
  geom_text(aes(
    label = paste0(round(runtime, 2), "s")),
    vjust = -0.5,
    position = position_dodge(width = 1)) +
  facet_wrap(~Case) +
  scale_fill_manual(values = method_colors) +
  scale_y_sqrt(labels = scales::label_number(suffix = "s")) +
  labs(
    title = "Runtime Comparison",
    x = "Dataset",
    y = "Runtime (seconds, sqrt scale)",
    fill = "Method"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )



# Save plot
ggsave("figs/fig_runtime.png", fig_runtime, width = 16, height = 8, dpi = 300)


# F1 ----



# Figure 3
fig3 <- fig3ab / fig3c / fig_runtime  +
  plot_annotation(tag_levels = "a")
# print(fig3)
ggsave(plot = fig3, filename = "figs/Figure3.png", width = 15, height = 17)




