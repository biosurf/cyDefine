library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)


# Define colors to match your previous plots
method_colors <- c(
  "cyDefine" = "#A8D8A8",
  "CyAnno" = "#FA8072",
  "Spectre" = "#87CEEB",
  "CyTOF Linear\nClassifier" = "#F4A460"
)


# Runtimes ----

df_runtimes <- readRDS("results/runtimes.rds")

# Runtime bar plot by dataset and method
fig_runtime <- df_runtimes %>%
  ggplot(aes(x = dataset, y = runtime, fill = method)) +
  geom_col(position = position_dodge(width = 1), alpha = 0.8) +
  geom_text(aes(
    label = paste0(round(runtime, 1), "s")),
    vjust = -0.5,
    position = position_dodge(width = 1)) +
  facet_wrap(~unassigned_label) +
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

df_f1 <- readRDS("results/F1.rds")

df_f1$method[df_f1$method == "CLC"] <- "CyTOF Linear\nClassifier"

## Figure 3ab - Heatmap ----
# Calculate average F1 scores by dataset and method
df_heatmap <- df_f1 |>
  mutate(method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))) |>
  group_by(dataset, method, unassigned) |>
  summarise(
    avg_f1 = median(f1_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(dataset, unassigned) |>
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
create_heatmap <- function(data, condition_filter, title) {
  data |>
    dplyr::filter(unassigned == condition_filter) |>
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
  df_heatmap,
  "w_unassigned",
  "Average F1 score\n(including unassigned cells)"
)

# Create heatmap for excluding unassigned cells
fig3b <- create_heatmap(
  df_heatmap,
  "wo_unassigned",
  "Average F1 score\n(excluding unassigned cells)"
) +
  ylab("")

# Create heatmap for excluding unassigned cells
fig3z <- create_heatmap(
  df_heatmap,
  "rm_unassigned",
  "Average F1 score\n(disregarding unassigned cells)"
) +
  ylab("")


# Combine heatmaps side by side
library(patchwork)
fig3ab <- fig3a + fig3b + fig3z +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "a")

# print(fig3ab)

# Save plots
ggsave("figs/Figure3ab.png", fig3ab, width = 10, height = 5, dpi = 300)


## Figure 3c Boxplot ----
fig3c <- df_f1 |>
  mutate(method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))) |>
  dplyr::filter(stringr::str_detect(file, "w_unassigned")) |>
  ggplot(aes(x = dataset, y = f1_score, fill = method)) +
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


# Figure 3
fig3 <- fig3ab / fig3c / fig_runtime  +
  plot_annotation(tag_levels = "a")
# print(fig3)
ggsave(plot = fig3, filename = "figs/Figure3.png", width = 15, height = 17)




