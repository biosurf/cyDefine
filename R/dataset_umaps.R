library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cyCombine)
library(cyDefine)


read_dataset <- function(dataset, platform = c("fcm", "cytof")) {
  platform <- match.arg(platform)
  path <- file.path("../ob-flow-datasets/prepared", platform, dataset)
  files <- list.files(path, pattern = "\\.csv\\.zst$", full.names = TRUE)
  df_list <- lapply(files, function(file) {
    tmp_csv <- tempfile(fileext = ".csv")
    # ret <- system(paste("/opt/homebrew/Caskroom/miniconda/base/bin/zstd -d", file, "-o", tmp_csv))
    ret <- system(paste("/services/tools/mamba-org/mamba/0.24.0/bin/zstd -d", file, "-o", tmp_csv))
    df <- readr::read_csv(tmp_csv, show_col_types = FALSE)
    df$filename <- basename(file)
    unlink(tmp_csv)
    return(df)
  })
  df <- do.call(dplyr::bind_rows, df_list)
  return(df)

}

p_theme <- function() {
  ggplot2::theme(title = element_text(size = 18), legend.text = element_text(size = 11), legend.justification = "left")
}

# ChikVirus ----

dataset <- "FR-FCM-Z238/ChikVirusPBMC_cytof"
platform <- "cytof"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("Xe131Di", "Cs133Di", "Ba138Di", "Ce140Di", "Pr141Di", "Ce142Di", "Pr141Di", "Lu176Di", "Os189Di", "Pt195Di", "CD45_acute", "CD45_conv", "filename") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T)
df$label[df$label == "Ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_chikV <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("ChikVirusPBMC", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_chikV <- p_chikV & p_theme()
# p_chikV
ggsave(p_chikV, filename = "figs/umap_chikV.png", width = 10, height = 6)

# Bodenmiller ----

dataset <- "BodenmillerXL/PBMC_cytof"
platform <- "cytof"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T)
df$label[df$label == "ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Bodenmiller <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("BodenmillerXL", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Bodenmiller <- p_Bodenmiller & p_theme()
# p_Bodenmiller
ggsave(p_Bodenmiller, filename = "figs/umap_Bodenmiller.png", width = 10, height = 6)


# Samusik ----

dataset <- "Samusik/MouseBoneMarrow_cytof"
platform <- "cytof"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T)
df$label[df$label == "ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Samusik <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("Samusik", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Samusik <- p_Samusik & p_theme()
# p_Samusik
ggsave(p_Samusik, filename = "figs/umap_Samusik.png", width = 10, height = 6)



# Levine32 ----

dataset <- "Levine/HumanBoneMarrow_cytof"
platform <- "cytof"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T)
df$label[df$label == "ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Levine32 <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("Levine32", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Levine32 <- p_Levine32 & p_theme()
# p_Levine32
ggsave(p_Levine32, filename = "figs/umap_Levine32.png", width = 10, height = 6)


# Z2KP - healthy ----

dataset <- "FR-FCM-Z2KP-healthy/PBMC_fcm"
platform <- "fcm"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename", "Time", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W" ) ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T, cofactor = 150)
df$label[df$label == "Ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Z2KPH <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("Z2KP - Healthy", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Z2KPH <- p_Z2KPH & p_theme()
# p_Z2KPH
ggsave(p_Z2KPH, filename = "figs/umap_Z2KPH.png", width = 10, height = 6)

# Z2KP - Covid ----

dataset <- "FR-FCM-Z2KP-covid/covidPBMC_fcm"
platform <- "fcm"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}


markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename", "Time", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W" ) ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T, cofactor = 150)
df$label[df$label == "Ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Z2KPC <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("Z2KP - Covid", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Z2KPC <- p_Z2KPC & p_theme()
# p_Z2KPC
ggsave(p_Z2KPC, filename = "figs/umap_Z2KPC.png", width = 10, height = 6)



# FlowCyt ----

dataset <- "FlowCyt/HumanBoneMarrow_fcm"
platform <- "fcm"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}

df <- tidyr::drop_na(df)

markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename", "FS INT", "SS INT") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T, cofactor = 150)
df$label[df$label == "ungated"] <- "unassigned"

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_FlowCyt <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("FlowCyt", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_FlowCyt <- p_FlowCyt & p_theme()
# p_FlowCyt
ggsave(p_FlowCyt, filename = "figs/umap_FlowCyt.png", width = 10, height = 6)


# Z3YR ----

dataset <- "FR-FCM-Z3YR/StimBlood_cytof"
platform <- "cytof"
raw_file <- file.path("../ob-flow-datasets/prepared", platform, dataset, "raw.RDS")

if (file.exists(raw_file)) {
  df <- readRDS(raw_file)
} else {
  df <- read_dataset(dataset = dataset, platform = platform)
  saveRDS(df, file = raw_file)
}

# df <- tidyr::drop_na(df)

markers <- cyCombine::get_markers(df)
markers <- markers[!markers %in% c("filename") ]

df <- cyCombine::transform_asinh(df, markers = markers, .keep = T, cofactor = 5)
df$label[df$label == "unlabeled"] <- "unassigned"

df$label <- df$label |> stringr::str_remove_all(".*/")

celltype_colors <- cyDefine::get_distinct_colors(sort(unique(df$label)), add_unassigned = T)

p_Z3YR <- cyDefine::plot_umap(
  df[df$label != "unassigned", c(markers, "label")],
  df[,c(markers, "label")],
  colors = celltype_colors,
  col = "label",
  title = c("Z3YR", ""),
  build_umap_on = "reference",
  down_sample = T,
  sample_n = 20000)
p_Z3YR <- p_Z3YR & p_theme()
# p_Z3YR
ggsave(p_Z3YR, filename = "figs/umap_Z3YR.png", width = 10, height = 6)



# Full UMAP ----
umap_full <- p_chikV / p_Bodenmiller / p_Samusik / p_Levine32 / p_Z3YR / p_Z2KPH / p_Z2KPC / p_FlowCyt


ggsave(umap_full, filename = "figs/umap_full.png", width = 20, height = 6*9)



ggsave(p_chikV / p_Bodenmiller / p_Samusik, filename = "figs/umap_part1.png", width = 20, height = 6*3)
ggsave(p_Samusik / p_Levine32 / p_Z3YR, filename = "figs/umap_part2.png", width = 20, height = 6*3)
ggsave(p_Z2KPH / p_Z2KPC / p_FlowCyt, filename = "figs/umap_part3.png", width = 20, height = 6*3)
