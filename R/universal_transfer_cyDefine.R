library(argparse)
library(dplyr)
# devtools::load_all("../tools/cyCombine")
# devtools::load_all("../cyDefine")
library(cyCombine)
library(cyDefine)
library(ggplot2)
library(patchwork)
library(aricode)

# Define argument parser
parser <- ArgumentParser(description="Run cyDefine with universal reference")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
# parser$add_argument("--unassigned", dest="unassigned", type="logical", help="train on unassigned?", default = FALSE)
# parser$add_argument("--unassigned_name", dest="unassigned_name", type="character", help="unassigned name", default = NULL)
parser$add_argument("--mtry", dest="mtry", type="numeric", help="mtry", default = NULL)
parser$add_argument("--seed", dest="seed", type="numeric", help="seed", default = 145)
parser$add_argument("--threads", dest="threads", type="numeric", help="threads", default = 4)
# parser$add_argument("--metadata", dest="metadata", type="character", help="input file #2")

# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
mtry <- opt$mtry
seed <- opt$seed
threads <- opt$threads
# unassigned <- opt$unassigned
# if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]

query <- bind_rows(input$reference, input$query)

if (file.exists("data/pbmc_reference.rds")) {
  reference <- readRDS("data/pbmc_reference.rds")
} else {
  reference <- cyDefine::get_reference("pbmc", store = TRUE, path = "data")
}

if (name == "Levine13") {
  map_specific_from <- NULL
  map_specific_to = NULL
} else if (name == "Levine32") {
  map_specific_from <- c("Flt3")
  map_specific_to = c("CD135")
} else if (name == "Samusik") {
  map_specific_from <- c("CD16_32", "cKit", "TCRb", "MHCII", "CD45.2")
  map_specific_to = c("CD16", "CD117", "abTCR", "HLA-DR", "CD45")
} else if (name == "POISED") {
  map_specific_from <- c("integrin")
  map_specific_to = c("Integrin-7")
} else {
  map_specific_from <- NULL
  map_specific_to = NULL
}


query_mapped <- map_marker_names(
  query,
  query_markers = input$markers,
  ref_markers = cyDefine::pbmc_markers,
  using_pbmc = TRUE,
  map_specific_from = map_specific_from,
  map_specific_to = map_specific_to)


if (is.null(mtry)) mtry <- ceiling(length(input$markers)/2)
# metadata <- opt$metadata

t <- system.time({
message("Classifying ", name, " using cyDefine")

## Run cyDefine
classified <- cyDefine(
  reference = reference,
  query = query_mapped,
  markers = cyCombine::get_markers(query_mapped),
  adapt_reference = TRUE,
  using_pbmc = TRUE,
  batch_correct = TRUE,
  norm_method = "scale",
  xdim = c(1,5),
  ydim = 5,
  identify_unassigned = TRUE,
  train_on_unassigned = FALSE,
  identify_type = "probability",
  probability_threshold = 0.8,
  unassigned_name = input$unassigned_name,
  # MAD_factor = 2.5,
  seed = seed,
  num.threads = threads,
  num.trees = 600,
  mtry = mtry
)
})

if (input$unassigned_name != "unassigned") {
  classified$query$predicted_celltype[classified$query$predicted_celltype == "unassigned"] <- input$unassigned_name
}

output <- data.frame(
  "celltype" = classified$query$celltype,
  "model_prediction" = classified$query$model_prediction,
  "predicted_celltype" = classified$query$predicted_celltype
)


# if (FALSE) {
p1 <- plot_umap(
  classified$reference,
  classified$query,
  markers = cyCombine::get_markers(query_mapped))
# p1
p2 <- plot_umap(
  classified$query,
  col = "celltype",
  title = paste(name, "-", "Celltype"),
  markers = cyCombine::get_markers(query_mapped)) +
  plot_umap(
    classified$query,
    col = "predicted_celltype",
    title = paste(name, "-", "Predicted celltype"),
    markers = cyCombine::get_markers(query_mapped))
# p2

# Store results
ggsave(paste0("figs/", name, "_universal.png"), p1)
ggsave(paste0("figs/", name, "_true-v-predicted.png"), p2, width = 15)
# }
suffix <- NULL#ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
saveRDS(output, paste0(dir, "/", name, "_cyDefine_universal", suffix, ".rds"))

cat("Runtime (elapsed): ", t["elapsed"], file = paste0(dir, "/", name, "_cyDefine_universal", suffix, "_runtime.txt"))


message("Runtime: ", t["elapsed"])

stop()
message("Computing ARI...")
output <- readRDS("data/POISED/POISED_cyDefine_universal.rds")
# output$model_prediction[output$model_prediction == "HSPC / cDC"] <- "cDC"
tb <- table(output$model_prediction, output$celltype)
mapping_vector <- apply(tb, 1, function(row) {
  colnames(tb)[which.max(row)]
})
mapping_vector[["Unknown"]] <- "Unknown"
mapping_vector[["mCD1"]] <- "cDC"
mapping_vector[["mDC2"]] <- "cDC"
mapping_vector[["Classical_mono"]] <- "CD14 Mono"
mapping_vector[["pDCs"]] <- "pDC"
mapping_vector[["NKT"]] <- "Unknown"


cellmapper <- tibble::tribble(
  ~ref, ~query,
  'Unknown', 'Unknown',
  'CD4 Proliferating / CD4 TEM / CD4 TCM', 'Teffector_EM',
  'CD4 Proliferating / CD4 TEM / CD4 TCM', 'Teffector_CM',
  'CD4 Naive', 'naive_T_effectors',
  'CD4 CTL', 'Unknown',
  'CD8 Proliferating / CD8 TCM / CD8 TEM', 'CD8_Effectory_memory',
  'CD8 Proliferating / CD8 TCM / CD8 TEM', 'CD8_Central_memory',
  'CD8 Naive', 'CD8_Naive',
  'Treg', 'Treg_Mem',
  'Treg', 'Treg_Naive',
  'MAIT / gdT', 'gdTCR_Mem',
  'MAIT / gdT', 'gdTCR_Naive',
  'pDC', 'pDCs',
  'cDC', 'mDC2',
  'cDC', 'mCD1',
  'dnT', 'Unknown',
  'Unknown', 'non_classical_mono',
  'Monocyte', 'Classical_mono',
  'Unknown', 'intermediate_mono',
  'Plasmablast', 'Unknown',
  'B intermediate / B memory', 'B_cells_Mem',
  'B naive', 'B_cells_Naive',
  'NK Proliferating / NK', 'NK_Cells',
  'Unknown', 'NKT',
  'ILC', 'Unknown',
  'Unknown', 'CD4_PeaReactive',
  'Unknown', 'CD8_PeaReactive'
)

output$celltype_mapped <- cellmapper$ref[match(output$celltype, cellmapper$query)]

output$predicted_mapped <- output$predicted_celltype
output$predicted_mapped[output$predicted_celltype %in% cellmapper$ref[cellmapper$query == "Unknown"]] <- "Unknown"

output$predicted_mapped[output$predicted_mapped %in% c("CD14 Mono", "CD16 Mono")] <- "Monocyte"
output$predicted_mapped[output$predicted_mapped == 'NK_CD56bright'] <- 'NK Proliferating / NK'
aricode::ARI(output$celltype_mapped , output$predicted_mapped)
cyDefine:::compute_f1(output$celltype_mapped , output$predicted_mapped) |> median()

# output$predicted_celltype[output$predicted_celltype == "HSPC / cDC"] <- "cDC"

output$ct <- mapping_vector[output$celltype]

o <- tibble::as_tibble(list("ct" = output$ct, "predicted_celltype" = output$predicted_celltype))
o <- dplyr::filter(o, predicted_celltype != "Unknown")
aricode::ARI(o$ct , o$predicted_celltype)
cyDefine:::compute_f1(o$ct , o$predicted_celltype) |> median()
