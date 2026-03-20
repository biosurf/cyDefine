library(argparse)
library(dplyr)
library(cyCombine)
# library(cyDefine)
devtools::load_all("../cyDefine")


# Define argument parser
parser <- ArgumentParser(description="Run cyDefine")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
# parser$add_argument("--data", dest="data", type="character", help="input file #1")
parser$add_argument("--unassigned", dest="unassigned", type="logical", help="train on unassigned?", default = FALSE)
parser$add_argument("--unassigned_name", dest="unassigned_name", type="character", help="unassigned name", default = NULL)
parser$add_argument("--adapt_reference", dest="adapt_reference", type="logical", help="adapt_reference?", default = FALSE)
parser$add_argument("--batch_correct", dest="batch_correct", type="logical", help="batch_correct?", default = FALSE)
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
adapt_reference <- opt$adapt_reference
batch_correct <- opt$batch_correct
mtry <- opt$mtry
seed <- opt$seed
threads <- opt$threads
unassigned <- opt$unassigned
if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]


if (is.null(mtry)) mtry <- floor(length(input$markers)/3)
# metadata <- opt$metadata

t <- system.time({
message("Classifying ", name, " using cyDefine")

## Run cyDefine
classified <- cyDefine(
  reference = input$reference,
  query = input$query,
  markers = input$markers,
  adapt_reference = adapt_reference,
  using_pbmc = FALSE,
  batch_correct = batch_correct,
  identify_unassigned = TRUE,
  identify_type = "probability",
  probability_threshold = 0.8,
  # MAD_factor = 2.5,
  train_on_unassigned = unassigned,
  unassigned_name = input$unassigned_name,
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

# Store results
suffix <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
saveRDS(output, paste0(dir, "/", name, "_cyDefine", suffix, ".rds"))

cat("Runtime (elapsed): ", t["elapsed"], file = paste0(dir, "/", name, "_cyDefine", suffix, "_runtime.txt"))


message("Runtime: ", t["elapsed"])
