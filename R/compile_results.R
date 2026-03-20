library(dplyr)
library(stringr)

# F1 ----
# Load all F1 score files recursively from data folder
f1_files <- list.files("data", pattern = "_F1\\.rds$", recursive = TRUE, full.names = TRUE)
f1_files <- f1_files[!stringr::str_detect(f1_files, "universal")]
# f1_files <- f1_files[stringr::str_detect(f1_files, "unassigned")]
# Create a list to store results
f1_results <- list()

# Load each file and extract information from filename
for (file in f1_files) {
	  # Read the RDS file
	  f1_data <- readRDS(file)
  if (length(f1_data) == 0) {
	      message(file)
      next}

    # Extract basename without extension
    basename <- tools::file_path_sans_ext(basename(file))

    # Parse filename to extract dataset, method, and other info
    # Assuming format: dataset_method_suffix_F1.rds
    parts <- strsplit(basename, "_")[[1]]

      # Remove "F1" from the end
      parts <- parts[-length(parts)]
      # Extract dataset (first part) and method (combine remaining parts)
      dataset <- parts[1]
        method <- parts[2]
        unassigned <- paste(parts[-1:-2], collapse = "_")

	  # Store results with metadata
	  f1_results[[basename]] <- data.frame(
					           celltype = names(f1_data),
						       f1_score = f1_data,
						       dataset = dataset,
						           method = method,
						           unassigned = unassigned,
							       file = basename, row.names = NULL
							     )
}

# Combine all results into one dataframe
df_f1 <- dplyr::bind_rows(f1_results)

# Clean up method names to match the plot
df <- df_f1 |>
  mutate(
	     method = case_when(
				      str_detect(method, "CLC") ~ "CyTOF Linear\nClassifier",
				            TRUE ~ method
				          ),
	     method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))
	       )

saveRDS(df_f1, "results/F1.rds")



# Runtime ----
#
# Load all F1 score files recursively from data folder
runtime_files <- list.files("data", pattern = "_runtime\\.txt$", recursive = TRUE, full.names = TRUE)
runtime_files <- runtime_files[!stringr::str_detect(runtime_files, "universal")]
# Create a list to store results
runtimes <- list()

# Load each file and extract information from filename
for (file in runtime_files) {
	  # Read the RDS file
	  runtime <- readLines(file, warn = FALSE)
  runtime <- as.numeric(stringr::str_remove(runtime, ".*:\\s"))

    if (is.na(runtime)) {
	        message(file)
      next}

    # Extract basename without extension
    basename <- tools::file_path_sans_ext(basename(file))

      # Parse filename to extract dataset, method, and other info
      # Assuming format: dataset_method_suffix_F1.rds
      parts <- strsplit(basename, "_")[[1]]

      # Remove "F1" from the end
      parts <- parts[-length(parts)]
        # Extract dataset (first part) and method (combine remaining parts)
        dataset <- parts[1]
        method <- parts[2]
	  unassigned <- paste(parts[-1:-2], collapse = "_")

	  # Store results with metadata
	  runtimes[[basename]] <- data.frame(
					         runtime = runtime,
						     dataset = dataset,
						     method = method,
						         unassigned = unassigned,
						         file = basename, row.names = NULL
							   )
}

# Combine all results into one dataframe
df_runtimes <- dplyr::bind_rows(runtimes)

# Clean up method names to match the plot
df_runtimes <- df_runtimes |>
  mutate(
	     method = case_when(
				      str_detect(method, "CLC") ~ "CyTOF Linear\nClassifier",
				            TRUE ~ method
				          ),
	     method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine")),
	         unassigned_label = case_when(
					            unassigned == "w_unassigned" ~ "Including unassigned",
						          unassigned == "wo_unassigned" ~ "Excluding unassigned",
						          TRUE ~ unassigned
							      ),
	     unassigned_label = factor(unassigned_label, levels = c("Including unassigned", "Excluding unassigned"))
	       )

saveRDS(df_runtimes, "results/runtimes.rds")
message("Results compiled.")
