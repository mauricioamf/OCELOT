# lasso_parallel.R

library(readr)
library(glmnet)
library(doParallel)
library(foreach)
library(doRNG)

set.seed(123)

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Starting parallel lasso\n")

# Setup path
path <- "~/OCELOT/manuscript/code/"

# Load data
cat("Loading expression data...\n")
expression_data <- read.csv(
  file.path(path, "Sce_count_matrix.csv"),
  row.names = 1, header = TRUE, sep = "\t"
)

expression_data[expression_data == "-" | is.na(expression_data)] <- 0
rownames(expression_data) <- trimws(rownames(expression_data))
cat("Expression data loaded. Dimensions:", dim(expression_data), "\n")

# Load GRN
cat("Loading TF list...\n")

TF_TG_GS <- read.delim(file.path(path, "Sce_GRN_yeastract.csv"), header = TRUE, sep = "\t")
TF_TG <- data.frame(
  TF = trimws(TF_TG_GS$TF),
  Target = trimws(TF_TG_GS$Target)
)

# Identify TFs and non-TFs
tf_genes <- unique(TF_TG$TF)
tfs <- intersect(tf_genes, rownames(expression_data))
non_tf_genes <- setdiff(rownames(expression_data), tfs)
cat(length(tfs), "TFs and", length(non_tf_genes), "target genes identified\n")

# Prepare TF matrix
tf_matrix <- as.matrix(expression_data[tfs, , drop = FALSE])
tf_matrix <- t(tf_matrix)
colnames(tf_matrix) <- make.names(colnames(tf_matrix), unique = TRUE)
cat("TF matrix prepared. Dimensions:", dim(tf_matrix), "\n")

# Setup checkpoint directory
checkpoint_dir <- file.path(path, "checkpoint_lasso")
if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir)
}

processed_genes <- gsub("\\.rds$", "", list.files(checkpoint_dir, pattern = "\\.rds$"))
remaining_genes <- setdiff(non_tf_genes, processed_genes)
cat("Remaining genes to process:", length(remaining_genes), "\n")

# Setup parallel backend
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores)) n_cores <- 1
cat("Initializing parallel backend with", n_cores, "cores\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary variables to workers
clusterExport(cl, varlist = c("expression_data", "tf_matrix"), envir = environment())

# Verify workers can access data
test_result <- clusterEvalQ(cl, {
  exists("expression_data") && exists("tf_matrix") && 
    !is.null(dim(expression_data)) && !is.null(dim(tf_matrix))
})
cat("Data accessibility test:", ifelse(all(unlist(test_result)), "PASSED", "FAILED"), "\n")

# Process genes in parallel
cat("Starting parallel processing of", length(non_tf_genes), "genes\n")
results <- foreach(gene = non_tf_genes, 
                   .packages = "glmnet", 
                   .errorhandling = "pass", .verbose = TRUE) %dopar% {
                     tryCatch({
                       cat("Processing gene:", gene, "\n")
                       response <- as.numeric(expression_data[gene, ])
                       
                       # Fit model
                       cv_model <- cv.glmnet(
                         x = tf_matrix,
                         y = response,
                         family = "gaussian",
                         nfolds = 10,
                         alpha = 1,
                         intercept = FALSE,
                         type.measure = "mse"
                       )
                       
                       # Extract coefficients
                       coef_vec <- as.vector(coef(cv_model, s = "lambda.min"))[-1]
                       
                       list(gene = gene, status = "success", coefficients = coef_vec)
                     }, error = function(e) {
                       msg <- paste("Error processing gene", gene, ":", e$message)
                       cat(msg, "\n")
                       list(gene = gene, status = "error", message = msg)
                     })
                   }

# Stop cluster and clean up
cat("Parallel processing completed. Stopping cluster...\n")
stopCluster(cl)

# Process results
cat("Processing results...\n")
success_count <- 0
skipped_count <- 0
error_count <- 0
result_list <- list()

for (res in results) {
  if (res$status == "success") {
    success_count <- success_count + 1
    result_list[[res$gene]] <- res$coefficients
  } else if (res$status == "skipped") {
    skipped_count <- skipped_count + 1
  } else {
    error_count <- error_count + 1
  }
}

cat("Results summary:\n")
cat("- Successfully processed:", success_count, "\n")
cat("- Skipped (constant response):", skipped_count, "\n")
cat("- Errors:", error_count, "\n")

if (success_count > 0) {
  # Convert to matrix
  coef_matrix <- do.call(rbind, result_list)
  rownames(coef_matrix) <- names(result_list)
  colnames(coef_matrix) <- colnames(tf_matrix)
  
  # Filter zero rows
  non_zero_rows <- rowSums(coef_matrix != 0) > 0
  coef_matrix <- coef_matrix[non_zero_rows, , drop = FALSE]
  
  # Save output
  output_file <- file.path(path, "Sce_GRN_lasso_raw.csv")
  cat("Saving results to", output_file, "\n")
  write.csv(coef_matrix, output_file)
} else {
  cat("No valid results to save\n")
}

# Final cleanup
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Script completed\n")
