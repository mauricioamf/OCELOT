# Load required libraries
library(caret)
library(reshape2)
library(data.table)
library(doParallel)
library(foreach)

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Starting parallel linear validation\n")

# Define the base path for data files
path <- "~/OCELOT/manuscript/code/"

# Load lasso continuous results
cat("Loading LASSO linear results...\n")
lasso_res_df <- read.csv(file.path(path, "Sce_GRN_lasso_raw.csv"), row.names = 1)

# Load filtered expression data
cat("Loading expression data...\n")
expression_data <- read.csv(file.path(path, "Sce_count_matrix.csv"), sep = "\t", row.names = 1)
expression_data[expression_data == "-" | is.na(expression_data)] <- 0  # Replace missing or invalid values with 0
rownames(expression_data) <- trimws(rownames(expression_data))

# Load gold standard TF-TG relationships
TF_TG_GS <- read.delim(file.path(path, "Sce_GRN_yeastract.csv"), header = TRUE, sep = "\t")
TF_TG <- data.frame(
  TF = trimws(TF_TG_GS$TF),
  Target = trimws(TF_TG_GS$Target)
)

# Identify TF and non-TF genes
tf_genes <- unique(TF_TG$TF)  
tfs <- intersect(tf_genes, rownames(expression_data))  
tfs <- trimws(tfs)  
non_tf_genes <- setdiff(rownames(expression_data), tfs)  

# Prepare TF expression data
tf_data <- t(expression_data[tfs, , drop = FALSE])  
tf_data <- as.data.frame(tf_data)
tf_data[] <- lapply(tf_data, as.numeric)  

# Define continuous validation thresholds
r2_threshold <- 0.7          # e.g., TFs must explain at least 70% of the variance
correlation_threshold <- 0.5 # e.g., Pearson correlation > 0.5

# Ensure lasso results are available
if (nrow(lasso_res_df) == 0) {
  stop("The lasso_res_df dataframe is empty or not properly initialized.\n")
}

# Setup parallel backend
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores)) n_cores <- 7
cat("Initializing parallel backend with", n_cores, "cores\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary variables to workers
clusterExport(cl, varlist = c("expression_data", "tf_data", "lasso_res_df"), envir = environment())

genes_to_validate <- rownames(lasso_res_df)
cat("Starting parallel validation of", length(genes_to_validate), "genes\n")

# Process genes in parallel
results <- foreach(gene = genes_to_validate, 
                   .packages = c("stats"), 
                   .errorhandling = "pass") %dopar% {
                     tryCatch({
                       coefficients <- as.numeric(lasso_res_df[gene, ])
                       selected_predictors <- which(coefficients != 0)  # Select predictors with non-zero coefficients
                       
                       if (length(selected_predictors) == 0) {
                         return(list(gene = gene, status = "skipped", reason = "no predictors"))
                       }
                       
                       selected_tf_data <- tf_data[, selected_predictors, drop = FALSE]
                       predictors <- colnames(tf_data)[selected_predictors]
                       
                       # Use continuous response data (No binarization)
                       response <- as.numeric(expression_data[gene, ])
                       
                       # Skip if response has zero variance (will cause regression errors)
                       if (var(response) == 0) {
                         return(list(gene = gene, status = "skipped", reason = "zero variance"))
                       }
                       
                       # Split data into training and testing sets
                       set.seed(123) # Maintains exact same test/train split for every gene
                       train_index <- sample(1:length(response), size = 0.7 * length(response))
                       train_x <- selected_tf_data[train_index, , drop = FALSE]
                       test_x <- selected_tf_data[-train_index, , drop = FALSE]
                       train_y <- response[train_index]
                       test_y <- response[-train_index]
                       
                       # Train standard linear regression model (Ordinary Least Squares)
                       model <- tryCatch(
                         lm(train_y ~ ., data = data.frame(train_x, train_y)),
                         error = function(e) return(NULL)
                       )
                       
                       if (is.null(model)) {
                         return(list(gene = gene, status = "skipped", reason = "lm error"))
                       }
                       
                       # Predict continuous values on test data
                       pred_values <- predict(model, newdata = data.frame(test_x))
                       
                       # Calculate R-squared on test set
                       ss_res <- sum((test_y - pred_values)^2)
                       ss_tot <- sum((test_y - mean(test_y))^2)
                       r_squared <- ifelse(ss_tot > 0, 1 - (ss_res / ss_tot), 0)
                       
                       # Calculate Pearson Correlation
                       correlation <- cor(test_y, pred_values, use = "complete.obs", method = "pearson")
                       if (is.na(correlation)) correlation <- 0
                       
                       # Create gene's result dataframe
                       gene_df <- data.frame(
                         Gene = rep(gene, length(predictors)),
                         Predictor = predictors,
                         Coefficient = coefficients[selected_predictors],
                         R_Squared = rep(r_squared, length(predictors)),
                         Correlation = rep(correlation, length(predictors)),
                         stringsAsFactors = FALSE
                       )
                       
                       list(gene = gene, status = "success", data = gene_df)
                       
                     }, error = function(e) {
                       list(gene = gene, status = "error", message = e$message)
                     })
                   }

# Stop cluster
cat("Parallel processing completed. Stopping cluster...\n")
stopCluster(cl)

# Process results
cat("Processing results...\n")
success_count <- 0
skipped_count <- 0
error_count <- 0
valid_dfs <- list()

for (res in results) {
  # Check if result is a proper list (sometimes foreach returns the error directly if .errorhandling="pass")
  if (is.list(res) && !is.null(res$status)) {
    if (res$status == "success") {
      success_count <- success_count + 1
      valid_dfs[[length(valid_dfs) + 1]] <- res$data
    } else if (res$status == "skipped") {
      skipped_count <- skipped_count + 1
    } else {
      error_count <- error_count + 1
      cat("Error in gene", res$gene, ":", res$message, "\n")
    }
  } else {
    error_count <- error_count + 1
  }
}

cat("Results summary:\n")
cat("- Successfully evaluated:", success_count, "\n")
cat("- Skipped (No predictors/Zero variance/LM Error):", skipped_count, "\n")
cat("- Errors:", error_count, "\n")

if (length(valid_dfs) > 0) {
  # Combine all successful dataframes
  validation_results <- do.call(rbind, valid_dfs)
  
  # Filter results by R-squared and Correlation thresholds
  filtered_result_df <- validation_results[
    validation_results$R_Squared >= r2_threshold & validation_results$Correlation >= correlation_threshold, 
  ]
  
  # Save filtered validation results
  write.csv(filtered_result_df, file.path(path, "Sce_GRN_lasso_values.csv"), row.names = FALSE)
  cat("Saved long-format validation values.\n")
  
  # Create coefficient matrix from filtered results
  if(nrow(filtered_result_df) > 0) {
    coef_matrix <- dcast(setDT(filtered_result_df), Gene ~ Predictor, value.var = "Coefficient", fill = 0)
    
    if(ncol(coef_matrix) > 1) {
      coef_matrix_mat <- as.matrix(coef_matrix[,-1])  # Convert to matrix
      rownames(coef_matrix_mat) <- coef_matrix$Gene
      
      # Save validated coefficient matrix
      write.csv(coef_matrix_mat, file.path(path, "Sce_GRN_lasso.csv"), row.names = TRUE)
      cat("Saved validated wide-format matrix.\n")
    }
  } else {
    cat("No edges passed the strict validation thresholds (R2 >=", r2_threshold, "and Cor >=", correlation_threshold, ").\n")
  }
  
} else {
  cat("No genes successfully generated validation metrics.\n")
}

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Script completed\n")