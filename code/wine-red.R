rm(list = ls())
set.seed(2022L)

#######################################################################
# Source and Library ##################################################
library(regmdc)
library(foreach)  # for parallelization
library(doParallel)  # for parallelization
library(earth)  # for mars fitting
library(caret)  # for cross-validation
library(dplyr)
#######################################################################

#######################################################################
# Core Utilization ####################################################
num_cores <- 5L
registerDoParallel(num_cores)
#######################################################################

#######################################################################
# The Number of Repetition ############################################
num_rep <- 25L
#######################################################################

#######################################################################
# Data Cleaning #######################################################
wine_data <- read.csv("../data/winequality-red.csv", sep = ";")

wine_data <- wine_data %>% 
  mutate_at(c(1:11), function(x) {(x - min(x)) / (max(x) - min(x))}) %>% 
  as.matrix()

y <- wine_data[, 12]
X <- wine_data[, -12]
n <- nrow(X)
d <- ncol(X)
s <- 2L

num_bins <- rep(100L, d)
for(col in (1L:ncol(X))) {
  if (length(unique(X[, col])) < 100L) {
    num_bins[col] <- NA
  }
}
#######################################################################
  
estimation_results <- foreach(rep = 1L:num_rep, .combine = 'rbind', .errorhandling = "remove") %dopar% {
  #######################################################################
  # Data Split ##########################################################
  set.seed(2022L + rep)
  
  training_ratio <- .8
  training <- sample(length(y), size = length(y) * training_ratio)
  
  X_training <- X[training, ]
  y_training <- y[training]
  X_test <- X[-training, ]
  y_test <- y[-training]
  #######################################################################
  
  #######################################################################
  # (1) The usual MARS ##################################################
  df <- data.frame(y = y_training, X_training)
  
  # Build a usual MARS model
  mars_model <- earth(
    y ~ ., 
    data = df,
    degree = s  # interaction restriction
  )
  
  # Compute the fitted values at the test set
  mars_fit_test <- as.vector(predict(mars_model, X_test))
  
  # Compute the prediction error with the test set
  mars_pred_error <- mean((y_test - mars_fit_test)**2)
  #######################################################################
  
  
  #######################################################################
  # (2) Our model #######################################################
  # Parameter searching
  V_set <- c(10, 30, 50)
  parameters <- expand.grid(V = V_set) 
  
  # Perform k-fold cross-validation
  k <- 5L
  folds <- createFolds(y_training, k = k, list = TRUE, returnTrain = FALSE)
  
  cv_results <- matrix(, nrow = length(V_set), ncol = 2L)
  cv_results[, 1L] <- V_set 
  
  for (para_index in (1L:nrow(parameters))) {
    V <- parameters[para_index, 1L]
    cv_pred_errors <- numeric(k)
    
    for (fold_index in (1L:k)) {
      X_cv_training <- X_training[-folds[[fold_index]], ]
      y_cv_training <- y_training[-folds[[fold_index]]]
      X_cv_test <- X_training[folds[[fold_index]], ]
      y_cv_test <- y_training[folds[[fold_index]]]
      
      tryCatch({
        # Build a model
        regmdc_model <- regmdc(X_cv_training, y_cv_training, s, method = "mars", 
                               V, number_of_bins = num_bins)
        
        # Compute the fitted values at the cv-test data
        our_fit_cv_test <- predict_regmdc(regmdc_model, X_cv_test)
        
        # Compute the prediction error with the cv-test data
        cv_pred_errors[fold_index] <- mean((y_cv_test - our_fit_cv_test)**2)
      }, error = function(err) {
        cv_pred_errors[fold_index] <- NA
      })
    }
    
    cv_results[para_index, 2L] <- mean(cv_pred_errors, na.rm = TRUE)
  }
  
  V_best <- cv_results[which.min(cv_results[, 2L]), 1L]
  
  tryCatch({
    # Build a model
    regmdc_model <- regmdc(X_training, y_training, s, method = "mars", V_best,
                           number_of_bins = num_bins)
  }, error = function(err) {
    cv_results_sorted <- cv_results[order(cv_results[, 2L], decreasing = FALSE), ]
    V_best <- cv_results_sorted[2L, 1L]
    regmdc_model <- regmdc(X_training, y_training, s, method = "mars", 
                           V_best, number_of_bins = num_bins)
  })
  
  # Compute the fitted values at the test set
  our_fit_test <- as.vector(predict_regmdc(regmdc_model, X_test))
  
  # Compute the prediction error with the test set
  our_pred_error <- mean((y_test - our_fit_test)**2)
  #######################################################################
  
  c(mars_pred_error = mars_pred_error, 
    our_pred_error = our_pred_error, 
    V_best = V_best)
}

avg <- apply(estimation_results[, -3L], MARGIN = 2L, FUN = mean)
se <- apply(estimation_results[, -3L], MARGIN = 2L, FUN = function(x) {
  sqrt(var(x) / length(x))
})

mars_pred_err <- quantile(estimation_results[, 1L], prob = seq(0, 1, by = 0.25)) 
our_pred_err <- quantile(estimation_results[, 2L], prob = seq(0, 1, by = 0.25))

summary_statistics <- rbind(avg, se, cbind(mars_pred_err, our_pred_err))
#######################################################################
# Output Results ######################################################
results <- list(
  estimation_results = estimation_results,
  summary_statistics = summary_statistics
)
#######################################################################

