rm(list = ls())
set.seed(2022L)

#######################################################################
# Source and Library ##################################################
library(regmdc)
library(foreach)  # for parallelization
library(doParallel)  # for parallelization
library(earth)  # for mars fitting
library(caret)  # for cross-validation
#######################################################################

#######################################################################
# Core Utilization ####################################################
num_cores <- 10L
registerDoParallel(num_cores)
#######################################################################

#######################################################################
# The Number of Repetitions ###########################################
num_reps <- 25L
# The Number of Samples ###############################################
n <- 100L
#n <- 200L
#n <- 400L
#n <- 800L
# Dimension ###########################################################
d <- 2L
# Sigma ###############################################################
sigma <- 1.0
# Interaction Restriction #############################################
s <- 2L
# The Number of Samples for Loss Approximation ########################
num_sample_loss <- 1000L
# The Number of Bins for Approximation ################################
num_bins <- 25L
#######################################################################

#######################################################################
# Target function #####################################################
#######################################################################
fstar <- function(x) {
  (
    5.0 * sin(4.0 / (sqrt((x[1])**2 + (x[2])**2) + 0.001)) + 7.5
  )
}
#######################################################################

estimation_results <- foreach(rep = 1L:num_reps, .combine = 'rbind', .errorhandling = "remove") %dopar% {
  #####################################################################
  # Design Points #####################################################
  set.seed(2022L)
  
  X_design <- sapply((1L:d), simplify = TRUE, function(col) {
    runif(n)
  })
  
  # Compute the values of f* at the design points
  theta <- apply(X_design, MARGIN = 1L, FUN = fstar)
  # Compute the values of y at the design points
  y <- theta + sigma * rnorm(n)
  
  # Samples for approximating loss ####################################
  X_loss <- sapply((1L:d), simplify = TRUE, function(col) {
    runif(num_sample_loss)
  })
  theta_loss <- apply(X_loss, MARGIN = 1L, FUN = fstar)
  #####################################################################
  
  #####################################################################
  # (1) The usual MARS ################################################
  df <- data.frame(y = y, X_design)
  
  # Build a usual MARS model
  mars_model <- earth(
    y ~ ., 
    data = df,
    degree = s  # interaction restriction
  )
  
  # Compute loss 
  mars_fit_loss <- predict(mars_model, X_loss)
  mars_loss <- mean((theta_loss - mars_fit_loss)**2)
  #####################################################################
  
  #####################################################################
  # (2) Our model #####################################################
  # Parameter searching
  V_set <- c(500, 1000, 5000)  # n = 100, 200
  #V_set <- c(5000, 10000, 50000)  # n = 400, 800
  parameters <- expand.grid(V = V_set) 
  
  # Perform k-fold cross-validation
  k <- 10L
  folds <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  
  cv_results <- matrix(, nrow = length(V_set), ncol = 2L)
  cv_results[, 1L] <- V_set 
  
  for (para_index in (1L:nrow(parameters))) {
    V <- parameters[para_index, 1L]
    cv_pred_errors <- numeric(k)
    
    for (fold_index in (1L:k)) {
      X_cv_training <- X_design[-folds[[fold_index]], ]
      y_cv_training <- y[-folds[[fold_index]]]
      X_cv_test <- X_design[folds[[fold_index]], ]
      y_cv_test <- y[folds[[fold_index]]]
      
      tryCatch({
        # Build a model
        regmdc_model <- regmdc(X_cv_training, y_cv_training, s, method = "mars", 
                               V, number_of_bins = num_bins)
        # regmdc_model <- regmdc(X_cv_training, y_cv_training, s, method = "mars", V)
        
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
    regmdc_model <- regmdc(X_design, y, s, method = "mars", V_best,
                           number_of_bins = num_bins)
    # regmdc_model <- regmdc(X_design, y, s, method = "mars", V_best)
  }, error = function(err) {
    cv_results_sorted <- cv_results[order(cv_results[, 2L], decreasing = FALSE), ]
    V_best <- cv_results_sorted[2L, 1L]
    regmdc_model <- regmdc(X_design, y, s, method = "mars", 
                           V_best, number_of_bins = num_bins)
    # regmdc_model <- regmdc(X_design, y, s, method = "mars", V_best)
  })
  
  # Compute loss 
  our_fit_loss <- predict_regmdc(regmdc_model, X_loss)
  our_loss <- mean((theta_loss - our_fit_loss)**2)
  ####################################################################### 
  
  c(mars_loss = mars_loss, our_loss = our_loss, V_best = V_best)
}

avg <- apply(estimation_results[, -3L], MARGIN = 2L, FUN = mean)
se <- apply(estimation_results[, -3L], MARGIN = 2L, FUN = function(x) {
  sqrt(var(x) / length(x))
})

mars_loss <- quantile(estimation_results[, 1L], prob = seq(0, 1, by = 0.25)) 
our_loss <- quantile(estimation_results[, 2L], prob = seq(0, 1, by = 0.25))

summary_statistics <- rbind(avg, se, cbind(mars_loss, our_loss))
#######################################################################
# Output Results ######################################################
results <- list(
  estimation_results = estimation_results,
  summary_statistics = summary_statistics
)
#######################################################################
