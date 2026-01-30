








### Load libraries
library(betareg)


# Helper functions for extracting model information
extract_formula <- function(mod) {
  if (inherits(mod,"betareg")) {
    return(formula(mod))
  } else {
    return(formula(mod))
  }
}

extract_model_type <- function(mod) {
  if (inherits(mod,"betareg")) {
    return("beta")
  } else {
    return("glm")
  }
}

extract_family <- function(mod) {
  if (inherits(mod,"betareg")) {
    return(NULL)
  } else {
    return(family(mod))
  }
}


# Calculate cross-fitted model predictions
crossfit_predict <- function(formula,data,model_type = "glm",
                             family = binomial(),newdata_list = NULL,
                             folds,predict_type = "response") {
  
  n_folds <- max(folds)
  
  # if newdata_list is NULL, just predict on original data
  if (is.null(newdata_list)) {
    newdata_list <- list(original = data)
  }
  
  # initialize output list
  n_scenarios <- length(newdata_list)
  predictions_list <- vector("list",n_scenarios)
  names(predictions_list) <- names(newdata_list)
  
  # initialize prediction vectors for each scenario
  for (i in 1:n_scenarios) {
    predictions_list[[i]] <- numeric(nrow(newdata_list[[i]]))
  }
  
  # get fold assignments for prediction data
  pred_folds <- if ("fold" %in% names(newdata_list[[1]])) {
    newdata_list[[1]]$fold
  } else {
    folds
  }
  
  # loop through each fold
  for (k in 1:n_folds) {
    # training data: all folds except k
    train_idx <- which(folds != k)
    train_data <- data[train_idx,]
    
    # test indices for prediction
    test_idx <- which(pred_folds == k)
    
    if (length(test_idx) == 0) next
    
    # fit model once on training data
    if (model_type == "glm") {
      fit <- glm(formula,data = train_data,family = family)
    } else if (model_type == "beta") {
      fit <- betareg::betareg(formula,data = train_data)
    } else {
      stop("Unsupported model_type. Use 'glm' or 'beta'.")
    }
    
    # use the fitted model to predict on all scenarios
    for (i in 1:n_scenarios) {
      test_data <- newdata_list[[i]][test_idx, , drop = FALSE]
      predictions_list[[i]][test_idx] <- predict(fit,newdata = test_data,type = predict_type)
    }
  }
  
  # if only one scenario, return vector instead of list
  if (n_scenarios == 1) {
    return(predictions_list[[1]])
  }
  
  return(predictions_list)
}


# Calculates a robust TMLE estimator for the ATE for a single outcome
robustcompATE <- function(outcome_mod,ps_mod,covariates,
                          exposure_name,outcome_name,
                          trunc_level,scale_factor = 1.0,
                          crossFit = FALSE,n_folds = 5,folds = NULL) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  Y <- covariates[[outcome_name]]
  
  # create fold assignments if using cross-fitting
  if (crossFit) {
    # use provided folds or create new ones
    if (is.null(folds)) {
      folds <- sample(rep(1:n_folds,length.out = n))
    }
    covariates$fold <- folds
    
    # extract formulas and model information from fitted models
    ps_formula <- extract_formula(ps_mod)
    outcome_formula <- extract_formula(outcome_mod)
    
    ps_type <- extract_model_type(ps_mod)
    outcome_type <- extract_model_type(outcome_mod)
    
    ps_family <- extract_family(ps_mod)
    outcome_family <- extract_family(outcome_mod)
  }
  
  ### calculate propensity score P(A=1|W)
  if (crossFit) {
    ps_pred <- crossfit_predict(formula = ps_formula,
                                data = covariates,
                                model_type = ps_type,
                                family = ps_family,
                                newdata_list = list(original = covariates),
                                folds = folds,
                                predict_type = "response")
  } else {
    ps_pred <- predict(ps_mod,newdata = covariates,type = "response")
  }
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  
  ### calculate clever covariates for TMLE
  # Observed clever covariate (depends on actual treatment received)
  H_A_Y <- ifelse(exposure == 1,
                  1.0 / ps_trunc,
                  -1.0 / (1.0 - ps_trunc))
  
  # Counterfactual clever covariates (for all observations)
  H_0_cf <- -1.0 / (1.0 - ps_trunc)
  H_1_cf <- 1.0 / ps_trunc
  
  ### get outcome predictions
  covariates_0 <- covariates
  covariates_0[[exposure_name]] <- 0
  
  covariates_1 <- covariates
  covariates_1[[exposure_name]] <- 1
  
  if (crossFit) {
    outcome_preds <- crossfit_predict(formula = outcome_formula,
                                      data = covariates,
                                      model_type = outcome_type,
                                      family = outcome_family,
                                      newdata_list = list(original = covariates,
                                                          a0 = covariates_0,
                                                          a1 = covariates_1),
                                      folds = folds,
                                      predict_type = "response")
    pred_y_A <- outcome_preds$original
    pred_y_0 <- outcome_preds$a0
    pred_y_1 <- outcome_preds$a1
  } else {
    pred_y_A <- predict(outcome_mod,newdata = covariates,type = "response")
    pred_y_0 <- predict(outcome_mod,newdata = covariates_0,type = "response")
    pred_y_1 <- predict(outcome_mod,newdata = covariates_1,type = "response")
  }
  
  ### estimate fluctuation parameter
  data_tmle <- data.frame(Y = Y,
                          pred_vals = pred_y_A,
                          H_A_Y = H_A_Y)
  eps <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                  data = data_tmle,family = "binomial"))
  
  ### update predictions using fluctuation parameter
  # Use counterfactual clever covariates for counterfactual predictions
  pred_star_y_A <- plogis(qlogis(pred_y_A) + eps * H_A_Y)
  pred_star_y_0 <- plogis(qlogis(pred_y_0) + eps * H_0_cf)
  pred_star_y_1 <- plogis(qlogis(pred_y_1) + eps * H_1_cf)
  
  ### calculate ATE using targeted estimates (on original scale)
  ate_star_unscaled <- pred_star_y_1 - pred_star_y_0
  ate_est_unscaled <- mean(ate_star_unscaled)
  ate_star_scaled <- ate_star_unscaled * scale_factor
  ate_est_scaled <- ate_est_unscaled * scale_factor
  
  ### calculate efficient influence function
  D_Y <- H_A_Y * (Y - pred_star_y_A)
  D_W <- ate_star_unscaled - ate_est_unscaled
  EIF <- (D_Y + D_W) * scale_factor
  
  ### calculate standard errors and confidence intervals
  var_ate <- var(EIF) / n
  se_ate <- sqrt(var_ate)
  ci_lower <- ate_est_scaled - qnorm(0.975) * se_ate
  ci_upper <- ate_est_scaled + qnorm(0.975) * se_ate
  
  ### return results
  return(list(ATE = ate_est_scaled,
              ITE = ate_star_scaled,
              IF = EIF,
              SE = se_ate,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Calculate CATE using robust TMLE estimates
# CATE = (E[Y|Z=1,W] - E[Y|Z=0,W]) / (E[A|Z=1,W] - E[A|Z=0,W])
calculate_cate_robust <- function(data,outcome_mod_y,outcome_mod_a,iv_model,
                                  iv_name = "Z",exposure_name = "rehabIRF",outcome_name_01,
                                  scale_factor = 1.0,trunc_level_num = 0.001,
                                  trunc_level_denom = 0.001,
                                  crossFit = FALSE,n_folds = 5) {
  
  n <- nrow(data)
  
  # create shared fold assignments if using cross-fitting
  if (crossFit) {
    folds <- sample(rep(1:n_folds,length.out = n))
  } else {
    folds <- NULL
  }
  
  ### calculate ATE for numerator (effect of Z on Y)
  results_numerator <- robustcompATE(outcome_mod = outcome_mod_y,
                                     ps_mod = iv_model,
                                     covariates = data,
                                     exposure_name = iv_name,
                                     outcome_name = outcome_name_01,
                                     trunc_level = trunc_level_num,
                                     scale_factor = scale_factor,
                                     crossFit = crossFit,
                                     n_folds = n_folds,
                                     folds = folds)
  
  ### calculate ATE for denominator (effect of Z on A)
  results_denominator <- robustcompATE(outcome_mod = outcome_mod_a,
                                       ps_mod = iv_model,
                                       covariates = data,
                                       exposure_name = iv_name,
                                       outcome_name = exposure_name,
                                       trunc_level = trunc_level_denom,
                                       scale_factor = 1.0,
                                       crossFit = crossFit,
                                       n_folds = n_folds,
                                       folds = folds)
  
  ### calculate CATE and approximate 95% confidence interval using delta method
  numerator <- results_numerator$ATE
  denominator <- results_denominator$ATE
  cate <- numerator / denominator
  
  ### delta method for variance of ratio using influence functions
  IF_numerator <- results_numerator$IF
  IF_denominator <- results_denominator$IF
  
  # IF for ratio: (IF_num - CATE * IF_denom) / denom
  IF_cate <- (IF_numerator - cate * IF_denominator) / denominator
  
  var_cate <- var(IF_cate) / n
  se_cate <- sqrt(var_cate)
  ci_lower <- cate - qnorm(0.975) * se_cate
  ci_upper <- cate + qnorm(0.975) * se_cate
  
  ### return results
  return(list(CATE = cate,
              SE = se_cate,
              CI = c(lower = ci_lower,upper = ci_upper),
              numerator = numerator,
              denominator = denominator,
              numerator_se = results_numerator$SE,
              denominator_se = results_denominator$SE,
              IF = IF_cate,
              IF_numerator = IF_numerator,
              IF_denominator = IF_denominator))
}


# Calculate robust ATE using doubly robust g-estimation with IV
calculate_ate_iv_robust <- function(outcome_mod_init,iv_model,data,
                                    exp_model_formula,delta_formula,
                                    iv_name = "Z",
                                    exposure_name = "rehabIRF",
                                    outcome_name_01,
                                    scale_factor = 1.0,trunc_level = 0.001) {
  
  # extract IV and exposure variables
  Z <- data[[iv_name]]
  D <- data[[exposure_name]]
  
  ### calculate P_0^D(X)
  iv_prob <- predict(iv_model,newdata = data,type = "response")
  iv_prob_trunc <- pmax(pmin(iv_prob,1.0 - trunc_level),trunc_level)
  f_Z_given_X <- ifelse(Z == 1,iv_prob_trunc,1 - iv_prob_trunc)
  iv_weights_0 <- ifelse(Z == 0,1.0 / (1.0 - iv_prob_trunc),0)
  iv_weights_1 <- ifelse(Z == 1,1.0 / iv_prob_trunc,0)
  data_z0 <- data[Z == 0,]
  weights_z0 <- iv_weights_0[Z == 0]
  environment(exp_model_formula) <- environment()
  exp_model_weighted <- glm(exp_model_formula,
                            data = data_z0,
                            weights = weights_z0,
                            family = binomial())
  P_D0 <- predict(exp_model_weighted,newdata = data,type = "response")
  
  ### solve estimating equation (14) for delta_D(X)
  data_z1 <- data[Z == 1,]
  D_z1 <- D[Z == 1]
  P_D0_z1 <- P_D0[Z == 1]
  weights_z1 <- iv_weights_1[Z == 1]
  outcome_delta <- D_z1 - P_D0_z1
  rhs_terms <- as.character(delta_formula)[2]
  delta_d_formula <- as.formula(paste("outcome_delta ~",rhs_terms))
  environment(delta_d_formula) <- environment()
  delta_d_mod <- lm(delta_d_formula,
                    data = data_z1,
                    weights = weights_z1)
  delta_D_X <- predict(delta_d_mod,newdata = data)
  
  ### calculate P_0^Y(X)
  P_Y_init <- predict(outcome_mod_init,newdata = data,type = "response")
  H_Z_Y <- ifelse(Z == 0,1.0 / (1.0 - iv_prob_trunc),0)
  Y_values <- data[[outcome_name_01]]
  data_fluct_Y <- data.frame(Y = Y_values,
                             pred_vals = P_Y_init,
                             H_Z_Y = H_Z_Y)
  fluct_mod_Y <- glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_Z_Y,
                     data = data_fluct_Y,
                     family = binomial())
  eps_Y <- coef(fluct_mod_Y)
  H_Z_Y_Z0 <- 1.0 / (1.0 - iv_prob_trunc)
  data_Z0 <- data
  data_Z0[[iv_name]] <- 0
  P_Y_init_Z0 <- predict(outcome_mod_init,
                         newdata = data_Z0,
                         type = "response")
  P_Y0 <- plogis(qlogis(P_Y_init_Z0) + eps_Y * H_Z_Y_Z0)
  
  ### now solve estimating equation (14) for delta_Y(X)
  Y_z1 <- Y_values[Z == 1]
  P_Y0_z1 <- P_Y0[Z == 1]
  X_matrix_z1 <- model.matrix(delta_formula,data = data_z1)
  Y_minus_PY0 <- Y_z1 - P_Y0_z1
  D_minus_PD0 <- D_z1 - P_D0_z1
  W_z1 <- weights_z1
  A_matrix <- t(X_matrix_z1 * (W_z1 * D_minus_PD0)) %*% X_matrix_z1
  b_vector <- t(X_matrix_z1 * W_z1) %*% Y_minus_PY0
  alpha_hat <- solve(A_matrix,b_vector)
  X_matrix_full <- model.matrix(delta_formula,data = data)
  delta_X <- as.vector(X_matrix_full %*% alpha_hat)
  
  ### calculate multiply robust estimator (15)
  Y_minus_PY0 <- Y_values - P_Y0
  D_minus_PD0 <- D - P_D0
  residual_Y <- Y_minus_PY0 - D_minus_PD0 * delta_X
  IV_weight_mr <- ifelse(Z == 1,iv_weights_1,iv_weights_0)
  ate_vals <- ((2 * Z - 1) * (1.0 / delta_D_X) * IV_weight_mr * residual_Y + delta_X) * scale_factor
  ate_mr <- mean(ate_vals)
  
  ### calculate influence function and 95% confidence intervals
  IF_component_1 <- (2 * Z - 1) * (1.0 / delta_D_X) * IV_weight_mr * residual_Y
  IF_component_2 <- delta_X
  IF_total <- (IF_component_1 + IF_component_2 - mean(IF_component_1 + IF_component_2)) * scale_factor
  n <- nrow(data)
  var_ate <- var(IF_total) / n
  se_ate <- sqrt(var_ate)
  ci_lower <- ate_mr - qnorm(0.975) * se_ate
  ci_upper <- ate_mr + qnorm(0.975) * se_ate
  
  ### return results
  return(list(ATE = ate_mr,
              ITE = ate_vals,
              SE = se_ate,
              CI = c(lower = ci_lower,upper = ci_upper),
              IF = IF_total,
              delta_D_X = delta_D_X,
              delta_X = delta_X))
}


# Generate simulated data with instrumental variable for causal inference
generate_stroke_data_iv <- function(n = 1000,
                                    iv_name = "Z",
                                    exposure_name = "rehabIRF",
                                    iv_strength = 0.3,
                                    fixed_covariates = NULL) {
  
  ### generate or use fixed baseline covariates
  if (is.null(fixed_covariates)) {
    # age
    age <- rnorm(n,mean = 70,sd = 12)
    age_scaled <- (age - mean(age)) / sd(age)
    
    # post-discharge disability as the unmeasured confounder
    post_discharge_disability <- rbeta(n,shape1 = 2,shape2 = 3) * 100
    post_discharge_disability_scaled <- (post_discharge_disability - mean(post_discharge_disability)) / 
      sd(post_discharge_disability)
    
    # stroke severity
    stroke_severity <- rpois(n,lambda = 8)
    stroke_severity <- pmin(stroke_severity,42)
    stroke_severity_scaled <- (stroke_severity - mean(stroke_severity)) / 
      sd(stroke_severity)
    
    # comorbidity score
    comorbidity <- rpois(n,lambda = 2)
    comorbidity <- pmin(comorbidity,10)
    comorbidity_scaled <- (comorbidity - mean(comorbidity)) / sd(comorbidity)
  } else {
    # use fixed covariates
    age <- fixed_covariates$age
    age_scaled <- fixed_covariates$age_scaled
    stroke_severity <- fixed_covariates$stroke_severity
    stroke_severity_scaled <- fixed_covariates$stroke_severity_scaled
    comorbidity <- fixed_covariates$comorbidity
    comorbidity_scaled <- fixed_covariates$comorbidity_scaled
    post_discharge_disability_scaled <- fixed_covariates$post_discharge_disability_scaled
  }
  
  ### generate instrumental variable
  # IV depends weakly on observed covariates X
  # IV independent of the unmeasured confounder U
  logit_iv <- 0.0 + 
    0.05 * age_scaled + 
    0.05 * stroke_severity_scaled - 
    0.05 * comorbidity_scaled
  
  iv_prob <- plogis(logit_iv)
  IV <- rbinom(n,size = 1,prob = iv_prob)
  
  ### generate treatment assignment using linear probability model
  # no Z-U interaction on the additive scale
  # higher post_discharge_disability -> higher treatment probability
  ps <- 0.20 + 
    iv_strength * IV +
    0.02 * age_scaled +
    0.05 * post_discharge_disability_scaled +
    0.02 * stroke_severity_scaled +
    -0.01 * comorbidity_scaled
  
  # check for invalid probabilities
  if (any(ps < 0 | ps > 1)) {
    stop(paste("Invalid probabilities in treatment model. Range:",
               round(min(ps),3),"to",round(max(ps),3),
               "\nAdjust coefficients or iv_strength."))
  }
  
  # generate treatment
  exposure <- rbinom(n,size = 1,prob = ps)
  
  ### generate outcomes conditional on treatment, observed covariates, and unmeasured confounder
  ## ADL/IADL
  # higher post_discharge_disability -> higher three-month disability
  # treatment (exposure) -> lower three-month disability
  mu_adl_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.5 * exposure +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  
  mu_adl <- plogis(mu_adl_logit)
  
  phi_adl <- 15
  shape1_adl <- mu_adl * phi_adl
  shape2_adl <- (1 - mu_adl) * phi_adl
  
  OUT3_ADL_IADL_01 <- rbeta(n,shape1 = shape1_adl,shape2 = shape2_adl)
  OUT3_ADL_IADL <- OUT3_ADL_IADL_01 * 3.0
  
  ### calculate true ATE of treatment A on outcome Y
  mu_adl_0_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 0 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  
  mu_adl_1_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 1 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  
  true_ate_adl <- mean(plogis(mu_adl_1_logit) - plogis(mu_adl_0_logit)) * 3.0
  
  ### compile dataset
  data <- data.frame(id = 1:n,
                     age = age,
                     stroke_severity = stroke_severity,
                     comorbidity = comorbidity,
                     OUT3_ADL_IADL = OUT3_ADL_IADL,
                     propensity_score = ps,
                     iv_prob = iv_prob)
  
  # add IV and exposure with specified names
  data[[iv_name]] <- IV
  data[[exposure_name]] <- exposure
  
  ### return results
  return(list(data = data,
              true_ate_adl = true_ate_adl,
              simulation_params = list(n = n,
                                       iv_name = iv_name,
                                       exposure_name = exposure_name,
                                       phi_adl = phi_adl,
                                       iv_strength = iv_strength,
                                       prop_IV1 = mean(IV),
                                       prop_treated = mean(exposure)),
              fixed_covariates = data.frame(age = age,
                                            age_scaled = age_scaled,
                                            stroke_severity = stroke_severity,
                                            stroke_severity_scaled = stroke_severity_scaled,
                                            comorbidity = comorbidity,
                                            comorbidity_scaled = comorbidity_scaled,
                                            post_discharge_disability_scaled = post_discharge_disability_scaled)))
}


# Run a single simulation iteration for IV analysis
# Note: expects fixed_covariates to exist in parent environment (set by run_simulation_study_iv)
run_single_simulation_iv <- function(n = 1000,
                                     iv_name = "Z",
                                     exposure_name = "rehabIRF",
                                     trunc_level = 0.001,
                                     iv_strength = 0.5,
                                     crossFit = FALSE,
                                     n_folds = 5) {
  
  # generate data with IV (using fixed_covariates from parent environment)
  sim_data <- generate_stroke_data_iv(n = n,
                                      iv_name = iv_name,
                                      exposure_name = exposure_name,
                                      iv_strength = iv_strength,
                                      fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  
  # re-scale outcomes to 0-1 for modeling
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  
  # apply transformation for beta regression to avoid boundary issues
  N <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N - 1) + 0.5) / N
  
  #########
  ### Fit models with hardcoded formulas
  #########
  
  # IV propensity score model P(Z|W)
  iv_mod <- glm(Z ~ age + stroke_severity + comorbidity,
                data = data, family = binomial())
  
  # Outcome model E[Y|Z,W] for CATE numerator
  outcome_mod_y <- betareg(OUT3_ADL_IADL_01 ~ Z + age + stroke_severity + comorbidity,
                           data = data)
  
  # Exposure model E[A|Z,W] for CATE denominator
  outcome_mod_a <- glm(rehabIRF ~ Z + age + stroke_severity + comorbidity,
                       data = data, family = binomial())
  
  # Initial outcome model for IV-based ATE (without Z)
  outcome_mod_init <- betareg(OUT3_ADL_IADL_01 ~ age + stroke_severity + comorbidity,
                              data = data)
  
  # Propensity score model for naive ATE P(A|W)
  ps_mod_naive <- glm(rehabIRF ~ age + stroke_severity + comorbidity,
                      data = data, family = binomial())
  
  # Outcome model E[Y|A,W] for naive ATE
  outcome_mod_naive <- betareg(OUT3_ADL_IADL_01 ~ rehabIRF + age + stroke_severity + comorbidity,
                               data = data)
  
  # Formulas for IV-based ATE calculation
  exp_model_formula <- rehabIRF ~ age + stroke_severity + comorbidity
  delta_formula <- ~ age + stroke_severity + comorbidity
  
  # Initialize CATE results
  cate_adl <- NA
  cate_se_adl <- NA
  cate_ci_lower_adl <- NA
  cate_ci_upper_adl <- NA
  cate_numerator_adl <- NA
  cate_denominator_adl <- NA
  cate_numerator_se_adl <- NA
  cate_denominator_se_adl <- NA
  cate_converged <- FALSE
  
  # Initialize IV-based ATE results
  ate_adl <- NA
  ate_se_adl <- NA
  ate_ci_lower_adl <- NA
  ate_ci_upper_adl <- NA
  ate_converged <- FALSE
  
  # Initialize Naive ATE results
  ate_naive_adl <- NA
  ate_naive_se_adl <- NA
  ate_naive_ci_lower_adl <- NA
  ate_naive_ci_upper_adl <- NA
  ate_naive_converged <- FALSE
  
  ### Calculate CATE for ADL
  tryCatch({
    results_cate <- calculate_cate_robust(data = data,
                                          outcome_mod_y = outcome_mod_y,
                                          outcome_mod_a = outcome_mod_a,
                                          iv_model = iv_mod,
                                          iv_name = iv_name,
                                          exposure_name = exposure_name,
                                          outcome_name_01 = "OUT3_ADL_IADL_01",
                                          scale_factor = 3.0,
                                          trunc_level_num = trunc_level,
                                          trunc_level_denom = trunc_level,
                                          crossFit = crossFit,
                                          n_folds = n_folds)
    
    cate_adl <- results_cate$CATE
    cate_se_adl <- results_cate$SE
    cate_ci_lower_adl <- results_cate$CI["lower"]
    cate_ci_upper_adl <- results_cate$CI["upper"]
    cate_numerator_adl <- results_cate$numerator
    cate_denominator_adl <- results_cate$denominator
    cate_numerator_se_adl <- results_cate$numerator_se
    cate_denominator_se_adl <- results_cate$denominator_se
    cate_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("CATE estimation failed:", e$message))
  })
  
  ### Calculate IV-based ATE for ADL
  tryCatch({
    results_ate <- calculate_ate_iv_robust(outcome_mod_init = outcome_mod_init,
                                           iv_model = iv_mod,
                                           data = data,
                                           exp_model_formula = exp_model_formula,
                                           delta_formula = delta_formula,
                                           iv_name = iv_name,
                                           exposure_name = exposure_name,
                                           outcome_name_01 = "OUT3_ADL_IADL_01",
                                           scale_factor = 3.0,
                                           trunc_level = trunc_level)
    
    ate_adl <- results_ate$ATE
    ate_se_adl <- results_ate$SE
    ate_ci_lower_adl <- results_ate$CI["lower"]
    ate_ci_upper_adl <- results_ate$CI["upper"]
    ate_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("IV-based ATE estimation failed:", e$message))
  })
  
  ### Calculate naive ATE for ADL (ignores unmeasured confounding)
  tryCatch({
    results_ate_naive <- robustcompATE(outcome_mod = outcome_mod_naive,
                                       ps_mod = ps_mod_naive,
                                       covariates = data,
                                       exposure_name = exposure_name,
                                       outcome_name = "OUT3_ADL_IADL_01",
                                       trunc_level = trunc_level,
                                       scale_factor = 3.0,
                                       crossFit = crossFit,
                                       n_folds = n_folds)
    
    ate_naive_adl <- results_ate_naive$ATE
    ate_naive_se_adl <- results_ate_naive$SE
    ate_naive_ci_lower_adl <- results_ate_naive$CI["lower"]
    ate_naive_ci_upper_adl <- results_ate_naive$CI["upper"]
    ate_naive_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("Naive ATE estimation failed:", e$message))
  })
  
  ### Calculate bias and coverage
  # CATE bias and coverage
  if (cate_converged && !is.na(cate_adl) && is.finite(cate_adl)) {
    cate_bias <- cate_adl - sim_data$true_ate_adl
    cate_covers <- (cate_ci_lower_adl <= sim_data$true_ate_adl) & 
      (sim_data$true_ate_adl <= cate_ci_upper_adl)
  } else {
    cate_bias <- NA
    cate_covers <- NA
  }
  
  # IV-based ATE bias and coverage
  if (ate_converged && !is.na(ate_adl) && is.finite(ate_adl)) {
    ate_bias <- ate_adl - sim_data$true_ate_adl
    ate_covers <- (ate_ci_lower_adl <= sim_data$true_ate_adl) & 
      (sim_data$true_ate_adl <= ate_ci_upper_adl)
  } else {
    ate_bias <- NA
    ate_covers <- NA
  }
  
  # Naive ATE bias and coverage
  if (ate_naive_converged && !is.na(ate_naive_adl) && is.finite(ate_naive_adl)) {
    ate_naive_bias <- ate_naive_adl - sim_data$true_ate_adl
    ate_naive_covers <- (ate_naive_ci_lower_adl <= sim_data$true_ate_adl) & 
      (sim_data$true_ate_adl <= ate_naive_ci_upper_adl)
  } else {
    ate_naive_bias <- NA
    ate_naive_covers <- NA
  }
  
  ### Return results
  data.frame(n = n,
             iv_name = iv_name,
             exposure_name = exposure_name,
             trunc_level = trunc_level,
             iv_strength = iv_strength,
             crossFit = crossFit,
             # CATE results for ADL
             cate_adl = cate_adl,
             cate_se_adl = cate_se_adl,
             cate_ci_lower_adl = cate_ci_lower_adl,
             cate_ci_upper_adl = cate_ci_upper_adl,
             cate_bias = cate_bias,
             cate_covers = cate_covers,
             cate_converged = cate_converged,
             cate_numerator_adl = cate_numerator_adl,
             cate_denominator_adl = cate_denominator_adl,
             cate_numerator_se_adl = cate_numerator_se_adl,
             cate_denominator_se_adl = cate_denominator_se_adl,
             # IV-based ATE results for ADL
             ate_adl = ate_adl,
             ate_se_adl = ate_se_adl,
             ate_ci_lower_adl = ate_ci_lower_adl,
             ate_ci_upper_adl = ate_ci_upper_adl,
             ate_bias = ate_bias,
             ate_covers = ate_covers,
             ate_converged = ate_converged,
             # Naive ATE results for ADL
             ate_naive_adl = ate_naive_adl,
             ate_naive_se_adl = ate_naive_se_adl,
             ate_naive_ci_lower_adl = ate_naive_ci_lower_adl,
             ate_naive_ci_upper_adl = ate_naive_ci_upper_adl,
             ate_naive_bias = ate_naive_bias,
             ate_naive_covers = ate_naive_covers,
             ate_naive_converged = ate_naive_converged,
             # True values
             true_ate_adl = sim_data$true_ate_adl,
             # Proportions
             prop_IV = mean(data[[iv_name]]),
             prop_exposure = mean(data[[exposure_name]]))
}


# Run simulation study for IV analysis
run_simulation_study_iv <- function(n_sims = 500,n = 1000,
                                    iv_name = "Z",
                                    exposure_name = "rehabIRF",
                                    trunc_level = 0.001,iv_strength = 0.5,
                                    crossFit = FALSE,n_folds = 5) {
  
  # Generate fixed covariates once and assign to parent environment
  initial_data <- generate_stroke_data_iv(n = n,
                                          iv_name = iv_name,
                                          exposure_name = exposure_name,
                                          iv_strength = iv_strength)
  # assign covariates to parent environment
  # considerably faster than passing as an argument for each simulation
  fixed_covariates <<- initial_data$fixed_covariates
  
  # initialize results list
  results_list <- vector("list",n_sims)
  
  # initialize progress bar
  pb <- txtProgressBar(min = 0,max = n_sims,style = 3)
  
  # run sims sequentially with progress bar
  for (i in 1:n_sims) {
    results_list[[i]] <- tryCatch({
      run_single_simulation_iv(n = n,
                               iv_name = iv_name,
                               exposure_name = exposure_name,
                               trunc_level = trunc_level,
                               iv_strength = iv_strength,
                               crossFit = crossFit,
                               n_folds = n_folds)
    },error = function(e) {
      cat(paste("\nError in simulation",i,":",e$message,"\n"))
      return(NULL)
    })
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  # remove NULL results from failed simulations
  results_list <- results_list[!sapply(results_list,is.null)]
  results <- do.call(rbind,results_list)
  
  return(results)
}


# Summarize IV simulation results for CATE, ATE, and naive ATE
summarize_simulation_iv <- function(results) {
  
  # naive ATE summary statistics
  ate_naive_summary <- list(mean_estimate = mean(results$ate_naive_adl,na.rm = TRUE),
                            mean_bias = mean(results$ate_naive_bias,na.rm = TRUE),
                            mean_se = mean(results$ate_naive_se_adl,na.rm = TRUE),
                            empirical_se = sd(results$ate_naive_adl,na.rm = TRUE),
                            coverage = mean(results$ate_naive_covers,na.rm = TRUE),
                            convergence_rate = mean(results$ate_naive_converged,na.rm = TRUE))
  
  # CATE summary statistics
  cate_summary <- list(mean_estimate = mean(results$cate_adl,na.rm = TRUE),
                       mean_bias = mean(results$cate_bias,na.rm = TRUE),
                       mean_se = mean(results$cate_se_adl,na.rm = TRUE),
                       empirical_se = sd(results$cate_adl,na.rm = TRUE),
                       coverage = mean(results$cate_covers,na.rm = TRUE),
                       convergence_rate = mean(results$cate_converged,na.rm = TRUE),
                       mean_numerator = mean(results$cate_numerator_adl,na.rm = TRUE),
                       mean_denominator = mean(results$cate_denominator_adl,na.rm = TRUE),
                       mean_numerator_se = mean(results$cate_numerator_se_adl,na.rm = TRUE),
                       mean_denominator_se = mean(results$cate_denominator_se_adl,na.rm = TRUE))
  
  # ATE (IV-based) summary statistics
  ate_summary <- list(mean_estimate = mean(results$ate_adl,na.rm = TRUE),
                      mean_bias = mean(results$ate_bias,na.rm = TRUE),
                      mean_se = mean(results$ate_se_adl,na.rm = TRUE),
                      empirical_se = sd(results$ate_adl,na.rm = TRUE),
                      coverage = mean(results$ate_covers,na.rm = TRUE),
                      convergence_rate = mean(results$ate_converged,na.rm = TRUE))
  
  # IV diagnostics
  iv_diagnostics <- list(mean_prop_IV = mean(results$prop_IV,na.rm = TRUE),
                         mean_prop_exposure = mean(results$prop_exposure,na.rm = TRUE),
                         iv_strength_setting = results$iv_strength[1],
                         iv_name = results$iv_name[1],
                         exposure_name = results$exposure_name[1],
                         # First stage strength from CATE denominator
                         mean_first_stage = mean(results$cate_denominator_adl,na.rm = TRUE))
  
  # true values should be constant across simulations
  true_values <- list(true_ate_adl = mean(results$true_ate_adl,na.rm = TRUE))
  
  summary_stats <- list(ATE_naive = ate_naive_summary,
                        CATE = cate_summary,
                        ATE = ate_summary,
                        IV_diagnostics = iv_diagnostics,
                        true_values = true_values,
                        n_simulations = nrow(results),
                        n = results$n[1])
  
  return(summary_stats)
}


# Print formatted summary of IV simulation results
print_simulation_summary_iv <- function(summary_stats) {
  
  cat("\n")
  cat("IV SIMULATION RESULTS\n")
  cat("=====================\n")
  cat(sprintf("N simulations: %d | n: %d | True ATE: %.4f\n",
              summary_stats$n_simulations,
              summary_stats$n,
              summary_stats$true_values$true_ate_adl))
  cat(sprintf("IV strength: %.2f | First stage (CATE denom): %.4f\n",
              summary_stats$IV_diagnostics$iv_strength_setting,
              summary_stats$IV_diagnostics$mean_first_stage))
  cat(sprintf("IV: '%s' | Exposure: '%s'\n",
              summary_stats$IV_diagnostics$iv_name,
              summary_stats$IV_diagnostics$exposure_name))
  cat(sprintf("P(IV=1): %.3f | P(Treated): %.3f\n\n",
              summary_stats$IV_diagnostics$mean_prop_IV,
              summary_stats$IV_diagnostics$mean_prop_exposure))
  
  # Create results table with model-based SE
  cat("Method      Est      Bias     Model_SE Emp_SE   Coverage Conv\n")
  cat("----------  -------  -------  -------- -------  -------- ----\n")
  
  # Naive ATE results
  cat(sprintf("ATE (Naive) %7.4f  %7.4f  %7.4f  %7.4f  %7.3f  %.2f\n",
              summary_stats$ATE_naive$mean_estimate,
              summary_stats$ATE_naive$mean_bias,
              summary_stats$ATE_naive$mean_se,
              summary_stats$ATE_naive$empirical_se,
              summary_stats$ATE_naive$coverage,
              summary_stats$ATE_naive$convergence_rate))
  
  # CATE results
  cat(sprintf("CATE        %7.4f  %7.4f  %7.4f  %7.4f  %7.3f  %.2f\n",
              summary_stats$CATE$mean_estimate,
              summary_stats$CATE$mean_bias,
              summary_stats$CATE$mean_se,
              summary_stats$CATE$empirical_se,
              summary_stats$CATE$coverage,
              summary_stats$CATE$convergence_rate))
  
  # ATE (IV) results
  cat(sprintf("ATE (IV)    %7.4f  %7.4f  %7.4f  %7.4f  %7.3f  %.2f\n",
              summary_stats$ATE$mean_estimate,
              summary_stats$ATE$mean_bias,
              summary_stats$ATE$mean_se,
              summary_stats$ATE$empirical_se,
              summary_stats$ATE$coverage,
              summary_stats$ATE$convergence_rate))
  
  # CATE components
  cat("\nCATE Components:\n")
  cat(sprintf("  Numerator (ITT on Y):   %.4f (SE: %.4f)\n",
              summary_stats$CATE$mean_numerator,
              summary_stats$CATE$mean_numerator_se))
  cat(sprintf("  Denominator (ITT on A): %.4f (SE: %.4f)\n",
              summary_stats$CATE$mean_denominator,
              summary_stats$CATE$mean_denominator_se))
  
  cat("\n")
}



# number of sims
N_SIMS <- 1000

# run simulation study with N = 1500
sim_results <- run_simulation_study_iv(n_sims = N_SIMS,n = 1500,
                                       trunc_level = 0.001,iv_strength = 0.4,
                                       crossFit = FALSE)

# summarize and print
sim_summary <- summarize_simulation_iv(sim_results)
print_simulation_summary_iv(sim_summary)





# number of sims
N_SIMS <- 100

# run simulation study with N = 1500
sim_results <- run_simulation_study_iv(n_sims = N_SIMS,n = 1500,
                                       trunc_level = 0.001,iv_strength = 0.4,
                                       crossFit = TRUE)

# summarize and print
sim_summary <- summarize_simulation_iv(sim_results)
print_simulation_summary_iv(sim_summary)









