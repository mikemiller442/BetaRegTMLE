








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


# Calculates a robust TMLE estimator for the ATT
robustcompATT <- function(outcome_mod,ps_mod,covariates,
                          exposure_name,outcome_name,
                          scale_factor = 1.0,trunc_level = 0.05,
                          max_iter = 100,tol = 1e-14) {
  
  ### calculate propensity scores and outcome predictions
  A <- covariates[[exposure_name]]
  P_A1 <- mean(A)
  pred_ps_raw <- predict(ps_mod,newdata = covariates,type = "response")
  pred_ps <- pmax(pmin(pred_ps_raw,1.0 - trunc_level),trunc_level)
  covariates_0 <- covariates
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates
  covariates_1[[exposure_name]] <- 1
  pred_y_A <- predict(outcome_mod,newdata = covariates,type = "response")
  pred_y_0 <- predict(outcome_mod,newdata = covariates_0,type = "response")
  pred_y_1 <- predict(outcome_mod,newdata = covariates_1,type = "response")
  outcome_values <- covariates[[outcome_name]]
  
  ### iterative TMLE
  for (iter in 1:max_iter) {
    
    # save previous estimates for convergence check
    psi_att_old <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    
    ### update outcome predictions
    H_0_Y <- -1.0 * pred_ps / (P_A1 * (1.0 - pred_ps))
    H_1_Y <- 1.0 / P_A1
    H_A_Y <- ifelse(A == 1,H_1_Y,H_0_Y)
    eps_y <- coef(glm(outcome_values ~ -1 + offset(qlogis(pred_y_A)) + H_A_Y,
                      family = "binomial"))
    pred_y_A <- plogis(qlogis(pred_y_A) + eps_y * H_A_Y)
    pred_y_0 <- plogis(qlogis(pred_y_0) + eps_y * H_0_Y)
    pred_y_1 <- plogis(qlogis(pred_y_1) + eps_y * H_1_Y)
    
    ### update propensity score
    psi_att <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    H_A_PS <- (1.0 / P_A1) * (pred_y_1 - pred_y_0 - psi_att)
    eps_ps <- coef(glm(A ~ -1 + offset(qlogis(pred_ps)) + H_A_PS,
                       family = "binomial"))
    pred_ps <- plogis(qlogis(pred_ps) + eps_ps * H_A_PS)
    
    # check convergence
    psi_att_new <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    if (abs(psi_att_new - psi_att_old) < tol) {
      break
    }
  }
  
  # calculate efficient influence function for ATT
  psi_att <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
  H_0_Y <- -1.0 * pred_ps / (P_A1 * (1.0 - pred_ps))
  H_1_Y <- 1.0 / P_A1
  H_A_Y <- ifelse(A == 1,H_1_Y,H_0_Y)
  D_Y <- H_A_Y * (outcome_values - pred_y_A)
  D_W <- (A / P_A1) * (pred_y_1 - pred_y_0 - psi_att)
  EIF <- D_Y + D_W
  mean_EIF <- mean(EIF)
  psi_att_scaled <- psi_att * scale_factor
  EIF_scaled <- EIF * scale_factor
  
  # calculate 95% confidence interval
  se <- sqrt(var(EIF_scaled) / length(EIF_scaled))
  ci_lower <- psi_att_scaled - qnorm(0.975) * se
  ci_upper <- psi_att_scaled + qnorm(0.975) * se
  
  ### return results
  return(list(ATT = psi_att_scaled,
              IF = EIF_scaled,
              SE = se,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Calculates a robust TMLE estimator for the ATE with missing at random outcomes
robustcompATE_MAR <- function(outcome_mod,ps_mod,missing_mod,
                              covariates_full,exposure_name,
                              outcome_observed,outcome_values,outcome_name,
                              trunc_level,scale_factor = 1.0) {
  
  n <- nrow(covariates_full)
  exposure <- covariates_full[[exposure_name]]
  
  ### calculate propensity score for exposure
  ps_pred <- predict(ps_mod,newdata = covariates_full,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  
  ### calculate probability of being observed
  prob_observed <- predict(missing_mod,newdata = covariates_full,type = "response")
  prob_observed_trunc <- pmax(pmin(prob_observed,1.0 - trunc_level),trunc_level)
  
  ### calculate observed data clever covariates for TMLE
  H_0_Y <- ifelse(exposure == 0 & outcome_observed == 1,
                  -1.0 / ((1.0 - ps_trunc) * prob_observed_trunc),
                  0)
  H_1_Y <- ifelse(exposure == 1 & outcome_observed == 1,
                  1.0 / (ps_trunc * prob_observed_trunc),
                  0)
  H_A_Y <- ifelse(exposure == 1,H_1_Y,H_0_Y)
  
  ### estimate fluctuation parameter using observed data only
  pred_y <- predict(outcome_mod,newdata = covariates_full,type = "response")
  obs_idx <- which(outcome_observed == 1)
  data_tmle <- data.frame(Y = outcome_values[obs_idx],
                          pred_vals = pred_y[obs_idx],
                          H_A_Y = H_A_Y[obs_idx])
  eps <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                  data = data_tmle,family = "binomial"))
  
  ### counterfactual clever covariates for all observations
  # need to compute probability of observation at COUNTERFACTUAL exposure levels
  covariates_0 <- covariates_full
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates_full
  covariates_1[[exposure_name]] <- 1
  
  prob_observed_a0 <- predict(missing_mod,newdata = covariates_0,type = "response")
  prob_observed_a0_trunc <- pmax(pmin(prob_observed_a0,1.0 - trunc_level),trunc_level)
  
  prob_observed_a1 <- predict(missing_mod,newdata = covariates_1,type = "response")
  prob_observed_a1_trunc <- pmax(pmin(prob_observed_a1,1.0 - trunc_level),trunc_level)
  
  H_0_cf <- -1.0 / ((1.0 - ps_trunc) * prob_observed_a0_trunc)
  H_1_cf <- 1.0 / (ps_trunc * prob_observed_a1_trunc)
  
  ### update predictions using fluctuation parameter
  covariates_0 <- covariates_full
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates_full
  covariates_1[[exposure_name]] <- 1
  
  pred_y_0 <- predict(outcome_mod,newdata = covariates_0,type = "response")
  pred_y_1 <- predict(outcome_mod,newdata = covariates_1,type = "response")
  pred_y_A <- predict(outcome_mod,newdata = covariates_full,type = "response")
  
  # use observed-data clever covariate for pred_star_y_A in first EIF component
  pred_star_y_A <- plogis(qlogis(pred_y_A) + eps * H_A_Y)
  # use counterfactual clever covariates for plug-in estimator
  pred_star_y_0 <- plogis(qlogis(pred_y_0) + eps * H_0_cf)
  pred_star_y_1 <- plogis(qlogis(pred_y_1) + eps * H_1_cf)
  
  ### calculate ATE using targeted estimates
  ate_star <- pred_star_y_1 - pred_star_y_0
  ate_star_scaled <- ate_star * scale_factor
  
  ### calculate scaled efficient influence function
  D_Y <- ifelse(outcome_observed == 1,
                H_A_Y * (outcome_values - pred_star_y_A),
                0)
  D_W <- ate_star - mean(ate_star)
  EIF <- (D_Y + D_W) * scale_factor
  
  ### calculate standard errors and confidence intervals
  var_ate <- var(EIF) / n
  se_ate <- sqrt(var_ate)
  ate_mean <- mean(ate_star_scaled)
  ci_lower <- ate_mean - qnorm(0.975) * se_ate
  ci_upper <- ate_mean + qnorm(0.975) * se_ate
  
  ### return results
  return(list(ATE = ate_mean,
              ITE = ate_star_scaled,
              IF = EIF,
              SE = se_ate,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Calculates a robust TMLE estimator for the ATT with missing outcomes
robustcompATT_MAR <- function(outcome_mod,ps_mod,missing_mod,
                              covariates_full,exposure_name,
                              outcome_observed,outcome_values,
                              outcome_name,trunc_level = 0.05,
                              scale_factor = 1.0,
                              max_iter = 100,tol = 1e-14) {
  
  ### calculate propensity scores and outcome predictions
  A <- covariates_full[[exposure_name]]
  prob_obs <- predict(missing_mod,type = "response")
  covariates_0 <- covariates_full
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates_full
  covariates_1[[exposure_name]] <- 1
  prob_obs_0 <- predict(missing_mod,newdata = covariates_0,type = "response")
  prob_obs_0_trunc <- pmax(pmin(prob_obs_0,1.0 - trunc_level),trunc_level)
  prob_obs_1 <- predict(missing_mod,newdata = covariates_1,type = "response")
  prob_obs_1_trunc <- pmax(pmin(prob_obs_1,1.0 - trunc_level),trunc_level)
  trunc_weights_missing_0 <- ifelse(outcome_observed == 1,1.0 / prob_obs_0_trunc,0)
  trunc_weights_missing_1 <- ifelse(outcome_observed == 1,1.0 / prob_obs_1_trunc,0)
  n_total <- nrow(covariates_full)
  P_A1 <- mean(A)
  pred_ps_raw <- predict(ps_mod,newdata = covariates_full,type = "response")
  pred_ps <- pmax(pmin(pred_ps_raw,1.0 - trunc_level),trunc_level)
  pred_y_A <- predict(outcome_mod,newdata = covariates_full,type = "response")
  pred_y_0 <- predict(outcome_mod,newdata = covariates_0,type = "response")
  pred_y_1 <- predict(outcome_mod,newdata = covariates_1,type = "response")
  
  # iterative TMLE
  for (iter in 1:max_iter) {
    
    # save previous estimate for convergence check
    psi_att_old <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    
    ### update outcome predictions
    H_0_Y <- -1.0 * (pred_ps * trunc_weights_missing_0) / (P_A1 * (1.0 - pred_ps))
    H_1_Y <- trunc_weights_missing_1 / P_A1
    H_A_Y <- ifelse(A == 1,H_1_Y,H_0_Y)
    obs_indices <- which(outcome_observed == 1)
    data_tmle_y <- data.frame(Y = outcome_values[obs_indices],
                              pred_vals = pred_y_A[obs_indices],
                              H_A_Y = H_A_Y[obs_indices])
    eps_y <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                      data = data_tmle_y,
                      family = "binomial"))
    pred_y_A <- plogis(qlogis(pred_y_A) + eps_y * H_A_Y)
    pred_y_0 <- plogis(qlogis(pred_y_0) + eps_y * H_0_Y)
    pred_y_1 <- plogis(qlogis(pred_y_1) + eps_y * H_1_Y)
    
    ### update propensity score
    psi_att <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    H_A_PS <- (1.0 / P_A1) * (pred_y_1 - pred_y_0 - psi_att)
    H_A_PS_observed <- H_A_PS * outcome_observed
    data_tmle_ps <- data.frame(A = A[obs_indices],
                               ps_vals = pred_ps[obs_indices],
                               H_A_PS = H_A_PS[obs_indices])
    eps_ps <- coef(glm(A ~ -1 + offset(qlogis(ps_vals)) + H_A_PS,
                       data = data_tmle_ps,
                       family = "binomial"))
    pred_ps <- plogis(qlogis(pred_ps) + eps_ps * H_A_PS)
    
    # check convergence
    psi_att_new <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
    if (abs(psi_att_new - psi_att_old) < tol) {
      break
    }
  }
  
  # calculate efficient influence function for ATT
  psi_att <- mean((pred_ps / P_A1) * (pred_y_1 - pred_y_0))
  H_0_Y <- -1.0 * (pred_ps * trunc_weights_missing_0) / (P_A1 * (1.0 - pred_ps))
  H_1_Y <- trunc_weights_missing_1 / P_A1
  H_A_Y <- ifelse(A == 1,H_1_Y,H_0_Y)
  D_Y <- outcome_observed * H_A_Y * (outcome_values - pred_y_A)
  D_Y[is.na(D_Y)] <- 0
  D_W <- (A / P_A1) * (pred_y_1 - pred_y_0 - psi_att)
  EIF <- D_Y + D_W
  mean_EIF <- mean(EIF)
  psi_att_scaled <- psi_att * scale_factor
  EIF_scaled <- EIF * scale_factor
  
  # calculate 95% confidence interval
  var_att <- var(EIF) / n_total
  se_att_scaled <- sqrt(var_att) * scale_factor
  ci_lower <- psi_att_scaled - qnorm(0.975) * se_att_scaled
  ci_upper <- psi_att_scaled + qnorm(0.975) * se_att_scaled
  
  # return results
  return(list(ATT = psi_att_scaled,
              IF = EIF_scaled,
              SE = se_att_scaled,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


### Generate simulated data
generate_stroke_data <- function(n = 1000,
                                 exposure_name = "rehabIRF",
                                 missing_outcome = FALSE,
                                 prob_observed_baseline = 0.8,
                                 fixed_covariates = NULL) {
  
  ### generate or use fixed baseline covariates
  if (is.null(fixed_covariates)) {
    # age
    age <- rnorm(n,mean = 70,sd = 12)
    age_scaled <- (age - mean(age)) / sd(age)
    
    # post-discharge disability
    # higher values = more disability
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
    age <- fixed_covariates$age
    age_scaled <- fixed_covariates$age_scaled
    post_discharge_disability <- fixed_covariates$post_discharge_disability
    post_discharge_disability_scaled <- fixed_covariates$post_discharge_disability_scaled
    stroke_severity <- fixed_covariates$stroke_severity
    stroke_severity_scaled <- fixed_covariates$stroke_severity_scaled
    comorbidity <- fixed_covariates$comorbidity
    comorbidity_scaled <- fixed_covariates$comorbidity_scaled
  }
  
  ### generate treatment assignment
  # propensity score depends on covariates
  # higher post_discharge_disability -> higher treatment probability
  logit_ps <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled 
  
  ps <- plogis(logit_ps)
  exposure <- rbinom(n,size = 1,prob = ps)
  
  ### generate outcomes using beta regression structure
  # true data generating mechanism for mean
  # higher post_discharge_disability -> higher follow-up disability
  # treatment (exposure) -> lower follow-up disability
  mu_adl_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.5 * exposure +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled +
    -0.15 * exposure * post_discharge_disability_scaled
  
  mu_adl <- plogis(mu_adl_logit)
  
  # beta distribution parameters
  phi_adl <- 15
  shape1_adl <- mu_adl * phi_adl
  shape2_adl <- (1.0 - mu_adl) * phi_adl
  
  # generate outcome on 0-1 scale
  OUT3_ADL_IADL_01 <- rbeta(n,shape1 = shape1_adl,shape2 = shape2_adl)
  
  # re-scale to 0-3 range
  OUT3_ADL_IADL <- OUT3_ADL_IADL_01 * 3.0
  
  ### generate missing outcome indicators if requested
  if (missing_outcome) {
    logit_observed_adl <- qlogis(prob_observed_baseline) + 
      0.3 * age_scaled + 
      0.2 * post_discharge_disability_scaled + 
      0.3 * stroke_severity_scaled + 
      0.15 * comorbidity_scaled +
      0.4 * exposure
    
    prob_observed_adl <- plogis(logit_observed_adl)
    observed_adl <- rbinom(n,size = 1,prob = prob_observed_adl)
    
    # set missing outcomes to NA
    OUT3_ADL_IADL_with_missing <- ifelse(observed_adl == 1,OUT3_ADL_IADL,NA)
  } else {
    # no missingness
    observed_adl <- rep(1,n)
    OUT3_ADL_IADL_with_missing <- OUT3_ADL_IADL
    prob_observed_adl <- rep(1,n)
  }
  
  ### calculate true ATEs for validation
  mu_adl_0_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 0 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled + -0.15 * 0 * post_discharge_disability_scaled
  
  mu_adl_1_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 1 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled + -0.15 * 1 * post_discharge_disability_scaled
  
  true_ate_adl <- mean(plogis(mu_adl_1_logit) - plogis(mu_adl_0_logit)) * 3.0
  
  ### calculate true ATTs for validation
  # ATT = E[E[Y | A = 1, X] - E[Y | A = 0, X] | A = 1]
  # using standardization: ATT = E[(Y(1) - Y(0)) * (f(1|X)] / E[f(1|X))]
  # where f(1|X) is the propensity score P(A=1|X)
  numerator_att <- mean((plogis(mu_adl_1_logit) - plogis(mu_adl_0_logit)) * ps)
  denominator_att <- mean(ps)
  true_att_adl <- (numerator_att / denominator_att) * 3.0
  
  ### compile dataset
  data <- data.frame(id = 1:n,
                     age = age,
                     post_discharge_disability = post_discharge_disability,
                     stroke_severity = stroke_severity,
                     comorbidity = comorbidity,
                     OUT3_ADL_IADL = OUT3_ADL_IADL_with_missing,
                     OUT3_ADL_IADL_complete = OUT3_ADL_IADL,
                     observed_adl = observed_adl,
                     prob_observed_adl = prob_observed_adl,
                     propensity_score = ps)
  
  # add exposure variable with the specified name
  data[[exposure_name]] <- exposure
  
  ### return results
  return(list(data = data,
              true_ate_adl = true_ate_adl,
              true_att_adl = true_att_adl,
              simulation_params = list(n = n,
                                       exposure_name = exposure_name,
                                       phi_adl = phi_adl,
                                       missing_outcome = missing_outcome,
                                       prob_observed_baseline = prob_observed_baseline,
                                       actual_missing_prop_adl = mean(observed_adl == 0)),
              fixed_covariates = data.frame(age = age,
                                            age_scaled = age_scaled,
                                            post_discharge_disability = post_discharge_disability,
                                            post_discharge_disability_scaled = post_discharge_disability_scaled,
                                            stroke_severity = stroke_severity,
                                            stroke_severity_scaled = stroke_severity_scaled,
                                            comorbidity = comorbidity,
                                            comorbidity_scaled = comorbidity_scaled)))
}


# Run a single simulation iteration for both ATE and ATT
run_single_simulation <- function(n = 1000,
                                  exposure_name = "rehabIRF",
                                  trunc_level = 0.05,
                                  missing_outcome = FALSE,
                                  prob_observed_baseline = 0.8,
                                  crossFit = FALSE,
                                  n_folds = 5) {
  
  # expects fixed_covariates to exist in parent environment
  sim_data <- generate_stroke_data(n = n,
                                   exposure_name = exposure_name,
                                   missing_outcome = missing_outcome,
                                   prob_observed_baseline = prob_observed_baseline,
                                   fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  
  # determine which outcomes are observed
  observed_adl <- data$observed_adl
  
  # re-scale outcomes to 0-1 for modeling
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  
  # apply transformation for beta regression to avoid boundary issues
  N_adl <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N_adl - 1) + 0.5) / N_adl
  
  ### estimate propensity score model (hardcoded formula)
  ps_model <- glm(rehabIRF ~ age + post_discharge_disability + stroke_severity + comorbidity,
                  data = data,family = binomial())
  
  if (missing_outcome) {
    # missing data model (hardcoded formula)
    missing_model_adl <- glm(observed_adl ~ rehabIRF + age + post_discharge_disability + 
                               stroke_severity + comorbidity,
                             data = data,family = binomial())
    
    ### estimate outcome models on observed data only
    data_obs_adl <- data[observed_adl == 1,]
    
    # outcome model (hardcoded formula with interaction term to match DGM)
    adl_mod <- betareg(OUT3_ADL_IADL_01 ~ age + rehabIRF * post_discharge_disability + 
                         stroke_severity + comorbidity,
                       data = data_obs_adl)
    
    ### apply TMLE for ATE with missing outcomes
    results_ate <- robustcompATE_MAR(outcome_mod = adl_mod,
                                     ps_mod = ps_model,
                                     missing_mod = missing_model_adl,
                                     covariates_full = data,
                                     exposure_name = exposure_name,
                                     outcome_observed = observed_adl,
                                     trunc_level = trunc_level,
                                     outcome_values = data$OUT3_ADL_IADL_01,
                                     outcome_name = "OUT3_ADL_IADL_01",
                                     scale_factor = 3.0)
    
    ### apply TMLE for ATT with missing outcomes
    results_att <- robustcompATT_MAR(outcome_mod = adl_mod,
                                     ps_mod = ps_model,
                                     missing_mod = missing_model_adl,
                                     covariates_full = data,
                                     exposure_name = exposure_name,
                                     outcome_observed = observed_adl,
                                     outcome_values = data$OUT3_ADL_IADL_01,
                                     outcome_name = "OUT3_ADL_IADL_01",
                                     trunc_level = trunc_level,
                                     scale_factor = 3.0,
                                     max_iter = 100,
                                     tol = 1e-14)
    
  } else {
    # outcome model (hardcoded formula with interaction term to match DGM)
    adl_mod <- betareg(OUT3_ADL_IADL_01 ~ age + rehabIRF * post_discharge_disability + 
                         stroke_severity + comorbidity,
                       data = data)
    
    ### apply TMLE for ATE
    results_ate <- robustcompATE(outcome_mod = adl_mod,
                                 ps_mod = ps_model,
                                 covariates = data,
                                 exposure_name = exposure_name,
                                 outcome_name = "OUT3_ADL_IADL_01",
                                 trunc_level = trunc_level,
                                 scale_factor = 3.0,
                                 crossFit = crossFit,
                                 n_folds = n_folds)
    
    ### apply TMLE for ATT
    results_att <- robustcompATT(outcome_mod = adl_mod,
                                 ps_mod = ps_model,
                                 covariates = data,
                                 exposure_name = exposure_name,
                                 outcome_name = "OUT3_ADL_IADL_01",
                                 trunc_level = trunc_level,
                                 scale_factor = 3.0)
  }
  
  ### extract results
  data.frame(n = n,
             exposure_name = exposure_name,
             trunc_level = trunc_level,
             missing_outcome = missing_outcome,
             crossFit = crossFit,
             # ATE results
             ate_adl = results_ate$ATE,
             ate_se_adl = results_ate$SE,
             ate_ci_lower_adl = results_ate$CI["lower"],
             ate_ci_upper_adl = results_ate$CI["upper"],
             # ATT results
             att_adl = results_att$ATT,
             att_se_adl = results_att$SE,
             att_ci_lower_adl = results_att$CI["lower"],
             att_ci_upper_adl = results_att$CI["upper"],
             # true values
             true_ate_adl = sim_data$true_ate_adl,
             true_att_adl = sim_data$true_att_adl)
}


# Run full simulation study for both ATE and ATT
run_simulation_study <- function(n_sims = 500,
                                 n = 1000,
                                 exposure_name = "rehabIRF",
                                 trunc_level = 0.05,
                                 missing_outcome = FALSE,
                                 prob_observed_baseline = 0.8,
                                 crossFit = FALSE,
                                 n_folds = 5) {
  
  # generate fixed covariates once and assign to parent environment
  initial_data <- generate_stroke_data(n = n,
                                       exposure_name = exposure_name,
                                       missing_outcome = FALSE)
  # assign covariates to parent environment
  # considerably faster than passing as an argument for each simulation
  fixed_covariates <<- initial_data$fixed_covariates
  
  results_list <- vector("list",n_sims)
  
  pb <- txtProgressBar(min = 0,max = n_sims,style = 3)
  
  for (i in 1:n_sims) {
    results_list[[i]] <- run_single_simulation(n = n,
                                               exposure_name = exposure_name,
                                               trunc_level = trunc_level,
                                               missing_outcome = missing_outcome,
                                               prob_observed_baseline = prob_observed_baseline,
                                               crossFit = crossFit,
                                               n_folds = n_folds)
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  results <- do.call(rbind,results_list)
  
  # clean up global variable
  rm(fixed_covariates,envir = .GlobalEnv)
  
  # calculate simulation performance metrics for ATE
  results$ate_bias_adl <- results$ate_adl - results$true_ate_adl
  results$ate_coverage_adl <- (results$ate_ci_lower_adl <= results$true_ate_adl) & 
    (results$ate_ci_upper_adl >= results$true_ate_adl)
  
  # calculate simulation performance metrics for ATT
  results$att_bias_adl <- results$att_adl - results$true_att_adl
  results$att_coverage_adl <- (results$att_ci_lower_adl <= results$true_att_adl) & 
    (results$att_ci_upper_adl >= results$true_att_adl)
  
  return(results)
}


# Summarize simulation results for both ATE and ATT
summarize_simulation <- function(results) {
  
  # ATE summary statistics
  ate_stats <- c(mean_estimate = mean(results$ate_adl),
                 true_effect = mean(results$true_ate_adl),
                 mean_bias = mean(results$ate_bias_adl),
                 mean_se = mean(results$ate_se_adl),
                 empirical_se = sd(results$ate_adl),
                 coverage = mean(results$ate_coverage_adl))
  
  # ATT summary statistics
  att_stats <- c(mean_estimate = mean(results$att_adl),
                 true_effect = mean(results$true_att_adl),
                 mean_bias = mean(results$att_bias_adl),
                 mean_se = mean(results$att_se_adl),
                 empirical_se = sd(results$att_adl),
                 coverage = mean(results$att_coverage_adl))
  
  # create summary table
  summary_table <- data.frame(Estimand = c("ATE","ATT"),
                              Mean_Estimate = round(c(ate_stats["mean_estimate"],att_stats["mean_estimate"]),4),
                              True_Effect = round(c(ate_stats["true_effect"],att_stats["true_effect"]),4),
                              Mean_Bias = round(c(ate_stats["mean_bias"],att_stats["mean_bias"]),4),
                              Mean_SE = round(c(ate_stats["mean_se"],att_stats["mean_se"]),4),
                              Empirical_SE = round(c(ate_stats["empirical_se"],att_stats["empirical_se"]),4),
                              Coverage = round(c(ate_stats["coverage"],att_stats["coverage"]),3))
  
  cat("\n=== Simulation Summary ===\n")
  cat(sprintf("N simulations: %d | Sample size: %d | Missing outcome: %s | Cross-fitting: %s\n\n",
              nrow(results),
              results$n[1],
              results$missing_outcome[1],
              results$crossFit[1]))
  print(summary_table,row.names = FALSE)
  
  return(summary_table)
}



### number of sims
N_SIMS <- 1500


# without cross-fitting
sim_results <- run_simulation_study(n_sims = N_SIMS,
                                    n = 1875,
                                    trunc_level = 0.005,
                                    missing_outcome = FALSE,
                                    crossFit = FALSE)
summary_results <- summarize_simulation(sim_results)

sim_results <- run_simulation_study(n_sims = N_SIMS,n = 1875,
                                    trunc_level = 0.005,
                                    missing_outcome = TRUE,
                                    crossFit = FALSE)
summary_results <- summarize_simulation(sim_results)



### number of sims
N_SIMS <- 100

# with cross-fitting
sim_results_cf <- run_simulation_study(n_sims = N_SIMS,
                                       n = 1875,
                                       trunc_level = 0.005,
                                       missing_outcome = FALSE,
                                       crossFit = TRUE,
                                       n_folds = 5)
sim_results_cf <- summarize_simulation(sim_results_cf)




# ### completely observed outcomes
# # N = 500
# sim_results <- run_simulation_study(n_sims = N_SIMS,n = 625,
#                                     trunc_level = 0.005,
#                                     missing_outcome = FALSE,
#                                     prob_observed_baseline = 0.8,
#                                     use_fixed_covariates = TRUE)
# summary_results <- summarize_simulation(sim_results)
# 
# # N = 1500
# sim_results <- run_simulation_study(n_sims = N_SIMS,n = 1875,
#                                     trunc_level = 0.005,
#                                     missing_outcome = FALSE,
#                                     prob_observed_baseline = 0.8,
#                                     use_fixed_covariates = TRUE)
# summary_results <- summarize_simulation(sim_results)
# 
# 
# 
# ### missing outcomes
# # N = 500
# sim_results <- run_simulation_study(n_sims = N_SIMS,n = 625,
#                                     trunc_level = 0.005,
#                                     missing_outcome = TRUE,
#                                     prob_observed_baseline = 0.8,
#                                     use_fixed_covariates = TRUE)
# summary_results <- summarize_simulation(sim_results)








