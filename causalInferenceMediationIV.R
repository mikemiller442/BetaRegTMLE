








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


# Calculates a robust TMLE estimator for the NDE
robustcompNDE <- function(outcome_mod,mediator_mod,ps_mod,covariates,
                          exposure_name,mediator_name,outcome_name,
                          trunc_level,scale_factor = 1.0,
                          use_regression_psi = FALSE) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  mediator <- covariates[[mediator_name]]
  outcome <- covariates[[outcome_name]]
  
  ### calculate propensity score and mediator densities
  ps_pred <- predict(ps_mod,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  covariates_m0 <- covariates
  covariates_m0[[exposure_name]] <- 0
  med_prob_0 <- predict(mediator_mod,newdata = covariates_m0,type = "response")
  covariates_m1 <- covariates
  covariates_m1[[exposure_name]] <- 1
  med_prob_1 <- predict(mediator_mod,newdata = covariates_m1,type = "response")
  med_density_obs <- ifelse(exposure == 0,
                            ifelse(mediator == 1,med_prob_0,1 - med_prob_0),
                            ifelse(mediator == 1,med_prob_1,1 - med_prob_1))
  med_density_0_obs <- ifelse(mediator == 1,med_prob_0,1 - med_prob_0)
  med_density_obs <- pmax(med_density_obs,trunc_level)
  med_density_0_obs <- pmax(med_density_0_obs,trunc_level)
  
  ### calculate targeted outcome regressions
  H_Y <- ifelse(exposure == 1,
                (med_density_0_obs / med_density_obs) / ps_trunc,
                -1 / (1 - ps_trunc))
  pred_y <- predict(outcome_mod,newdata = covariates,type = "response")
  data_tmle_y <- data.frame(Y = outcome,
                            pred_vals = pred_y,
                            H_Y = H_Y)
  
  eps_y <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_Y,
                    data = data_tmle_y,family = "binomial"))
  pred_y_star <- plogis(qlogis(pred_y) + eps_y * H_Y)
  
  ### calculate pseudo-outcome obs_diff (needed for both approaches)
  covariates_y0_obs <- covariates
  covariates_y0_obs[[exposure_name]] <- 0
  pred_y0_obs <- predict(outcome_mod,newdata = covariates_y0_obs,type = "response")
  H_Y_0_obs <- -1 / (1 - ps_trunc)
  pred_y0_obs_star <- plogis(qlogis(pred_y0_obs) + eps_y * H_Y_0_obs)
  
  covariates_y1_obs <- covariates
  covariates_y1_obs[[exposure_name]] <- 1
  pred_y1_obs <- predict(outcome_mod,newdata = covariates_y1_obs,type = "response")
  H_Y_1_obs <- ifelse(mediator == 1,
                      med_prob_0 / (med_prob_1 * ps_trunc),
                      (1 - med_prob_0) / ((1 - med_prob_1) * ps_trunc))
  pred_y1_obs_star <- plogis(qlogis(pred_y1_obs) + eps_y * H_Y_1_obs)
  obs_diff <- pred_y1_obs_star - pred_y0_obs_star
  
  ### scale obs_diff for TMLE targeting
  min_obs_diff <- min(obs_diff)
  max_obs_diff <- max(obs_diff)
  range_obs_diff <- max_obs_diff - min_obs_diff
  if (range_obs_diff < 1e-10) {
    range_obs_diff <- 1
  }
  obs_diff_scaled <- (obs_diff - min_obs_diff) / range_obs_diff
  obs_diff_scaled <- (obs_diff_scaled * (n - 1) + 0.5) / n
  
  if (use_regression_psi) {
    ### regression-based approach for NDE
    # among A=0,regress obs_diff_scaled on covariates with IPW
    data_a0 <- data.frame(obs_diff = obs_diff_scaled[exposure == 0],
                          age = covariates$age[exposure == 0],
                          post_discharge_disability = covariates$post_discharge_disability[exposure == 0],
                          stroke_severity = covariates$stroke_severity[exposure == 0],
                          comorbidity = covariates$comorbidity[exposure == 0])
    
    # weights = 1/(1-ps_trunc) for those with A=0
    weights_a0 <- 1 / (1 - ps_trunc[exposure == 0])
    
    # fit weighted linear model
    nde_reg <- lm(obs_diff ~ age + post_discharge_disability + stroke_severity + comorbidity,
                  data = data_a0,weights = weights_a0)
    
    # predict for all observations
    pred_data <- data.frame(age = covariates$age,
                            post_discharge_disability = covariates$post_discharge_disability,
                            stroke_severity = covariates$stroke_severity,
                            comorbidity = covariates$comorbidity)
    
    nde_cond_scaled <- predict(nde_reg,newdata = pred_data)
    
    # ensure predictions are in valid range for logistic targeting
    nde_cond_scaled <- pmax(pmin(nde_cond_scaled,1 - 1e-6),1e-6)
    
  } else {
    ### mediator density weighting for NDE
    # calculate NDE parameter using targeted outcome regression
    covariates_y0m0 <- covariates
    covariates_y0m0[[exposure_name]] <- 0
    covariates_y0m0[[mediator_name]] <- 0
    pred_y0m0 <- predict(outcome_mod,newdata = covariates_y0m0,type = "response")
    H_Y_0m0 <- -1 / (1 - ps_trunc)
    pred_y0m0_star <- plogis(qlogis(pred_y0m0) + eps_y * H_Y_0m0)
    covariates_y0m1 <- covariates
    covariates_y0m1[[exposure_name]] <- 0
    covariates_y0m1[[mediator_name]] <- 1
    pred_y0m1 <- predict(outcome_mod,newdata = covariates_y0m1,type = "response")
    H_Y_0m1 <- -1 / (1 - ps_trunc)
    pred_y0m1_star <- plogis(qlogis(pred_y0m1) + eps_y * H_Y_0m1)
    E_y0_m0 <- pred_y0m0_star * (1 - med_prob_0) + pred_y0m1_star * med_prob_0
    covariates_y1m0 <- covariates
    covariates_y1m0[[exposure_name]] <- 1
    covariates_y1m0[[mediator_name]] <- 0
    pred_y1m0 <- predict(outcome_mod,newdata = covariates_y1m0,type = "response")
    H_Y_1m0 <- (1 - med_prob_0) / ((1 - med_prob_1) * ps_trunc)
    pred_y1m0_star <- plogis(qlogis(pred_y1m0) + eps_y * H_Y_1m0)
    covariates_y1m1 <- covariates
    covariates_y1m1[[exposure_name]] <- 1
    covariates_y1m1[[mediator_name]] <- 1
    pred_y1m1 <- predict(outcome_mod,newdata = covariates_y1m1,type = "response")
    H_Y_1m1 <- med_prob_0 / (med_prob_1 * ps_trunc)
    pred_y1m1_star <- plogis(qlogis(pred_y1m1) + eps_y * H_Y_1m1)
    E_y1_m0 <- pred_y1m0_star * (1 - med_prob_0) + pred_y1m1_star * med_prob_0
    nde_cond <- E_y1_m0 - E_y0_m0
    
    # scale nde_cond to match obs_diff_scaled
    nde_cond_scaled <- (nde_cond - min_obs_diff) / range_obs_diff
    nde_cond_scaled <- (nde_cond_scaled * (n - 1) + 0.5) / n
  }
  
  ### TMLE targeting step (same for both approaches)
  H_Z <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
  
  data_tmle_nde <- data.frame(pseudo_outcome = obs_diff_scaled,
                              pred_vals = nde_cond_scaled,
                              H_Z = H_Z)
  # filter to only those with A=0
  data_tmle_nde_subset <- data_tmle_nde[exposure == 0,]
  # calculate fluctuation parameter
  eps_nde <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z,
                      data = data_tmle_nde_subset,
                      family = "binomial"))
  
  # update the NDE estimates using fluctuation parameter
  nde_cond_star_scaled <- plogis(qlogis(nde_cond_scaled) + eps_nde * H_Z)
  # transform back to original scale
  nde_cond_star <- ((nde_cond_star_scaled * n - 0.5) / (n - 1)) * range_obs_diff + min_obs_diff
  # scale to final outcome scale
  nde_vec <- nde_cond_star * scale_factor
  nde_est <- mean(nde_vec)
  
  ### calculate the efficient influence function
  D_Y <- H_Y * (outcome - pred_y_star)
  mean_obs_diff_A0 <- mean(obs_diff[exposure == 0])
  D_M <- ifelse(exposure == 0,
                (1 / (1 - ps_trunc)) * (obs_diff - mean_obs_diff_A0),
                0)
  obs_diff_star <- plogis(qlogis(obs_diff_scaled) + eps_nde * H_Z)
  obs_diff_star <- ((obs_diff_star * n - 0.5) / (n - 1)) * range_obs_diff + min_obs_diff
  D_Z <- H_Z * (obs_diff_star - nde_cond_star)
  D_W <- nde_vec - nde_est
  EIF <- (D_Y + D_M + D_Z + D_W) * scale_factor
  
  ### calculate standard errors and confidence intervals
  var_nde <- var(EIF) / n
  se_nde <- sqrt(var_nde)
  ci_lower <- nde_est - qnorm(0.975) * se_nde
  ci_upper <- nde_est + qnorm(0.975) * se_nde
  
  ### return results
  return(list(NDE = nde_est,
              NDE_individual = nde_vec,
              IF = EIF,
              SE = se_nde,
              CI = c(lower = ci_lower,upper = ci_upper)))
}



# Calculates a robust TMLE estimator for the NIE
robustcompNIE <- function(outcome_mod,mediator_mod,ps_mod,covariates,
                          exposure_name,mediator_name,outcome_name,
                          trunc_level,scale_factor = 1.0,
                          use_regression_psi = FALSE) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  mediator <- covariates[[mediator_name]]
  outcome <- covariates[[outcome_name]]
  
  ### calculate propensity score and mediator densities
  ps_pred <- predict(ps_mod,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  covariates_m0 <- covariates
  covariates_m0[[exposure_name]] <- 0
  med_prob_0 <- predict(mediator_mod,newdata = covariates_m0,type = "response")
  covariates_m1 <- covariates
  covariates_m1[[exposure_name]] <- 1
  med_prob_1 <- predict(mediator_mod,newdata = covariates_m1,type = "response")
  med_prob_0 <- pmax(pmin(med_prob_0,1.0 - trunc_level),trunc_level)
  med_prob_1 <- pmax(pmin(med_prob_1,1.0 - trunc_level),trunc_level)
  
  ### calculate targeted outcome regressions
  Q_Z_0 <- ifelse(mediator == 1,med_prob_0,1 - med_prob_0)
  Q_Z_1 <- ifelse(mediator == 1,med_prob_1,1 - med_prob_1)
  H_Y <- ifelse(exposure == 1,
                (1 - Q_Z_0 / Q_Z_1) / ps_trunc,
                0)
  pred_y <- predict(outcome_mod,newdata = covariates,type = "response")
  data_tmle_y <- data.frame(Y = outcome,
                            pred_vals = pred_y,
                            H_Y = H_Y)
  eps_y <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_Y,
                    data = data_tmle_y,family = "binomial"))
  pred_y_star <- plogis(qlogis(pred_y) + eps_y * H_Y)
  
  ### calculate targeted outcome at A=1 with observed M
  covariates_a1_obs <- covariates
  covariates_a1_obs[[exposure_name]] <- 1
  pred_y_a1_obs <- predict(outcome_mod,newdata = covariates_a1_obs,type = "response")
  
  Q_Z_0_obs <- ifelse(mediator == 1,med_prob_0,1 - med_prob_0)
  Q_Z_1_obs <- ifelse(mediator == 1,med_prob_1,1 - med_prob_1)
  H_Y_a1_obs <- (1 - Q_Z_0_obs / Q_Z_1_obs) / ps_trunc
  
  Q_bar_Y_1Z_star <- plogis(qlogis(pred_y_a1_obs) + eps_y * H_Y_a1_obs)
  
  if (use_regression_psi) {
    ### regression based approach for psi_NIE
    # fit regression among A=1 to get E[Q_bar | W,A=1]
    data_a1 <- data.frame(Q_bar = Q_bar_Y_1Z_star[exposure == 1],
                          age = covariates$age[exposure == 1],
                          post_discharge_disability = covariates$post_discharge_disability[exposure == 1],
                          stroke_severity = covariates$stroke_severity[exposure == 1],
                          comorbidity = covariates$comorbidity[exposure == 1])
    psi_reg_a1 <- betareg(Q_bar ~ age + post_discharge_disability + stroke_severity + comorbidity,
                          data = data_a1)
    
    # fit regression among A=0 to get E[Q_bar | W,A=0]
    data_a0 <- data.frame(Q_bar = Q_bar_Y_1Z_star[exposure == 0],
                          age = covariates$age[exposure == 0],
                          post_discharge_disability = covariates$post_discharge_disability[exposure == 0],
                          stroke_severity = covariates$stroke_severity[exposure == 0],
                          comorbidity = covariates$comorbidity[exposure == 0])
    psi_reg_a0 <- betareg(Q_bar ~ age + post_discharge_disability + stroke_severity + comorbidity,
                          data = data_a0)
    
    # predict for all observations
    pred_data <- data.frame(age = covariates$age,
                            post_discharge_disability = covariates$post_discharge_disability,
                            stroke_severity = covariates$stroke_severity,
                            comorbidity = covariates$comorbidity)
    psi_nie_z_1_reg <- predict(psi_reg_a1,newdata = pred_data,type = "response")
    psi_nie_z_0_reg <- predict(psi_reg_a0,newdata = pred_data,type = "response")
    
    # TMLE targeting for psi_nie_z_1 (using A=1 observations)
    H_Z_1 <- ifelse(exposure == 1,1 / ps_trunc,0)
    data_tmle_psi1 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star,
                                 pred_vals = psi_nie_z_1_reg,
                                 H_Z_1 = H_Z_1)
    data_tmle_psi1$pred_vals <- pmax(pmin(data_tmle_psi1$pred_vals,1 - 1e-6),1e-6)
    data_tmle_psi1_subset <- data_tmle_psi1[exposure == 1,]
    
    eps_psi_1 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_1,
                          data = data_tmle_psi1_subset,family = "binomial"))
    
    psi_nie_z_1_star <- plogis(qlogis(psi_nie_z_1_reg) + eps_psi_1 * H_Z_1)
    
    # TMLE targeting for psi_nie_z_0 (using A=0 observations)
    H_Z_0 <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
    data_tmle_psi0 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star,
                                 pred_vals = psi_nie_z_0_reg,
                                 H_Z_0 = H_Z_0)
    data_tmle_psi0$pred_vals <- pmax(pmin(data_tmle_psi0$pred_vals,1 - 1e-6),1e-6)
    data_tmle_psi0_subset <- data_tmle_psi0[exposure == 0,]
    
    eps_psi_0 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_0,
                          data = data_tmle_psi0_subset,family = "binomial"))
    
    psi_nie_z_0_star <- plogis(qlogis(psi_nie_z_0_reg) + eps_psi_0 * H_Z_0)
  } else {
    ### mediator density weighting for psi_NIE
    # calculate NIE parameter using targeted outcome regression
    covariates_y1m0 <- covariates
    covariates_y1m0[[exposure_name]] <- 1
    covariates_y1m0[[mediator_name]] <- 0
    pred_y1m0 <- predict(outcome_mod,newdata = covariates_y1m0,type = "response")
    covariates_y1m1 <- covariates
    covariates_y1m1[[exposure_name]] <- 1
    covariates_y1m1[[mediator_name]] <- 1
    pred_y1m1 <- predict(outcome_mod,newdata = covariates_y1m1,type = "response")
    Q_Z_0_m0 <- 1 - med_prob_0
    Q_Z_1_m0 <- 1 - med_prob_1
    H_Y_1m0 <- (1 - Q_Z_0_m0 / Q_Z_1_m0) / ps_trunc
    Q_Z_0_m1 <- med_prob_0
    Q_Z_1_m1 <- med_prob_1
    H_Y_1m1 <- (1 - Q_Z_0_m1 / Q_Z_1_m1) / ps_trunc
    pred_y1m0_star <- plogis(qlogis(pred_y1m0) + eps_y * H_Y_1m0)
    pred_y1m1_star <- plogis(qlogis(pred_y1m1) + eps_y * H_Y_1m1)
    psi_nie_z_0 <- pred_y1m0_star * (1 - med_prob_0) + pred_y1m1_star * med_prob_0
    psi_nie_z_1 <- pred_y1m0_star * (1 - med_prob_1) + pred_y1m1_star * med_prob_1
    
    ### calculate targeted estimate of NIE
    H_Z_1 <- ifelse(exposure == 1,1 / ps_trunc,0)
    H_Z_0 <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
    # calculate regression model for A=1
    data_tmle_a1 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star,
                               pred_vals = psi_nie_z_1,
                               H_Z_1 = H_Z_1)
    data_tmle_a1_subset <- data_tmle_a1[exposure == 1,]
    eps_nie_1 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_1,
                          data = data_tmle_a1_subset,
                          family = "binomial"))
    # calculate regression model for A=0
    data_tmle_a0 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star,
                               pred_vals = psi_nie_z_0,
                               H_Z_0 = H_Z_0)
    data_tmle_a0_subset <- data_tmle_a0[exposure == 0,]
    eps_nie_0 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_0,
                          data = data_tmle_a0_subset,
                          family = "binomial"))
    
    # calculate targeted estimates for all observations
    psi_nie_z_1_star <- plogis(qlogis(psi_nie_z_1) + eps_nie_1 * H_Z_1)
    psi_nie_z_0_star <- plogis(qlogis(psi_nie_z_0) + eps_nie_0 * H_Z_0)
  }
  
  ### calculate final NIE estimate
  nie_vec <- (psi_nie_z_1_star - psi_nie_z_0_star) * scale_factor
  nie_est <- mean(nie_vec)
  
  ### calculate the efficient influence function
  D_Y <- H_Y * (outcome - pred_y_star)
  psi_nie_z_A_star <- ifelse(exposure == 1,psi_nie_z_1_star,psi_nie_z_0_star)
  residual <- Q_bar_Y_1Z_star - psi_nie_z_A_star
  H_Z <- (2 * exposure - 1) / ifelse(exposure == 1,ps_trunc,1 - ps_trunc)
  D_M <- H_Z * residual
  D_W <- nie_vec - nie_est
  EIF <- (D_Y + D_M + D_W) * scale_factor
  
  ### calculate standard errors and confidence intervals
  var_nie <- var(EIF) / n
  se_nie <- sqrt(var_nie)
  ci_lower <- nie_est - qnorm(0.975) * se_nie
  ci_upper <- nie_est + qnorm(0.975) * se_nie
  
  ### return results
  return(list(NIE = nie_est,
              NIE_individual = nie_vec,
              IF = EIF,
              SE = se_nie,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Calculate cross-fitted model predictions
crossfit_predict <- function(formula,data,model_type = "glm",
                             family = binomial(),newdata_list = NULL,
                             folds,predict_type = "response") {
  
  n_folds <- max(folds)
  
  # if newdata_list is NULL,just predict on original data
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


# Calculates a robust TMLE estimator for the conditional ATE: E[Y|Z=1,A=1,W] - E[Y|Z=0,A=1,W]
robustcompATE_cond <- function(outcome_mod,ps_mod,cond_mod,
                               covariates_full,exposure_name,
                               condition_var,outcome_values,
                               outcome_name,trunc_level,
                               condition_value,
                               scale_factor = 1.0,
                               crossFit = FALSE,n_folds = 5,folds = NULL) {
  
  n <- nrow(covariates_full)
  exposure <- covariates_full[[exposure_name]]
  condition <- covariates_full[[condition_var]]
  Y <- outcome_values
  
  # helper function for safe truncation
  safe_trunc <- function(x,level = trunc_level) {
    x <- pmax(pmin(x,1.0 - level),level)
    return(x)
  }
  
  # create fold assignments if using cross-fitting
  if (crossFit) {
    # use provided folds or create new ones
    if (is.null(folds)) {
      folds <- sample(rep(1:n_folds,length.out = n))
    }
    covariates_full$fold <- folds
    
    # extract formulas and model information from fitted models
    ps_formula <- extract_formula(ps_mod)
    cond_formula <- extract_formula(cond_mod)
    outcome_formula <- extract_formula(outcome_mod)
    
    ps_type <- extract_model_type(ps_mod)
    cond_type <- extract_model_type(cond_mod)
    outcome_type <- extract_model_type(outcome_mod)
    
    ps_family <- extract_family(ps_mod)
    cond_family <- extract_family(cond_mod)
    outcome_family <- extract_family(outcome_mod)
  }
  
  ### calculate propensity score for exposure P(Z=1|W)
  if (crossFit) {
    ps_pred <- crossfit_predict(formula = ps_formula,
                                data = covariates_full,
                                model_type = ps_type,
                                family = ps_family,
                                newdata_list = list(original = covariates_full),
                                folds = folds,
                                predict_type = "response")
  } else {
    ps_pred <- predict(ps_mod,newdata = covariates_full,type = "response")
  }
  ps_trunc <- safe_trunc(ps_pred)
  
  ### calculate probability of conditioning event P(condition_var=1|W)
  if (crossFit) {
    prob_cond_1 <- crossfit_predict(formula = cond_formula,
                                    data = covariates_full,
                                    model_type = cond_type,
                                    family = cond_family,
                                    newdata_list = list(original = covariates_full),
                                    folds = folds,
                                    predict_type = "response")
  } else {
    prob_cond_1 <- predict(cond_mod,newdata = covariates_full,type = "response")
  }
  
  # P(condition_var = condition_value | W)
  if (condition_value == 1) {
    prob_cond <- prob_cond_1
  } else {
    prob_cond <- 1 - prob_cond_1
  }
  prob_cond_trunc <- safe_trunc(prob_cond)
  
  ### calculate clever covariates for TMLE conditional on condition_var = condition_value
  # H includes indicator I(condition == condition_value) weighted by 1/P(condition == condition_value|W)
  H_0_Y <- ifelse(exposure == 0 & condition == condition_value,
                  -1.0 / ((1.0 - ps_trunc) * prob_cond_trunc),
                  0)
  H_1_Y <- ifelse(exposure == 1 & condition == condition_value,
                  1.0 / (ps_trunc * prob_cond_trunc),
                  0)
  H_A_Y <- ifelse(exposure == 1,H_1_Y,H_0_Y)
  
  ### estimate fluctuation parameter using data where condition == condition_value only
  # Create counterfactual datasets
  covariates_0_cond <- covariates_full
  covariates_0_cond[[condition_var]] <- condition_value
  covariates_0_cond[[exposure_name]] <- 0
  
  covariates_1_cond <- covariates_full
  covariates_1_cond[[condition_var]] <- condition_value
  covariates_1_cond[[exposure_name]] <- 1
  
  covariates_A_cond <- covariates_full
  covariates_A_cond[[condition_var]] <- condition_value
  
  if (crossFit) {
    outcome_preds <- crossfit_predict(formula = outcome_formula,
                                      data = covariates_full,
                                      model_type = outcome_type,
                                      family = outcome_family,
                                      newdata_list = list(original = covariates_full,
                                                          z0_cond = covariates_0_cond,
                                                          z1_cond = covariates_1_cond,
                                                          zA_cond = covariates_A_cond),
                                      folds = folds,
                                      predict_type = "response")
    pred_y <- outcome_preds$original
    pred_y_0 <- outcome_preds$z0_cond
    pred_y_1 <- outcome_preds$z1_cond
    pred_y_A <- outcome_preds$zA_cond
  } else {
    pred_y <- predict(outcome_mod,newdata = covariates_full,type = "response")
    pred_y_0 <- predict(outcome_mod,newdata = covariates_0_cond,type = "response")
    pred_y_1 <- predict(outcome_mod,newdata = covariates_1_cond,type = "response")
    pred_y_A <- predict(outcome_mod,newdata = covariates_A_cond,type = "response")
  }
  
  cond_idx <- which(condition == condition_value)
  
  data_tmle <- data.frame(Y = Y[cond_idx],
                          pred_vals = pred_y[cond_idx],
                          H_A_Y = H_A_Y[cond_idx])
  eps <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                  data = data_tmle,family = "binomial"))
  
  ### update predictions using fluctuation parameter
  # counterfactual clever covariates for all observations
  H_0_cf <- -1.0 / ((1.0 - ps_trunc) * prob_cond_trunc)
  H_1_cf <- 1.0 / (ps_trunc * prob_cond_trunc)
  
  # use counterfactual clever covariates for the updates
  pred_star_y_A <- plogis(qlogis(pred_y_A) + eps * H_A_Y)
  pred_star_y_0 <- plogis(qlogis(pred_y_0) + eps * H_0_cf)
  pred_star_y_1 <- plogis(qlogis(pred_y_1) + eps * H_1_cf)
  
  ### calculate conditional ATE using targeted estimates
  ate_cond_star <- pred_star_y_1 - pred_star_y_0
  ate_cond_star_scaled <- ate_cond_star * scale_factor
  
  ### calculate scaled efficient influence function
  D_Y <- H_A_Y * (Y - pred_star_y_A)
  D_W <- ate_cond_star - mean(ate_cond_star)
  EIF <- (D_Y + D_W) * scale_factor
  
  ### calculate standard errors and confidence intervals
  var_ate <- var(EIF) / n
  se_ate <- sqrt(var_ate)
  ate_mean <- mean(ate_cond_star_scaled)
  ci_lower <- ate_mean - qnorm(0.975) * se_ate
  ci_upper <- ate_mean + qnorm(0.975) * se_ate
  
  ### return results
  return(list(ATE = ate_mean,
              ITE = ate_cond_star_scaled,
              IF = EIF,
              SE = se_ate,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Calculates theta(a_prime, a_star) using one-step estimator
# Uses plug-in estimates plus EIF correction
# Supports cross-fitting for nuisance parameter estimation
robustTheta_reduced_form <- function(outcome_mod,mediator_mod,iv_mod,ps_mod,covariates,
                                     exposure_name,mediator_name,iv_name,outcome_name,
                                     a_prime,a_star,
                                     trunc_level,scale_factor = 1.0,
                                     crossFit = FALSE,n_folds = 5,folds = NULL) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  mediator <- covariates[[mediator_name]]
  L <- covariates[[iv_name]]
  outcome <- covariates[[outcome_name]]
  
  # helper function for safe truncation
  safe_trunc <- function(x,level = trunc_level) {
    x <- pmax(pmin(x,1.0 - level),level)
    return(x)
  }
  
  # create fold assignments if using cross-fitting
  if (crossFit) {
    # use provided folds or create new ones
    if (is.null(folds)) {
      folds <- sample(rep(1:n_folds,length.out = n))
    }
    covariates$fold <- folds
    
    # extract formulas and model information from fitted models
    ps_formula <- extract_formula(ps_mod)
    iv_formula <- extract_formula(iv_mod)
    mediator_formula <- extract_formula(mediator_mod)
    outcome_formula <- extract_formula(outcome_mod)
    
    ps_type <- extract_model_type(ps_mod)
    iv_type <- extract_model_type(iv_mod)
    mediator_type <- extract_model_type(mediator_mod)
    outcome_type <- extract_model_type(outcome_mod)
    
    ps_family <- extract_family(ps_mod)
    iv_family <- extract_family(iv_mod)
    mediator_family <- extract_family(mediator_mod)
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
    ps_pred <- predict(ps_mod,type = "response")
  }
  ps_trunc <- safe_trunc(ps_pred)
  
  # P(A = a_prime | W) and P(A = a_star | W)
  ps_a_prime <- ifelse(a_prime == 1,ps_trunc,1 - ps_trunc)
  ps_a_star <- ifelse(a_star == 1,ps_trunc,1 - ps_trunc)
  
  ### calculate P(L=1|W) for instrument density
  if (crossFit) {
    pL1_given_W <- crossfit_predict(formula = iv_formula,
                                    data = covariates,
                                    model_type = iv_type,
                                    family = iv_family,
                                    newdata_list = list(original = covariates),
                                    folds = folds,
                                    predict_type = "response")
  } else {
    pL1_given_W <- predict(iv_mod,newdata = covariates,type = "response")
  }
  pL1_given_W <- safe_trunc(pL1_given_W)
  
  ### calculate mediator density predictions P(M=1|A,W)
  # create counterfactual datasets
  covariates_a0 <- covariates
  covariates_a0[[exposure_name]] <- 0
  
  covariates_a1 <- covariates
  covariates_a1[[exposure_name]] <- 1
  
  if (crossFit) {
    mediator_preds <- crossfit_predict(formula = mediator_formula,
                                       data = covariates,
                                       model_type = mediator_type,
                                       family = mediator_family,
                                       newdata_list = list(a0 = covariates_a0,
                                                           a1 = covariates_a1),
                                       folds = folds,
                                       predict_type = "response")
    pM1_given_a0 <- mediator_preds$a0
    pM1_given_a1 <- mediator_preds$a1
  } else {
    pM1_given_a0 <- predict(mediator_mod,newdata = covariates_a0,type = "response")
    pM1_given_a1 <- predict(mediator_mod,newdata = covariates_a1,type = "response")
  }
  
  # P(M=1|A=a_star, W) - this is γ(1|Z=z*,W) in the EIF
  gamma <- ifelse(a_star == 1,pM1_given_a1,pM1_given_a0)
  
  ### calculate outcome predictions Q(A,L,W)
  # create counterfactual datasets for outcome
  # Q(A=a_prime, L=0, W) = μ(0, z', w)
  covariates_a_prime_L0 <- covariates
  covariates_a_prime_L0[[exposure_name]] <- a_prime
  covariates_a_prime_L0[[iv_name]] <- 0
  
  # Q(A=a_prime, L=1, W) = μ(1, z', w)
  covariates_a_prime_L1 <- covariates
  covariates_a_prime_L1[[exposure_name]] <- a_prime
  covariates_a_prime_L1[[iv_name]] <- 1
  
  if (crossFit) {
    outcome_preds <- crossfit_predict(formula = outcome_formula,
                                      data = covariates,
                                      model_type = outcome_type,
                                      family = outcome_family,
                                      newdata_list = list(a_prime_L0 = covariates_a_prime_L0,
                                                          a_prime_L1 = covariates_a_prime_L1),
                                      folds = folds,
                                      predict_type = "response")
    Q_L0 <- outcome_preds$a_prime_L0
    Q_L1 <- outcome_preds$a_prime_L1
  } else {
    Q_L0 <- predict(outcome_mod,newdata = covariates_a_prime_L0,type = "response")
    Q_L1 <- predict(outcome_mod,newdata = covariates_a_prime_L1,type = "response")
  }
  
  ### calculate plug-in estimate of theta
  # f_theta(z',z*)(w) = sum_m mu(m,z',w) * gamma(m|z*,w) = gamma * Q_L1 + (1-gamma) * Q_L0
  theta_individual <- gamma * Q_L1 + (1 - gamma) * Q_L0
  
  ### calculate EIF components
  # D_Y = I(z=z') * gamma(L|z*,w) / [P(z'|w) * P(L|w)] * (Y - mu(l,z',w))
  # Uses outcome model evaluated at observed INSTRUMENT L
  Q_at_L <- ifelse(L == 1, Q_L1, Q_L0)  # mu(l, z', w)
  pL_at_obs <- ifelse(L == 1, pL1_given_W, 1 - pL1_given_W)  # P(L=l|w)
  gamma_at_L <- ifelse(L == 1, gamma, 1 - gamma)  # gamma(L=l|z*,w)
  density_ratio <- gamma_at_L / pL_at_obs
  
  H_Y <- ifelse(exposure == a_prime,
                density_ratio / ps_a_prime,
                0)
  D_Y <- H_Y * (outcome - Q_at_L)
  
  # D_M = I(z=z*) / P(z*|w) * (mu(m,z',w) - f_theta(z',z*)(w))
  # Uses outcome model evaluated at observed MEDIATOR M (substituted for L)
  Q_at_M <- ifelse(mediator == 1, Q_L1, Q_L0)  # mu(m, z', w)
  
  H_M <- ifelse(exposure == a_star,
                1 / ps_a_star,
                0)
  D_M <- H_M * (Q_at_M - theta_individual)
  
  ### one-step estimator
  theta_est_unscaled <- mean(D_Y + D_M + theta_individual)
  theta_est <- theta_est_unscaled * scale_factor
  
  # D_W: centering term (for variance calculation)
  D_W <- theta_individual - mean(theta_individual)
  
  # full EIF (centered)
  EIF <- (D_Y + D_M + D_W) * scale_factor
  
  ### calculate standard errors
  var_theta <- var(EIF,na.rm = TRUE) / n
  se_theta <- sqrt(var_theta)
  ci_lower <- theta_est - qnorm(0.975) * se_theta
  ci_upper <- theta_est + qnorm(0.975) * se_theta
  
  return(list(theta = theta_est,
              theta_individual = theta_individual * scale_factor,
              SE = se_theta,
              CI = c(lower = ci_lower,upper = ci_upper),
              IF = EIF,
              D_Y = D_Y * scale_factor,
              D_M = D_M * scale_factor,
              D_W = D_W * scale_factor))
}


# Calculates the complier interventional indirect and direct effects
# IIE = (theta(1,1) - theta(1,0)) / first_stage
# IDE = (theta(1,1) - theta(0,0)) - IIE = TE_reduced - IIE
# IDE_alt = (theta(1,0) - theta(0,0)) / first_stage
robustcompCMA_IV <- function(outcome_mod, mediator_mod, mediator_mod_with_iv,
                             iv_mod, ps_mod, covariates,
                             exposure_name, mediator_name, iv_name, outcome_name,
                             trunc_level_num, trunc_level_denom, scale_factor = 1.0,
                             crossFit = FALSE, n_folds = 5) {
  
  n <- nrow(covariates)
  
  # create shared fold assignments if using cross-fitting
  if (crossFit) {
    folds <- sample(rep(1:n_folds, length.out = n))
  } else {
    folds <- NULL
  }
  
  ### calculate theta(1,1)
  theta_11 <- robustTheta_reduced_form(outcome_mod = outcome_mod,
                                       mediator_mod = mediator_mod,
                                       iv_mod = iv_mod,
                                       ps_mod = ps_mod,
                                       covariates = covariates,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       iv_name = iv_name,
                                       outcome_name = outcome_name,
                                       a_prime = 1,
                                       a_star = 1,
                                       trunc_level = trunc_level_num,
                                       scale_factor = scale_factor,
                                       crossFit = crossFit,
                                       n_folds = n_folds,
                                       folds = folds)
  
  ### calculate theta(1,0)
  theta_10 <- robustTheta_reduced_form(outcome_mod = outcome_mod,
                                       mediator_mod = mediator_mod,
                                       iv_mod = iv_mod,
                                       ps_mod = ps_mod,
                                       covariates = covariates,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       iv_name = iv_name,
                                       outcome_name = outcome_name,
                                       a_prime = 1,
                                       a_star = 0,
                                       trunc_level = trunc_level_num,
                                       scale_factor = scale_factor,
                                       crossFit = crossFit,
                                       n_folds = n_folds,
                                       folds = folds)
  
  ### calculate theta(0,0)
  theta_00 <- robustTheta_reduced_form(outcome_mod = outcome_mod,
                                       mediator_mod = mediator_mod,
                                       iv_mod = iv_mod,
                                       ps_mod = ps_mod,
                                       covariates = covariates,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       iv_name = iv_name,
                                       outcome_name = outcome_name,
                                       a_prime = 0,
                                       a_star = 0,
                                       trunc_level = trunc_level_num,
                                       scale_factor = scale_factor,
                                       crossFit = crossFit,
                                       n_folds = n_folds,
                                       folds = folds)
  
  ### calculate first stage: E[M|Z=1,A=1,W] - E[M|Z=0,A=1,W]
  first_stage <- robustcompATE_cond(outcome_mod = mediator_mod_with_iv,
                                    ps_mod = iv_mod,
                                    cond_mod = ps_mod,
                                    covariates_full = covariates,
                                    exposure_name = iv_name,
                                    condition_var = exposure_name,
                                    outcome_values = covariates[[mediator_name]],
                                    outcome_name = mediator_name,
                                    trunc_level = trunc_level_denom,
                                    condition_value = 1,
                                    scale_factor = 1.0,
                                    crossFit = crossFit,
                                    n_folds = n_folds,
                                    folds = folds)
  
  ### Calculate reduced form components
  # IIE_reduced = theta(1,1) - theta(1,0)
  iie_reduced <- theta_11$theta - theta_10$theta
  IF_iie_reduced <- theta_11$IF - theta_10$IF
  
  # TE_reduced = theta(1,1) - theta(0,0)
  te_reduced <- theta_11$theta - theta_00$theta
  IF_te_reduced <- theta_11$IF - theta_00$IF
  
  # IDE_reduced = theta(1,0) - theta(0,0) (for alternative IDE estimator)
  ide_reduced <- theta_10$theta - theta_00$theta
  IF_ide_reduced <- theta_10$IF - theta_00$IF
  
  ### Standard errors for reduced form estimates
  var_iie_reduced <- var(IF_iie_reduced) / n
  se_iie_reduced <- sqrt(var_iie_reduced)
  
  var_te_reduced <- var(IF_te_reduced) / n
  se_te_reduced <- sqrt(var_te_reduced)
  
  var_ide_reduced <- var(IF_ide_reduced) / n
  se_ide_reduced <- sqrt(var_ide_reduced)
  
  ### Calculate IIE = IIE_reduced / first_stage
  FS <- first_stage$ATE
  IF_FS <- first_stage$IF
  iie_iv <- iie_reduced / FS
  
  # Delta method IF for IIE
  IF_iie <- (IF_iie_reduced - iie_iv * IF_FS) / FS
  
  # Standard error and CI for IIE
  var_iie <- var(IF_iie) / n
  se_iie <- sqrt(var_iie)
  ci_lower_iie <- iie_iv - qnorm(0.975) * se_iie
  ci_upper_iie <- iie_iv + qnorm(0.975) * se_iie
  
  ### Calculate IDE = TE_reduced - IIE (original estimator)
  ide_iv <- te_reduced - iie_iv
  
  # Delta method IF for IDE = IF_TE_reduced - IF_IIE
  IF_ide <- IF_te_reduced - IF_iie
  
  # Standard error and CI for IDE
  var_ide <- var(IF_ide) / n
  se_ide <- sqrt(var_ide)
  ci_lower_ide <- ide_iv - qnorm(0.975) * se_ide
  ci_upper_ide <- ide_iv + qnorm(0.975) * se_ide
  
  ### Calculate IDE_alt = IDE_reduced / first_stage (alternative estimator)
  ide_alt_iv <- ide_reduced / FS
  
  # Delta method IF for IDE_alt
  IF_ide_alt <- (IF_ide_reduced - ide_alt_iv * IF_FS) / FS
  
  # Standard error and CI for IDE_alt
  var_ide_alt <- var(IF_ide_alt) / n
  se_ide_alt <- sqrt(var_ide_alt)
  ci_lower_ide_alt <- ide_alt_iv - qnorm(0.975) * se_ide_alt
  ci_upper_ide_alt <- ide_alt_iv + qnorm(0.975) * se_ide_alt
  
  return(list(
    # Theta estimates
    theta_11 = theta_11$theta,
    theta_10 = theta_10$theta,
    theta_00 = theta_00$theta,
    theta_11_SE = theta_11$SE,
    theta_10_SE = theta_10$SE,
    theta_00_SE = theta_00$SE,
    theta_11_CI = theta_11$CI,
    theta_10_CI = theta_10$CI,
    theta_00_CI = theta_00$CI,
    IF_theta_11 = theta_11$IF,
    IF_theta_10 = theta_10$IF,
    IF_theta_00 = theta_00$IF,
    # Reduced form estimates
    IIE_reduced = iie_reduced,
    IIE_reduced_SE = se_iie_reduced,
    TE_reduced = te_reduced,
    TE_reduced_SE = se_te_reduced,
    IDE_reduced = ide_reduced,
    IDE_reduced_SE = se_ide_reduced,
    IF_IIE_reduced = IF_iie_reduced,
    IF_TE_reduced = IF_te_reduced,
    IF_IDE_reduced = IF_ide_reduced,
    # First stage
    first_stage = FS,
    first_stage_SE = first_stage$SE,
    first_stage_CI = first_stage$CI,
    IF_first_stage = IF_FS,
    # IIE (IV-adjusted)
    IIE = iie_iv,
    IIE_SE = se_iie,
    IIE_CI = c(lower = ci_lower_iie, upper = ci_upper_iie),
    IF_IIE = IF_iie,
    # IDE (IV-adjusted, original: TE_reduced - IIE)
    IDE = ide_iv,
    IDE_SE = se_ide,
    IDE_CI = c(lower = ci_lower_ide, upper = ci_upper_ide),
    IF_IDE = IF_ide,
    # IDE_alt (IV-adjusted, alternative: IDE_reduced / FS)
    IDE_alt = ide_alt_iv,
    IDE_alt_SE = se_ide_alt,
    IDE_alt_CI = c(lower = ci_lower_ide_alt, upper = ci_upper_ide_alt),
    IF_IDE_alt = IF_ide_alt
  ))
}


### Generate simulated data with IV-mediation structure
generate_stroke_data_iv_mediation <- function(n = 1000,
                                              exposure_name = "insured",
                                              mediator_name = "rehabIRF",
                                              iv_name = "Z",
                                              iv_strength = 0.10,
                                              iv_strength_uninsured_ratio = 0.3,
                                              fixed_covariates = NULL) {
  
  ### generate or use fixed baseline covariates
  if (is.null(fixed_covariates)) {
    age <- rnorm(n,mean = 70,sd = 12)
    age_scaled <- (age - mean(age)) / sd(age)
    
    post_discharge_disability <- rbeta(n,shape1 = 2,shape2 = 3) * 100
    post_discharge_disability_scaled <- (post_discharge_disability - mean(post_discharge_disability)) / 
      sd(post_discharge_disability)
    
    stroke_severity <- rpois(n,lambda = 8)
    stroke_severity <- pmin(stroke_severity,42)
    stroke_severity_scaled <- (stroke_severity - mean(stroke_severity)) / sd(stroke_severity)
    
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
  
  ### generate exposure
  logit_ps_exposure <- 0.85 +
    0.05 * age_scaled + 
    0.06 * stroke_severity_scaled + 
    -0.04 * comorbidity_scaled
  
  ps_exposure <- plogis(logit_ps_exposure)
  insured <- rbinom(n,size = 1,prob = ps_exposure)
  
  ### generate instrumental variable
  logit_iv <- 0.0 + 
    0.05 * age_scaled + 
    0.05 * stroke_severity_scaled - 
    0.05 * comorbidity_scaled
  
  iv_prob <- plogis(logit_iv)
  IV <- rbinom(n,size = 1,prob = iv_prob)
  
  ### calculate IV strength components
  iv_strength_uninsured <- iv_strength * iv_strength_uninsured_ratio
  iv_interaction_effect <- iv_strength - iv_strength_uninsured
  
  mediator_intercept_base <- 0.32 - 0.4 * iv_strength_uninsured
  
  ps_mediator <- mediator_intercept_base +
    iv_strength_uninsured * IV +
    iv_interaction_effect * IV * insured +
    0.10 * insured +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.005 * insured * age_scaled +
    0.003 * insured * stroke_severity_scaled +
    -0.003 * insured * comorbidity_scaled
  
  if (any(ps_mediator < 0 | ps_mediator > 1)) {
    stop(paste("Invalid probabilities in mediator model. Range:",
               round(min(ps_mediator),3),"to",round(max(ps_mediator),3),
               "\nAdjust coefficients or iv_strength."))
  }
  
  rehabIRF <- rbinom(n,size = 1,prob = ps_mediator)
  
  true_iv_ame_on_mediator <- mean(iv_strength_uninsured * (1 - insured) + iv_strength * insured)
  
  ps_m_z1_a1 <- mediator_intercept_base +
    iv_strength_uninsured * 1 +
    iv_interaction_effect * 1 * 1 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 1 +
    0.005 * 1 * age_scaled +
    0.003 * 1 * stroke_severity_scaled +
    -0.003 * 1 * comorbidity_scaled
  
  ps_m_z0_a1 <- mediator_intercept_base +
    iv_strength_uninsured * 0 +
    iv_interaction_effect * 0 * 1 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 1 +
    0.005 * 1 * age_scaled +
    0.003 * 1 * stroke_severity_scaled +
    -0.003 * 1 * comorbidity_scaled
  
  ps_m_z1_a0 <- mediator_intercept_base +
    iv_strength_uninsured * 1 +
    iv_interaction_effect * 1 * 0 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 0 +
    0.005 * 0 * age_scaled +
    0.003 * 0 * stroke_severity_scaled +
    -0.003 * 0 * comorbidity_scaled
  
  ps_m_z0_a0 <- mediator_intercept_base +
    iv_strength_uninsured * 0 +
    iv_interaction_effect * 0 * 0 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 0 +
    0.005 * 0 * age_scaled +
    0.003 * 0 * stroke_severity_scaled +
    -0.003 * 0 * comorbidity_scaled
  
  true_first_stage <- mean(ps_m_z1_a1 - ps_m_z0_a1)
  true_first_stage_uninsured <- mean(ps_m_z1_a0 - ps_m_z0_a0)
  true_first_stage_insured <- mean(ps_m_z1_a1 - ps_m_z0_a1)
  
  ### generate outcome using logit link for mean
  logit_mu_adl <- 0.20 + 
    0.80 * post_discharge_disability_scaled +
    -0.25 * insured +
    -0.50 * rehabIRF +
    -0.15 * insured * rehabIRF +
    0.12 * stroke_severity_scaled + 
    0.08 * age_scaled + 
    0.04 * comorbidity_scaled
  
  mu_adl <- plogis(logit_mu_adl)
  
  phi_adl <- 15
  shape1_adl <- mu_adl * phi_adl
  shape2_adl <- (1.0 - mu_adl) * phi_adl
  
  OUT3_ADL_IADL_01 <- rbeta(n,shape1 = shape1_adl,shape2 = shape2_adl)
  OUT3_ADL_IADL <- OUT3_ADL_IADL_01 * 3.0
  
  ### calculate true causal effects using logit link
  logit_mu_00 <- 0.20 + 
    0.80 * post_discharge_disability_scaled +
    -0.25 * 0 + -0.50 * 0 + -0.15 * 0 * 0 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_00 <- plogis(logit_mu_00) * 3.0
  
  logit_mu_01 <- 0.20 + 
    0.80 * post_discharge_disability_scaled +
    -0.25 * 0 + -0.50 * 1 + -0.15 * 0 * 1 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_01 <- plogis(logit_mu_01) * 3.0
  
  logit_mu_10 <- 0.20 + 
    0.80 * post_discharge_disability_scaled +
    -0.25 * 1 + -0.50 * 0 + -0.15 * 1 * 0 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_10 <- plogis(logit_mu_10) * 3.0
  
  logit_mu_11 <- 0.20 + 
    0.80 * post_discharge_disability_scaled +
    -0.25 * 1 + -0.50 * 1 + -0.15 * 1 * 1 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_11 <- plogis(logit_mu_11) * 3.0
  
  ### calculate mediator propensity scores under each exposure level
  ps_m0_natural <- mediator_intercept_base +
    iv_strength_uninsured * IV +
    iv_interaction_effect * IV * 0 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 0 +
    0.005 * 0 * age_scaled +
    0.003 * 0 * stroke_severity_scaled +
    -0.003 * 0 * comorbidity_scaled
  
  ps_m1_natural <- mediator_intercept_base +
    iv_strength_uninsured * IV +
    iv_interaction_effect * IV * 1 +
    0.01 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.10 * 1 +
    0.005 * 1 * age_scaled +
    0.003 * 1 * stroke_severity_scaled +
    -0.003 * 1 * comorbidity_scaled
  
  if (any(ps_m0_natural < 0 | ps_m0_natural > 1)) {
    stop(paste("Invalid natural mediator probabilities (A=0). Range:",
               round(min(ps_m0_natural),3),"to",round(max(ps_m0_natural),3),
               "\nAdjust coefficients in mediator model."))
  }
  
  if (any(ps_m1_natural < 0 | ps_m1_natural > 1)) {
    stop(paste("Invalid natural mediator probabilities (A=1). Range:",
               round(min(ps_m1_natural),3),"to",round(max(ps_m1_natural),3),
               "\nAdjust coefficients in mediator model."))
  }
  
  # total effect
  true_te <- mean((mu_11 * ps_m1_natural + mu_10 * (1 - ps_m1_natural)) - 
                    (mu_01 * ps_m0_natural + mu_00 * (1 - ps_m0_natural)))
  
  # natural direct effect
  true_nde <- mean((mu_11 * ps_m0_natural + mu_10 * (1 - ps_m0_natural)) - 
                     (mu_01 * ps_m0_natural + mu_00 * (1 - ps_m0_natural)))
  
  # natural indirect effect
  true_nie <- mean((mu_11 * ps_m1_natural + mu_10 * (1 - ps_m1_natural)) - 
                     (mu_11 * ps_m0_natural + mu_10 * (1 - ps_m0_natural)))
  
  true_prop_mediated <- true_nie / true_te
  
  ### compile dataset
  data <- data.frame(id = 1:n,
                     age = age,
                     stroke_severity = stroke_severity,
                     comorbidity = comorbidity,
                     OUT3_ADL_IADL = OUT3_ADL_IADL,
                     propensity_score_exposure = ps_exposure,
                     propensity_score_mediator = ps_mediator,
                     propensity_score_mediator_natural_a0 = ps_m0_natural,
                     propensity_score_mediator_natural_a1 = ps_m1_natural,
                     iv_prob = iv_prob)
  
  data[[exposure_name]] <- insured
  data[[mediator_name]] <- rehabIRF
  data[[iv_name]] <- IV
  
  return(list(data = data,
              true_te = true_te,
              true_nde = true_nde,
              true_nie = true_nie,
              true_prop_mediated = true_prop_mediated,
              true_iv_ame_on_mediator = true_iv_ame_on_mediator,
              true_first_stage = true_first_stage,
              true_first_stage_insured = true_first_stage_insured,
              true_first_stage_uninsured = true_first_stage_uninsured,
              iv_strength_insured = iv_strength,
              iv_strength_uninsured = iv_strength_uninsured,
              simulation_params = list(n = n,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       iv_name = iv_name,
                                       phi_adl = phi_adl,
                                       iv_strength = iv_strength,
                                       iv_strength_uninsured_ratio = iv_strength_uninsured_ratio,
                                       iv_strength_uninsured = iv_strength_uninsured,
                                       iv_interaction_effect = iv_interaction_effect,
                                       mediator_intercept = mediator_intercept_base,
                                       prop_insured = mean(insured),
                                       prop_rehabIRF = mean(rehabIRF),
                                       prop_IV1 = mean(IV),
                                       mean_ps_mediator = mean(ps_mediator),
                                       range_ps_mediator = range(ps_mediator),
                                       mean_ps_mediator_natural_a0 = mean(ps_m0_natural),
                                       mean_ps_mediator_natural_a1 = mean(ps_m1_natural)),
              fixed_covariates = data.frame(age = age,
                                            age_scaled = age_scaled,
                                            post_discharge_disability = post_discharge_disability,
                                            post_discharge_disability_scaled = post_discharge_disability_scaled,
                                            stroke_severity = stroke_severity,
                                            stroke_severity_scaled = stroke_severity_scaled,
                                            comorbidity = comorbidity,
                                            comorbidity_scaled = comorbidity_scaled)))
}


# Run a single simulation iteration
run_single_simulation_iv_mediation <- function(n = 1000,
                                               exposure_name = "insured",
                                               mediator_name = "rehabIRF",
                                               iv_name = "Z",
                                               crossFit = FALSE,
                                               n_folds = 5,
                                               iv_strength = 0.3,
                                               trunc_level = 0.001) {
  
  sim_data <- generate_stroke_data_iv_mediation(n = n,
                                                exposure_name = exposure_name,
                                                mediator_name = mediator_name,
                                                iv_name = iv_name,
                                                iv_strength = iv_strength,
                                                fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  N_adl <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N_adl - 1) + 0.5) / N_adl
  
  # Initialize IV-based results
  iie_adl <- NA; iie_se_adl <- NA; iie_ci_lower_adl <- NA; iie_ci_upper_adl <- NA
  ide_adl <- NA; ide_se_adl <- NA; ide_ci_lower_adl <- NA; ide_ci_upper_adl <- NA
  ide_alt_adl <- NA; ide_alt_se_adl <- NA; ide_alt_ci_lower_adl <- NA; ide_alt_ci_upper_adl <- NA
  iie_reduced <- NA; iie_reduced_se <- NA
  te_reduced <- NA; te_reduced_se <- NA
  ide_reduced <- NA; ide_reduced_se <- NA
  first_stage_est <- NA; first_stage_se <- NA; first_stage_ci_lower <- NA; first_stage_ci_upper <- NA
  theta_11 <- NA; theta_11_se <- NA
  theta_10 <- NA; theta_10_se <- NA
  theta_00 <- NA; theta_00_se <- NA
  cma_converged <- FALSE
  
  # Initialize naive results
  naive_nie_adl <- NA; naive_nie_se_adl <- NA; naive_nie_ci_lower_adl <- NA; naive_nie_ci_upper_adl <- NA
  naive_nie_converged <- FALSE
  naive_nde_adl <- NA; naive_nde_se_adl <- NA; naive_nde_ci_lower_adl <- NA; naive_nde_ci_upper_adl <- NA
  naive_nde_converged <- FALSE
  
  ### Fit models
  ps_mod <- glm(insured ~ age+stroke_severity+comorbidity,
                data = data, family = binomial())
  
  mediator_mod <- glm(rehabIRF ~ age+stroke_severity+comorbidity
                      +insured*(age+stroke_severity+comorbidity),
                      data = data, family = binomial())
  
  mediator_mod_with_iv <- glm(rehabIRF ~ age+stroke_severity+comorbidity
                              +insured*(Z+age+stroke_severity+comorbidity),
                              data = data, family = binomial())
  
  iv_mod <- glm(Z ~ age+stroke_severity+comorbidity,
                data = data, family = binomial())
  
  outcome_mod <- betareg(OUT3_ADL_IADL_01 ~ age+stroke_severity+comorbidity
                         +insured*(Z+age+stroke_severity+comorbidity),
                         data = data)
  
  outcome_mod_with_med <- betareg(OUT3_ADL_IADL_01 ~ age+stroke_severity+comorbidity
                                  +insured*(rehabIRF+age+stroke_severity+comorbidity),
                                  data = data)
  
  ### Calculate IV-based IIE and IDE using combined function
  tryCatch({
    results_cma <- robustcompCMA_IV(outcome_mod = outcome_mod,
                                    mediator_mod = mediator_mod,
                                    mediator_mod_with_iv = mediator_mod_with_iv,
                                    iv_mod = iv_mod,
                                    ps_mod = ps_mod,
                                    covariates = data,
                                    exposure_name = exposure_name,
                                    mediator_name = mediator_name,
                                    iv_name = iv_name,
                                    outcome_name = "OUT3_ADL_IADL_01",
                                    trunc_level_num = trunc_level,
                                    trunc_level_denom = trunc_level,
                                    scale_factor = 3.0,
                                    crossFit = crossFit,
                                    n_folds = n_folds)
    
    # Extract IIE results
    iie_adl <- results_cma$IIE
    iie_se_adl <- results_cma$IIE_SE
    iie_ci_lower_adl <- results_cma$IIE_CI["lower"]
    iie_ci_upper_adl <- results_cma$IIE_CI["upper"]
    iie_reduced <- results_cma$IIE_reduced
    iie_reduced_se <- results_cma$IIE_reduced_SE
    
    # Extract IDE results
    ide_adl <- results_cma$IDE
    ide_se_adl <- results_cma$IDE_SE
    ide_ci_lower_adl <- results_cma$IDE_CI["lower"]
    ide_ci_upper_adl <- results_cma$IDE_CI["upper"]
    te_reduced <- results_cma$TE_reduced
    te_reduced_se <- results_cma$TE_reduced_SE
    
    # Extract IDE_alt results
    ide_alt_adl <- results_cma$IDE_alt
    ide_alt_se_adl <- results_cma$IDE_alt_SE
    ide_alt_ci_lower_adl <- results_cma$IDE_alt_CI["lower"]
    ide_alt_ci_upper_adl <- results_cma$IDE_alt_CI["upper"]
    ide_reduced <- results_cma$IDE_reduced
    ide_reduced_se <- results_cma$IDE_reduced_SE
    
    # Extract first stage results
    first_stage_est <- results_cma$first_stage
    first_stage_se <- results_cma$first_stage_SE
    first_stage_ci_lower <- results_cma$first_stage_CI["lower"]
    first_stage_ci_upper <- results_cma$first_stage_CI["upper"]
    
    # Extract theta results
    theta_11 <- results_cma$theta_11
    theta_11_se <- results_cma$theta_11_SE
    theta_10 <- results_cma$theta_10
    theta_10_se <- results_cma$theta_10_SE
    theta_00 <- results_cma$theta_00
    theta_00_se <- results_cma$theta_00_SE
    
    cma_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("CMA estimation failed:", e$message))
  })
  
  ### Calculate naive NIE
  tryCatch({
    results_naive_nie <- robustcompNIE(outcome_mod = outcome_mod_with_med,
                                       mediator_mod = mediator_mod,
                                       ps_mod = ps_mod,
                                       covariates = data,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       outcome_name = "OUT3_ADL_IADL_01",
                                       trunc_level = trunc_level,
                                       scale_factor = 3.0)
    
    naive_nie_adl <- results_naive_nie$NIE
    naive_nie_se_adl <- results_naive_nie$SE
    naive_nie_ci_lower_adl <- results_naive_nie$CI["lower"]
    naive_nie_ci_upper_adl <- results_naive_nie$CI["upper"]
    naive_nie_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("Naive NIE estimation failed:", e$message))
  })
  
  ### Calculate naive NDE
  tryCatch({
    results_naive_nde <- robustcompNDE(outcome_mod = outcome_mod_with_med,
                                       mediator_mod = mediator_mod,
                                       ps_mod = ps_mod,
                                       covariates = data,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       outcome_name = "OUT3_ADL_IADL_01",
                                       trunc_level = trunc_level,
                                       scale_factor = 3.0)
    
    naive_nde_adl <- results_naive_nde$NDE
    naive_nde_se_adl <- results_naive_nde$SE
    naive_nde_ci_lower_adl <- results_naive_nde$CI["lower"]
    naive_nde_ci_upper_adl <- results_naive_nde$CI["upper"]
    naive_nde_converged <- TRUE
    
  }, error = function(e) {
    warning(paste("Naive NDE estimation failed:", e$message))
  })
  
  # Calculate bias and coverage for IIE
  iie_bias <- ifelse(cma_converged && !is.na(iie_adl) && is.finite(iie_adl),
                     iie_adl - sim_data$true_nie, NA)
  iie_covers <- ifelse(cma_converged && !is.na(iie_adl) && is.finite(iie_adl),
                       (iie_ci_lower_adl <= sim_data$true_nie) & (sim_data$true_nie <= iie_ci_upper_adl), NA)
  
  # Calculate bias and coverage for IDE
  ide_bias <- ifelse(cma_converged && !is.na(ide_adl) && is.finite(ide_adl),
                     ide_adl - sim_data$true_nde, NA)
  ide_covers <- ifelse(cma_converged && !is.na(ide_adl) && is.finite(ide_adl),
                       (ide_ci_lower_adl <= sim_data$true_nde) & (sim_data$true_nde <= ide_ci_upper_adl), NA)
  
  # Calculate bias and coverage for IDE_alt
  ide_alt_bias <- ifelse(cma_converged && !is.na(ide_alt_adl) && is.finite(ide_alt_adl),
                         ide_alt_adl - sim_data$true_nde, NA)
  ide_alt_covers <- ifelse(cma_converged && !is.na(ide_alt_adl) && is.finite(ide_alt_adl),
                           (ide_alt_ci_lower_adl <= sim_data$true_nde) & (sim_data$true_nde <= ide_alt_ci_upper_adl), NA)
  
  # Calculate bias and coverage for first stage
  first_stage_bias <- ifelse(cma_converged && !is.na(first_stage_est) && is.finite(first_stage_est),
                             first_stage_est - sim_data$true_first_stage, NA)
  first_stage_covers <- ifelse(cma_converged && !is.na(first_stage_est) && is.finite(first_stage_est),
                               (first_stage_ci_lower <= sim_data$true_first_stage) & 
                                 (sim_data$true_first_stage <= first_stage_ci_upper), NA)
  
  # Calculate bias and coverage for naive NIE
  naive_nie_bias <- ifelse(naive_nie_converged && !is.na(naive_nie_adl) && is.finite(naive_nie_adl),
                           naive_nie_adl - sim_data$true_nie, NA)
  naive_nie_covers <- ifelse(naive_nie_converged && !is.na(naive_nie_adl) && is.finite(naive_nie_adl),
                             (naive_nie_ci_lower_adl <= sim_data$true_nie) & 
                               (sim_data$true_nie <= naive_nie_ci_upper_adl), NA)
  
  # Calculate bias and coverage for naive NDE
  naive_nde_bias <- ifelse(naive_nde_converged && !is.na(naive_nde_adl) && is.finite(naive_nde_adl),
                           naive_nde_adl - sim_data$true_nde, NA)
  naive_nde_covers <- ifelse(naive_nde_converged && !is.na(naive_nde_adl) && is.finite(naive_nde_adl),
                             (naive_nde_ci_lower_adl <= sim_data$true_nde) & 
                               (sim_data$true_nde <= naive_nde_ci_upper_adl), NA)
  
  data.frame(n = n, iv_strength = iv_strength, trunc_level = trunc_level, crossFit = crossFit,
             # IIE results
             iie_adl = iie_adl, iie_se_adl = iie_se_adl,
             iie_ci_lower_adl = iie_ci_lower_adl, iie_ci_upper_adl = iie_ci_upper_adl,
             iie_bias = iie_bias, iie_covers = iie_covers,
             iie_reduced = iie_reduced, iie_reduced_se = iie_reduced_se,
             # IDE results (original: TE_reduced - IIE)
             ide_adl = ide_adl, ide_se_adl = ide_se_adl,
             ide_ci_lower_adl = ide_ci_lower_adl, ide_ci_upper_adl = ide_ci_upper_adl,
             ide_bias = ide_bias, ide_covers = ide_covers,
             te_reduced = te_reduced, te_reduced_se = te_reduced_se,
             # IDE_alt results (alternative: IDE_reduced / FS)
             ide_alt_adl = ide_alt_adl, ide_alt_se_adl = ide_alt_se_adl,
             ide_alt_ci_lower_adl = ide_alt_ci_lower_adl, ide_alt_ci_upper_adl = ide_alt_ci_upper_adl,
             ide_alt_bias = ide_alt_bias, ide_alt_covers = ide_alt_covers,
             ide_reduced = ide_reduced, ide_reduced_se = ide_reduced_se,
             # Theta estimates
             theta_11 = theta_11, theta_11_se = theta_11_se,
             theta_10 = theta_10, theta_10_se = theta_10_se,
             theta_00 = theta_00, theta_00_se = theta_00_se,
             # First stage
             first_stage_est = first_stage_est, first_stage_se = first_stage_se,
             first_stage_ci_lower = first_stage_ci_lower, first_stage_ci_upper = first_stage_ci_upper,
             first_stage_bias = first_stage_bias, first_stage_covers = first_stage_covers,
             # Convergence
             cma_converged = cma_converged,
             # Naive estimators
             naive_nie_adl = naive_nie_adl, naive_nie_se_adl = naive_nie_se_adl,
             naive_nie_ci_lower_adl = naive_nie_ci_lower_adl, naive_nie_ci_upper_adl = naive_nie_ci_upper_adl,
             naive_nie_bias = naive_nie_bias, naive_nie_covers = naive_nie_covers,
             naive_nie_converged = naive_nie_converged,
             naive_nde_adl = naive_nde_adl, naive_nde_se_adl = naive_nde_se_adl,
             naive_nde_ci_lower_adl = naive_nde_ci_lower_adl, naive_nde_ci_upper_adl = naive_nde_ci_upper_adl,
             naive_nde_bias = naive_nde_bias, naive_nde_covers = naive_nde_covers,
             naive_nde_converged = naive_nde_converged,
             # True values
             true_iie_adl = sim_data$true_nie, true_ide_adl = sim_data$true_nde,
             true_te_adl = sim_data$true_te, true_iv_ame_on_mediator = sim_data$true_iv_ame_on_mediator,
             true_first_stage = sim_data$true_first_stage,
             prop_insured = mean(data[[exposure_name]]),
             prop_rehabIRF = mean(data[[mediator_name]]),
             prop_IV = mean(data[[iv_name]]))
}


run_simulation_study_iv_mediation <- function(n_sims = 100, n = 1000,
                                              exposure_name = "insured",
                                              mediator_name = "rehabIRF",
                                              iv_name = "Z",
                                              crossFit = FALSE, n_folds = 5,
                                              iv_strength = 0.3, trunc_level = 0.001) {
  
  initial_data <- generate_stroke_data_iv_mediation(n = n,
                                                    exposure_name = exposure_name,
                                                    mediator_name = mediator_name,
                                                    iv_name = iv_name,
                                                    iv_strength = iv_strength)
  fixed_covariates <<- initial_data$fixed_covariates
  
  results_list <- vector("list", n_sims)
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  
  for (i in 1:n_sims) {
    results_list[[i]] <- run_single_simulation_iv_mediation(n = n,
                                                            exposure_name = exposure_name,
                                                            mediator_name = mediator_name,
                                                            iv_name = iv_name,
                                                            crossFit = crossFit,
                                                            n_folds = n_folds,
                                                            iv_strength = iv_strength,
                                                            trunc_level = trunc_level)
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  results_df <- do.call(rbind, results_list)
  rm(fixed_covariates, envir = .GlobalEnv)
  
  return(results_df)
}


summarize_simulation_iv_mediation <- function(results_df) {
  
  # IIE summary
  iie_summary <- list(mean_estimate = mean(results_df$iie_adl, na.rm = TRUE),
                      mean_bias = mean(results_df$iie_bias, na.rm = TRUE),
                      mean_se = mean(results_df$iie_se_adl, na.rm = TRUE),
                      empirical_se = sd(results_df$iie_adl, na.rm = TRUE),
                      coverage = mean(results_df$iie_covers, na.rm = TRUE),
                      rmse = sqrt(mean(results_df$iie_bias^2, na.rm = TRUE)),
                      convergence_rate = sum(results_df$cma_converged, na.rm = TRUE) / nrow(results_df))
  
  # IDE summary (original: TE_reduced - IIE)
  ide_summary <- list(mean_estimate = mean(results_df$ide_adl, na.rm = TRUE),
                      mean_bias = mean(results_df$ide_bias, na.rm = TRUE),
                      mean_se = mean(results_df$ide_se_adl, na.rm = TRUE),
                      empirical_se = sd(results_df$ide_adl, na.rm = TRUE),
                      coverage = mean(results_df$ide_covers, na.rm = TRUE),
                      rmse = sqrt(mean(results_df$ide_bias^2, na.rm = TRUE)),
                      convergence_rate = sum(results_df$cma_converged, na.rm = TRUE) / nrow(results_df))
  
  # IDE_alt summary (alternative: IDE_reduced / FS)
  ide_alt_summary <- list(mean_estimate = mean(results_df$ide_alt_adl, na.rm = TRUE),
                          mean_bias = mean(results_df$ide_alt_bias, na.rm = TRUE),
                          mean_se = mean(results_df$ide_alt_se_adl, na.rm = TRUE),
                          empirical_se = sd(results_df$ide_alt_adl, na.rm = TRUE),
                          coverage = mean(results_df$ide_alt_covers, na.rm = TRUE),
                          rmse = sqrt(mean(results_df$ide_alt_bias^2, na.rm = TRUE)),
                          convergence_rate = sum(results_df$cma_converged, na.rm = TRUE) / nrow(results_df))
  
  # Theta summaries (single set - used for both IIE and IDE)
  theta_11_summary <- list(mean_estimate = mean(results_df$theta_11, na.rm = TRUE),
                           mean_se = mean(results_df$theta_11_se, na.rm = TRUE),
                           empirical_se = sd(results_df$theta_11, na.rm = TRUE))
  
  theta_10_summary <- list(mean_estimate = mean(results_df$theta_10, na.rm = TRUE),
                           mean_se = mean(results_df$theta_10_se, na.rm = TRUE),
                           empirical_se = sd(results_df$theta_10, na.rm = TRUE))
  
  theta_00_summary <- list(mean_estimate = mean(results_df$theta_00, na.rm = TRUE),
                           mean_se = mean(results_df$theta_00_se, na.rm = TRUE),
                           empirical_se = sd(results_df$theta_00, na.rm = TRUE))
  
  # Reduced form summaries
  iie_reduced_summary <- list(mean_estimate = mean(results_df$iie_reduced, na.rm = TRUE),
                              mean_se = mean(results_df$iie_reduced_se, na.rm = TRUE),
                              empirical_se = sd(results_df$iie_reduced, na.rm = TRUE))
  
  te_reduced_summary <- list(mean_estimate = mean(results_df$te_reduced, na.rm = TRUE),
                             mean_se = mean(results_df$te_reduced_se, na.rm = TRUE),
                             empirical_se = sd(results_df$te_reduced, na.rm = TRUE))
  
  ide_reduced_summary <- list(mean_estimate = mean(results_df$ide_reduced, na.rm = TRUE),
                              mean_se = mean(results_df$ide_reduced_se, na.rm = TRUE),
                              empirical_se = sd(results_df$ide_reduced, na.rm = TRUE))
  
  # First stage summary
  first_stage_summary <- list(mean_estimate = mean(results_df$first_stage_est, na.rm = TRUE),
                              mean_bias = mean(results_df$first_stage_bias, na.rm = TRUE),
                              mean_se = mean(results_df$first_stage_se, na.rm = TRUE),
                              empirical_se = sd(results_df$first_stage_est, na.rm = TRUE),
                              coverage = mean(results_df$first_stage_covers, na.rm = TRUE),
                              rmse = sqrt(mean(results_df$first_stage_bias^2, na.rm = TRUE)))
  
  # Naive estimator summaries
  naive_nie_summary <- list(mean_estimate = mean(results_df$naive_nie_adl, na.rm = TRUE),
                            mean_bias = mean(results_df$naive_nie_bias, na.rm = TRUE),
                            mean_se = mean(results_df$naive_nie_se_adl, na.rm = TRUE),
                            empirical_se = sd(results_df$naive_nie_adl, na.rm = TRUE),
                            coverage = mean(results_df$naive_nie_covers, na.rm = TRUE),
                            rmse = sqrt(mean(results_df$naive_nie_bias^2, na.rm = TRUE)),
                            convergence_rate = sum(results_df$naive_nie_converged, na.rm = TRUE) / nrow(results_df))
  
  naive_nde_summary <- list(mean_estimate = mean(results_df$naive_nde_adl, na.rm = TRUE),
                            mean_bias = mean(results_df$naive_nde_bias, na.rm = TRUE),
                            mean_se = mean(results_df$naive_nde_se_adl, na.rm = TRUE),
                            empirical_se = sd(results_df$naive_nde_adl, na.rm = TRUE),
                            coverage = mean(results_df$naive_nde_covers, na.rm = TRUE),
                            rmse = sqrt(mean(results_df$naive_nde_bias^2, na.rm = TRUE)),
                            convergence_rate = sum(results_df$naive_nde_converged, na.rm = TRUE) / nrow(results_df))
  
  # IV diagnostics
  iv_diagnostics <- list(iv_strength_setting = results_df$iv_strength[1],
                         mean_prop_IV1 = mean(results_df$prop_IV, na.rm = TRUE),
                         mean_prop_exposure = mean(results_df$prop_insured, na.rm = TRUE),
                         mean_prop_mediator = mean(results_df$prop_rehabIRF, na.rm = TRUE),
                         mean_true_iv_ame = mean(results_df$true_iv_ame_on_mediator, na.rm = TRUE))
  
  # True values
  true_values <- list(true_iie = mean(results_df$true_iie_adl, na.rm = TRUE),
                      true_ide = mean(results_df$true_ide_adl, na.rm = TRUE),
                      true_te = mean(results_df$true_te_adl, na.rm = TRUE),
                      true_first_stage = mean(results_df$true_first_stage, na.rm = TRUE),
                      prop_mediated = mean(results_df$true_iie_adl, na.rm = TRUE) / 
                        mean(results_df$true_te_adl, na.rm = TRUE))
  
  list(IIE = iie_summary, 
       IDE = ide_summary,
       IDE_alt = ide_alt_summary,
       theta_11 = theta_11_summary, 
       theta_10 = theta_10_summary,
       theta_00 = theta_00_summary,
       IIE_reduced = iie_reduced_summary,
       TE_reduced = te_reduced_summary,
       IDE_reduced = ide_reduced_summary,
       First_Stage = first_stage_summary,
       Naive_NIE = naive_nie_summary, 
       Naive_NDE = naive_nde_summary,
       IV_diagnostics = iv_diagnostics, 
       true_values = true_values,
       n_simulations = nrow(results_df), 
       n = results_df$n[1])
}


print_simulation_summary_iv_mediation <- function(summary_stats) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("                      IV-MEDIATION SIMULATION RESULTS\n")
  cat("================================================================================\n")
  cat(sprintf("N simulations: %d | Sample size: %d | IV strength: %.2f\n",
              summary_stats$n_simulations, summary_stats$n,
              summary_stats$IV_diagnostics$iv_strength_setting))
  cat(sprintf("P(IV=1): %.3f | P(Exposure=1): %.3f | P(Mediator=1): %.3f\n",
              summary_stats$IV_diagnostics$mean_prop_IV1,
              summary_stats$IV_diagnostics$mean_prop_exposure,
              summary_stats$IV_diagnostics$mean_prop_mediator))
  cat(sprintf("Mean IV effect on Mediator: %.4f\n",
              summary_stats$IV_diagnostics$mean_true_iv_ame))
  
  cat("\n--------------------------------------------------------------------------------\n")
  cat("TRUE CAUSAL EFFECTS\n")
  cat("--------------------------------------------------------------------------------\n")
  cat(sprintf("True IIE (NIE):      %.4f\n", summary_stats$true_values$true_iie))
  cat(sprintf("True IDE (NDE):      %.4f\n", summary_stats$true_values$true_ide))
  cat(sprintf("True TE:             %.4f\n", summary_stats$true_values$true_te))
  cat(sprintf("True First Stage:    %.4f\n", summary_stats$true_values$true_first_stage))
  cat(sprintf("Proportion Mediated: %.3f\n", summary_stats$true_values$prop_mediated))
  
  cat("\n--------------------------------------------------------------------------------\n")
  cat("THETA PARAMETER ESTIMATES\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Parameter    Est       Model_SE  Emp_SE    SE_Ratio\n")
  cat("------------ --------  --------  --------  --------\n")
  cat(sprintf("theta(1,1)   %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$theta_11$mean_estimate,
              summary_stats$theta_11$mean_se,
              summary_stats$theta_11$empirical_se,
              summary_stats$theta_11$mean_se / summary_stats$theta_11$empirical_se))
  cat(sprintf("theta(1,0)   %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$theta_10$mean_estimate,
              summary_stats$theta_10$mean_se,
              summary_stats$theta_10$empirical_se,
              summary_stats$theta_10$mean_se / summary_stats$theta_10$empirical_se))
  cat(sprintf("theta(0,0)   %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$theta_00$mean_estimate,
              summary_stats$theta_00$mean_se,
              summary_stats$theta_00$empirical_se,
              summary_stats$theta_00$mean_se / summary_stats$theta_00$empirical_se))
  
  cat("\n--------------------------------------------------------------------------------\n")
  cat("REDUCED FORM AND FIRST STAGE ESTIMATES\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Component              True      Est       Bias      Model_SE  Emp_SE    Coverage\n")
  cat("---------------------- --------  --------  --------  --------  --------  --------\n")
  cat(sprintf("IIE_reduced (θ11-θ10)  NA        %.4f    NA        %.4f    %.4f    NA\n",
              summary_stats$IIE_reduced$mean_estimate,
              summary_stats$IIE_reduced$mean_se,
              summary_stats$IIE_reduced$empirical_se))
  cat(sprintf("TE_reduced (θ11-θ00)   NA        %.4f    NA        %.4f    %.4f    NA\n",
              summary_stats$TE_reduced$mean_estimate,
              summary_stats$TE_reduced$mean_se,
              summary_stats$TE_reduced$empirical_se))
  cat(sprintf("IDE_reduced (θ10-θ00)  NA        %.4f    NA        %.4f    %.4f    NA\n",
              summary_stats$IDE_reduced$mean_estimate,
              summary_stats$IDE_reduced$mean_se,
              summary_stats$IDE_reduced$empirical_se))
  cat(sprintf("First Stage            %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_first_stage,
              summary_stats$First_Stage$mean_estimate,
              summary_stats$First_Stage$mean_bias,
              summary_stats$First_Stage$mean_se,
              summary_stats$First_Stage$empirical_se,
              summary_stats$First_Stage$coverage))
  
  cat("\n--------------------------------------------------------------------------------\n")
  cat("IV-BASED ESTIMATORS (accounting for unmeasured confounding)\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("IIE  = (θ11-θ10) / FS\n")
  cat("IDE  = (θ11-θ00) - IIE  [TE_reduced - IIE]\n")
  cat("IDE* = (θ10-θ00) / FS   [alternative estimator]\n")
  cat("\n")
  cat("Estimand    True      Est       Bias      Model_SE  Emp_SE    Coverage\n")
  cat("----------- --------  --------  --------  --------  --------  --------\n")
  cat(sprintf("IIE         %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_iie,
              summary_stats$IIE$mean_estimate,
              summary_stats$IIE$mean_bias,
              summary_stats$IIE$mean_se,
              summary_stats$IIE$empirical_se,
              summary_stats$IIE$coverage))
  cat(sprintf("IDE         %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_ide,
              summary_stats$IDE$mean_estimate,
              summary_stats$IDE$mean_bias,
              summary_stats$IDE$mean_se,
              summary_stats$IDE$empirical_se,
              summary_stats$IDE$coverage))
  cat(sprintf("IDE*        %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_ide,
              summary_stats$IDE_alt$mean_estimate,
              summary_stats$IDE_alt$mean_bias,
              summary_stats$IDE_alt$mean_se,
              summary_stats$IDE_alt$empirical_se,
              summary_stats$IDE_alt$coverage))
  
  cat("\n--------------------------------------------------------------------------------\n")
  cat("NAIVE ESTIMATORS (ignoring unmeasured confounding)\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("Estimand    True      Est       Bias      Model_SE  Emp_SE    Coverage\n")
  cat("----------- --------  --------  --------  --------  --------  --------\n")
  cat(sprintf("Naive NIE   %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_iie,
              summary_stats$Naive_NIE$mean_estimate,
              summary_stats$Naive_NIE$mean_bias,
              summary_stats$Naive_NIE$mean_se,
              summary_stats$Naive_NIE$empirical_se,
              summary_stats$Naive_NIE$coverage))
  cat(sprintf("Naive NDE   %.4f    %.4f    %.4f    %.4f    %.4f    %.3f\n",
              summary_stats$true_values$true_ide,
              summary_stats$Naive_NDE$mean_estimate,
              summary_stats$Naive_NDE$mean_bias,
              summary_stats$Naive_NDE$mean_se,
              summary_stats$Naive_NDE$empirical_se,
              summary_stats$Naive_NDE$coverage))
  cat("================================================================================\n")
}


# Run simulation
N_SIMS <- 100

sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 6000,
                                               crossFit = FALSE,
                                               iv_strength = 0.30)

sim_summary <- summarize_simulation_iv_mediation(sim_study)
print_simulation_summary_iv_mediation(sim_summary)





# number of sims
N_SIMS <- 10

# N=5000
sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 6000,
                                               crossFit = TRUE,
                                               iv_strength = 0.30)

# display results
sim_summary <- summarize_simulation_iv_mediation(sim_study)
print_simulation_summary_iv_mediation(sim_summary)





