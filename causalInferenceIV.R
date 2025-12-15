








### Load libraries
library(betareg)



# Calculates a robust TMLE estimator for the ATE for a single outcome
robustcompATE <- function(outcome_mod,ps_mod,covariates,
                          exposure_name,outcome_name,
                          trunc_level,scale_factor = 1.0) {
  
  ### calculate clever covariates for TMLE
  ps_pred <- predict(ps_mod,newdata = covariates,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  exposure <- covariates[[exposure_name]]
  trunc_weights_0 <- ifelse(exposure == 0,-1.0 / (1.0 - ps_trunc),0)
  trunc_weights_1 <- ifelse(exposure == 1,1.0 / ps_trunc,0)
  H_0_Y <- -1.0 * trunc_weights_0
  H_1_Y <- trunc_weights_1
  H_A_Y <- ifelse(exposure == 1,H_1_Y,H_0_Y)
  
  ### estimate fluctuation parameter
  pred_y <- predict(outcome_mod,newdata = covariates,type = "response")
  outcome_values <- covariates[[outcome_name]]
  data_tmle <- data.frame(Y = outcome_values,
                          pred_vals = pred_y,
                          H_A_Y = H_A_Y)
  eps <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                  data = data_tmle,family = "binomial"))
  
  ### update predictions using fluctuation parameter
  pred_y_A <- predict(outcome_mod,covariates,type = "response")
  covariates_0 <- covariates
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates
  covariates_1[[exposure_name]] <- 1
  pred_y_0 <- predict(outcome_mod,covariates_0,type = "response")
  pred_y_1 <- predict(outcome_mod,covariates_1,type = "response")
  pred_star_y_A <- plogis(qlogis(pred_y_A) + eps * H_A_Y)
  pred_star_y_0 <- plogis(qlogis(pred_y_0) + eps * H_0_Y)
  pred_star_y_1 <- plogis(qlogis(pred_y_1) + eps * H_1_Y)
  
  ### calculate ATE using targeted estimates
  ate_star <- pred_star_y_1 - pred_star_y_0
  ate_star_scaled <- ate_star * scale_factor
  
  ### calculate efficient influence function
  D_Y <- H_A_Y * (outcome_values - pred_star_y_A)
  D_W <- ate_star - mean(ate_star)
  EIF <- D_Y + D_W
  
  ### calculate standard errors and confidence intervals
  n <- length(EIF)
  var_ate <- var(EIF) / n
  se_ate_scaled <- sqrt(var_ate) * scale_factor
  ate_mean <- mean(ate_star_scaled)
  ci_lower <- ate_mean - qnorm(0.975) * se_ate_scaled
  ci_upper <- ate_mean + qnorm(0.975) * se_ate_scaled
  
  ### return results
  return(list(ATE = ate_mean,
              ITE = ate_star_scaled,
              IF = EIF,
              SE = se_ate_scaled,
              CI = c(lower = ci_lower,upper = ci_upper)))
}



# Calculate ATE using robust TMLE estimates but ignoring unmeasured confounding
calculate_ate_naive_robust <- function(data,exposure_name = "rehabIRF",
                                       outcome_name_01,
                                       scale_factor = 1.0,trunc_level = 0.05) {
  
  ps_formula <- as.formula(paste(exposure_name,"~ age + stroke_severity + comorbidity"))
  ps_model <- glm(ps_formula,data = data,family = binomial())
  outcome_formula <- as.formula(paste(outcome_name_01,"~",exposure_name,
                                      "+ age + stroke_severity + comorbidity"))
  outcome_mod <- betareg(outcome_formula,data = data)
  results_ate <- robustcompATE(outcome_mod = outcome_mod,
                               ps_mod = ps_model,
                               covariates = data,
                               exposure_name = exposure_name,
                               outcome_name = outcome_name_01,
                               trunc_level = trunc_level,
                               scale_factor = scale_factor)
  
  ### return results
  return(list(ATE = results_ate$ATE,
              ITE = results_ate$ITE,
              SE = results_ate$SE,
              CI = results_ate$CI,
              IF = results_ate$IF))
}



# Calculate LATE using robust TMLE estimates
calculate_late_robust <- function(data,iv_name = "Z",
                                  exposure_name = "rehabIRF",
                                  outcome_name_01,
                                  scale_factor = 1.0,trunc_level = 0.05) {
  
  ### calculate ATE for numerator
  iv_formula <- as.formula(paste(iv_name,"~ age + stroke_severity + comorbidity"))
  iv_model <- glm(iv_formula,data = data,family = binomial())
  outcome_formula_y <- as.formula(paste(outcome_name_01,"~",iv_name,
                                        "+ age + stroke_severity + comorbidity"))
  outcome_mod_y <- betareg(outcome_formula_y,data = data)
  results_numerator <- robustcompATE(outcome_mod = outcome_mod_y,
                                     ps_mod = iv_model,
                                     covariates = data,
                                     exposure_name = iv_name,
                                     outcome_name = outcome_name_01,
                                     trunc_level = trunc_level,
                                     scale_factor = scale_factor)
  
  ### calculate ATE for denominator
  outcome_formula_a <- as.formula(paste(exposure_name,"~",iv_name,
                                        "+ age + stroke_severity + comorbidity"))
  outcome_mod_a <- glm(outcome_formula_a,data = data,family = binomial())
  results_denominator <- robustcompATE(outcome_mod = outcome_mod_a,
                                       ps_mod = iv_model,
                                       covariates = data,
                                       exposure_name = iv_name,
                                       outcome_name = exposure_name,
                                       trunc_level = trunc_level,
                                       scale_factor = 1.0)
  
  ### calculate LATE and approximate 95% confidence interval using delta method
  numerator <- results_numerator$ATE
  denominator <- results_denominator$ATE
  late <- numerator / denominator
  var_late <- ((1/denominator)^2 * results_numerator$SE^2
               + (numerator/denominator^2)^2 * results_denominator$SE^2)
  se_late <- sqrt(var_late)
  ci_lower <- late - qnorm(0.975) * se_late
  ci_upper <- late + qnorm(0.975) * se_late
  
  ### return results
  return(list(LATE = late,
              SE = se_late,
              CI = c(lower = ci_lower,upper = ci_upper),
              numerator = numerator,
              denominator = denominator,
              numerator_se = results_numerator$SE,
              denominator_se = results_denominator$SE))
}



# Calculate robust ATE using doubly robust g-estimation with IV
calculate_ate_iv_robust <- function(outcome_mod_init,iv_model,data,
                                    exp_model_formula,delta_formula,
                                    iv_name = "Z",
                                    exposure_name = "rehabIRF",
                                    outcome_name_01,
                                    scale_factor = 1.0,trunc_level = 0.05) {
  
  # Extract IV and exposure variables
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
    # Use fixed covariates
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
  ps <- 0.28 + 
    iv_strength * IV +
    0.03 * age_scaled +
    0.08 * post_discharge_disability_scaled +
    0.03 * stroke_severity_scaled +
    -0.02 * comorbidity_scaled
  
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
    0.1 * comorbidity_scaled +
    -0.15 * exposure * post_discharge_disability_scaled
  
  mu_adl <- plogis(mu_adl_logit)
  
  phi_adl <- 15
  shape1_adl <- mu_adl * phi_adl
  shape2_adl <- (1 - mu_adl) * phi_adl
  
  OUT3_ADL_IADL_01 <- rbeta(n,shape1 = shape1_adl,shape2 = shape2_adl)
  OUT3_ADL_IADL <- OUT3_ADL_IADL_01 * 3.0
  
  ### calculate true ATE of treatment A on outcome Y
  mu_adl_0_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 0 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled + -0.15 * 0 * post_discharge_disability_scaled
  
  mu_adl_1_logit <- -1.0 + 1.0 * post_discharge_disability_scaled + -0.5 * 1 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 
    0.1 * comorbidity_scaled + -0.15 * 1 * post_discharge_disability_scaled
  
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
run_single_simulation_iv <- function(n = 1000,
                                     iv_name = "Z",
                                     exposure_name = "rehabIRF",
                                     trunc_level = 0.05,
                                     iv_strength = 0.5,
                                     fixed_covariates = NULL,
                                     return_vectors = TRUE) {
  
  # generate data with IV (using fixed covariates if provided)
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
  
  ### Fit models needed for ATE calculation
  # IV propensity score model
  iv_formula <- as.formula(paste(iv_name,"~ age + stroke_severity + comorbidity"))
  iv_model <- glm(iv_formula,data = data,family = binomial())
  
  # Initial outcome model
  outcome_model <- betareg(OUT3_ADL_IADL_01 ~ age + stroke_severity + comorbidity,
                           data = data)
  
  # Define formulas for ATE calculation
  exp_model_formula <- as.formula(paste(exposure_name,"~ age + stroke_severity + comorbidity"))
  delta_formula <- ~ age + stroke_severity + comorbidity
  
  ### calculate naive ATE for ADL
  ### ignores unmeasured confounding
  ate_adl_naive <- calculate_ate_naive_robust(data = data,
                                              exposure_name = exposure_name,
                                              outcome_name_01 = "OUT3_ADL_IADL_01",
                                              scale_factor = 3.0,
                                              trunc_level = trunc_level)
  
  ### calculate LATE for ADL
  late_adl <- calculate_late_robust(data = data,
                                    iv_name = iv_name,
                                    exposure_name = exposure_name,
                                    outcome_name_01 = "OUT3_ADL_IADL_01",
                                    scale_factor = 3.0,
                                    trunc_level = trunc_level)
  
  ### calculate robust ATE
  ate_adl <- calculate_ate_iv_robust(outcome_mod_init = outcome_model,
                                     iv_model = iv_model,
                                     data = data,
                                     exp_model_formula = exp_model_formula,
                                     delta_formula = delta_formula,
                                     iv_name = iv_name,
                                     exposure_name = exposure_name,
                                     outcome_name_01 = "OUT3_ADL_IADL_01",
                                     scale_factor = 3.0,
                                     trunc_level = trunc_level)
  
  ### extract results
  if (return_vectors) {
    # return full results including vectors
    results_df <- data.frame(n = n,
                             iv_name = iv_name,
                             exposure_name = exposure_name,
                             trunc_level = trunc_level,
                             iv_strength = iv_strength,
                             # LATE results for ADL
                             late_adl = late_adl$LATE,
                             late_se_adl = late_adl$SE,
                             late_ci_lower_adl = late_adl$CI["lower"],
                             late_ci_upper_adl = late_adl$CI["upper"],
                             late_numerator_adl = late_adl$numerator,
                             late_denominator_adl = late_adl$denominator,
                             late_numerator_se_adl = late_adl$numerator_se,
                             late_denominator_se_adl = late_adl$denominator_se,
                             # ATE results for ADL (IV-based)
                             ate_adl = ate_adl$ATE,
                             ate_se_adl = ate_adl$SE,
                             ate_ci_lower_adl = ate_adl$CI["lower"],
                             ate_ci_upper_adl = ate_adl$CI["upper"],
                             # Naive ATE results for ADL
                             ate_naive_adl = ate_adl_naive$ATE,
                             ate_naive_se_adl = ate_adl_naive$SE,
                             ate_naive_ci_lower_adl = ate_adl_naive$CI["lower"],
                             ate_naive_ci_upper_adl = ate_adl_naive$CI["upper"],
                             # true values
                             true_ate_adl = sim_data$true_ate_adl)
    
    # add vector results as separate list elements
    results_df$ate_vals_adl <- list(ate_adl$ITE)
    results_df$delta_D_X <- list(ate_adl$delta_D_X)
    results_df$ate_naive_vals_adl <- list(ate_adl_naive$ITE)
    
  } else {
    # return only scalar results to combine across simulations
    results_df <- data.frame(n = n,
                             iv_name = iv_name,
                             exposure_name = exposure_name,
                             trunc_level = trunc_level,
                             iv_strength = iv_strength,
                             # LATE results for ADL
                             late_adl = late_adl$LATE,
                             late_se_adl = late_adl$SE,
                             late_ci_lower_adl = late_adl$CI["lower"],
                             late_ci_upper_adl = late_adl$CI["upper"],
                             late_numerator_adl = late_adl$numerator,
                             late_denominator_adl = late_adl$denominator,
                             late_numerator_se_adl = late_adl$numerator_se,
                             late_denominator_se_adl = late_adl$denominator_se,
                             # ATE results for ADL (IV-based)
                             ate_adl = ate_adl$ATE,
                             ate_se_adl = ate_adl$SE,
                             ate_ci_lower_adl = ate_adl$CI["lower"],
                             ate_ci_upper_adl = ate_adl$CI["upper"],
                             # Naive ATE results for ADL
                             ate_naive_adl = ate_adl_naive$ATE,
                             ate_naive_se_adl = ate_adl_naive$SE,
                             ate_naive_ci_lower_adl = ate_adl_naive$CI["lower"],
                             ate_naive_ci_upper_adl = ate_adl_naive$CI["upper"],
                             # true values
                             true_ate_adl = sim_data$true_ate_adl)
    
    # optionally add summary statistics of the vectors
    results_df$ate_vals_adl_mean <- mean(ate_adl$ITE)
    results_df$ate_vals_adl_sd <- sd(ate_adl$ITE)
    results_df$ate_vals_adl_min <- min(ate_adl$ITE)
    results_df$ate_vals_adl_max <- max(ate_adl$ITE)
    results_df$delta_D_X_mean <- mean(ate_adl$delta_D_X)
    results_df$delta_D_X_sd <- sd(ate_adl$delta_D_X)
    # naive ATE summary statistics
    results_df$ate_naive_vals_adl_mean <- mean(ate_adl_naive$ITE)
    results_df$ate_naive_vals_adl_sd <- sd(ate_adl_naive$ITE)
    results_df$ate_naive_vals_adl_min <- min(ate_adl_naive$ITE)
    results_df$ate_naive_vals_adl_max <- max(ate_adl_naive$ITE)
  }
  
  return(results_df)
}



# Run simulation study for IV analysis
run_simulation_study_iv <- function(n_sims = 500,n = 1000,
                                    iv_name = "Z",
                                    exposure_name = "rehabIRF",
                                    trunc_level = 0.05,iv_strength = 0.5,
                                    use_fixed_covariates = TRUE) {
  
  # generate fixed covariates once if requested
  fixed_covariates <- NULL
  if (use_fixed_covariates) {
    # generate covariates once
    initial_data <- generate_stroke_data_iv(n = n,
                                            iv_name = iv_name,
                                            exposure_name = exposure_name,
                                            iv_strength = iv_strength)
    fixed_covariates <- initial_data$fixed_covariates
  }
  
  # run sims sequentially
  results_list <- lapply(1:n_sims,function(i) {
    tryCatch({
      run_single_simulation_iv(n = n,
                               iv_name = iv_name,
                               exposure_name = exposure_name,
                               trunc_level = trunc_level,
                               iv_strength = iv_strength,
                               fixed_covariates = fixed_covariates,
                               return_vectors = FALSE)
    },error = function(e) {
      cat(paste("Error in simulation",i,":",e$message,"\n"))
      return(NULL)
    })
  })
  
  results_list <- results_list[!sapply(results_list,is.null)]
  results <- do.call(rbind,results_list)
  
  # calculate simulation performance metrics for LATE
  results$late_bias_adl <- results$late_adl - results$true_ate_adl
  
  # calculate LATE coverage (now available from function output)
  results$late_coverage_adl <- (results$late_ci_lower_adl <= results$true_ate_adl) & 
    (results$late_ci_upper_adl >= results$true_ate_adl)
  
  # calculate simulation performance metrics for ATE (IV-based)
  results$ate_bias_adl <- results$ate_adl - results$true_ate_adl
  
  # calculate coverage for ATE (IV-based) - now using CI from function
  results$ate_coverage_adl <- (results$ate_ci_lower_adl <= results$true_ate_adl) & 
    (results$ate_ci_upper_adl >= results$true_ate_adl)
  
  # calculate simulation performance metrics for naive ATE
  results$ate_naive_bias_adl <- results$ate_naive_adl - results$true_ate_adl
  
  # calculate coverage for naive ATE - now using CI from function
  results$ate_naive_coverage_adl <- (results$ate_naive_ci_lower_adl <= results$true_ate_adl) & 
    (results$ate_naive_ci_upper_adl >= results$true_ate_adl)
  
  return(results)
}


# Summarize IV simulation results for LATE, ATE, and naive ATE
summarize_simulation_iv <- function(results) {
  
  # naive ATE summary statistics
  ate_naive_summary <- list(mean_estimate = mean(results$ate_naive_adl),
                            mean_bias = mean(results$ate_naive_bias_adl),
                            mean_se = mean(results$ate_naive_se_adl),
                            empirical_se = sd(results$ate_naive_adl),
                            coverage = mean(results$ate_naive_coverage_adl))
  
  # LATE summary statistics
  late_summary <- list(mean_estimate = mean(results$late_adl),
                       mean_bias = mean(results$late_bias_adl),
                       mean_se = mean(results$late_se_adl),
                       empirical_se = sd(results$late_adl),
                       coverage = mean(results$late_coverage_adl),
                       mean_numerator = mean(results$late_numerator_adl),
                       mean_denominator = mean(results$late_denominator_adl),
                       mean_numerator_se = mean(results$late_numerator_se_adl),
                       mean_denominator_se = mean(results$late_denominator_se_adl))
  
  # ATE summary statistics (IV-based)
  ate_summary <- list(mean_estimate = mean(results$ate_adl),
                      mean_bias = mean(results$ate_bias_adl),
                      mean_se = mean(results$ate_se_adl),
                      empirical_se = sd(results$ate_adl),
                      coverage = mean(results$ate_coverage_adl),
                      mean_delta_D = mean(results$mean_delta_D_adl),
                      mean_delta_X = mean(results$mean_delta_X_adl))
  
  # IV diagnostics
  iv_diagnostics <- list(mean_first_stage_strength = mean(results$first_stage_strength),
                         sd_first_stage_strength = sd(results$first_stage_strength),
                         mean_prop_IV1 = mean(results$prop_IV1),
                         mean_prop_treated = mean(results$prop_treated),
                         iv_strength_setting = results$iv_strength[1],
                         iv_name = results$iv_name[1],
                         exposure_name = results$exposure_name[1])
  
  # true values should be constant across simulations
  true_values <- list(true_ate_adl = mean(results$true_ate_adl))
  
  summary_stats <- list(ATE_naive = ate_naive_summary,
                        LATE = late_summary,
                        ATE = ate_summary,
                        IV_diagnostics = iv_diagnostics,
                        true_values = true_values,
                        n_simulations = nrow(results))
  
  return(summary_stats)
}

# Print formatted summary of IV simulation results
print_simulation_summary_iv <- function(summary_stats) {
  
  cat("\n")
  cat("IV SIMULATION RESULTS\n")
  cat("=====================\n")
  cat(sprintf("N simulations: %d | True ATE: %.4f | IV strength: %.2f | First stage: %.4f\n",
              summary_stats$n_simulations,
              summary_stats$true_values$true_ate_adl,
              summary_stats$IV_diagnostics$iv_strength_setting,
              summary_stats$IV_diagnostics$mean_first_stage_strength))
  cat(sprintf("IV: '%s' | Exposure: '%s'\n",
              summary_stats$IV_diagnostics$iv_name,
              summary_stats$IV_diagnostics$exposure_name))
  cat(sprintf("P(IV=1): %.3f | P(Treated): %.3f\n\n",
              summary_stats$IV_diagnostics$mean_prop_IV1,
              summary_stats$IV_diagnostics$mean_prop_treated))
  
  # Create results table with model-based SE
  cat("Method      Est      Bias     Model_SE Emp_SE   Coverage\n")
  cat("----------  -------  -------  -------- -------  --------\n")
  
  # Naive ATE results
  cat(sprintf("ATE (Naive) %.4f   %.4f   %.4f   %.4f   %.3f\n",
              summary_stats$ATE_naive$mean_estimate,
              summary_stats$ATE_naive$mean_bias,
              summary_stats$ATE_naive$mean_se,
              summary_stats$ATE_naive$empirical_se,
              summary_stats$ATE_naive$coverage))
  
  # LATE results (now includes coverage)
  cat(sprintf("LATE        %.4f   %.4f   %.4f   %.4f   %.3f\n",
              summary_stats$LATE$mean_estimate,
              summary_stats$LATE$mean_bias,
              summary_stats$LATE$mean_se,
              summary_stats$LATE$empirical_se,
              summary_stats$LATE$coverage))
  
  # ATE results (IV-based)
  cat(sprintf("ATE (IV)    %.4f   %.4f   %.4f   %.4f   %.3f\n",
              summary_stats$ATE$mean_estimate,
              summary_stats$ATE$mean_bias,
              summary_stats$ATE$mean_se,
              summary_stats$ATE$empirical_se,
              summary_stats$ATE$coverage))
  
  cat("\n")
}







# number of sims
N_SIMS <- 800
set.seed(80924)

# run simulation study with N = 500
sim_results <- run_simulation_study_iv(n_sims = N_SIMS,n = 500,
                                       trunc_level = 0.005,iv_strength = 0.35,
                                       use_fixed_covariates = TRUE)

# summarize and print
sim_summary <- summarize_simulation_iv(sim_results)
print_simulation_summary_iv(sim_summary)

# run simulation study with N = 1500
sim_results <- run_simulation_study_iv(n_sims = N_SIMS,n = 1500,
                                       trunc_level = 0.005,iv_strength = 0.35,
                                       use_fixed_covariates = TRUE)

# summarize and print
sim_summary <- summarize_simulation_iv(sim_results)
print_simulation_summary_iv(sim_summary)

# run simulation study with N = 5000
sim_results <- run_simulation_study_iv(n_sims = N_SIMS,n = 5000,
                                       trunc_level = 0.005,iv_strength = 0.35,
                                       use_fixed_covariates = TRUE)

# summarize and print
sim_summary <- summarize_simulation_iv(sim_results)
print_simulation_summary_iv(sim_summary)













