








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



# Calculates the numerator for the complier subgroup indirect effect
robustcompNIE_reduced_form <- function(outcome_mod,mediator_mod_no_iv,iv_mod,ps_mod,covariates,
                                       exposure_name,mediator_name,iv_name,outcome_name,
                                       trunc_level,scale_factor = 1.0) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  mediator <- covariates[[mediator_name]]
  iv <- covariates[[iv_name]]
  outcome <- covariates[[outcome_name]]
  
  ### calculate propensity score and mediator densities
  ps_pred <- predict(ps_mod,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  covariates_a0 <- covariates
  covariates_a0[[exposure_name]] <- 0
  med_prob_a0 <- predict(mediator_mod_no_iv,newdata = covariates_a0,type = "response")
  covariates_a1 <- covariates
  covariates_a1[[exposure_name]] <- 1
  med_prob_a1 <- predict(mediator_mod_no_iv,newdata = covariates_a1,type = "response")
  
  med_prob_a0 <- pmax(pmin(med_prob_a0,1.0 - trunc_level),trunc_level)
  med_prob_a1 <- pmax(pmin(med_prob_a1,1.0 - trunc_level),trunc_level)
  
  ### evaluate mediator density at observed instrument values
  Q_M_at_Z_a0 <- ifelse(iv == 1,med_prob_a0,1 - med_prob_a0)
  Q_M_at_Z_a1 <- ifelse(iv == 1,med_prob_a1,1 - med_prob_a1)
  
  ### target outcome regressions
  H_Y <- ifelse(exposure == 1,
                (1 - Q_M_at_Z_a0 / Q_M_at_Z_a1) / ps_trunc,
                0)
  pred_y <- predict(outcome_mod,newdata = covariates,type = "response")
  data_tmle_y <- data.frame(Y = outcome,
                            pred_vals = pred_y,
                            H_Y = H_Y)
  eps_y <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_Y,
                    data = data_tmle_y,family = "binomial"))
  pred_y_star <- plogis(qlogis(pred_y) + eps_y * H_Y)
  
  ### calculate outcome predictions for counterfactual exposure values
  covariates_a1_z0_pred <- covariates
  covariates_a1_z0_pred[[exposure_name]] <- 1
  covariates_a1_z0_pred[[iv_name]] <- 0
  pred_y_a1_z0 <- predict(outcome_mod,newdata = covariates_a1_z0_pred,type = "response")
  covariates_a1_z1_pred <- covariates
  covariates_a1_z1_pred[[exposure_name]] <- 1
  covariates_a1_z1_pred[[iv_name]] <- 1
  pred_y_a1_z1 <- predict(outcome_mod,newdata = covariates_a1_z1_pred,type = "response")
  
  # target updates for Z=0 and Z=1
  H_Y_a1_z0 <- (1 - (1 - med_prob_a0) / (1 - med_prob_a1)) / ps_trunc
  H_Y_a1_z1 <- (1 - med_prob_a0 / med_prob_a1) / ps_trunc
  pred_y_a1_z0_star <- plogis(qlogis(pred_y_a1_z0) + eps_y * H_Y_a1_z0)
  pred_y_a1_z1_star <- plogis(qlogis(pred_y_a1_z1) + eps_y * H_Y_a1_z1)
  
  ### weight by mediator densities
  psi_nie_z_0 <- pred_y_a1_z0_star * (1 - med_prob_a0) + pred_y_a1_z1_star * med_prob_a0
  psi_nie_z_1 <- pred_y_a1_z0_star * (1 - med_prob_a1) + pred_y_a1_z1_star * med_prob_a1
  
  ### target psi_NIE,Z
  H_Z_1 <- ifelse(exposure == 1,1 / ps_trunc,0)
  H_Z_0 <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
  
  # get Q_bar for observed data
  covariates_set_a1 <- covariates
  covariates_set_a1[[exposure_name]] <- 1
  pred_y_a1_obs <- predict(outcome_mod,newdata = covariates_set_a1,type = "response")
  H_Y_a1_obs <- ifelse(iv == 1,H_Y_a1_z1,H_Y_a1_z0)
  Q_bar_Y_1Z_star <- plogis(qlogis(pred_y_a1_obs) + eps_y * H_Y_a1_obs)
  
  # scale predictions for targeting
  all_preds <- c(psi_nie_z_0,psi_nie_z_1,Q_bar_Y_1Z_star)
  min_all <- min(all_preds)
  max_all <- max(all_preds)
  range_all <- max_all - min_all
  if (range_all < 1e-10) range_all <- 1
  
  scale_pred <- function(x) {
    x_scaled <- (x - min_all) / range_all
    x_scaled <- (x_scaled * (n - 1) + 0.5) / n
    x_scaled <- pmax(pmin(x_scaled,1 - 1e-10),1e-10)
    return(x_scaled)
  }
  
  psi_nie_z_0_scaled <- scale_pred(psi_nie_z_0)
  psi_nie_z_1_scaled <- scale_pred(psi_nie_z_1)
  Q_bar_Y_1Z_star_scaled <- scale_pred(Q_bar_Y_1Z_star)
  
  # target for A=1
  data_tmle_a1 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star_scaled,
                             pred_vals = psi_nie_z_1_scaled,
                             H_Z_1 = H_Z_1)
  data_tmle_a1_subset <- data_tmle_a1[exposure == 1,]
  eps_nie_1 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_1,
                        data = data_tmle_a1_subset,family = "binomial"))
  
  # target for A=0
  data_tmle_a0 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star_scaled,
                             pred_vals = psi_nie_z_0_scaled,
                             H_Z_0 = H_Z_0)
  data_tmle_a0_subset <- data_tmle_a0[exposure == 0,]
  eps_nie_0 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_0,
                        data = data_tmle_a0_subset,family = "binomial"))
  
  # calculate targeted estimates
  psi_nie_z_1_star_scaled <- plogis(qlogis(psi_nie_z_1_scaled) + eps_nie_1 * H_Z_1)
  psi_nie_z_0_star_scaled <- plogis(qlogis(psi_nie_z_0_scaled) + eps_nie_0 * H_Z_0)
  
  # transform back
  unscale_pred <- function(x_scaled) {
    (x_scaled * n - 0.5) / (n - 1) * range_all + min_all
  }
  
  psi_nie_z_1_star <- unscale_pred(psi_nie_z_1_star_scaled)
  psi_nie_z_0_star <- unscale_pred(psi_nie_z_0_star_scaled)
  nie_vec <- (psi_nie_z_1_star - psi_nie_z_0_star) * scale_factor
  nie_est <- mean(nie_vec)
  
  ### calculate the efficient influence function
  D_Y <- H_Y * (outcome - pred_y_star)
  psi_nie_z_A_star <- ifelse(exposure == 1,psi_nie_z_1_star,psi_nie_z_0_star)
  residual <- Q_bar_Y_1Z_star - psi_nie_z_A_star
  H_Z <- (2 * exposure - 1) / ifelse(exposure == 1,ps_trunc,1 - ps_trunc)
  D_Z <- H_Z * residual
  D_W <- nie_vec - nie_est
  EIF <- (D_Y + D_Z + D_W) * scale_factor
  
  ### calculate standard errors
  var_nie <- var(EIF) / n
  se_nie <- sqrt(var_nie)
  ci_lower <- nie_est - qnorm(0.975) * se_nie
  ci_upper <- nie_est + qnorm(0.975) * se_nie
  
  return(list(NIE = nie_est,
              NIE_individual = nie_vec,
              IF = EIF,
              SE = se_nie,
              CI = c(lower = ci_lower,upper = ci_upper)))
}


# Robust estimator for the natural indirect effect in the complier subgroup
robustcompNIE_IV <- function(outcome_mod,mediator_mod_no_iv,composite_mod,iv_mod,
                             ps_mod,covariates,
                             exposure_name,mediator_name,iv_name,outcome_name,
                             composite_outcome_name,
                             trunc_level,scale_factor = 1.0) {
  
  n <- nrow(covariates)
  
  # indirect effect using the instrument in the outcome expectations
  reduced_form <- robustcompNIE_reduced_form(outcome_mod = outcome_mod,
                                             mediator_mod_no_iv = mediator_mod_no_iv,
                                             iv_mod = iv_mod,
                                             ps_mod = ps_mod,
                                             covariates = covariates,
                                             exposure_name = exposure_name,
                                             mediator_name = mediator_name,
                                             iv_name = iv_name,
                                             outcome_name = outcome_name,
                                             trunc_level = trunc_level,
                                             scale_factor = scale_factor)
  
  # effect of the instrument on the composite of exposure = 1 and mediator = 1
  first_stage <- robustcompATE(outcome_mod = composite_mod,
                               ps_mod = iv_mod,
                               covariates = covariates,
                               exposure_name = iv_name,
                               outcome_name = composite_outcome_name,
                               trunc_level = trunc_level,
                               scale_factor = 1.0)
  
  # ratio estimator for complier subgroup indirect effect
  nie_iv <- reduced_form$NIE / first_stage$ATE
  
  # calculate standard errors using the delta method
  var_ratio <- (nie_iv^2) * ((reduced_form$SE^2) / (reduced_form$NIE^2) +
                               (first_stage$SE^2) / (first_stage$ATE^2))
  se_nie_iv <- sqrt(var_ratio)
  ci_lower <- nie_iv - qnorm(0.975) * se_nie_iv
  ci_upper <- nie_iv + qnorm(0.975) * se_nie_iv
  
  return(list(NIE = nie_iv,
              NIE_reduced = reduced_form$NIE,
              first_stage = first_stage$ATE,
              SE = se_nie_iv,
              SE_reduced = reduced_form$SE,
              SE_first_stage = first_stage$SE,
              CI = c(lower = ci_lower,upper = ci_upper)))
}



### Generate simulated data with IV-mediation structure
generate_stroke_data_iv_mediation <- function(n = 1000,
                                              exposure_name = "insured",
                                              mediator_name = "rehabIRF",
                                              iv_name = "Z",
                                              iv_strength = 0.3,
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
  
  ### generate exposure (insured)
  logit_ps_exposure <- 0.85 +
    0.05 * age_scaled + 
    0.06 * stroke_severity_scaled + 
    -0.04 * comorbidity_scaled
  
  ps_exposure <- plogis(logit_ps_exposure)
  insured <- rbinom(n,size = 1,prob = ps_exposure)
  
  ### generate instrumental variable (Z)
  logit_iv <- 0.0 + 
    0.05 * age_scaled + 
    0.05 * stroke_severity_scaled - 
    0.05 * comorbidity_scaled
  
  iv_prob <- plogis(logit_iv)
  IV <- rbinom(n,size = 1,prob = iv_prob)
  
  ### generate mediator (rehabIRF) using logistic model
  logit_ps_mediator <- -0.8 +
    iv_strength * IV +
    0.6 * insured +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled
  
  ps_mediator <- plogis(logit_ps_mediator)
  rehabIRF <- rbinom(n,size = 1,prob = ps_mediator)
  
  ### calculate true average marginal effect of Z on M
  # P(M=1|W,A,Z=1) - P(M=1|W,A,Z=0) averaged over the population
  logit_ps_mediator_z1 <- -0.8 +
    iv_strength * 1 +
    0.6 * insured +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled
  
  logit_ps_mediator_z0 <- -0.8 +
    iv_strength * 0 +
    0.6 * insured +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled
  
  ps_mediator_z1 <- plogis(logit_ps_mediator_z1)
  ps_mediator_z0 <- plogis(logit_ps_mediator_z0)
  
  true_iv_ame_on_mediator <- mean(ps_mediator_z1 - ps_mediator_z0)
  
  ### generate outcome (OUT3_ADL_IADL)
  mu_adl_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * insured +
    -0.5 * rehabIRF +
    -0.25 * insured * rehabIRF +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  
  mu_adl <- plogis(mu_adl_logit)
  
  phi_adl <- 15
  shape1_adl <- mu_adl * phi_adl
  shape2_adl <- (1.0 - mu_adl) * phi_adl
  
  OUT3_ADL_IADL_01 <- rbeta(n,shape1 = shape1_adl,shape2 = shape2_adl)
  OUT3_ADL_IADL <- OUT3_ADL_IADL_01 * 3.0
  
  ### calculate true causal effects
  
  # Y(0,0)
  mu_00_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 0 + -0.5 * 0 + -0.25 * 0 * 0 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 0.1 * comorbidity_scaled
  mu_00 <- plogis(mu_00_logit) * 3.0
  
  # Y(0,1)
  mu_01_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 0 + -0.5 * 1 + -0.25 * 0 * 1 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 0.1 * comorbidity_scaled
  mu_01 <- plogis(mu_01_logit) * 3.0
  
  # Y(1,0)
  mu_10_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 1 + -0.5 * 0 + -0.25 * 1 * 0 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 0.1 * comorbidity_scaled
  mu_10 <- plogis(mu_10_logit) * 3.0
  
  # Y(1,1)
  mu_11_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 1 + -0.5 * 1 + -0.25 * 1 * 1 +
    0.4 * stroke_severity_scaled + 0.2 * age_scaled + 0.1 * comorbidity_scaled
  mu_11 <- plogis(mu_11_logit) * 3.0
  
  # Calculate mediator propensity scores under each exposure level WITHOUT IV
  logit_ps_m0_natural <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled +
    0.6 * 0
  ps_m0_natural <- plogis(logit_ps_m0_natural)
  
  logit_ps_m1_natural <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled +
    0.6 * 1
  ps_m1_natural <- plogis(logit_ps_m1_natural)
  
  # Total Effect
  true_te <- mean((mu_11 * ps_m1_natural + mu_10 * (1 - ps_m1_natural)) - 
                    (mu_01 * ps_m0_natural + mu_00 * (1 - ps_m0_natural)))
  
  # Natural Direct Effect
  true_nde <- mean((mu_11 * ps_m0_natural + mu_10 * (1 - ps_m0_natural)) - 
                     (mu_01 * ps_m0_natural + mu_00 * (1 - ps_m0_natural)))
  
  # Natural Indirect Effect
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
                     iv_prob = iv_prob)
  
  data[[exposure_name]] <- insured
  data[[mediator_name]] <- rehabIRF
  data[[iv_name]] <- IV
  
  return(list(
    data = data,
    true_te = true_te,
    true_nde = true_nde,
    true_nie = true_nie,
    true_prop_mediated = true_prop_mediated,
    true_iv_ame_on_mediator = true_iv_ame_on_mediator,
    simulation_params = list(n = n,
                             exposure_name = exposure_name,
                             mediator_name = mediator_name,
                             iv_name = iv_name,
                             phi_adl = phi_adl,
                             iv_strength = iv_strength,
                             prop_insured = mean(insured),
                             prop_rehabIRF = mean(rehabIRF),
                             prop_IV1 = mean(IV),
                             mean_ps_mediator = mean(ps_mediator),
                             range_ps_mediator = range(ps_mediator)),
    fixed_covariates = data.frame(age = age,
                                  age_scaled = age_scaled,
                                  post_discharge_disability = post_discharge_disability,
                                  post_discharge_disability_scaled = post_discharge_disability_scaled,
                                  stroke_severity = stroke_severity,
                                  stroke_severity_scaled = stroke_severity_scaled,
                                  comorbidity = comorbidity,
                                  comorbidity_scaled = comorbidity_scaled)))
}





# Updated simulation function with composite outcome model
run_single_simulation_iv_mediation <- function(n = 1000,
                                               exposure_name = "insured",
                                               mediator_name = "rehabIRF",
                                               iv_name = "Z",
                                               iv_strength = 0.3,
                                               trunc_level = 0.01,
                                               fixed_covariates = NULL) {
  
  sim_data <- generate_stroke_data_iv_mediation(n = n,
                                                exposure_name = exposure_name,
                                                mediator_name = mediator_name,
                                                iv_name = iv_name,
                                                iv_strength = iv_strength,
                                                fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  
  # Re-scale outcomes
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  N_adl <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N_adl - 1) + 0.5) / N_adl
  
  # Create composite outcome: M * A
  data$composite_MA <- data[[mediator_name]] * data[[exposure_name]]
  
  ### Fit models
  # Propensity score for A
  ps_formula <- as.formula(paste(exposure_name,"~ age + stroke_severity + comorbidity"))
  ps_mod <- glm(ps_formula,data = data,family = binomial())
  
  # Mediator model P(M|W,A) - WITHOUT instrument (for reduced form NIE)
  mediator_formula_no_iv <- as.formula(paste(mediator_name,
                                             "~ age + stroke_severity + comorbidity +",
                                             exposure_name))
  mediator_mod_no_iv <- glm(mediator_formula_no_iv,data = data,family = binomial())
  
  # Composite model P(M*A=1|W,Z) - for first stage
  composite_formula <- as.formula(paste("composite_MA ~ age + stroke_severity + comorbidity +",
                                        iv_name))
  composite_mod <- glm(composite_formula,data = data,family = binomial())
  
  # IV model P(Z|W)
  iv_formula <- as.formula(paste(iv_name,"~ age + stroke_severity + comorbidity"))
  iv_mod <- glm(iv_formula,data = data,family = binomial())
  
  # Outcome model E[Y|W,A,Z]
  outcome_formula <- as.formula(paste(
    "OUT3_ADL_IADL_01 ~ age + stroke_severity + comorbidity +",
    exposure_name,"+",iv_name,"+",exposure_name,"*",iv_name
  ))
  outcome_mod <- betareg(outcome_formula,data = data)
  
  ### Calculate NIE using IV approach
  nie_adl <- NA
  nie_se_adl <- NA
  nie_ci_lower_adl <- NA
  nie_ci_upper_adl <- NA
  nie_reduced <- NA
  first_stage_est <- NA
  nie_converged <- FALSE
  
  tryCatch({
    results_nie <- robustcompNIE_IV(outcome_mod = outcome_mod,
                                    mediator_mod_no_iv = mediator_mod_no_iv,
                                    composite_mod = composite_mod,
                                    iv_mod = iv_mod,
                                    ps_mod = ps_mod,
                                    covariates = data,
                                    exposure_name = exposure_name,
                                    mediator_name = mediator_name,
                                    iv_name = iv_name,
                                    outcome_name = "OUT3_ADL_IADL_01",
                                    composite_outcome_name = "composite_MA",
                                    trunc_level = trunc_level,
                                    scale_factor = 3.0)
    
    nie_adl <- results_nie$NIE
    nie_se_adl <- results_nie$SE
    nie_ci_lower_adl <- results_nie$CI["lower"]
    nie_ci_upper_adl <- results_nie$CI["upper"]
    nie_reduced <- results_nie$NIE_reduced
    first_stage_est <- results_nie$first_stage
    nie_converged <- TRUE
  },error = function(e) {
    warning(paste("NIE estimation failed:",e$message))
  })
  
  # Calculate bias and coverage
  if (nie_converged && !is.na(nie_adl) && is.finite(nie_adl)) {
    nie_bias <- nie_adl - sim_data$true_nie
    nie_covers <- (nie_ci_lower_adl <= sim_data$true_nie) & 
      (sim_data$true_nie <= nie_ci_upper_adl)
  } else {
    nie_bias <- NA
    nie_covers <- NA
  }
  
  # IV diagnostics - first stage on composite
  first_stage_lm <- lm(as.formula(paste("composite_MA ~ age + stroke_severity + comorbidity +",
                                        iv_name)),
                       data = data)
  fs_summary <- summary(first_stage_lm)
  iv_idx <- which(names(coef(first_stage_lm)) == iv_name)
  iv_coef <- coef(first_stage_lm)[iv_idx]
  iv_tstat <- fs_summary$coefficients[iv_idx,"t value"]
  iv_fstat <- iv_tstat^2
  
  data.frame(n = n,
             iv_strength = iv_strength,
             trunc_level = trunc_level,
             nie_adl = nie_adl,
             nie_se_adl = nie_se_adl,
             nie_ci_lower_adl = nie_ci_lower_adl,
             nie_ci_upper_adl = nie_ci_upper_adl,
             nie_bias = nie_bias,
             nie_covers = nie_covers,
             nie_converged = nie_converged,
             nie_reduced = nie_reduced,
             first_stage_est = first_stage_est,
             true_nie_adl = sim_data$true_nie,
             true_nde_adl = sim_data$true_nde,
             true_te_adl = sim_data$true_te,
             true_iv_ame_on_mediator = sim_data$true_iv_ame_on_mediator,
             prop_insured = mean(data[[exposure_name]]),
             prop_rehabIRF = mean(data[[mediator_name]]),
             prop_IV = mean(data[[iv_name]]),
             prop_composite = mean(data$composite_MA),
             iv_first_stage_coef = iv_coef,
             iv_first_stage_fstat = iv_fstat)
}


# Full simulation study
run_simulation_study_iv_mediation <- function(n_sims = 100,
                                              n = 1000,
                                              iv_strength = 0.3,
                                              trunc_level = 0.05) {
  
  results_list <- vector("list",n_sims)
  
  pb <- txtProgressBar(min = 0,max = n_sims,style = 3)
  
  for (i in 1:n_sims) {
    results_list[[i]] <- run_single_simulation_iv_mediation(
      n = n,
      iv_strength = iv_strength,
      trunc_level = trunc_level
    )
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  results_df <- do.call(rbind,results_list)
  
  # Summary statistics
  summary_stats <- data.frame(n = n,
                              n_sims = n_sims,
                              iv_strength = iv_strength,
                              n_converged_nie = sum(results_df$nie_converged,na.rm = TRUE),
                              nie_mean = mean(results_df$nie_adl,na.rm = TRUE),
                              nie_true = mean(results_df$true_nie_adl,na.rm = TRUE),
                              nie_bias = mean(results_df$nie_bias,na.rm = TRUE),
                              nie_empirical_se = sd(results_df$nie_adl,na.rm = TRUE),
                              nie_mean_se = mean(results_df$nie_se_adl,na.rm = TRUE),
                              nie_coverage = mean(results_df$nie_covers,na.rm = TRUE),
                              nie_rmse = sqrt(mean(results_df$nie_bias^2,na.rm = TRUE)),
                              mean_reduced_form = mean(results_df$nie_reduced,na.rm = TRUE),
                              sd_reduced_form = sd(results_df$nie_reduced,na.rm = TRUE),
                              mean_first_stage = mean(results_df$first_stage_est,na.rm = TRUE),
                              sd_first_stage = sd(results_df$first_stage_est,na.rm = TRUE),
                              mean_iv_fstat = mean(results_df$iv_first_stage_fstat,na.rm = TRUE),
                              mean_iv_coef = mean(results_df$iv_first_stage_coef,na.rm = TRUE),
                              mean_true_iv_ame = mean(results_df$true_iv_ame_on_mediator,na.rm = TRUE),
                              mean_prop_exposure = mean(results_df$prop_insured,na.rm = TRUE),
                              mean_prop_mediator = mean(results_df$prop_rehabIRF,na.rm = TRUE),
                              mean_prop_iv = mean(results_df$prop_IV,na.rm = TRUE),
                              mean_prop_composite = mean(results_df$prop_composite,na.rm = TRUE))
  
  return(list(results = results_df,summary = summary_stats))
}


# Run simulation
# set.seed(80924)
cat("\n=== Running IV-Mediation Simulation Study ===\n")
sim_study_iv <- run_simulation_study_iv_mediation(n_sims = 1000,n = 10000,iv_strength = 0.4)

# Display results
cat("\n=== IV-Mediation Performance Summary ===\n")

cat("\nData Characteristics:\n")
cat("  Mean P(Exposure=1):   ",round(sim_study_iv$summary$mean_prop_exposure,3),"\n")
cat("  Mean P(Mediator=1):   ",round(sim_study_iv$summary$mean_prop_mediator,3),"\n")
cat("  Mean P(Instrument=1): ",round(sim_study_iv$summary$mean_prop_iv,3),"\n")
cat("  Mean P(Composite=1):  ",round(sim_study_iv$summary$mean_prop_composite,3),"\n")

cat("\nReduced Form (Numerator):\n")
cat("  Mean estimate:        ",round(sim_study_iv$summary$mean_reduced_form,4),"\n")
cat("  Empirical SE:         ",round(sim_study_iv$summary$sd_reduced_form,4),"\n")

cat("\nFirst Stage (Denominator):\n")
cat("  Mean estimate:        ",round(sim_study_iv$summary$mean_first_stage,4),"\n")
cat("  Empirical SE:         ",round(sim_study_iv$summary$sd_first_stage,4),"\n")

cat("\nNatural Indirect Effect (NIE) - IV Approach:\n")
cat("  True NIE:             ",round(sim_study_iv$summary$nie_true,4),"\n")
cat("  Mean estimate:        ",round(sim_study_iv$summary$nie_mean,4),"\n")
cat("  Bias:                 ",round(sim_study_iv$summary$nie_bias,4),"\n")
cat("  Empirical SE:         ",round(sim_study_iv$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study_iv$summary$nie_mean_se,4),"\n")
cat("  Coverage:             ",round(sim_study_iv$summary$nie_coverage,3),"\n")
cat("  Convergence rate:     ",round(sim_study_iv$summary$n_converged_nie / sim_study_iv$summary$n_sims,3),"\n")

cat("\nIV Diagnostics:\n")
cat("  Mean first-stage F-stat:      ",round(sim_study_iv$summary$mean_iv_fstat,2),"\n")
cat("  Mean IV coefficient (LPM):    ",round(sim_study_iv$summary$mean_iv_coef,4),"\n")
cat("  Mean true IV AME on mediator: ",round(sim_study_iv$summary$mean_true_iv_ame,4),"\n")












