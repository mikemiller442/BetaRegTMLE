








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



# Calculates a robust TMLE estimator for the NDE
robustcompNDE <- function(outcome_mod,mediator_mod,ps_mod,covariates,
                          exposure_name,mediator_name,outcome_name,
                          trunc_level,scale_factor = 1.0) {
  
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
  
  ### calculate NDE parameter using targeted outcome regression
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
  
  ### calculate targeted estimate of NDE using regression-based estimator
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
  H_Z <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
  min_obs_diff <- min(obs_diff)
  max_obs_diff <- max(obs_diff)
  range_obs_diff <- max_obs_diff - min_obs_diff
  # avoid division by zero
  if (range_obs_diff < 1e-10) {
    range_obs_diff <- 1
  }
  obs_diff_scaled <- (obs_diff - min_obs_diff) / range_obs_diff
  # apply transformation to avoid boundary
  obs_diff_scaled <- (obs_diff_scaled * (n - 1) + 0.5) / n
  # have to also transform the g-computation estimates
  nde_cond_scaled <- (nde_cond - min_obs_diff) / range_obs_diff
  nde_cond_scaled <- (nde_cond_scaled * (n - 1) + 0.5) / n
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
  nde_cond_star <- nde_cond_star_scaled * range_obs_diff + min_obs_diff
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
  obs_diff_star <- obs_diff_star * range_obs_diff + min_obs_diff
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
                          trunc_level,scale_factor = 1.0) {
  
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
  
  ### calculate NIE parameter using targeted outcome regression
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
  
  ### calculate targeted estimate of NIE using regression-based estimator
  covariates_a1_obs <- covariates
  covariates_a1_obs[[exposure_name]] <- 1
  pred_y_a1_obs <- predict(outcome_mod,newdata = covariates_a1_obs,type = "response")
  
  Q_Z_0_obs <- ifelse(mediator == 1,med_prob_0,1 - med_prob_0)
  Q_Z_1_obs <- ifelse(mediator == 1,med_prob_1,1 - med_prob_1)
  H_Y_a1_obs <- (1 - Q_Z_0_obs / Q_Z_1_obs) / ps_trunc
  
  Q_bar_Y_1Z_star <- plogis(qlogis(pred_y_a1_obs) + eps_y * H_Y_a1_obs)
  
  ### calculate targeted estimate of NIE
  H_Z_1 <- ifelse(exposure == 1,1 / ps_trunc,0)
  H_Z_0 <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
  # find common range for scaling
  all_preds <- c(psi_nie_z_0,psi_nie_z_1,Q_bar_Y_1Z_star)
  min_all <- min(all_preds)
  max_all <- max(all_preds)
  range_all <- max_all - min_all
  if (range_all < 1e-10) {
    range_all <- 1
  }
  scale_pred <- function(x) {
    x_scaled <- (x - min_all) / range_all
    x_scaled <- (x_scaled * (n - 1) + 0.5) / n
    x_scaled <- pmax(pmin(x_scaled,1 - 1e-10),1e-10)
    return(x_scaled)
  }
  psi_nie_z_0_scaled <- scale_pred(psi_nie_z_0)
  psi_nie_z_1_scaled <- scale_pred(psi_nie_z_1)
  Q_bar_Y_1Z_star_scaled <- scale_pred(Q_bar_Y_1Z_star)
  # calculate regression model for A=1
  data_tmle_a1 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star_scaled,
                             pred_vals = psi_nie_z_1_scaled,
                             H_Z_1 = H_Z_1)
  data_tmle_a1_subset <- data_tmle_a1[exposure == 1,]
  eps_nie_1 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_1,
                        data = data_tmle_a1_subset,
                        family = "binomial"))
  # calculate regression model for A=0
  data_tmle_a0 <- data.frame(pseudo_outcome = Q_bar_Y_1Z_star_scaled,
                             pred_vals = psi_nie_z_0_scaled,
                             H_Z_0 = H_Z_0)
  data_tmle_a0_subset <- data_tmle_a0[exposure == 0,]
  eps_nie_0 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_0,
                        data = data_tmle_a0_subset,
                        family = "binomial"))
  
  # calculate targeted estimates for all observations
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



### Generate simulated data with mediation structure
generate_stroke_data_mediation <- function(n = 1000,
                                           exposure_name = "insured",
                                           mediator_name = "rehabIRF",
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
  
  ### generate exposure assignment (insured)
  # weak associations with covariates, ~70% insured
  logit_ps_exposure <- 0.85 +
    0.05 * age_scaled + 
    0.08 * post_discharge_disability_scaled +
    0.06 * stroke_severity_scaled + 
    -0.04 * comorbidity_scaled 
  
  ps_exposure <- plogis(logit_ps_exposure)
  insured <- rbinom(n,size = 1,prob = ps_exposure)
  
  ### generate mediator assignment (rehabIRF)
  # propensity score depends on covariates AND exposure (insured)
  # higher post_discharge_disability -> higher treatment probability
  # insured -> higher probability of rehabIRF
  logit_ps_mediator <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled +
    0.6 * insured
  
  ps_mediator <- plogis(logit_ps_mediator)
  rehabIRF <- rbinom(n,size = 1,prob = ps_mediator)
  
  ### generate outcomes using beta regression structure
  # true data generating mechanism for mean
  # Key mediation structure:
  # - direct effect of insured
  # - indirect effect through rehabIRF
  # - interaction between insured and rehabIRF
  # - NO three-way interaction with post_discharge_disability
  mu_adl_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * insured +
    -0.5 * rehabIRF +
    -0.25 * insured * rehabIRF +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  
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
      0.4 * insured
    
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
  
  ### calculate true causal effects for validation
  
  # Total Effect (TE): E[Y(1, M(1)) - Y(0, M(0))]
  # Natural Direct Effect (NDE): E[Y(1, M(0)) - Y(0, M(0))]
  # Natural Indirect Effect (NIE): E[Y(1, M(1)) - Y(1, M(0))]
  
  # For each individual, calculate potential outcomes under all combinations
  # Y(a, m) where a = insured, m = rehabIRF
  
  # Y(0, 0): not insured, no rehab
  mu_00_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 0 + 
    -0.5 * 0 + 
    -0.25 * 0 * 0 +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  mu_00 <- plogis(mu_00_logit) * 3.0
  
  # Y(0, 1): not insured, with rehab
  mu_01_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 0 + 
    -0.5 * 1 + 
    -0.25 * 0 * 1 +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  mu_01 <- plogis(mu_01_logit) * 3.0
  
  # Y(1, 0): insured, no rehab
  mu_10_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 1 + 
    -0.5 * 0 + 
    -0.25 * 1 * 0 +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  mu_10 <- plogis(mu_10_logit) * 3.0
  
  # Y(1, 1): insured, with rehab
  mu_11_logit <- -1.0 + 
    1.0 * post_discharge_disability_scaled +
    -0.3 * 1 + 
    -0.5 * 1 + 
    -0.25 * 1 * 1 +
    0.4 * stroke_severity_scaled + 
    0.2 * age_scaled + 
    0.1 * comorbidity_scaled
  mu_11 <- plogis(mu_11_logit) * 3.0
  
  # For NIE and NDE, we need E[M(a)] = P(M=1|A=a)
  # These are the propensity scores for mediator under each exposure level
  
  # P(M=1|A=0, X)
  logit_ps_m0 <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled +
    0.6 * 0
  ps_m0 <- plogis(logit_ps_m0)
  
  # P(M=1|A=1, X)
  logit_ps_m1 <- -0.8 +
    0.15 * age_scaled + 
    0.8 * post_discharge_disability_scaled +
    0.15 * stroke_severity_scaled + 
    -0.1 * comorbidity_scaled +
    0.6 * 1
  ps_m1 <- plogis(logit_ps_m1)
  
  # Total Effect: E[Y(1, M(1)) - Y(0, M(0))]
  # = E[Y(1,1)*P(M=1|A=1) + Y(1,0)*P(M=0|A=1)] - E[Y(0,1)*P(M=1|A=0) + Y(0,0)*P(M=0|A=0)]
  true_te <- mean((mu_11 * ps_m1 + mu_10 * (1 - ps_m1)) - 
                    (mu_01 * ps_m0 + mu_00 * (1 - ps_m0)))
  
  # Natural Direct Effect: E[Y(1, M(0)) - Y(0, M(0))]
  # = E[Y(1,1)*P(M=1|A=0) + Y(1,0)*P(M=0|A=0)] - E[Y(0,1)*P(M=1|A=0) + Y(0,0)*P(M=0|A=0)]
  true_nde <- mean((mu_11 * ps_m0 + mu_10 * (1 - ps_m0)) - 
                     (mu_01 * ps_m0 + mu_00 * (1 - ps_m0)))
  
  # Natural Indirect Effect: E[Y(1, M(1)) - Y(1, M(0))]
  # = E[Y(1,1)*P(M=1|A=1) + Y(1,0)*P(M=0|A=1)] - E[Y(1,1)*P(M=1|A=0) + Y(1,0)*P(M=0|A=0)]
  true_nie <- mean((mu_11 * ps_m1 + mu_10 * (1 - ps_m1)) - 
                     (mu_11 * ps_m0 + mu_10 * (1 - ps_m0)))
  
  # Proportion mediated
  true_prop_mediated <- true_nie / true_te
  
  ### compile dataset
  data <- data.frame(id = 1:n,
                     age = age,
                     post_discharge_disability = post_discharge_disability,
                     stroke_severity = stroke_severity,
                     comorbidity = comorbidity,
                     insured = insured,
                     rehabIRF = rehabIRF,
                     OUT3_ADL_IADL = OUT3_ADL_IADL_with_missing,
                     OUT3_ADL_IADL_complete = OUT3_ADL_IADL,
                     observed_adl = observed_adl,
                     prob_observed_adl = prob_observed_adl,
                     propensity_score_exposure = ps_exposure,
                     propensity_score_mediator = ps_mediator)
  
  ### return results
  return(list(data = data,
              true_te = true_te,
              true_nde = true_nde,
              true_nie = true_nie,
              true_prop_mediated = true_prop_mediated,
              simulation_params = list(n = n,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       phi_adl = phi_adl,
                                       missing_outcome = missing_outcome,
                                       prob_observed_baseline = prob_observed_baseline,
                                       actual_missing_prop_adl = mean(observed_adl == 0),
                                       prop_insured = mean(insured),
                                       prop_rehabIRF = mean(rehabIRF)),
              fixed_covariates = data.frame(age = age,
                                            age_scaled = age_scaled,
                                            post_discharge_disability = post_discharge_disability,
                                            post_discharge_disability_scaled = post_discharge_disability_scaled,
                                            stroke_severity = stroke_severity,
                                            stroke_severity_scaled = stroke_severity_scaled,
                                            comorbidity = comorbidity,
                                            comorbidity_scaled = comorbidity_scaled)))
}

### Example usage
set.seed(80924)
sim_data <- generate_stroke_data_mediation(n = 1000)

# Check key properties
with(sim_data,{
  cat("Sample size:",simulation_params$n,"\n")
  cat("Proportion insured:",round(simulation_params$prop_insured,3),"\n")
  cat("Proportion in rehabIRF:",round(simulation_params$prop_rehabIRF,3),"\n")
  cat("\nTrue causal effects:\n")
  cat("  Total Effect:",round(true_te,4),"\n")
  cat("  Natural Direct Effect:",round(true_nde,4),"\n")
  cat("  Natural Indirect Effect:",round(true_nie,4),"\n")
  cat("  Proportion Mediated:",round(true_prop_mediated,4),"\n")
})

# Quick check of associations
summary(glm(insured ~ age + post_discharge_disability + stroke_severity + comorbidity,
            data = sim_data$data,family = binomial))
summary(glm(rehabIRF ~ age + post_discharge_disability + stroke_severity + comorbidity + insured,
            data = sim_data$data,family = binomial))











# Add NIE calculation to the simulation function
run_single_simulation_mediation_full <- function(n = 1000,
                                                 exposure_name = "insured",
                                                 mediator_name = "rehabIRF",
                                                 trunc_level = 0.05,
                                                 missing_outcome = FALSE,
                                                 prob_observed_baseline = 0.8,
                                                 fixed_covariates = NULL) {
  
  # Generate data
  sim_data <- generate_stroke_data_mediation(n = n,
                                             exposure_name = exposure_name,
                                             mediator_name = mediator_name,
                                             missing_outcome = missing_outcome,
                                             prob_observed_baseline = prob_observed_baseline,
                                             fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  observed_adl <- data$observed_adl
  
  # Re-scale outcomes to 0-1 for modeling
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  N_adl <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N_adl - 1) + 0.5) / N_adl
  
  # Fit models
  ps_exposure_formula <- as.formula(paste(exposure_name,
                                          "~ age + post_discharge_disability + stroke_severity + comorbidity"))
  ps_exposure_model <- glm(ps_exposure_formula,data = data,family = binomial())
  
  mediator_formula <- as.formula(paste(mediator_name,
                                       "~ age + post_discharge_disability + stroke_severity + comorbidity +",
                                       exposure_name))
  mediator_model <- glm(mediator_formula,data = data,family = binomial())
  
  outcome_formula <- as.formula(paste(
    "OUT3_ADL_IADL_01 ~ age + post_discharge_disability + stroke_severity + comorbidity +",
    exposure_name,"+",mediator_name,"+",
    exposure_name,"*",mediator_name
  ))
  adl_mod <- betareg(outcome_formula,data = data)
  
  # Calculate NDE
  tryCatch({
    results_nde <- robustcompNDE(outcome_mod = adl_mod,
                                 mediator_mod = mediator_model,
                                 ps_mod = ps_exposure_model,
                                 covariates = data,
                                 exposure_name = exposure_name,
                                 mediator_name = mediator_name,
                                 outcome_name = "OUT3_ADL_IADL_01",
                                 trunc_level = trunc_level,
                                 scale_factor = 3.0)
    
    nde_adl <- results_nde$NDE
    nde_se_adl <- results_nde$SE
    nde_ci_lower_adl <- results_nde$CI["lower"]
    nde_ci_upper_adl <- results_nde$CI["upper"]
    nde_converged <- TRUE
  },error = function(e) {
    warning(paste("NDE estimation failed:",e$message))
    nde_adl <<- NA
    nde_se_adl <<- NA
    nde_ci_lower_adl <<- NA
    nde_ci_upper_adl <<- NA
    nde_converged <<- FALSE
  })
  
  # Calculate NIE
  tryCatch({
    results_nie <- robustcompNIE(outcome_mod = adl_mod,
                                 mediator_mod = mediator_model,
                                 ps_mod = ps_exposure_model,
                                 covariates = data,
                                 exposure_name = exposure_name,
                                 mediator_name = mediator_name,
                                 outcome_name = "OUT3_ADL_IADL_01",
                                 trunc_level = trunc_level,
                                 scale_factor = 3.0)
    
    nie_adl <- results_nie$NIE
    nie_se_adl <- results_nie$SE
    nie_ci_lower_adl <- results_nie$CI["lower"]
    nie_ci_upper_adl <- results_nie$CI["upper"]
    nie_converged <- TRUE
  },error = function(e) {
    warning(paste("NIE estimation failed:",e$message))
    nie_adl <<- NA
    nie_se_adl <<- NA
    nie_ci_lower_adl <<- NA
    nie_ci_upper_adl <<- NA
    nie_converged <<- FALSE
  })
  
  # Calculate biases and coverage
  if (nde_converged) {
    nde_bias <- nde_adl - sim_data$true_nde
    nde_covers <- (nde_ci_lower_adl <= sim_data$true_nde) & 
      (sim_data$true_nde <= nde_ci_upper_adl)
  } else {
    nde_bias <- NA
    nde_covers <- NA
  }
  
  if (nie_converged) {
    nie_bias <- nie_adl - sim_data$true_nie
    nie_covers <- (nie_ci_lower_adl <= sim_data$true_nie) & 
      (sim_data$true_nie <= nie_ci_upper_adl)
  } else {
    nie_bias <- NA
    nie_covers <- NA
  }
  
  # Extract results
  results_df <- data.frame(n = n,
                           exposure_name = exposure_name,
                           mediator_name = mediator_name,
                           trunc_level = trunc_level,
                           missing_outcome = missing_outcome,
                           # NDE results
                           nde_adl = nde_adl,
                           nde_se_adl = nde_se_adl,
                           nde_ci_lower_adl = nde_ci_lower_adl,
                           nde_ci_upper_adl = nde_ci_upper_adl,
                           nde_bias = nde_bias,
                           nde_covers = nde_covers,
                           nde_converged = nde_converged,
                           # NIE results
                           nie_adl = nie_adl,
                           nie_se_adl = nie_se_adl,
                           nie_ci_lower_adl = nie_ci_lower_adl,
                           nie_ci_upper_adl = nie_ci_upper_adl,
                           nie_bias = nie_bias,
                           nie_covers = nie_covers,
                           nie_converged = nie_converged,
                           # True values
                           true_nde_adl = sim_data$true_nde,
                           true_nie_adl = sim_data$true_nie,
                           true_te_adl = sim_data$true_te,
                           true_prop_mediated = sim_data$true_prop_mediated,
                           # Sample characteristics
                           prop_insured = mean(data[[exposure_name]]),
                           prop_rehabIRF = mean(data[[mediator_name]]))
  
  return(results_df)
}

# Test it with full mediation analysis
test_sim_full <- run_single_simulation_mediation_full(n = 1000)

cat("\n=== Full Mediation Analysis Results ===\n")
cat("\nTrue effects:\n")
cat("  Total Effect:",round(test_sim_full$true_te_adl,4),"\n")
cat("  Natural Direct Effect:",round(test_sim_full$true_nde_adl,4),"\n")
cat("  Natural Indirect Effect:",round(test_sim_full$true_nie_adl,4),"\n")
cat("  Proportion Mediated:",round(test_sim_full$true_prop_mediated,4),"\n")

cat("\nEstimated effects:\n")
cat("  NDE:",round(test_sim_full$nde_adl,4),
    " (Bias:",round(test_sim_full$nde_bias,4),")\n")
cat("  NIE:",round(test_sim_full$nie_adl,4),
    " (Bias:",round(test_sim_full$nie_bias,4),")\n")
cat("  Sum (NDE+NIE):",round(test_sim_full$nde_adl + test_sim_full$nie_adl,4),"\n")

### Run multiple simulations to assess performance for both NDE and NIE
run_simulation_study_mediation_full <- function(n_sims = 100,
                                                n = 1000,
                                                trunc_level = 0.05) {
  
  results_list <- vector("list",n_sims)
  
  # Progress bar
  pb <- txtProgressBar(min = 0,max = n_sims,style = 3)
  
  for (i in 1:n_sims) {
    results_list[[i]] <- run_single_simulation_mediation_full(n = n,
                                                              trunc_level = trunc_level)
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  # Combine results
  results_df <- do.call(rbind,results_list)
  
  # Calculate summary statistics
  summary_stats <- data.frame(n = n,
                              n_sims = n_sims,
                              n_converged_nde = sum(results_df$nde_converged,na.rm = TRUE),
                              n_converged_nie = sum(results_df$nie_converged,na.rm = TRUE),
                              # NDE performance
                              nde_mean = mean(results_df$nde_adl,na.rm = TRUE),
                              nde_true = unique(results_df$true_nde_adl)[1],
                              nde_bias = mean(results_df$nde_bias,na.rm = TRUE),
                              nde_empirical_se = sd(results_df$nde_adl,na.rm = TRUE),
                              nde_mean_se = mean(results_df$nde_se_adl,na.rm = TRUE),
                              nde_coverage = mean(results_df$nde_covers,na.rm = TRUE),
                              nde_rmse = sqrt(mean(results_df$nde_bias^2,na.rm = TRUE)),
                              # NIE performance
                              nie_mean = mean(results_df$nie_adl,na.rm = TRUE),
                              nie_true = unique(results_df$true_nie_adl)[1],
                              nie_bias = mean(results_df$nie_bias,na.rm = TRUE),
                              nie_empirical_se = sd(results_df$nie_adl,na.rm = TRUE),
                              nie_mean_se = mean(results_df$nie_se_adl,na.rm = TRUE),
                              nie_coverage = mean(results_df$nie_covers,na.rm = TRUE),
                              nie_rmse = sqrt(mean(results_df$nie_bias^2,na.rm = TRUE)),
                              # Total effect from sum
                              te_from_sum_mean = mean(results_df$nde_adl + results_df$nie_adl,na.rm = TRUE),
                              te_true = unique(results_df$true_te_adl)[1])
  
  return(list(results = results_df,
              summary = summary_stats))
}

### Example: Run simulation study
cat("\n=== Running Full Mediation Simulation Study ===\n")
sim_study_full <- run_simulation_study_mediation_full(n_sims = 2000,n = 7000)

cat("\nSimulation Study Summary:\n")
print(sim_study_full$summary)

# Print detailed summary table
cat("\n=== Detailed Performance Summary ===\n")
cat("\nNatural Direct Effect (NDE):\n")
cat("  True value:        ",round(sim_study_full$summary$nde_true,4),"\n")
cat("  Mean estimate:     ",round(sim_study_full$summary$nde_mean,4),"\n")
cat("  Bias:              ",round(sim_study_full$summary$nde_bias,4),"\n")
cat("  Empirical SE:      ",round(sim_study_full$summary$nde_empirical_se,4),"\n")
cat("  Mean estimated SE: ",round(sim_study_full$summary$nde_mean_se,4),"\n")
cat("  Coverage:          ",round(sim_study_full$summary$nde_coverage,3),"\n")
cat("  Convergence rate:  ",round(sim_study_full$summary$n_converged_nde / sim_study_full$summary$n_sims,3),"\n")

cat("\nNatural Indirect Effect (NIE):\n")
cat("  True value:        ",round(sim_study_full$summary$nie_true,4),"\n")
cat("  Mean estimate:     ",round(sim_study_full$summary$nie_mean,4),"\n")
cat("  Bias:              ",round(sim_study_full$summary$nie_bias,4),"\n")
cat("  Empirical SE:      ",round(sim_study_full$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE: ",round(sim_study_full$summary$nie_mean_se,4),"\n")
cat("  Coverage:          ",round(sim_study_full$summary$nie_coverage,3),"\n")
cat("  Convergence rate:  ",round(sim_study_full$summary$n_converged_nie / sim_study_full$summary$n_sims,3),"\n")

cat("\nTotal Effect Check:\n")
cat("  True TE:           ",round(sim_study_full$summary$te_true,4),"\n")
cat("  Mean(NDE + NIE):   ",round(sim_study_full$summary$te_from_sum_mean,4),"\n")
cat("  Difference:        ",round(sim_study_full$summary$te_from_sum_mean - sim_study_full$summary$te_true,4),"\n")












