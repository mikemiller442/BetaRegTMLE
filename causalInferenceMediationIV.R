








### Load libraries
library(betareg)


# Calculates a robust TMLE estimator for the ATE with missing at random outcomes
robustcompATE_MAR <- function(outcome_mod,ps_mod,missing_mod,
                              covariates_full,exposure_name,
                              outcome_observed,outcome_values,outcome_name,
                              trunc_level,scale_factor = 1.0) {
  
  n <- nrow(covariates_full)
  exposure <- covariates_full[[exposure_name]]
  
  ### calculate propensity score for exposure (instrument)
  ps_pred <- predict(ps_mod,newdata = covariates_full,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  
  ### calculate probability of being observed (exposure = 1)
  prob_observed <- predict(missing_mod,newdata = covariates_full,type = "response")
  prob_observed_trunc <- pmax(pmin(prob_observed,1.0 - trunc_level),trunc_level)
  
  ### calculate clever covariates for TMLE
  trunc_weights_0 <- ifelse(exposure == 0 & outcome_observed == 1,
                            -1.0 / ((1.0 - ps_trunc) * prob_observed_trunc),0)
  trunc_weights_1 <- ifelse(exposure == 1 & outcome_observed == 1,
                            1.0 / (ps_trunc * prob_observed_trunc),0)
  H_0_Y <- -1.0 * ifelse(outcome_observed == 1,1 / ((1.0 - ps_trunc) * prob_observed_trunc),0)
  H_1_Y <- ifelse(outcome_observed == 1,1 / (ps_trunc * prob_observed_trunc),0)
  H_A_Y <- ifelse(exposure == 1,H_1_Y,H_0_Y)
  
  ### estimate fluctuation parameter using observed data only
  pred_y <- predict(outcome_mod,newdata = covariates_full,type = "response")
  obs_idx <- which(outcome_observed == 1)
  
  data_tmle <- data.frame(Y = outcome_values[obs_idx],
                          pred_vals = pred_y[obs_idx],
                          H_A_Y = H_A_Y[obs_idx])
  eps <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_A_Y,
                  data = data_tmle,family = "binomial"))
  
  ### update predictions using fluctuation parameter
  covariates_0 <- covariates_full
  covariates_0[[exposure_name]] <- 0
  covariates_1 <- covariates_full
  covariates_1[[exposure_name]] <- 1
  
  pred_y_0 <- predict(outcome_mod,newdata = covariates_0,type = "response")
  pred_y_1 <- predict(outcome_mod,newdata = covariates_1,type = "response")
  pred_y_A <- predict(outcome_mod,newdata = covariates_full,type = "response")
  
  pred_star_y_A <- plogis(qlogis(pred_y_A) + eps * H_A_Y)
  pred_star_y_0 <- plogis(qlogis(pred_y_0) + eps * H_0_Y)
  pred_star_y_1 <- plogis(qlogis(pred_y_1) + eps * H_1_Y)
  
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


# Calculates the reduced form for the complier subgroup indirect effect
robustcompNIE_reduced_form <- function(outcome_mod,mediator_mod,iv_mod,ps_mod,covariates,
                                       exposure_name,mediator_name,iv_name,outcome_name,
                                       trunc_level,scale_factor = 1.0) {
  
  n <- nrow(covariates)
  exposure <- covariates[[exposure_name]]
  mediator <- covariates[[mediator_name]]
  iv <- covariates[[iv_name]]
  outcome <- covariates[[outcome_name]]
  
  ### calculate propensity score P(A=1|W)
  ps_pred <- predict(ps_mod,type = "response")
  ps_trunc <- pmax(pmin(ps_pred,1.0 - trunc_level),trunc_level)
  
  ### calculate P(Z=1|W) for instrument density
  iv_prob <- predict(iv_mod,newdata = covariates,type = "response")
  iv_prob_trunc <- pmax(pmin(iv_prob,1.0 - trunc_level),trunc_level)
  
  # P(Z=z|W) for observed Z values
  pZ_at_obs <- ifelse(iv == 1,iv_prob_trunc,1 - iv_prob_trunc)
  
  ### calculate P(M|A,W) using mediator model
  # predict P(M=1|A=0,W)
  covariates_a0 <- covariates
  covariates_a0[[exposure_name]] <- 0
  pM1_given_a0 <- predict(mediator_mod,newdata = covariates_a0,type = "response")
  
  # predict P(M=1|A=1,W)
  covariates_a1 <- covariates
  covariates_a1[[exposure_name]] <- 1
  pM1_given_a1 <- predict(mediator_mod,newdata = covariates_a1,type = "response")
  
  # P(M=m|A=0,W) evaluated at observed Z
  pM_at_Z_given_a0 <- ifelse(iv == 1,pM1_given_a0,1 - pM1_given_a0)
  # P(M=m|A=1,W) evaluated at observed Z
  pM_at_Z_given_a1 <- ifelse(iv == 1,pM1_given_a1,1 - pM1_given_a1)
  
  ### calculate density ratio for first component
  density_ratio <- pM_at_Z_given_a0 / pZ_at_obs
  H_Y <- ifelse(exposure == 1,
                (1 - density_ratio) / ps_trunc,
                0)
  
  # get outcome predictions using the original model (conditions on Z)
  pred_y <- predict(outcome_mod,newdata = covariates,type = "response")
  
  # fit fluctuation model using observed data
  data_tmle_y <- data.frame(Y = outcome,
                            pred_vals = pred_y,
                            H_Y = H_Y)
  eps_y <- coef(glm(Y ~ -1 + offset(qlogis(pred_vals)) + H_Y,
                    data = data_tmle_y,family = "binomial"))
  
  
  pred_y_star <- plogis(qlogis(pred_y) + eps_y * H_Y)
  
  ### calculate Q_bar_Y(W, A=1, Z) for Z=0 and Z=1
  covariates_a1_z0 <- covariates
  covariates_a1_z0[[exposure_name]] <- 1
  covariates_a1_z0[[iv_name]] <- 0
  pred_y_a1_z0 <- predict(outcome_mod,newdata = covariates_a1_z0,type = "response")
  
  covariates_a1_z1 <- covariates
  covariates_a1_z1[[exposure_name]] <- 1
  covariates_a1_z1[[iv_name]] <- 1
  pred_y_a1_z1 <- predict(outcome_mod,newdata = covariates_a1_z1,type = "response")
  
  # calculate clever covariates for Q_bar at Z=0 and Z=1
  # P(M=0|A=0,W) / P(Z=0|W) for Z=0
  density_ratio_z0 <- (1 - pM1_given_a0) / (1 - iv_prob_trunc)
  H_Y_z0 <- (1 - density_ratio_z0) / ps_trunc
  
  # P(M=1|A=0,W) / P(Z=1|W) for Z=1
  density_ratio_z1 <- pM1_given_a0 / iv_prob_trunc
  H_Y_z1 <- (1 - density_ratio_z1) / ps_trunc
  
  # update Q_bar predictions
  Q_bar_Y_a1_z0_star <- plogis(qlogis(pred_y_a1_z0) + eps_y * H_Y_z0)
  Q_bar_Y_a1_z1_star <- plogis(qlogis(pred_y_a1_z1) + eps_y * H_Y_z1)
  
  ### calculate psi_NIE using mediator density
  # psi_NIE,Z(w, a) = E_M[Q_bar_Y(w, 1, M) | W=w, A=a]
  # For binary M: = P(M=1|W,A=a) * Q_bar_Y(w,1,Z=1) + P(M=0|W,A=a) * Q_bar_Y(w,1,Z=0)
  
  # psi_NIE(W, A=1) = P(M=1|W,A=1) * Q_bar*(W,1,Z=1) + P(M=0|W,A=1) * Q_bar*(W,1,Z=0)
  psi_nie_a1 <- pM1_given_a1 * Q_bar_Y_a1_z1_star + (1 - pM1_given_a1) * Q_bar_Y_a1_z0_star
  # psi_NIE(W, A=0) = P(M=1|W,A=0) * Q_bar*(W,1,Z=1) + P(M=0|W,A=0) * Q_bar*(W,1,Z=0)
  psi_nie_a0 <- pM1_given_a0 * Q_bar_Y_a1_z1_star + (1 - pM1_given_a0) * Q_bar_Y_a1_z0_star
  
  ### target psi_NIE
  # target psi_NIE(W, A=1) and psi_NIE(W, A=0) using observations with A=1 and A=0 respectively
  # substitute M for Z in outcome model and calculate calculate Q_bar at observed M
  covariates_a1_m_obs <- covariates
  covariates_a1_m_obs[[exposure_name]] <- 1
  covariates_a1_m_obs[[iv_name]] <- mediator
  pred_y_a1_m_obs <- predict(outcome_mod,newdata = covariates_a1_m_obs,type = "response")
  
  # clever covariate for Z=M
  density_ratio_at_m <- ifelse(mediator == 1,
                               pM1_given_a0 / iv_prob_trunc,
                               (1 - pM1_given_a0) / (1 - iv_prob_trunc))
  H_Y_at_m <- (1 - density_ratio_at_m) / ps_trunc
  Q_bar_Y_a1_m_star <- plogis(qlogis(pred_y_a1_m_obs) + eps_y * H_Y_at_m)
  
  # clever covariates for targeting psi_NIE
  H_Z_1 <- ifelse(exposure == 1,1 / ps_trunc,0)
  H_Z_0 <- ifelse(exposure == 0,1 / (1 - ps_trunc),0)
  
  # scale psi_nie for beta regression targeting
  min_psi <- min(c(psi_nie_a1,psi_nie_a0))
  max_psi <- max(c(psi_nie_a1,psi_nie_a0))
  range_psi <- max_psi - min_psi
  if (range_psi < 1e-10) range_psi <- 1
  
  scale_psi <- function(x) {
    x_scaled <- (x - min_psi) / range_psi
    (x_scaled * (n - 1) + 0.5) / n
  }
  
  unscale_psi <- function(x_scaled) {
    x_unscaled <- (x_scaled * n - 0.5) / (n - 1)
    x_unscaled * range_psi + min_psi
  }
  
  psi_nie_a1_scaled <- scale_psi(psi_nie_a1)
  psi_nie_a1_scaled <- pmax(pmin(psi_nie_a1_scaled,1 - 1e-6),1e-6)
  
  psi_nie_a0_scaled <- scale_psi(psi_nie_a0)
  psi_nie_a0_scaled <- pmax(pmin(psi_nie_a0_scaled,1 - 1e-6),1e-6)
  
  Q_bar_Y_a1_m_scaled <- scale_psi(Q_bar_Y_a1_m_star)
  Q_bar_Y_a1_m_scaled <- pmax(pmin(Q_bar_Y_a1_m_scaled,1 - 1e-6),1e-6)
  
  # target for A=1 using Q_bar_Y(W,1,M) as pseudo-outcome
  data_tmle_a1 <- data.frame(pseudo_outcome = Q_bar_Y_a1_m_scaled,
                             pred_vals = psi_nie_a1_scaled,
                             H_Z_1 = H_Z_1)
  data_tmle_a1_subset <- data_tmle_a1[exposure == 1,]
  eps_nie_1 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_1,
                        data = data_tmle_a1_subset,family = "binomial"))
  
  # target for A=0 using Q_bar_Y(W,1,M) as pseudo-outcome
  data_tmle_a0 <- data.frame(pseudo_outcome = Q_bar_Y_a1_m_scaled,
                             pred_vals = psi_nie_a0_scaled,
                             H_Z_0 = H_Z_0)
  data_tmle_a0_subset <- data_tmle_a0[exposure == 0,]
  eps_nie_0 <- coef(glm(pseudo_outcome ~ -1 + offset(qlogis(pred_vals)) + H_Z_0,
                        data = data_tmle_a0_subset,family = "binomial"))
  
  # calculate targeted estimates
  psi_nie_a1_star_scaled <- plogis(qlogis(psi_nie_a1_scaled) + eps_nie_1 * H_Z_1)
  psi_nie_a0_star_scaled <- plogis(qlogis(psi_nie_a0_scaled) + eps_nie_0 * H_Z_0)
  
  # transform back to original scale
  psi_nie_a1_star <- unscale_psi(psi_nie_a1_star_scaled)
  psi_nie_a0_star <- unscale_psi(psi_nie_a0_star_scaled)
  Q_bar_Y_a1_m_unscaled <- unscale_psi(Q_bar_Y_a1_m_scaled)
  
  ### calculate reduced form NIE
  nie_vec <- (psi_nie_a1_star - psi_nie_a0_star) * scale_factor
  nie_est <- mean(nie_vec)
  
  ### calculate the efficient influence function
  D_Y <- H_Y * (outcome - pred_y_star)
  psi_nie_A_star <- ifelse(exposure == 1,psi_nie_a1_star,psi_nie_a0_star)
  residual <- Q_bar_Y_a1_m_unscaled - psi_nie_A_star
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
              CI = c(lower = ci_lower,upper = ci_upper),
              psi_nie_a1 = mean(psi_nie_a1_star),
              psi_nie_a0 = mean(psi_nie_a0_star),
              eps_y = eps_y,
              eps_nie_1 = eps_nie_1,
              eps_nie_0 = eps_nie_0))
}


# Robust estimator for the natural indirect effect in the complier subgroup
robustcompNIE_IV <- function(outcome_mod,mediator_mod,mediator_mod_with_iv,iv_mod,
                             ps_mod,covariates,
                             exposure_name,mediator_name,iv_name,outcome_name,
                             trunc_level_num,trunc_level_denom,scale_factor = 1.0) {
  
  n <- nrow(covariates)
  
  ### NIE in the reduced form
  reduced_form <- robustcompNIE_reduced_form(outcome_mod = outcome_mod,
                                             mediator_mod = mediator_mod,
                                             iv_mod = iv_mod,
                                             ps_mod = ps_mod,
                                             covariates = covariates,
                                             exposure_name = exposure_name,
                                             mediator_name = mediator_name,
                                             iv_name = iv_name,
                                             outcome_name = outcome_name,
                                             trunc_level = trunc_level_num,
                                             scale_factor = scale_factor)
  
  
  ### calculate conditional ATE E[M|Z=1,A=1,W] - E[M|Z=0,A=1,W]
  covariates_exposure_1 <- covariates
  covariates_exposure_1[[exposure_name]] <- 1
  first_stage <- robustcompATE_MAR(outcome_mod = mediator_mod_with_iv,
                                   ps_mod = iv_mod,
                                   missing_mod = ps_mod,
                                   covariates_full = covariates_exposure_1,
                                   exposure_name = iv_name,
                                   outcome_observed = covariates[[exposure_name]],
                                   outcome_values = covariates[[mediator_name]],
                                   outcome_name = mediator_name,
                                   trunc_level = trunc_level_denom,
                                   scale_factor = 1.0)
  
  ### ratio estimator
  nie_iv <- reduced_form$NIE / first_stage$ATE
  
  ### delta method for influence function of ratio
  IF_reduced <- reduced_form$IF
  IF_first_stage <- first_stage$IF
  IF_nie <- (IF_reduced - nie_iv * IF_first_stage) / first_stage$ATE
  
  ### calculate standard errors using delta method influence function
  var_nie_iv <- var(IF_nie) / n
  se_nie_iv <- sqrt(var_nie_iv)
  ci_lower <- nie_iv - qnorm(0.975) * se_nie_iv
  ci_upper <- nie_iv + qnorm(0.975) * se_nie_iv
  
  return(list(NIE = nie_iv,
              NIE_reduced = reduced_form$NIE,
              first_stage = first_stage$ATE,
              SE = se_nie_iv,
              SE_reduced = reduced_form$SE,
              SE_first_stage = first_stage$SE,
              CI = c(lower = ci_lower,upper = ci_upper),
              CI_reduced = reduced_form$CI,
              CI_first_stage = first_stage$CI,
              IF = IF_nie,
              IF_reduced = IF_reduced,
              IF_first_stage = IF_first_stage))
}


### Generate simulated data with IV-mediation structure
generate_stroke_data_iv_mediation <- function(n = 1000,
                                              exposure_name = "insured",
                                              mediator_name = "rehabIRF",
                                              iv_name = "Z",
                                              iv_strength = 0.10,
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
  
  ### generate mediator (rehabIRF) using linear probability model
  mediator_intercept_base <- 0.35 - 0.4 * iv_strength
  
  ps_mediator <- mediator_intercept_base +
    iv_strength * IV +
    0.20 * insured +
    0.01 * age_scaled + 
    0.03 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.005 * insured * age_scaled +
    0.003 * insured * stroke_severity_scaled +
    -0.003 * insured * comorbidity_scaled
  
  # check for invalid probabilities
  if (any(ps_mediator < 0 | ps_mediator > 1)) {
    stop(paste("Invalid probabilities in mediator model. Range:",
               round(min(ps_mediator),3),"to",round(max(ps_mediator),3),
               "\nAdjust coefficients or iv_strength."))
  }
  
  rehabIRF <- rbinom(n,size = 1,prob = ps_mediator)
  
  ### calculate true average marginal effect of Z on M
  true_iv_ame_on_mediator <- iv_strength
  
  ### generate outcome (OUT3_ADL_IADL) using LOGIT link for mean
  logit_mu_adl <- 0.20 + 
    0.25 * post_discharge_disability_scaled +
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
  
  ### calculate true causal effects using LOGIT link
  # Y(0,0)
  logit_mu_00 <- 0.20 + 
    0.25 * post_discharge_disability_scaled +
    -0.25 * 0 + -0.50 * 0 + -0.15 * 0 * 0 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_00 <- plogis(logit_mu_00) * 3.0
  
  # Y(0,1)
  logit_mu_01 <- 0.20 + 
    0.25 * post_discharge_disability_scaled +
    -0.25 * 0 + -0.50 * 1 + -0.15 * 0 * 1 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_01 <- plogis(logit_mu_01) * 3.0
  
  # Y(1,0)
  logit_mu_10 <- 0.20 + 
    0.25 * post_discharge_disability_scaled +
    -0.25 * 1 + -0.50 * 0 + -0.15 * 1 * 0 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_10 <- plogis(logit_mu_10) * 3.0
  
  # Y(1,1)
  logit_mu_11 <- 0.20 + 
    0.25 * post_discharge_disability_scaled +
    -0.25 * 1 + -0.50 * 1 + -0.15 * 1 * 1 +
    0.12 * stroke_severity_scaled + 0.08 * age_scaled + 0.04 * comorbidity_scaled
  mu_11 <- plogis(logit_mu_11) * 3.0
  
  # Calculate mediator propensity scores under each exposure level INCLUDING IV (LPM)
  ps_m0_natural <- mediator_intercept_base +
    iv_strength * IV +
    0.01 * age_scaled + 
    0.03 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.20 * 0 +
    0.005 * 0 * age_scaled +
    0.003 * 0 * stroke_severity_scaled +
    -0.003 * 0 * comorbidity_scaled
  
  ps_m1_natural <- mediator_intercept_base +
    iv_strength * IV +
    0.01 * age_scaled + 
    0.03 * post_discharge_disability_scaled +
    0.01 * stroke_severity_scaled + 
    -0.005 * comorbidity_scaled +
    0.20 * 1 +
    0.005 * 1 * age_scaled +
    0.003 * 1 * stroke_severity_scaled +
    -0.003 * 1 * comorbidity_scaled
  
  # check for invalid probabilities
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
              simulation_params = list(n = n,
                                       exposure_name = exposure_name,
                                       mediator_name = mediator_name,
                                       iv_name = iv_name,
                                       phi_adl = phi_adl,
                                       iv_strength = iv_strength,
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
                                               iv_strength = 0.3,
                                               trunc_level = 0.001,
                                               fixed_covariates = NULL) {
  
  sim_data <- generate_stroke_data_iv_mediation(n = n,
                                                exposure_name = exposure_name,
                                                mediator_name = mediator_name,
                                                iv_name = iv_name,
                                                iv_strength = iv_strength,
                                                fixed_covariates = fixed_covariates)
  
  data <- sim_data$data
  
  # re-scale outcomes
  data$OUT3_ADL_IADL_01 <- data$OUT3_ADL_IADL / 3.0
  N_adl <- nrow(data)
  data$OUT3_ADL_IADL_01 <- (data$OUT3_ADL_IADL_01 * (N_adl - 1) + 0.5) / N_adl
  
  ### fit models
  # propensity score for A: P(A=1|W)
  ps_formula <- as.formula(paste(exposure_name,"~ age + stroke_severity + comorbidity"))
  ps_mod <- glm(ps_formula,data = data,family = binomial())
  
  # mediator model P(M|A,W) - marginal mediator density conditional on exposure
  mediator_formula <- as.formula(paste(mediator_name,
                                       "~ age + stroke_severity + comorbidity +",
                                       exposure_name,"+",
                                       exposure_name,"* age +",
                                       exposure_name,"* stroke_severity +",
                                       exposure_name,"* comorbidity"))
  mediator_mod <- glm(mediator_formula,data = data,family = binomial())
  
  # mediator model P(M|A,W,Z) - fit on entire dataset with exposure interactions
  mediator_formula_with_iv <- as.formula(paste(mediator_name,
                                               "~ age + stroke_severity + comorbidity +",
                                               iv_name,"+",
                                               exposure_name,"+",
                                               exposure_name,"* age +",
                                               exposure_name,"* stroke_severity +",
                                               exposure_name,"* comorbidity"))
  mediator_mod_with_iv <- glm(mediator_formula_with_iv,
                              data = data,
                              family = binomial())
  
  # IV model P(Z|W)
  iv_formula <- as.formula(paste(iv_name,"~ age + stroke_severity + comorbidity"))
  iv_mod <- glm(iv_formula,data = data,family = binomial())
  
  # outcome model E[Y|W,A,Z]
  outcome_formula <- as.formula(paste(
    "OUT3_ADL_IADL_01 ~ age + stroke_severity + comorbidity +",
    exposure_name,"+",iv_name,"+",
    exposure_name,"*",iv_name,"+",
    exposure_name,"* age +",
    exposure_name,"* stroke_severity +",
    exposure_name,"* comorbidity"
  ))
  outcome_mod <- betareg(outcome_formula,data = data)
  
  ### calculate NIE using IV approach
  nie_adl <- NA
  nie_se_adl <- NA
  nie_ci_lower_adl <- NA
  nie_ci_upper_adl <- NA
  nie_reduced <- NA
  nie_reduced_se <- NA
  first_stage_est <- NA
  first_stage_se <- NA
  nie_converged <- FALSE
  
  tryCatch({
    results_nie <- robustcompNIE_IV(outcome_mod = outcome_mod,
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
                                    scale_factor = 3.0)
    
    nie_adl <- results_nie$NIE
    nie_se_adl <- results_nie$SE
    nie_ci_lower_adl <- results_nie$CI["lower"]
    nie_ci_upper_adl <- results_nie$CI["upper"]
    nie_reduced <- results_nie$NIE_reduced
    nie_reduced_se <- results_nie$SE_reduced
    first_stage_est <- results_nie$first_stage
    first_stage_se <- results_nie$SE_first_stage
    nie_converged <- TRUE
  },error = function(e) {
    warning(paste("NIE estimation failed:",e$message))
  })
  
  # calculate bias and coverage
  if (nie_converged && !is.na(nie_adl) && is.finite(nie_adl)) {
    nie_bias <- nie_adl - sim_data$true_nie
    nie_covers <- (nie_ci_lower_adl <= sim_data$true_nie) & 
      (sim_data$true_nie <= nie_ci_upper_adl)
  } else {
    nie_bias <- NA
    nie_covers <- NA
  }
  
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
             nie_reduced_se = nie_reduced_se,
             first_stage_est = first_stage_est,
             first_stage_se = first_stage_se,
             true_nie_adl = sim_data$true_nie,
             true_nde_adl = sim_data$true_nde,
             true_te_adl = sim_data$true_te,
             true_iv_ame_on_mediator = sim_data$true_iv_ame_on_mediator,
             prop_insured = mean(data[[exposure_name]]),
             prop_rehabIRF = mean(data[[mediator_name]]),
             prop_IV = mean(data[[iv_name]]))
}


# full simulation study
run_simulation_study_iv_mediation <- function(n_sims = 100,
                                              n = 1000,
                                              exposure_name = "insured",
                                              mediator_name = "rehabIRF",
                                              iv_name = "Z",
                                              iv_strength = 0.3,
                                              trunc_level = 0.001,
                                              use_fixed_covariates = TRUE) {
  
  # generate fixed covariates once if requested
  fixed_covariates <- NULL
  if (use_fixed_covariates) {
    initial_data <- generate_stroke_data_iv_mediation(n = n,
                                                      exposure_name = exposure_name,
                                                      mediator_name = mediator_name,
                                                      iv_name = iv_name,
                                                      iv_strength = iv_strength)
    fixed_covariates <- initial_data$fixed_covariates
  }
  
  results_list <- vector("list",n_sims)
  
  pb <- txtProgressBar(min = 0,max = n_sims,style = 3)
  
  for (i in 1:n_sims) {
    results_list[[i]] <- run_single_simulation_iv_mediation(
      n = n,
      exposure_name = exposure_name,
      mediator_name = mediator_name,
      iv_name = iv_name,
      iv_strength = iv_strength,
      trunc_level = trunc_level,
      fixed_covariates = fixed_covariates
    )
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  results_df <- do.call(rbind,results_list)
  
  # summary statistics
  summary_stats <- data.frame(n = n,
                              n_sims = n_sims,
                              iv_strength = iv_strength,
                              n_converged_nie = sum(results_df$nie_converged,na.rm = TRUE),
                              nie_mean = mean(results_df$nie_adl,na.rm = TRUE),
                              nie_true = mean(results_df$true_nie_adl,na.rm = TRUE),
                              nde_true = mean(results_df$true_nde_adl,na.rm = TRUE),
                              te_true = mean(results_df$true_te_adl,na.rm = TRUE),
                              nie_bias = mean(results_df$nie_bias,na.rm = TRUE),
                              nie_empirical_se = sd(results_df$nie_adl,na.rm = TRUE),
                              nie_mean_se = mean(results_df$nie_se_adl,na.rm = TRUE),
                              nie_coverage = mean(results_df$nie_covers,na.rm = TRUE),
                              nie_rmse = sqrt(mean(results_df$nie_bias^2,na.rm = TRUE)),
                              mean_reduced_form = mean(results_df$nie_reduced,na.rm = TRUE),
                              empirical_se_reduced_form = sd(results_df$nie_reduced,na.rm = TRUE),
                              mean_se_reduced_form = mean(results_df$nie_reduced_se,na.rm = TRUE),
                              mean_first_stage = mean(results_df$first_stage_est,na.rm = TRUE),
                              empirical_se_first_stage = sd(results_df$first_stage_est,na.rm = TRUE),
                              mean_se_first_stage = mean(results_df$first_stage_se,na.rm = TRUE),
                              mean_true_iv_ame = mean(results_df$true_iv_ame_on_mediator,na.rm = TRUE),
                              mean_prop_exposure = mean(results_df$prop_insured,na.rm = TRUE),
                              mean_prop_mediator = mean(results_df$prop_rehabIRF,na.rm = TRUE),
                              mean_prop_iv = mean(results_df$prop_IV,na.rm = TRUE))
  
  return(list(results = results_df,summary = summary_stats))
}


# number of sims
N_SIMS <- 1500

### use weak instrument
# N=1000
sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 1000,
                                               iv_strength = 0.15,
                                               use_fixed_covariates = TRUE)


# display results
cat("\n=== IV-Mediation Performance Summary ===\n")

cat("\nData Characteristics:\n")
cat("  Mean P(Exposure=1):   ",round(sim_study$summary$mean_prop_exposure,3),"\n")
cat("  Mean P(Mediator=1):   ",round(sim_study$summary$mean_prop_mediator,3),"\n")
cat("  Mean P(Instrument=1): ",round(sim_study$summary$mean_prop_iv,3),"\n")

cat("\nTrue Causal Effects:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  True NDE:             ",round(sim_study$summary$nde_true,4),"\n")
cat("  True TE:              ",round(sim_study$summary$te_true,4),"\n")
cat("  Proportion Mediated:  ",round(sim_study$summary$nie_true / sim_study$summary$te_true,4),"\n")

cat("\nReduced Form (Numerator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_reduced_form,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_reduced_form,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_reduced_form,4),"\n")

cat("\nFirst Stage (Denominator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_first_stage,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_first_stage,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_first_stage,4),"\n")

cat("\nNatural Indirect Effect (NIE) - IV Approach:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  Mean estimate:        ",round(sim_study$summary$nie_mean,4),"\n")
cat("  Bias:                 ",round(sim_study$summary$nie_bias,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$nie_mean_se,4),"\n")
cat("  Coverage:             ",round(sim_study$summary$nie_coverage,3),"\n")
cat("  Convergence rate:     ",round(sim_study$summary$n_converged_nie / sim_study$summary$n_sims,3),"\n")




# N=8000
sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 8000,
                                               iv_strength = 0.05,
                                               use_fixed_covariates = FALSE)

print(sim_study$summary)

# display results
cat("\n=== IV-Mediation Performance Summary ===\n")

cat("\nData Characteristics:\n")
cat("  Mean P(Exposure=1):   ",round(sim_study$summary$mean_prop_exposure,3),"\n")
cat("  Mean P(Mediator=1):   ",round(sim_study$summary$mean_prop_mediator,3),"\n")
cat("  Mean P(Instrument=1): ",round(sim_study$summary$mean_prop_iv,3),"\n")

cat("\nTrue Causal Effects:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  True NDE:             ",round(sim_study$summary$nde_true,4),"\n")
cat("  True TE:              ",round(sim_study$summary$te_true,4),"\n")
cat("  Proportion Mediated:  ",round(sim_study$summary$nie_true / sim_study$summary$te_true,4),"\n")

cat("\nReduced Form (Numerator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_reduced_form,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_reduced_form,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_reduced_form,4),"\n")

cat("\nFirst Stage (Denominator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_first_stage,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_first_stage,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_first_stage,4),"\n")

cat("\nNatural Indirect Effect (NIE) - IV Approach:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  Mean estimate:        ",round(sim_study$summary$nie_mean,4),"\n")
cat("  Bias:                 ",round(sim_study$summary$nie_bias,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$nie_mean_se,4),"\n")
cat("  Coverage:             ",round(sim_study$summary$nie_coverage,3),"\n")
cat("  Convergence rate:     ",round(sim_study$summary$n_converged_nie / sim_study$summary$n_sims,3),"\n")







### use stronger instrument
# N=1000
sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 1000,
                                               iv_strength = 0.5,
                                               use_fixed_covariates = TRUE)


# display results
cat("\n=== IV-Mediation Performance Summary ===\n")

cat("\nData Characteristics:\n")
cat("  Mean P(Exposure=1):   ",round(sim_study$summary$mean_prop_exposure,3),"\n")
cat("  Mean P(Mediator=1):   ",round(sim_study$summary$mean_prop_mediator,3),"\n")
cat("  Mean P(Instrument=1): ",round(sim_study$summary$mean_prop_iv,3),"\n")

cat("\nTrue Causal Effects:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  True NDE:             ",round(sim_study$summary$nde_true,4),"\n")
cat("  True TE:              ",round(sim_study$summary$te_true,4),"\n")
cat("  Proportion Mediated:  ",round(sim_study$summary$nie_true / sim_study$summary$te_true,4),"\n")

cat("\nReduced Form (Numerator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_reduced_form,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_reduced_form,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_reduced_form,4),"\n")

cat("\nFirst Stage (Denominator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_first_stage,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_first_stage,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_first_stage,4),"\n")

cat("\nNatural Indirect Effect (NIE) - IV Approach:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  Mean estimate:        ",round(sim_study$summary$nie_mean,4),"\n")
cat("  Bias:                 ",round(sim_study$summary$nie_bias,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$nie_mean_se,4),"\n")
cat("  Coverage:             ",round(sim_study$summary$nie_coverage,3),"\n")
cat("  Convergence rate:     ",round(sim_study$summary$n_converged_nie / sim_study$summary$n_sims,3),"\n")




# N=8000
sim_study <- run_simulation_study_iv_mediation(n_sims = N_SIMS,
                                               n = 8000,
                                               iv_strength = 0.45,
                                               use_fixed_covariates = FALSE)

print(sim_study$summary)

# display results
cat("\n=== IV-Mediation Performance Summary ===\n")

cat("\nData Characteristics:\n")
cat("  Mean P(Exposure=1):   ",round(sim_study$summary$mean_prop_exposure,3),"\n")
cat("  Mean P(Mediator=1):   ",round(sim_study$summary$mean_prop_mediator,3),"\n")
cat("  Mean P(Instrument=1): ",round(sim_study$summary$mean_prop_iv,3),"\n")

cat("\nTrue Causal Effects:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  True NDE:             ",round(sim_study$summary$nde_true,4),"\n")
cat("  True TE:              ",round(sim_study$summary$te_true,4),"\n")
cat("  Proportion Mediated:  ",round(sim_study$summary$nie_true / sim_study$summary$te_true,4),"\n")

cat("\nReduced Form (Numerator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_reduced_form,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_reduced_form,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_reduced_form,4),"\n")

cat("\nFirst Stage (Denominator):\n")
cat("  Mean estimate:        ",round(sim_study$summary$mean_first_stage,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$empirical_se_first_stage,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$mean_se_first_stage,4),"\n")

cat("\nNatural Indirect Effect (NIE) - IV Approach:\n")
cat("  True NIE:             ",round(sim_study$summary$nie_true,4),"\n")
cat("  Mean estimate:        ",round(sim_study$summary$nie_mean,4),"\n")
cat("  Bias:                 ",round(sim_study$summary$nie_bias,4),"\n")
cat("  Empirical SE:         ",round(sim_study$summary$nie_empirical_se,4),"\n")
cat("  Mean estimated SE:    ",round(sim_study$summary$nie_mean_se,4),"\n")
cat("  Coverage:             ",round(sim_study$summary$nie_coverage,3),"\n")
cat("  Convergence rate:     ",round(sim_study$summary$n_converged_nie / sim_study$summary$n_sims,3),"\n")












