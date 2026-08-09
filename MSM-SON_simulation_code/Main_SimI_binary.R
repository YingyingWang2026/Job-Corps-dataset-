# Enconding: UTF-8
# Simulation I: Estimation and model-selection performance under M_1 and M_2
# main.R — for binary outcome Y_2

library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(numDeriv)
library(CompQuadForm)
library(MASS)

source("EM1_binary.R") 
source("EM2_binary.R") 


set.seed(2026)


N_SAMPLES   <- c(1000,2000,5000) # Sample size
N_BOOTSTRAP <- 200 # Number of simulations
SEED_BASE   <- 2026
N_CORES     <- max(1, detectCores() - 1)


SAVE_DIR <- "result_for_SimI_binary" 
if (!dir.exists(SAVE_DIR)) dir.create(SAVE_DIR, recursive = TRUE)


# Data-generating mechanism M_1 

generate_data_M_1 <- function(n = 1000) {
  X      <- rnorm(n)
  Y_1_true <- rbinom(n, 1, plogis(-0.5 + 1.0 * X))
  Y_2_true <- rbinom(n, 1, plogis(-1.0 + 0.8*X + 1.5*Y_1_true + 0.5*X*Y_1_true))
  R_1 <- rbinom(n, 1, plogis(-0.8 + 0.5*X + 0.3*Y_1_true))
  R_2 <- rbinom(n, 1, plogis(-1.0 + 1.2 * R_1 + 0.4 * Y_2_true)) 
  data.frame(
    X, Y_1_obs = ifelse(R_1 == 1, Y_1_true, NA), R_1,
    Y_2_obs = ifelse(R_2 == 1, Y_2_true, NA), R_2, Y_1_true, Y_2_true
  )
}


# Data-generating mechanism M_2 

generate_data_M_2 <- function(n = 1000) {
  X      <- rnorm(n)
  Y_1_true <- rbinom(n, 1, plogis(-0.5 + 1.0 * X))
  Y_2_true <- rbinom(n, 1, plogis(-1.0 + 0.8*X + 1.5*Y_1_true + 0.5*X*Y_1_true))
  R_1 <- rbinom(n, 1, plogis(-0.8 + 0.5*X + 0.3*Y_1_true))
  R_2 <- rbinom(n, 1, plogis(-1.0 + 2.5 * Y_1_true + 0.4 * Y_2_true))
  data.frame(
    X, Y_1_obs = ifelse(R_1 == 1, Y_1_true, NA), R_1,
    Y_2_obs = ifelse(R_2 == 1, Y_2_true, NA), R_2, Y_1_true, Y_2_true
  )
}

# Calculate Tau_SC.
compute_Tau_SC <- function(params, data) {
  phi <- params$phi; xv <- data$X
  EY1 <- mean(plogis(phi[1] + phi[2]*xv + phi[3]*1 + phi[4]*xv*1))
  EY0 <- mean(plogis(phi[1] + phi[2]*xv + phi[3]*0 + phi[4]*xv*0))
  list(E_Y2_1 = EY1, E_Y2_0 = EY0, Tau_SC = EY1 - EY0)
}

compute_true_Tau_SC <- function(data) {
  b0 <- -1.0; b1 <- 0.8; b2 <- 1.5; b3 <- 0.5; xv <- data$X
  EY1 <- mean(plogis(b0 + b1*xv + b2 + b3*xv))
  EY0 <- mean(plogis(b0 + b1*xv))
  list(E_Y2_1 = EY1, E_Y2_0 = EY0, Tau_SC = EY1 - EY0)
}

#estimate single
compute_sandwich <- function(theta_hat, data, loglik_fun, indloglik_fun) {
  H <- numDeriv::hessian(func = loglik_fun, x = theta_hat, data = data)
  score <- numDeriv::jacobian(func = indloglik_fun, x = theta_hat, data = data)
  B <- crossprod(score)
  Ainv <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  V <- Ainv %*% B %*% Ainv
  (V + t(V)) / 2
}

compute_Tau_SC_from_theta <- function(theta, data, model) {
  params <- if (model == "M_1") theta_to_params_M_1(theta) else theta_to_params_M_2(theta)
  compute_Tau_SC(params, data)$Tau_SC
}

run_single <- function(n, data_gen_func, em_func, true_param_func, model) {
  data <- data_gen_func(n)
  result <- em_func(data)
  tp <- true_param_func(data)
  ea <- compute_Tau_SC(result$params, data)
  ta <- compute_true_Tau_SC(data)
  
  if (model == "M_1") {
    theta_hat <- params_to_theta_M_1(result$params)
    V_theta <- compute_sandwich(theta_hat, data, loglik_theta_M_1, individual_loglik_theta_M_1)
  } else if (model == "M_2") {
    theta_hat <- params_to_theta_M_2(result$params)
    V_theta <- compute_sandwich(theta_hat, data, loglik_theta_M_2, individual_loglik_theta_M_2)
  } else {
    stop("model must be either 'M_1' or 'M_2'")
  }
  
  grad_Tau_SC <- numDeriv::grad(
    func = function(theta) compute_Tau_SC_from_theta(theta, data, model),
    x = theta_hat
  )
  Tau_SC_var <- as.numeric(t(grad_Tau_SC) %*% V_theta %*% grad_Tau_SC)
  Tau_SC_var <- max(Tau_SC_var, 0)
  Tau_SC_se <- sqrt(Tau_SC_var)
  
  theta_true <- c(tp$alpha, tp$phi, tp$gamma, tp$rho)
  
  list(
    converged = result$converged,
    Tau_SC = ea$Tau_SC, E_Y2_1 = ea$E_Y2_1, E_Y2_0 = ea$E_Y2_0,
    true_Tau_SC = ta$Tau_SC, Tau_SC_var = Tau_SC_var, Tau_SC_se = Tau_SC_se,
    mae_alpha = mean(abs(result$params$alpha - tp$alpha)),
    mae_phi = mean(abs(result$params$phi - tp$phi)),
    mae_gamma = mean(abs(result$params$gamma - tp$gamma)),
    mae_rho = mean(abs(result$params$rho - tp$rho)),
    loglik = tail(result$loglik_history, 1),
    theta_hat = theta_hat,
    theta_true = theta_true
  )
}
# estimate Bootstrap simulation 
run_experiment <- function(n_samples, n_boot, data_gen_func, em_func,
                           true_param_func, label, n_cores, save_dir, model) {
  
  all_res <- list()
  
  export_funcs <- c(
    "generate_data_M_1", "generate_data_M_2",
    "initialize_parameters_M_1", "E_step_M_1", "M_step_M_1",
    "compute_observed_loglik_M_1", "EM_algorithm_M_1", "compute_true_parameters_M_1",
    "compute_individual_loglik_M_1",
    "params_to_theta_M_1", "theta_to_params_M_1",
    "loglik_theta_M_1", "individual_loglik_theta_M_1",
    "initialize_parameters_M_2", "E_step_M_2", "M_step_M_2",
    "compute_observed_loglik_M_2", "EM_algorithm_M_2", "compute_true_parameters_M_2",
    "compute_individual_loglik_M_2",
    "params_to_theta_M_2", "theta_to_params_M_2",
    "loglik_theta_M_2", "individual_loglik_theta_M_2",
    "compute_Tau_SC", "compute_true_Tau_SC", "compute_sandwich", "compute_Tau_SC_from_theta", "run_single",
    "SEED_BASE"
  )
  
  
  theta_names <- c("alpha_0", "alpha_1",
                   "phi_0", "phi_1", "phi_2", "phi_3",
                   "gamma_0", "gamma_1", "gamma_2",
                   "rho_0", "rho_1", "rho_2")
  
  for (n in n_samples) {
    cat(sprintf("\n  sample n = %d ...\n", n))
    t0 <- Sys.time()
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    clusterExport(cl, export_funcs, envir = environment())
    
    
    raw <- foreach(b = 1:n_boot, .packages = c("tidyverse", "numDeriv", "MASS")) %dopar% {
      set.seed(SEED_BASE + b)
      run_single(
        n = n,
        data_gen_func = data_gen_func,
        em_func = em_func,
        true_param_func = true_param_func,
        model = model
      )
    }
    
    stopCluster(cl)
    
    # Tau_SC results.
    Tau_SC_df <- do.call(rbind, lapply(seq_along(raw), function(b) {
      r <- raw[[b]]
      data.frame(
        iter = b, sample_size = n, converged = r$converged,
        Tau_SC = r$Tau_SC, E_Y2_1 = r$E_Y2_1, E_Y2_0 = r$E_Y2_0,
        true_Tau_SC = r$true_Tau_SC,
        bias = r$Tau_SC - r$true_Tau_SC,
        Tau_SC_var = r$Tau_SC_var, Tau_SC_se = r$Tau_SC_se,
        mae_alpha = r$mae_alpha, mae_phi = r$mae_phi,
        mae_gamma = r$mae_gamma, mae_rho = r$mae_rho,
        total_loglik = r$loglik
      )
    }))
    
    # Parameter estimation results.
    param_df <- do.call(rbind, lapply(seq_along(raw), function(b) {
      r <- raw[[b]]
      df <- data.frame(iter = b, sample_size = n)
      est <- as.data.frame(t(r$theta_hat))
      colnames(est) <- theta_names
      true_vals <- as.data.frame(t(r$theta_true))
      colnames(true_vals) <- paste0(theta_names, "_true")
      cbind(df, est, true_vals)
    }))
    
    # Delta-method 95% confidence interval and coverage.
    Tau_SC_df$CI_lower <- Tau_SC_df$Tau_SC - 1.96 * Tau_SC_df$Tau_SC_se
    Tau_SC_df$CI_upper <- Tau_SC_df$Tau_SC + 1.96 * Tau_SC_df$Tau_SC_se
    Tau_SC_df$coverage <- as.numeric(
      Tau_SC_df$true_Tau_SC >= Tau_SC_df$CI_lower & Tau_SC_df$true_Tau_SC <= Tau_SC_df$CI_upper
    )
    
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    
    cov_rate <- mean(Tau_SC_df$coverage, na.rm = TRUE)
    cat(sprintf(
      "    Completed! %.1fs | Coverage %.1f%% | Tau_SC %.4f | Bias %.4f | MC_SD %.4f | Mean_SE %.4f | Mean_Var %.6f\n",
      dt, cov_rate * 100, mean(Tau_SC_df$Tau_SC), mean(Tau_SC_df$bias),
      sd(Tau_SC_df$Tau_SC), mean(Tau_SC_df$Tau_SC_se), mean(Tau_SC_df$Tau_SC_var)
    ))
    cat(sprintf("    MAE: alpha=%.4f phi=%.4f gamma=%.4f rho=%.4f\n",
                mean(Tau_SC_df$mae_alpha), mean(Tau_SC_df$mae_phi),
                mean(Tau_SC_df$mae_gamma), mean(Tau_SC_df$mae_rho)))
    
    all_res[[as.character(n)]] <- list(
      Tau_SC_df = Tau_SC_df,
      param_df = param_df,
      time     = dt
    )
  }
  
  
  summary_df <- data.frame(
    sample_size = n_samples,
    Coverage  = sapply(all_res, function(x) mean(x$Tau_SC_df$coverage, na.rm = TRUE)),
    Tau_SC_mean  = sapply(all_res, function(x) mean(x$Tau_SC_df$Tau_SC)),
    Tau_SC_bias  = sapply(all_res, function(x) mean(x$Tau_SC_df$bias)),
    Tau_SC_MSE   = sapply(all_res, function(x) mean(x$Tau_SC_df$bias^2)),
    MC_SD     = sapply(all_res, function(x) sd(x$Tau_SC_df$Tau_SC)),
    Mean_SE   = sapply(all_res, function(x) mean(x$Tau_SC_df$Tau_SC_se)),
    Mean_Var  = sapply(all_res, function(x) mean(x$Tau_SC_df$Tau_SC_var))
  )
  cat("\nSummary:\n"); print(round(summary_df, 4))
  write.csv(summary_df,
            file.path(save_dir, sprintf("%s_summary.csv", label)),
            row.names = FALSE)
  
  all_res
}



# Step 1: candidate model M_1 fitted to data generated under M_1.
cat("\n\n*** Step 1: EM_M_1 on data_M_1 (correct specification) ***\n")
res1 <- run_experiment(N_SAMPLES, N_BOOTSTRAP, generate_data_M_1, EM_algorithm_M_1,
                       compute_true_parameters_M_1, "M_1_on_dataM_1", N_CORES, SAVE_DIR,
                       model = "M_1")


# Step 2: candidate model M_2 fitted to data generated under M_2.
cat("\n\n*** Step 2: EM_M_2 on data_M_2 (correct specification) ***\n")
res2 <- run_experiment(N_SAMPLES, N_BOOTSTRAP, generate_data_M_2, EM_algorithm_M_2,
                       compute_true_parameters_M_2, "M_2_on_dataM_2", N_CORES, SAVE_DIR,
                       model = "M_2")


# Step 3: candidate model M_1 fitted to data generated under M_2.
cat("\n\n*** Step 3: EM_M_1 on data_M_2 (M_1 misspecified) ***\n")
res3 <- run_experiment(N_SAMPLES, N_BOOTSTRAP, generate_data_M_2, EM_algorithm_M_1,
                       compute_true_parameters_M_1, "M_1_on_dataM_2", N_CORES, SAVE_DIR,
                       model = "M_1")


# Step 4: candidate model M_2 fitted to data generated under M_1.
cat("\n\n*** Step 4: EM_M_2 on data_M_1 (M_2 misspecified) ***\n")
res4 <- run_experiment(N_SAMPLES, N_BOOTSTRAP, generate_data_M_1, EM_algorithm_M_2,
                       compute_true_parameters_M_2, "M_2_on_dataM_1", N_CORES, SAVE_DIR,
                       model = "M_2")

# Vuong Two-Step Test base function

compute_scores <- function(theta, data, indiv_loglik_func, eps = 1e-5) {
  n <- nrow(data); p <- length(theta)
  S <- matrix(0, n, p)
  ll0 <- indiv_loglik_func(theta, data)
  for (j in 1:p) {
    theta_p <- theta; theta_p[j] <- theta_p[j] + eps
    S[, j] <- (indiv_loglik_func(theta_p, data) - ll0) / eps
  }
  S
}

safe_inverse <- function(Y_1) {
  Y_1 <- (Y_1 + t(Y_1)) / 2
  tryCatch(solve(Y_1), error = function(e) ginv(Y_1))
}

imhof_pval_from_W <- function(t_val, W, tol = 1e-10) {
  eig <- eigen(W, symmetric = FALSE, only.values = TRUE)$values
  max_abs_imag <- max(abs(Im(eig)))
  lam <- Re(eig)
  lam <- lam[is.finite(lam)]
  lam <- lam[abs(lam) > tol]^2
  
  if (length(lam) == 0) {
    return(list(pval = ifelse(t_val <= tol, 1, 0),
                n_lambda = 0L,
                max_abs_imag = max_abs_imag,
                ifault = 0L))
  }
  
  res <- CompQuadForm::imhof(t_val, lambda = lam)
  list(
    pval = pmax(0, pmin(1, res$Qq)),
    n_lambda = length(lam),
    max_abs_imag = max_abs_imag,
    ifault = ifelse(is.null(res$ifault), 0L, as.integer(res$ifault))
  )
}


# Single Vuong-test replication.
run_vuong_single <- function(sim_idx, n, seed_base, alpha,
                             data_gen_func,
                             theta_M_1, theta_M_2,
                             indiv_ll_M_1_func, indiv_ll_M_2_func,
                             total_ll_M_1_func, total_ll_M_2_func) {
  tryCatch({
    set.seed(seed_base + sim_idx)
    data <- data_gen_func(n)
    
    ll_M_1 <- indiv_ll_M_1_func(theta_M_1, data)
    ll_M_2 <- indiv_ll_M_2_func(theta_M_2, data)
    
    d <- ll_M_1 - ll_M_2
    LR_n <- sum(d)
    
    omega2_n <- mean((d - mean(d))^2)
    omega_n  <- sqrt(omega2_n)
    Tn <- n * omega2_n
    
    A_M_1 <- hessian(function(th) total_ll_M_1_func(th, data), theta_M_1) / n
    A_M_2 <- hessian(function(th) total_ll_M_2_func(th, data), theta_M_2) / n
    A_M_1 <- (A_M_1 + t(A_M_1)) / 2
    A_M_2 <- (A_M_2 + t(A_M_2)) / 2
    
    S_M_1 <- compute_scores(theta_M_1, data, indiv_ll_M_1_func)
    S_M_2 <- compute_scores(theta_M_2, data, indiv_ll_M_2_func)
    
    B_M_1  <- crossprod(S_M_1) / n
    B_M_2  <- crossprod(S_M_2) / n
    B_M_12 <- crossprod(S_M_1, S_M_2) / n
    
    A_M_1_inv <- safe_inverse(A_M_1)
    A_M_2_inv <- safe_inverse(A_M_2)
    
    W <- rbind(
      cbind(-B_M_1 %*% A_M_1_inv,  -B_M_12 %*% A_M_2_inv),
      cbind(t(B_M_12) %*% A_M_1_inv,  B_M_2 %*% A_M_2_inv)
    )
    
    
    imhof_out <- imhof_pval_from_W(Tn, W)
    pval <- imhof_out$pval
    
    if (pval > alpha) {
      decision <- "indistinguishable"
      Vn <- NA
      var_reject <- FALSE
    } else {
      Vn <- LR_n / (sqrt(n) * omega_n)
      var_reject <- TRUE
      
      if      (Vn >  1.96) decision <- "M_1_wins"
      else if (Vn < -1.96) decision <- "M_2_wins"
      else                 decision <- "M_1_eq_M_2"
    }
    
    data.frame(
      sim = sim_idx,
      pval = pval,
      decision = decision,
      LR_n = LR_n,
      omega_n = omega_n,
      Tn = Tn,
      Vn= Vn,
      variance_test_reject = var_reject,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    data.frame(
      sim = sim_idx,
      pval = NA,
      decision = "ERROR",
      LR_n = NA,
      omega_n = NA,
      Tn = NA,
      Vn = NA,
      variance_test_reject = NA,
      stringsAsFactors = FALSE
    )
  })
}



# Vuong batch simulations

run_vuong_batch <- function(n_samples, res_M_1, res_M_2, data_gen_func,
                            indiv_ll_M_1_func, indiv_ll_M_2_func,
                            total_ll_M_1_func, total_ll_M_2_func,
                            label, n_cores, save_dir,
                            seed_base = SEED_BASE,
                            alpha = 0.05) {
  
  cat("\n", strrep("=", 70), "\n  ", label,
      "\n", strrep("=", 70), "\n", sep="")
  model_selection_summary <- list()
  
  export_funcs <- c(
    "generate_data_M_1", "generate_data_M_2",
    "compute_individual_loglik_M_1", "compute_individual_loglik_M_2",
    "compute_observed_loglik_M_1",   "compute_observed_loglik_M_2",
    "theta_to_params_M_1", "theta_to_params_M_2",
    "individual_loglik_theta_M_1", "individual_loglik_theta_M_2",
    "loglik_theta_M_1", "loglik_theta_M_2",
    "compute_scores",  "safe_inverse",
    "imhof_pval_from_W", "run_vuong_single",
    "SEED_BASE"
  )
  
  for (n in n_samples) {
    cat(sprintf("\n  --- n = %d ---\n", n))
    t0 <- Sys.time()
    
    
    boot_M_1 <- res_M_1[[as.character(n)]]$param_df
    boot_M_2 <- res_M_2[[as.character(n)]]$param_df
    n_sim  <- min(nrow(boot_M_1), nrow(boot_M_2))
    
    
    theta_names <- c("alpha_0", "alpha_1",
                     "phi_0", "phi_1", "phi_2", "phi_3",
                     "gamma_0", "gamma_1", "gamma_2",
                     "rho_0", "rho_1", "rho_2")
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    clusterExport(cl, export_funcs, envir = environment())
    clusterEvalQ(cl, { library(numDeriv); library(CompQuadForm); library(MASS) })
    
    vuong_results <- foreach(b = 1:n_sim, .combine = "rbind") %dopar% {
      th_M_1 <- as.numeric(boot_M_1[b, theta_names])
      th_M_2 <- as.numeric(boot_M_2[b, theta_names])
      
      run_vuong_single(
        sim_idx = b,
        n = n,
        seed_base = seed_base,
        alpha = alpha,
        data_gen_func = data_gen_func,
        theta_M_1 = th_M_1,
        theta_M_2 = th_M_2,
        indiv_ll_M_1_func = indiv_ll_M_1_func,
        indiv_ll_M_2_func = indiv_ll_M_2_func,
        total_ll_M_1_func = total_ll_M_1_func,
        total_ll_M_2_func = total_ll_M_2_func
      )
    }
    
    stopCluster(cl)
    
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    valid <- vuong_results[vuong_results$decision != "ERROR", ]
    n_valid <- nrow(valid)
    n_reject <- sum(valid$variance_test_reject, na.rm = TRUE)
    n_err <- sum(vuong_results$decision == "ERROR")
    
    cat(sprintf("    Completed! %.1fs, valid replications: %d/%d\n", elapsed, n_valid, n_sim))
    cat(sprintf("    First-stage variance test: rejected %d (%.1f%%), not rejected %d (%.1f%%)\n",
                n_reject, 100*n_reject/n_valid,
                n_valid - n_reject, 100*(n_valid - n_reject)/n_valid))
    cat("    Final decision distribution:\n")
    dt <- table(valid$decision)
    for (d in names(dt)) {
      cat(sprintf("      %s: %3d (%5.1f%%)\n", d, dt[d], 100*dt[d]/n_valid))
    }
    if (n_err > 0) cat(sprintf("      ERROR: %d\n", n_err))
    decision_count <- function(value) sum(valid$decision == value, na.rm = TRUE)
    decision_rate <- function(value) if (n_valid > 0) decision_count(value) / n_valid else NA_real_
    model_selection_summary[[as.character(n)]] <- data.frame(
      sample_size = n, n_sim = n_sim, n_valid = n_valid, n_error = n_err,
      variance_reject_n = n_reject,
      variance_reject_rate = if (n_valid > 0) n_reject / n_valid else NA_real_,
      variance_not_reject_n = n_valid - n_reject,
      variance_not_reject_rate = if (n_valid > 0) (n_valid - n_reject) / n_valid else NA_real_,
      M_1_wins_n = decision_count("M_1_wins"), M_1_wins_rate = decision_rate("M_1_wins"),
      M_2_wins_n = decision_count("M_2_wins"), M_2_wins_rate = decision_rate("M_2_wins"),
      M_1_eq_M_2_n = decision_count("M_1_eq_M_2"), M_1_eq_M_2_rate = decision_rate("M_1_eq_M_2"),
      indistinguishable_n = decision_count("indistinguishable"),
      indistinguishable_rate = decision_rate("indistinguishable")
    )
  }
  model_selection_summary <- do.call(rbind, model_selection_summary)
  rownames(model_selection_summary) <- NULL
  write.csv(model_selection_summary,
            file.path(save_dir, sprintf("%s_summary.csv", label)),
            row.names = FALSE)
  model_selection_summary
}


# Vuong test output: true DGP = M_1, with M_2 misspecified.

cat("  Vuong test — DGP = M_1 (M_2 misspecified)\n") 
run_vuong_batch(
  n_samples       = N_SAMPLES,
  res_M_1         = res1,   
  res_M_2         = res4,  
  data_gen_func   = generate_data_M_1,
  indiv_ll_M_1_func = individual_loglik_theta_M_1,
  indiv_ll_M_2_func = individual_loglik_theta_M_2,
  total_ll_M_1_func = loglik_theta_M_1,
  total_ll_M_2_func = loglik_theta_M_2,
  label    = "vuong_dataM_1_M_2_misspecified",
  n_cores  = N_CORES,
  save_dir = SAVE_DIR,
  alpha = 0.05
)


# Vuong test output: true DGP = M_2, with M_1 misspecified. 
cat("  Vuong test — DGP = M_2 (M_1 misspecified)\n") 
run_vuong_batch(
  n_samples       = N_SAMPLES,
  res_M_1         = res3,   
  res_M_2         = res2,   
  data_gen_func   = generate_data_M_2,
  indiv_ll_M_1_func = individual_loglik_theta_M_1,
  indiv_ll_M_2_func = individual_loglik_theta_M_2,
  total_ll_M_1_func = loglik_theta_M_1,
  total_ll_M_2_func = loglik_theta_M_2,
  label    = "vuong_dataM_2_M_1_misspecified",
  n_cores  = N_CORES,
  save_dir = SAVE_DIR,
  alpha = 0.05
)
