# Enconding: UTF-8
# Simulation III: Estimation and model-selection performance under M_{1,2}
# main.R — for continuous outcome Y_2

library(foreach)
library(doSNOW)
library(MASS)
library(numDeriv)
library(CompQuadForm)

source("EM1_continuous.R")
source("EM2_continuous.R")

set.seed(2026)

N_SIM_VALUES <- c(100, 200, 500) # Number of simulations
N_SIZE      <- 1000 # Sample size 
N_CORES <- max(1, parallel::detectCores() - 1)
SEED_BASE   <- 2026
ALPHA       <- 0.05
M_PFI       <- 100
SAVE_DIR <- "result_for_SimIII_continuous"
dir.create(SAVE_DIR, recursive = TRUE, showWarnings = FALSE)

#  Data-generating mechanism  M_{1,2}
generate_data_P <- function(n = 1000) {
  X <- rnorm(n)
  Y_1_true <- rbinom(n, 1, plogis(-0.5 + 1.0 * X))
  mu_Y <- -1.0 + 0.8 * X + 1.5 * Y_1_true + 0.5 * X * Y_1_true
  Y_2_true <- rnorm(n, mean = mu_Y, sd = 1.0)
  R_1 <- rbinom(n, 1, plogis(-0.8 + 0.5 * X + 0.3 * Y_1_true))
  R_2 <- rbinom(n, 1, plogis(-1.0 + 2.5 * Y_1_true + 1.2 * R_1 + 0.4 * Y_2_true))
  
  data.frame(
    X,
    Y_1_obs = ifelse(R_1 == 1, Y_1_true, NA),
    R_1,
    Y_2_obs = ifelse(R_2 == 1, Y_2_true, NA),
    R_2,
    Y_1_true,
    Y_2_true
  )
}
# base_fuction

compute_Tau_SC <- function(params, data) {
  phi <- params$phi; xv <- data$X
  EY1 <- mean(phi$beta[1] + phi$beta[2]*xv + phi$beta[3]*1 + phi$beta[4]*xv*1)
  EY0 <- mean(phi$beta[1] + phi$beta[2]*xv + phi$beta[3]*0 + phi$beta[4]*xv*0)
  list(E_Y2_1 = EY1, E_Y2_0 = EY0, Tau_SC = EY1 - EY0)
}

compute_true_Tau_SC <- function(data) {
  xv <- data$X
  EY1 <- mean(-1.0 + 0.8*xv + 1.5 + 0.5*xv)
  EY0 <- mean(-1.0 + 0.8*xv)
  list(E_Y2_1 = EY1, E_Y2_0 = EY0, Tau_SC = EY1 - EY0)
}

compute_sandwich <- function(theta_hat,
                             data,
                             loglik_fun,
                             indloglik_fun) {
  
  theta_hat <- as.numeric(theta_hat)
  hessian_matrix <- numDeriv::hessian(func = loglik_fun, x = theta_hat, data = data)
  score <- numDeriv::jacobian( func = indloglik_fun,x = theta_hat,data = data)
  hessian_matrix <- as.matrix(hessian_matrix)
  score <- as.matrix(score)
  
  if (ncol(score) != length(theta_hat)) {
    stop(
      sprintf(
        "Score matrix has %d columns, but theta_hat has length %d.",
        ncol(score),
        length(theta_hat)
      )
    )
  }
  B <- crossprod(score)
  Ainv <- tryCatch(
    solve(-hessian_matrix),
    error = function(e) MASS::ginv(-hessian_matrix)
  )
  V <- Ainv %*% B %*% t(Ainv)
  V <- (V + t(V)) / 2
  V
}

empty_single_result <- function(error_message = NA_character_) {
  list(
    success = FALSE,
    converged = FALSE,
    valid = FALSE,
    error_message = error_message,
    Tau_SC = NA_real_,
    E_Y2_1 = NA_real_,
    E_Y2_0 = NA_real_,
    true_Tau_SC = NA_real_,
    Tau_SC_var = NA_real_,
    Tau_SC_se = NA_real_,
    loglik = NA_real_,
    theta_hat = rep(NA_real_, 13),
    theta_true = rep(NA_real_, 13)
  )
}

run_single <- function(n, data_gen_func, em_func, true_param_func, M, model_type) {
  
  tryCatch({
    
    data <- data_gen_func(n)
    result <- em_func(data, M = M)
    tp <- true_param_func(data)
    ea <- compute_Tau_SC(result$params, data)
    ta <- compute_true_Tau_SC(data)
    
    theta_hat <- c(result$params$alpha, result$params$phi$beta, result$params$phi$sigma, result$params$gamma, result$params$rho)
    
    if (model_type == "M_1") {
      V_theta <- compute_sandwich(theta_hat, data, loglik_theta_M_1, individual_loglik_theta_M_1)
    } else if (model_type == "M_2") {
      V_theta <- compute_sandwich(theta_hat, data, loglik_theta_M_2, individual_loglik_theta_M_2)
    } else {
      stop("model_type must be either 'M_1' or 'M_2'.")
    }
    
    xbar <- mean(data$X)
    grad <- c(0, 0, 0, 0, 1, xbar, 0, 0, 0, 0, 0, 0, 0)
    Tau_SC_var <- as.numeric(t(grad) %*% V_theta %*% grad)
    
    if (is.finite(Tau_SC_var) && Tau_SC_var < 0 && Tau_SC_var > -1e-10) Tau_SC_var <- 0
    
    Tau_SC_se <- if (is.finite(Tau_SC_var) && Tau_SC_var >= 0) sqrt(Tau_SC_var) else NA_real_
    theta_true <- c(tp$alpha, tp$phi$beta, tp$phi$sigma, tp$gamma, tp$rho)
    final_loglik <- if (!is.null(result$loglik_history) && length(result$loglik_history) > 0) tail(result$loglik_history, 1) else NA_real_
    valid <- is.finite(ea$Tau_SC) && is.finite(ta$Tau_SC) && is.finite(Tau_SC_var) && Tau_SC_var >= 0 && is.finite(Tau_SC_se)
    
    list(
      success = TRUE,
      converged = isTRUE(result$converged),
      valid = valid,
      error_message = NA_character_,
      Tau_SC = ea$Tau_SC,
      E_Y2_1 = ea$E_Y2_1,
      E_Y2_0 = ea$E_Y2_0,
      true_Tau_SC = ta$Tau_SC,
      Tau_SC_var = Tau_SC_var,
      Tau_SC_se = Tau_SC_se,
      loglik = final_loglik,
      theta_hat = theta_hat,
      theta_true = theta_true
    )
    
  }, error = function(e) {
    
    empty_single_result(conditionMessage(e))
    
  })
}
run_experiment <- function(n_samples, n_sim, data_gen_func, em_func, true_param_func, label, model_type, n_cores, M, save_dir) {
  all_res <- list()
  
  export_funcs <- c(
    "generate_data_P", "compute_sandwich",
    "initialize_parameters_M_1", "E_step_M_1", "M_step_M_1",
    "compute_observed_loglik_M_1", "EM_algorithm_M_1", "compute_true_parameters_M_1",
    "compute_individual_loglik_M_1", "params_to_theta_M_1", "theta_to_params_M_1",
    "loglik_theta_M_1", "individual_loglik_theta_M_1",
    "initialize_parameters_M_2", "E_step_M_2", "M_step_M_2",
    "compute_observed_loglik_M_2", "EM_algorithm_M_2", "compute_true_parameters_M_2",
    "compute_individual_loglik_M_2", "params_to_theta_M_2", "theta_to_params_M_2",
    "loglik_theta_M_2", "individual_loglik_theta_M_2",
    "compute_Tau_SC", "compute_true_Tau_SC", "empty_single_result", "run_single", "SEED_BASE"
  )
  
  theta_names <- c(
    "alpha_0", "alpha_1",
    "beta_0", "beta_1", "beta_2", "beta_3",
    "sigma",
    "gamma_0", "gamma_1", "gamma_2",
    "rho_0", "rho_1", "rho_2"
  )
  
  for (n in n_samples) {
    
    cat(sprintf("\nSample size n = %d, number of simulations = %d, M = %d\n", n, n_sim, M))
    
    t0 <- Sys.time()
    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    
    parallel::clusterExport(cl, export_funcs, envir = .GlobalEnv)
    
    parallel::clusterEvalQ(cl, {
      library(MASS)
      library(numDeriv)
      NULL
    })
    
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
    progress_fun <- function(k) utils::setTxtProgressBar(pb, k)
    snow_options <- list(progress = progress_fun)
    
    raw <- tryCatch({
      
      foreach(
        b = seq_len(n_sim),
        .packages = c("MASS", "numDeriv"),
        .options.snow = snow_options,
        .errorhandling = "pass"
      ) %dopar% {
        
        set.seed(SEED_BASE + b)
        run_single(n, data_gen_func, em_func, true_param_func, M, model_type)
        
      }
      
    }, finally = {
      
      close(pb)
      parallel::stopCluster(cl)
      
    })
    
    raw <- lapply(raw, function(r) {
      
      if (inherits(r, "error")) return(empty_single_result(conditionMessage(r)))
      
      r
      
    })
    
    ace_df <- do.call(rbind, lapply(seq_along(raw), function(b) {
      
      r <- raw[[b]]
      
      data.frame(
        iter = b,
        sample_size = n,
        success = r$success,
        converged = r$converged,
        valid = r$valid,
        error_message = r$error_message,
        Tau_SC = r$Tau_SC,
        E_Y2_1 = r$E_Y2_1,
        E_Y2_0 = r$E_Y2_0,
        true_Tau_SC = r$true_Tau_SC,
        bias = r$Tau_SC - r$true_Tau_SC,
        Tau_SC_var = r$Tau_SC_var,
        Tau_SC_se = r$Tau_SC_se,
        total_loglik = r$loglik
      )
      
    }))
    
    param_df <- do.call(rbind, lapply(seq_along(raw), function(b) {
      
      r <- raw[[b]]
      df <- data.frame(iter = b, sample_size = n, success = r$success, converged = r$converged, valid = r$valid, total_loglik = r$loglik)
      est <- as.data.frame(t(r$theta_hat))
      true_vals <- as.data.frame(t(r$theta_true))
      
      colnames(est) <- theta_names
      colnames(true_vals) <- paste0(theta_names, "_true")
      
      cbind(df, est, true_vals)
      
    }))
    
    ace_df$CI_lower <- ace_df$Tau_SC - 1.96 * ace_df$Tau_SC_se
    ace_df$CI_upper <- ace_df$Tau_SC + 1.96 * ace_df$Tau_SC_se
    
    ace_df$coverage <- ifelse(
      is.finite(ace_df$true_Tau_SC) & is.finite(ace_df$CI_lower) & is.finite(ace_df$CI_upper),
      as.numeric(ace_df$true_Tau_SC >= ace_df$CI_lower & ace_df$true_Tau_SC <= ace_df$CI_upper),
      NA_real_
    )
    
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    
    cat(sprintf("\nCompleted in %.1f seconds | Successful %d/%d | Converged %d/%d | Valid %d/%d\n",
                dt,
                sum(ace_df$success, na.rm = TRUE), n_sim,
                sum(ace_df$converged, na.rm = TRUE), n_sim,
                sum(ace_df$valid, na.rm = TRUE), n_sim))
    
    cat(sprintf("Coverage %.1f%% | Mean Tau_SC %.4f | Bias %.4f | MC_SD %.4f | Mean_SE %.4f | Mean_Var %.6f\n",
                mean(ace_df$coverage, na.rm = TRUE) * 100,
                mean(ace_df$Tau_SC, na.rm = TRUE),
                mean(ace_df$bias, na.rm = TRUE),
                sd(ace_df$Tau_SC, na.rm = TRUE),
                mean(ace_df$Tau_SC_se, na.rm = TRUE),
                mean(ace_df$Tau_SC_var, na.rm = TRUE)))
    
    all_res[[as.character(n)]] <- list(
      ace_df = ace_df,
      param_df = param_df,
      raw = raw,
      time = dt
    )
  }
  
  summary_df <- do.call(rbind, lapply(names(all_res), function(n_name) {
    
    df <- all_res[[n_name]]$ace_df
    
    data.frame(
      sample_size = as.numeric(n_name),
      Coverage = mean(df$coverage, na.rm = TRUE),
      Tau_SC_mean = mean(df$Tau_SC, na.rm = TRUE),
      Tau_SC_bias = mean(df$bias, na.rm = TRUE),
      Tau_SC_MSE = mean(df$bias^2, na.rm = TRUE),
      MC_SD = sd(df$Tau_SC, na.rm = TRUE),
      Mean_SE = mean(df$Tau_SC_se, na.rm = TRUE),
      Mean_Var = mean(df$Tau_SC_var, na.rm = TRUE)
    )
    
  }))
  
  rownames(summary_df) <- NULL
  
  cat("\nSummary results:\n")
  print(round(summary_df, 4))

  list(all_res = all_res, summary = summary_df)
}

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

safe_inverse <- function(matrix_to_invert) {
  matrix_to_invert <- (matrix_to_invert + t(matrix_to_invert)) / 2
  tryCatch(
    solve(matrix_to_invert),
    error = function(e) MASS::ginv(matrix_to_invert)
  )
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

run_vuong_type12_from_results <- function(sim_idx, n, result_M_1, result_M_2, data_gen_func, indiv_ll_M_1_func, indiv_ll_M_2_func, total_ll_M_1_func, total_ll_M_2_func) {
  
  tryCatch({
    
    if (!isTRUE(result_M_1$success) || !isTRUE(result_M_2$success)) {
      return(data.frame(sim = sim_idx, success_M_1 = isTRUE(result_M_1$success), success_M_2 = isTRUE(result_M_2$success), converged_M_1 = isTRUE(result_M_1$converged), converged_M_2 = isTRUE(result_M_2$converged), pval = NA_real_, decision = "ERROR_IN_EM", LR_n = NA_real_, omega2_n = NA_real_, omega_n = NA_real_, Tn = NA_real_, Vn = NA_real_, variance_test_reject = NA, n_lambda = NA_integer_, max_abs_imag = NA_real_, ifault = NA_integer_, error_message = paste(result_M_1$error_message, result_M_2$error_message, sep = " | "), stringsAsFactors = FALSE))
    }
    
    if (!isTRUE(result_M_1$converged) || !isTRUE(result_M_2$converged)) {
      return(data.frame(sim = sim_idx, success_M_1 = TRUE, success_M_2 = TRUE, converged_M_1 = isTRUE(result_M_1$converged), converged_M_2 = isTRUE(result_M_2$converged), pval = NA_real_, decision = "NOT_CONVERGED", LR_n = NA_real_, omega2_n = NA_real_, omega_n = NA_real_, Tn = NA_real_, Vn = NA_real_, variance_test_reject = NA, n_lambda = NA_integer_, max_abs_imag = NA_real_, ifault = NA_integer_, error_message = NA_character_, stringsAsFactors = FALSE))
    }
    
    set.seed(SEED_BASE + sim_idx)
    data <- data_gen_func(n)
    
    theta_M_1 <- as.numeric(result_M_1$theta_hat)
    theta_M_2 <- as.numeric(result_M_2$theta_hat)
    ll_M_1 <- as.numeric(indiv_ll_M_1_func(theta_M_1, data))
    ll_M_2 <- as.numeric(indiv_ll_M_2_func(theta_M_2, data))
    
    d <- ll_M_1 - ll_M_2
    n_obs <- length(d)
    LR_n <- sum(d)
    omega2_n <- mean((d - mean(d))^2)
    omega_n <- sqrt(omega2_n)
    Tn <- n_obs * omega2_n
    
    A_M_1 <- numDeriv::hessian(function(th) total_ll_M_1_func(th, data), theta_M_1) / n_obs
    A_M_2 <- numDeriv::hessian(function(th) total_ll_M_2_func(th, data), theta_M_2) / n_obs
    
    A_M_1 <- (A_M_1 + t(A_M_1)) / 2
    A_M_2 <- (A_M_2 + t(A_M_2)) / 2
    
    S_M_1 <- compute_scores(theta_M_1, data, indiv_ll_M_1_func)
    S_M_2 <- compute_scores(theta_M_2, data, indiv_ll_M_2_func)
    
    B_M_1 <- crossprod(S_M_1) / n_obs
    B_M_2 <- crossprod(S_M_2) / n_obs
    B_M_1G <- crossprod(S_M_1, S_M_2) / n_obs
    
    A_M_1_inv <- safe_inverse(A_M_1)
    A_M_2_inv <- safe_inverse(A_M_2)
    
    W <- rbind(
      cbind(-B_M_1 %*% A_M_1_inv, -B_M_1G %*% A_M_2_inv),
      cbind(t(B_M_1G) %*% A_M_1_inv, B_M_2 %*% A_M_2_inv)
    )
    
    imhof_out <- imhof_pval_from_W(Tn, W)
    pval <- imhof_out$pval
  
    if (pval > ALPHA) {
      decision <- "Cannot_distinguish"
      Vn <- NA_real_
      variance_test_reject <- FALSE
    } else {
      variance_test_reject <- TRUE
      
      if (!is.finite(omega_n) || omega_n <= 1e-12) {
        Vn <- NA_real_
        decision <- "ZERO_VARIANCE"
      } else {
        Vn <- LR_n / (sqrt(n_obs) * omega_n)
        zcrit <- qnorm(1 - ALPHA / 2)
        
        if (Vn > zcrit) {
          decision <- "M_1_wins"
        } else if (Vn < -zcrit) {
          decision <- "M_2_wins"
        } else {
          decision <- "H_eq_M_2"
        }
      }
    }
    
    data.frame(
      sim = sim_idx,
      success_M_1 = TRUE,
      success_M_2 = TRUE,
      converged_M_1 = TRUE,
      converged_M_2 = TRUE,
      pval = pval,
      decision = decision,
      LR_n = LR_n,
      omega2_n = omega2_n,
      omega_n = omega_n,
      Tn = Tn,
      Vn = Vn,
      variance_test_reject = variance_test_reject,
      n_lambda = imhof_out$n_lambda,
      max_abs_imag = imhof_out$max_abs_imag,
      ifault = imhof_out$ifault,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    
    data.frame(
      sim = sim_idx,
      success_M_1 = isTRUE(result_M_1$success),
      success_M_2 = isTRUE(result_M_2$success),
      converged_M_1 = isTRUE(result_M_1$converged),
      converged_M_2 = isTRUE(result_M_2$converged),
      pval = NA_real_,
      decision = "ERROR",
      LR_n = NA_real_,
      omega2_n = NA_real_,
      omega_n = NA_real_,
      Tn = NA_real_,
      Vn = NA_real_,
      variance_test_reject = NA,
      n_lambda = NA_integer_,
      max_abs_imag = NA_real_,
      ifault = NA_integer_,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    
  })
}
# Run the simulation sequentially with N_SIM = 100, 200, and 500, and combine the results.
M_1_summary_rows <- list()
M_2_summary_rows <- list()
model_selection_rows <- list()
simulation_results <- list()

for (N_SIM in N_SIM_VALUES) {
  cat(sprintf("\nConfiguration: N_SIM=%d, n=%d, alpha=%.2f, cores=%d\n", N_SIM, N_SIZE, ALPHA, N_CORES))

  res_M_1 <- run_experiment(
    N_SIZE, N_SIM, generate_data_P, EM_algorithm_M_1,
    compute_true_parameters_M_1, "M_1_on_dataP", "M_1", N_CORES, M_PFI, SAVE_DIR
  )
  res_M_2 <- run_experiment(
    N_SIZE, N_SIM, generate_data_P, EM_algorithm_M_2,
    compute_true_parameters_M_2, "M_2_on_dataP", "M_2", N_CORES, M_PFI, SAVE_DIR
  )

  M_1_summary_rows[[as.character(N_SIM)]] <- cbind(n_sim = N_SIM, res_M_1$summary)
  M_2_summary_rows[[as.character(N_SIM)]] <- cbind(n_sim = N_SIM, res_M_2$summary)
  sample_key <- as.character(N_SIZE)
  M_1_raw <- res_M_1$all_res[[sample_key]]$raw
  M_2_raw <- res_M_2$all_res[[sample_key]]$raw
  N_PAIR <- min(N_SIM, length(M_1_raw), length(M_2_raw))
  t_start_vuong <- Sys.time()

cl <- parallel::makeCluster(N_CORES)
doSNOW::registerDoSNOW(cl)

parallel::clusterEvalQ(cl, {
  library(numDeriv)
  library(CompQuadForm)
  library(MASS)
  NULL
})

export_vars_vuong <- c(
  "ALPHA",
  "SEED_BASE",
  "generate_data_P",
  "compute_scores",
  "safe_inverse",
  "imhof_pval_from_W",
  "loglik_theta_M_1",
  "loglik_theta_M_2",
  "individual_loglik_theta_M_1",
  "individual_loglik_theta_M_2",
  "compute_observed_loglik_M_1",
  "compute_observed_loglik_M_2",
  "compute_individual_loglik_M_1",
  "compute_individual_loglik_M_2",
  "theta_to_params_M_1",
  "theta_to_params_M_2",
  "run_vuong_type12_from_results"
)

missing_exports <- export_vars_vuong[!vapply(export_vars_vuong, exists, logical(1), envir = .GlobalEnv)]

if (length(missing_exports) > 0) {
  parallel::stopCluster(cl)
  stop("Missing objects: ", paste(missing_exports, collapse = ", "))
}

parallel::clusterExport(cl, export_vars_vuong, envir = .GlobalEnv)

pb_vuong <- utils::txtProgressBar(min = 0, max = N_PAIR, style = 3)
progress_vuong <- function(k) utils::setTxtProgressBar(pb_vuong, k)
opts_vuong <- list(progress = progress_vuong)

results_raw <- tryCatch({
  
  foreach(
    b = seq_len(N_PAIR),
    .packages = c("numDeriv", "CompQuadForm", "MASS"),
    .options.snow = opts_vuong,
    .errorhandling = "pass"
  ) %dopar% {
    
    run_vuong_type12_from_results(
      sim_idx = b,
      n = N_SIZE,
      result_M_1 = M_1_raw[[b]],
      result_M_2 = M_2_raw[[b]],
      data_gen_func = generate_data_P,
      indiv_ll_M_1_func = individual_loglik_theta_M_1,
      indiv_ll_M_2_func = individual_loglik_theta_M_2,
      total_ll_M_1_func = loglik_theta_M_1,
      total_ll_M_2_func = loglik_theta_M_2
    )
    
  }
  
}, finally = {
  
  close(pb_vuong)
  parallel::stopCluster(cl)
  
})

cat("\n")

t_end_vuong <- Sys.time()

results_raw <- lapply(seq_along(results_raw), function(b) {
  r <- results_raw[[b]]
  
  if (inherits(r, "error")) {
    return(data.frame(
      sim = b,
      success_M_1 = NA,
      success_M_2 = NA,
      converged_M_1 = NA,
      converged_M_2 = NA,
      pval = NA_real_,
      decision = "PARALLEL_ERROR",
      LR_n = NA_real_,
      omega2_n = NA_real_,
      omega_n = NA_real_,
      Tn = NA_real_,
      Vn = NA_real_,
      variance_test_reject = NA,
      n_lambda = NA_integer_,
      max_abs_imag = NA_real_,
      ifault = NA_integer_,
      error_message = conditionMessage(r),
      stringsAsFactors = FALSE
    ))
  }
  
  r
})

results_df <- do.call(rbind, results_raw)
rownames(results_df) <- NULL

stage1_flag <- results_df$success_M_1 %in% TRUE & results_df$success_M_2 %in% TRUE &
  results_df$converged_M_1 %in% TRUE & results_df$converged_M_2 %in% TRUE &
  is.finite(results_df$pval) & !is.na(results_df$variance_test_reject)
stage1 <- results_df[stage1_flag, , drop = FALSE]
stage2 <- stage1[stage1$variance_test_reject %in% TRUE & is.finite(stage1$Vn), , drop = FALSE]
safe_rate <- function(count, denominator) if (denominator > 0) count / denominator else NA_real_
variance_reject_n <- sum(stage1$variance_test_reject %in% TRUE)
variance_not_reject_n <- sum(stage1$variance_test_reject %in% FALSE)
M_1_wins_n <- sum(stage2$decision == "M_1_wins", na.rm = TRUE)
M_2_wins_n <- sum(stage2$decision == "M_2_wins", na.rm = TRUE)
M_1_eq_M_2_n <- sum(stage2$decision == "M_1_eq_M_2", na.rm = TRUE)

model_selection_rows[[as.character(N_SIM)]] <- data.frame(
  n_sim = N_SIM, sample_size = N_SIZE, n_total = nrow(results_df),
  n_excluded = nrow(results_df) - nrow(stage1), stage1_valid_n = nrow(stage1),
  variance_reject_n = variance_reject_n,
  variance_reject_rate = safe_rate(variance_reject_n, nrow(stage1)),
  variance_not_reject_n = variance_not_reject_n,
  variance_not_reject_rate = safe_rate(variance_not_reject_n, nrow(stage1)),
  stage2_valid_n = nrow(stage2),
  M_1_wins_n = M_1_wins_n, M_1_wins_rate = safe_rate(M_1_wins_n, nrow(stage2)),
  M_2_wins_n = M_2_wins_n, M_2_wins_rate = safe_rate(M_2_wins_n, nrow(stage2)),
  M_1_eq_M_2_n = M_1_eq_M_2_n,
  M_1_eq_M_2_rate = safe_rate(M_1_eq_M_2_n, nrow(stage2))
)
simulation_results[[as.character(N_SIM)]] <- list(
  M_1 = res_M_1, M_2 = res_M_2, model_selection = results_df,
  elapsed_vuong_minutes = as.numeric(difftime(t_end_vuong, t_start_vuong, units = "mins"))
)
}

M_1_summary <- do.call(rbind, M_1_summary_rows)
M_2_summary <- do.call(rbind, M_2_summary_rows)
model_selection_summary <- do.call(rbind, model_selection_rows)
rownames(M_1_summary) <- rownames(M_2_summary) <- rownames(model_selection_summary) <- NULL

write.csv(M_1_summary, file.path(SAVE_DIR, "M_1_on_dataP_summary.csv"), row.names = FALSE)
write.csv(M_2_summary, file.path(SAVE_DIR, "M_2_on_dataP_summary.csv"), row.names = FALSE)
write.csv(model_selection_summary, file.path(SAVE_DIR, "model_selection_dataP_summary.csv"), row.names = FALSE)
