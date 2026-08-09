# Enconding: UTF-8
# Simulation II: Stage 1 model-selection performance under Data-generating mechanism M_{\phi}
# main.R — for continuous outcome Y_2


library(parallel)
library(foreach)
library(doParallel)
library(numDeriv)
library(CompQuadForm)
library(MASS)
library(dplyr)
library(ggplot2)
library(scales)


source("EM1_continuous.R")
source("EM2_continuous.R")

set.seed(2026)

N_SIM       <- 1000      # Number of simulations
N_SIZE      <- 1000      # Sample size
N_CORES     <- max(1, detectCores() - 1)
SEED_BASE   <- 2026
M_PFI       <- 100
SAVE_DIR    <- "result_for_SimII"
dir.create(SAVE_DIR, recursive = TRUE, showWarnings = FALSE)

## Data-generating mechanism M_{\phi}
generate_data_N <- function(n = 1000) {
  X      <- rnorm(n)
  Y_1_true <- rbinom(n, 1, plogis(-0.5 + 1.0 * X))
  mu_Y   <- -1.0 + 0.8 * X + 1.5 * Y_1_true + 0.5 * X * Y_1_true
  Y_2_true <- rnorm(n, mean = mu_Y, sd = 1.0)
  R_1    <- rbinom(n, 1, plogis(0.8 + 0.3 * X + 1.2 * Y_1_true))
  R_2    <- rbinom(n, 1, plogis(0.8 + 0.3 * Y_2_true))
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

compute_scores_local <- function(theta, data, indiv_loglik_func, eps = 1e-5) {
  n_obs <- nrow(data)
  p     <- length(theta)
  S     <- matrix(0, n_obs, p)
  ll0   <- indiv_loglik_func(theta, data)
  
  for (j in seq_len(p)) {
    th_p <- theta
    th_p[j] <- th_p[j] + eps
    S[, j] <- (indiv_loglik_func(th_p, data) - ll0) / eps
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

run_vuong_type1_single <- function(sim_idx, n, seed_base, M,
                                   data_gen_func,
                                   em_M_1_func, em_M_2_func,
                                   indiv_ll_M_1_func, indiv_ll_M_2_func,
                                   total_ll_M_1_func, total_ll_M_2_func) {
  tryCatch({
    set.seed(seed_base + sim_idx)
    data <- data_gen_func(n)
    
    fit_M_1 <- em_M_1_func(data, M = M)
    fit_M_2 <- em_M_2_func(data, M = M)
    
    if (!fit_M_1$converged || !fit_M_2$converged) {
      return(data.frame(
        sim = sim_idx,
        converged = FALSE,
        Tn = NA_real_,
        omega2 = NA_real_,
        pval = NA_real_,
        reject_05 = NA_integer_,
        n_lambda = NA_integer_,
        max_abs_imag = NA_real_,
        ifault = NA_integer_
      ))
    }
    
    theta_M_1 <- params_to_theta_M_1(fit_M_1$params)
    theta_M_2 <- params_to_theta_M_2(fit_M_2$params)
    
    ll_M_1 <- indiv_ll_M_1_func(theta_M_1, data)
    ll_M_2 <- indiv_ll_M_2_func(theta_M_2, data)
    d    <- ll_M_1 - ll_M_2
    
    omega2_n <- mean((d - mean(d))^2)
    Tn       <- n * omega2_n
    
    A_M_1 <- hessian(function(th) total_ll_M_1_func(th, data), theta_M_1) / n
    A_M_2 <- hessian(function(th) total_ll_M_2_func(th, data), theta_M_2) / n
    A_M_1 <- (A_M_1 + t(A_M_1)) / 2
    A_M_2 <- (A_M_2 + t(A_M_2)) / 2
    
    S_M_1 <- compute_scores_local(theta_M_1, data, indiv_ll_M_1_func)
    S_M_2 <- compute_scores_local(theta_M_2, data, indiv_ll_M_2_func)
    
    B_M_1  <- crossprod(S_M_1) / n
    B_M_2  <- crossprod(S_M_2) / n
    B_M_12 <- crossprod(S_M_1, S_M_2) / n
    
    A_M_1_inv <- safe_inverse(A_M_1)
    A_M_2_inv <- safe_inverse(A_M_2)
    
    W <- rbind(
      cbind(-B_M_1 %*% A_M_1_inv, -B_M_12 %*% A_M_2_inv),
      cbind(t(B_M_12) %*% A_M_1_inv, B_M_2 %*% A_M_2_inv)
    )
    
    imhof_out <- imhof_pval_from_W(Tn, W)
    
    data.frame(
      sim = sim_idx,
      converged = TRUE,
      Tn = Tn,
      omega2 = omega2_n,
      pval = imhof_out$pval,
      reject_05 = as.integer(imhof_out$pval < 0.05),
      n_lambda = imhof_out$n_lambda,
      max_abs_imag = imhof_out$max_abs_imag,
      ifault = imhof_out$ifault
    )
  }, error = function(e) {
    data.frame(
      sim = sim_idx,
      converged = FALSE,
      Tn = NA_real_,
      omega2 = NA_real_,
      pval = NA_real_,
      reject_05 = NA_integer_,
      n_lambda = NA_integer_,
      max_abs_imag = NA_real_,
      ifault = NA_integer_
    )
  })
}


export_vars <- c(
  "generate_data_N",
  "compute_scores_local", "safe_inverse", "imhof_pval_from_W",
  "initialize_parameters_M_1", "E_step_M_1", "M_step_M_1",
  "compute_observed_loglik_M_1", "EM_algorithm_M_1",
  "compute_individual_loglik_M_1",
  "params_to_theta_M_1", "theta_to_params_M_1",
  "loglik_theta_M_1", "individual_loglik_theta_M_1",
  "initialize_parameters_M_2", "E_step_M_2", "M_step_M_2",
  "compute_observed_loglik_M_2", "EM_algorithm_M_2",
  "compute_individual_loglik_M_2",
  "params_to_theta_M_2", "theta_to_params_M_2",
  "loglik_theta_M_2", "individual_loglik_theta_M_2",
  "run_vuong_type1_single"
)

cl <- makeCluster(N_CORES)
registerDoParallel(cl)
clusterExport(cl, export_vars, envir = environment())
clusterEvalQ(cl, {
  library(numDeriv)
  library(CompQuadForm)
  library(MASS)
})

results_raw <- foreach(
  b = seq_len(N_SIM),
  .packages = c("numDeriv", "CompQuadForm", "MASS")
) %dopar% {
  run_vuong_type1_single(
    sim_idx = b,
    n = N_SIZE,
    seed_base = SEED_BASE,
    M = M_PFI,
    data_gen_func = generate_data_N,
    em_M_1_func = EM_algorithm_M_1,
    em_M_2_func = EM_algorithm_M_2,
    indiv_ll_M_1_func = individual_loglik_theta_M_1,
    indiv_ll_M_2_func = individual_loglik_theta_M_2,
    total_ll_M_1_func = loglik_theta_M_1,
    total_ll_M_2_func = loglik_theta_M_2
  )
}

stopCluster(cl)

# result for the cumulative calculations.
results_df <- do.call(rbind, results_raw) %>%
  arrange(sim) %>%
  dplyr::filter(converged == TRUE) %>%
  group_by(sim) %>%
  slice(1) %>%              
  ungroup() %>%
  mutate(
    cum_reject  = cumsum(reject_05),
    cum_valid   = row_number(),
    cum_type1   = cum_reject / cum_valid
  )


n_valid <- nrow(results_df)
n_failed <- N_SIM - n_valid

mc_se <- function(p, m) sqrt(p * (1 - p) / pmax(m, 1))
alpha_vec <- c(0.05, 0.10)

results_multi <- lapply(alpha_vec, function(a) {
  results_df %>%
    mutate(
      reject_a = as.integer(pval < a),
      cum_reject = cumsum(reject_a),
      cum_type1 = cum_reject / row_number(),
      mc_se = mc_se(cum_type1, row_number()),
      ci_lo = pmax(0, cum_type1 - 1.96 * mc_se),
      ci_hi = pmin(1, cum_type1 + 1.96 * mc_se),
      alpha = a,
      alpha_label = paste0("alpha == ", a)
    ) %>%
    dplyr::select(sim, cum_valid, cum_type1, mc_se, ci_lo, ci_hi, alpha, alpha_label)
}) %>%
  bind_rows() %>%
  mutate(alpha_label = factor(alpha_label,
                              levels = paste0("alpha == ", sort(alpha_vec))))
summary_table <- results_multi %>%
  group_by(alpha) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  dplyr::select(alpha, cum_valid, cum_type1, ci_lo, ci_hi) %>%
  mutate(across(c(cum_type1, ci_lo, ci_hi), ~ round(.x, 4)))

# Figure:Stage 1 rejection rates at two nominal levels
p1 <- ggplot(results_multi,
             aes(x = cum_valid, y = cum_type1,
                 color = alpha_label, group = alpha_label)) +
  geom_line(linewidth = 0.75) +
  geom_hline(
    data = data.frame(alpha = alpha_vec,
                      alpha_label = paste0("alpha == ", alpha_vec)),
    aes(yintercept = alpha, color = alpha_label),
    linetype = "dashed", linewidth = 0.5
  ) +
  facet_wrap(~ alpha_label, scales = "free_y", ncol = 2,
             labeller = label_parsed) +
  xlim(200, max(results_multi$cum_valid)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_color_manual(
    values = c("alpha == 0.05" = "#d6604d",
               "alpha == 0.1" = "#4dac26")
  ) +
  labs(
    title = "Cumulative Type I Error Rate",  
    x = "Number of Simulations",
    y = "Rejection Rate"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
ggsave(file.path(SAVE_DIR, "fig1_cumulative_type1.pdf"), p1,
       width = 10, height = 4, create.dir = TRUE)


