
##### EM2.R — Model M_2: P(R_2 | Y_1, Y_2) for binary outcome Y_2

initialize_parameters_M_2 <- function(data) { 
 alpha: P(Y_1|X)
  d1 <- data[!is.na(data$Y_1_obs), ]
  if (nrow(d1) > 10) {
    alpha <- coef(glm(Y_1_obs ~ X, data = d1, family = binomial()))
  } else {
    alpha <- c("(Intercept)" = 0, "X" = 0.5)
  }
  
 phi: P(Y_2|X,Y_1)
  d2 <- data[complete.cases(data$Y_1_obs, data$Y_2_obs), ]
  if (nrow(d2) > 10) {
    phi <- coef(glm(Y_2_obs ~ X * Y_1_obs, data = d2, family = binomial()))
    names(phi) <- c("(Intercept)", "X", "Y_1", "X:Y_1")
  } else {
    phi <- c("(Intercept)" = 0, "X" = 0.3, "Y_1" = 0.5, "X:Y_1" = 0.2)
  }
  
 gamma: P(R_1|X,Y_1)
  gamma <- tryCatch(
    coef(glm(R_1 ~ X + Y_1_true, data = data, family = binomial())),
    error = function(e) c("(Intercept)" = 0.5, "X" = 0.2, "Y_1" = 0.2)
  )
  
 rho: P(R_2 | Y_1, Y_2) 
  rho <- tryCatch(
    coef(glm(R_2 ~ Y_1_true + Y_2_true, data = data, family = binomial())),
    error = function(e) c("(Intercept)" = 0.5, "Y_1" = 0.2, "Y_2" = 0.2)
  )
  
 imputed_samples
  imp <- vector("list", nrow(data))
  for (i in seq_len(nrow(data))) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    ao <- data$Y_1_obs[i]; yo <- data$Y_2_obs[i]
    if      (r1 == 1 && r2 == 1) imp[[i]] <- list(ct="c1", X=x, Y_1=ao, R_1=r1, Y_2=yo, R_2=r2)
    else if (r1 == 1 && r2 == 0) imp[[i]] <- list(ct="c2", X=x, Y_1=ao, R_1=r1, R_2=r2)
    else if (r1 == 0 && r2 == 1) imp[[i]] <- list(ct="c3", X=x, R_1=r1, Y_2=yo, R_2=r2)
    else                         imp[[i]] <- list(ct="c4", X=x, R_1=r1, R_2=r2)
  }
  
  list(alpha = alpha, phi = phi, gamma = gamma, rho = rho, imp = imp)
}


# E step
E_step_M_2 <- function(data, params) {
  a_ <- params$alpha; p_ <- params$phi
  g_ <- params$gamma; r_ <- params$rho
  imp <- params$imp
  wd <- list()
  
  for (i in seq_along(imp)) {
    s <- imp[[i]]; x <- s$X; r1 <- s$R_1; r2 <- s$R_2
    
    if (s$ct == "c1") {
      wd[[length(wd) + 1]] <- list(X=x, Y_1=s$Y_1, R_1=r1, Y_2=s$Y_2, R_2=r2, w=1)
      
    } else if (s$ct == "c2") {
      a <- s$Y_1; ws <- numeric(2)
      for (y in 0:1) {
        py  <- plogis(p_[1] + p_[2]*x + p_[3]*a + p_[4]*x*a)
        pry <- plogis(r_[1] + r_[2]*a + r_[3]*y)            
        ws[y+1] <- ifelse(y == 1, py, 1 - py) * (1 - pry)
      }
      ws <- ws / sum(ws)
      for (y in 0:1)
        wd[[length(wd)+1]] <- list(X=x, Y_1=a, R_1=r1, Y_2=y, R_2=r2, w=ws[y+1])
      
    } else if (s$ct == "c3") {
      y <- s$Y_2; ws <- numeric(2)
      for (a in 0:1) {
        pa  <- plogis(a_[1] + a_[2]*x);  pa  <- ifelse(a == 1, pa, 1 - pa)
        pra <- plogis(g_[1] + g_[2]*x + g_[3]*a); pra0 <- 1 - pra
        py  <- plogis(p_[1] + p_[2]*x + p_[3]*a + p_[4]*x*a)
        py  <- ifelse(y == 1, py, 1 - py)
        pry <- plogis(r_[1] + r_[2]*a + r_[3]*y)             
        ws[a+1] <- pa * pra0 * py * pry
      }
      ws <- ws / sum(ws)
      for (a in 0:1)
        wd[[length(wd)+1]] <- list(X=x, Y_1=a, R_1=r1, Y_2=y, R_2=r2, w=ws[a+1])
      
    } else {
      ws <- matrix(0, 2, 2)
      for (a in 0:1) for (y in 0:1) {
        pa  <- plogis(a_[1] + a_[2]*x);  pa  <- ifelse(a == 1, pa, 1 - pa)
        pra <- plogis(g_[1] + g_[2]*x + g_[3]*a); pra0 <- 1 - pra
        py  <- plogis(p_[1] + p_[2]*x + p_[3]*a + p_[4]*x*a)
        py  <- ifelse(y == 1, py, 1 - py)
        pry <- plogis(r_[1] + r_[2]*a + r_[3]*y); pry0 <- 1 - pry  
        ws[a+1, y+1] <- pa * pra0 * py * pry0
      }
      ws <- ws / sum(ws)
      for (a in 0:1) for (y in 0:1)
        wd[[length(wd)+1]] <- list(X=x, Y_1=a, R_1=r1, Y_2=y, R_2=r2, w=ws[a+1,y+1])
    }
  }
  wd
}


# M step 
M_step_M_2 <- function(wd, imp) {
  df <- data.frame(
    X   = sapply(wd, `[[`, "X"),
    Y_1   = sapply(wd, `[[`, "Y_1"),
    R_1 = sapply(wd, `[[`, "R_1"),
    Y_2   = sapply(wd, `[[`, "Y_2"),
    R_2 = sapply(wd, `[[`, "R_2"),
    w   = sapply(wd, `[[`, "w")
  )
  
  alpha <- coef(glm(Y_1   ~ X,           data = df, family = binomial(), weights = w))
  phi   <- coef(glm(Y_2   ~ X + Y_1 + X:Y_1, data = df, family = binomial(), weights = w))
  gamma <- coef(glm(R_1 ~ X + Y_1,       data = df, family = binomial(), weights = w))
  rho   <- coef(glm(R_2 ~ Y_1 + Y_2,       data = df, family = binomial(), weights = w)) # Y_1 and Y_2 are included in model M_2.
  
  list(alpha = alpha, phi = phi, gamma = gamma, rho = rho, imp = imp)
}

# Observed-data log-likelihood 
compute_observed_loglik_M_2 <- function(data, params) {
  a_ <- params$alpha; p_ <- params$phi
  g_ <- params$gamma; r_ <- params$rho
  ll <- 0
  
  for (i in 1:nrow(data)) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    ao <- data$Y_1_obs[i]; yo <- data$Y_2_obs[i]
    pi <- 0
    for (a in 0:1) { if (!is.na(ao) && a != ao) next
      for (y in 0:1) { if (!is.na(yo) && y != yo) next
        pa  <- plogis(a_[1] + a_[2]*x);  pa  <- ifelse(a == 1, pa, 1 - pa)
        pra <- plogis(g_[1] + g_[2]*x + g_[3]*a)
        pra <- ifelse(r1 == 1, pra, 1 - pra)
        py  <- plogis(p_[1] + p_[2]*x + p_[3]*a + p_[4]*x*a)
        py  <- ifelse(y == 1, py, 1 - py)
        pry <- plogis(r_[1] + r_[2]*a + r_[3]*y)           
        pry <- ifelse(r2 == 1, pry, 1 - pry)
        pi  <- pi + pa * pra * py * pry
      }
    }
    ll <- ll + log(pi + 1e-10)
  }
  ll
}


# Individual observed-data log-likelihood 
compute_individual_loglik_M_2 <- function(data, params) {
  a_ <- params$alpha; p_ <- params$phi
  g_ <- params$gamma; r_ <- params$rho
  n <- nrow(data); ill <- numeric(n)
  
  for (i in 1:n) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    ao <- data$Y_1_obs[i]; yo <- data$Y_2_obs[i]
    pi <- 0
    for (a in 0:1) { if (!is.na(ao) && a != ao) next
      for (y in 0:1) { if (!is.na(yo) && y != yo) next
        pa  <- plogis(a_[1] + a_[2]*x);  pa  <- ifelse(a == 1, pa, 1 - pa)
        pra <- plogis(g_[1] + g_[2]*x + g_[3]*a)
        pra <- ifelse(r1 == 1, pra, 1 - pra)
        py  <- plogis(p_[1] + p_[2]*x + p_[3]*a + p_[4]*x*a)
        py  <- ifelse(y == 1, py, 1 - py)
        pry <- plogis(r_[1] + r_[2]*a + r_[3]*y)             
        pry <- ifelse(r2 == 1, pry, 1 - pry)
        pi  <- pi + pa * pra * py * pry
      }
    }
    ill[i] <- log(pi + 1e-10)
  }
  ill
}


# EM algorithm for model M_2. 
EM_algorithm_M_2 <- function(data, max_iter = 500, tol = 1e-4) {
  params <- initialize_parameters_M_2(data)
  lh <- numeric(max_iter)
  
  for (it in 1:max_iter) {
    wd <- E_step_M_2(data, params)
    params <- M_step_M_2(wd, params$imp)
    lh[it] <- compute_observed_loglik_M_2(data, params)
    if (it > 1 && abs(lh[it] - lh[it-1]) < tol) {
      lh <- lh[1:it]
      return(list(params = params, loglik_history = lh, converged = TRUE))
    }
  }
  list(params = params, loglik_history = lh, converged = FALSE)
}


# True parameters
compute_true_parameters_M_2 <- function(data) {
  list(
    alpha = coef(glm(Y_1_true ~ X,              data = data, family = binomial())),
    phi   = coef(glm(Y_2_true ~ X * Y_1_true,     data = data, family = binomial())),
    gamma = coef(glm(R_1    ~ X + Y_1_true,     data = data, family = binomial())),
    rho   = coef(glm(R_2    ~ Y_1_true + Y_2_true, data = data, family = binomial()))
  )
}


# Convert parameters. 
params_to_theta_M_2 <- function(params) {
  c(params$alpha, params$phi, params$gamma, params$rho)
}

theta_to_params_M_2 <- function(theta) {
  list(alpha = theta[1:2], phi = theta[3:6], gamma = theta[7:9], rho = theta[10:12])
}

loglik_theta_M_2 <- function(theta, data) {
  compute_observed_loglik_M_2(data, theta_to_params_M_2(theta))
}

individual_loglik_theta_M_2 <- function(theta, data) {
  compute_individual_loglik_M_2(data, theta_to_params_M_2(theta))
}
