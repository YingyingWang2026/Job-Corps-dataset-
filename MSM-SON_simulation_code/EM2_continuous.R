
##### EM2.R — Model M_2: P(R_2 | Y_1, Y_2) for continuous outcome Y_2

initialize_parameters_M_2 <- function(data, M = 100) {
  # alpha: P(Y_1|X)
  d1 <- data[!is.na(data$Y_1_obs), ]
  if (nrow(d1) > 10) {
    alpha <- tryCatch(coef(glm(Y_1_obs ~ X, data = d1, family = binomial())),
                      error = function(e) c("(Intercept)" = 0, "X" = 0.5))
  } else {
    alpha <- c("(Intercept)" = 0, "X" = 0.5)
  }
  
  # phi: Y_2|X,Y_1 
  d2 <- data[complete.cases(data$Y_1_obs, data$Y_2_obs), ]
  if (nrow(d2) > 10) {
    phi <- tryCatch({
      m <- lm(Y_2_obs ~ X * Y_1_obs, data = d2)
      list(beta = coef(m), sigma = summary(m)$sigma)
    }, error = function(e) {
      list(beta = c("(Intercept)"=-1, "X"=0.8, "Y_1_obs"=1.5, "X:Y_1_obs"=0.5), sigma=1.0)
    })
    names(phi$beta) <- c("(Intercept)", "X", "Y_1_obs", "X:Y_1_obs")
  } else {
    phi <- list(beta = c("(Intercept)"=-1, "X"=0.8, "Y_1_obs"=1.5, "X:Y_1_obs"=0.5), sigma=1.0)
  }
  
  # gamma: P(R_1|X,Y_1)
  gamma <- tryCatch({
    da <- data
    if (sum(is.na(data$Y_1_obs)) > 0) {
      py1 <- plogis(alpha[1] + alpha[2]*data$X)
      da$Y_1_imp <- ifelse(is.na(data$Y_1_obs), as.numeric(py1 > 0.5), data$Y_1_obs)
    } else { da$Y_1_imp <- data$Y_1_obs }
    g <- coef(glm(R_1 ~ X + Y_1_imp, data = da, family = binomial()))
    names(g) <- c("(Intercept)", "X", "Y_1_obs")
    g
  }, error = function(e) c("(Intercept)"=0.8, "X"=0.3, "Y_1_obs"=-4))
  
  # rho: P(R_2 | Y_1, Y_2) 
  rho <- tryCatch({
    dr <- data
    if (sum(is.na(data$Y_1_obs)) > 0) {
      py1 <- plogis(alpha[1] + alpha[2]*data$X)
      dr$Y_1_imp <- ifelse(is.na(data$Y_1_obs), as.numeric(py1 > 0.5), data$Y_1_obs)
    } else { dr$Y_1_imp <- data$Y_1_obs }
    if (sum(is.na(data$Y_2_obs)) > 0) {
      my <- mean(data$Y_2_obs, na.rm = TRUE)
      dr$Y_2_imp <- ifelse(is.na(data$Y_2_obs), ifelse(!is.na(my), my, 0), data$Y_2_obs)
    } else { dr$Y_2_imp <- data$Y_2_obs }
    r <- coef(glm(R_2 ~ Y_1_imp + Y_2_imp, data = dr, family = binomial()))
    names(r) <- c("(Intercept)", "Y_1_obs", "Y_2_obs")
    r
  }, error = function(e) c("(Intercept)"=0.8, "Y_1_obs"=2.5, "Y_2_obs"=0.3))
  
  # impute
  imp <- vector("list", nrow(data))
  for (i in seq_len(nrow(data))) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    y1o <- data$Y_1_obs[i]; y2o <- data$Y_2_obs[i]
    
    if (r1 == 1 && r2 == 1) {
      imp[[i]] <- list(ct="c1", X=x, Y_1=y1o, R_1=r1, Y_2=y2o, R_2=r2)
    } else if (r1 == 1 && r2 == 0) {
      mu_p <- phi$beta[1] + phi$beta[2]*x + phi$beta[3]*y1o + phi$beta[4]*x*y1o
      y2s <- rnorm(M, mean = mu_p, sd = phi$sigma)
      imp[[i]] <- list(ct="c2", X=x, Y_1=y1o, R_1=r1, R_2=r2,
                       Y_2_samples=y2s, h_mu2=mu_p, h_sigma2=phi$sigma)
    } else if (r1 == 0 && r2 == 1) {
      imp[[i]] <- list(ct="c3", X=x, R_1=r1, Y_2=y2o, R_2=r2)
    } else {
      pY1_1 <- plogis(alpha[1] + alpha[2]*x)
      pR1_0_Y1_1 <- 1 - plogis(gamma[1] + gamma[2]*x + gamma[3]*1)
      pR1_0_Y1_0 <- 1 - plogis(gamma[1] + gamma[2]*x + gamma[3]*0)
      n1 <- pY1_1 * pR1_0_Y1_1; n0 <- (1 - pY1_1) * pR1_0_Y1_0
      prob_y1 <- c(n0/(n0+n1), n1/(n0+n1))
      My10 <- round(M * prob_y1[1]); My11 <- M - My10
      
      sl <- list()
      for (y1 in 0:1) {
        My1 <- if (y1 == 0) My10 else My11
        if (My1 > 0) {
          mu_p <- phi$beta[1] + phi$beta[2]*x + phi$beta[3]*y1 + phi$beta[4]*x*y1
          y2s <- rnorm(My1, mean = mu_p, sd = phi$sigma)
          sl[[y1+1]] <- list(Y_1=y1, Y_2_samples=y2s, h_mu2=mu_p, h_sigma2=phi$sigma,
                            q_Y1_prob=prob_y1[y1+1])
        }
      }
      imp[[i]] <- list(ct="c4", X=x, R_1=r1, R_2=r2, samples_list=sl)
    }
  }
  
  list(alpha=alpha, phi=phi, gamma=gamma, rho=rho, imp=imp)
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
      wd[[length(wd)+1]] <- list(X=x, Y_1=s$Y_1, R_1=r1, Y_2=s$Y_2, R_2=r2, w=1)
      
    } else if (s$ct == "c2") {
      y1 <- s$Y_1; y2s <- s$Y_2_samples; hmu2 <- s$h_mu2; hsd2 <- s$h_sigma2
      Mm <- length(y2s); ws <- numeric(Mm)
      for (ss in 1:Mm) {
        y2_s <- y2s[ss]
        mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1+p_$beta[4]*x*y1
        fy2  <- dnorm(y2_s, mean=mu_y2, sd=p_$sigma)
        pr2_0 <- 1 - plogis(r_[1]+r_[2]*y1+r_[3]*y2_s)           
        hd  <- dnorm(y2_s, mean=hmu2, sd=hsd2)
        ws[ss] <- fy2 * pr2_0 / (hd + 1e-10)
      }
      ws <- ws/sum(ws)
      for (ss in 1:Mm)
        wd[[length(wd)+1]] <- list(X=x, Y_1=y1, R_1=r1, Y_2=y2s[ss], R_2=r2, w=ws[ss])
      
    } else if (s$ct == "c3") {
      y2 <- s$Y_2; ws <- numeric(2)
      for (y1 in 0:1) {
        py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1==1, py1, 1-py1)
        pr10 <- 1 - plogis(g_[1]+g_[2]*x+g_[3]*y1)
        mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1+p_$beta[4]*x*y1
        fy2 <- dnorm(y2, mean=mu_y2, sd=p_$sigma)
        pr2_1 <- plogis(r_[1]+r_[2]*y1+r_[3]*y2)                 
        ws[y1+1] <- py1 * pr10 * fy2 * pr2_1
      }
      ws <- ws/sum(ws)
      for (y1 in 0:1)
        wd[[length(wd)+1]] <- list(X=x, Y_1=y1, R_1=r1, Y_2=y2, R_2=r2, w=ws[y1+1])
      
    } else {
      sl <- s$samples_list
      all_w <- numeric(); all_s <- list()
      for (ag in sl) {
        y1 <- ag$Y_1; y2s <- ag$Y_2_samples; hmu2 <- ag$h_mu2; hsd2 <- ag$h_sigma2
        qY1 <- ag$q_Y1_prob; My1 <- length(y2s); wa <- numeric(My1)
        for (ss in 1:My1) {
          y2_s <- y2s[ss]
          py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1==1, py1, 1-py1)
          pr10 <- 1 - plogis(g_[1]+g_[2]*x+g_[3]*y1)
          mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1+p_$beta[4]*x*y1
          fy2 <- dnorm(y2_s, mean=mu_y2, sd=p_$sigma)
          pr2_0 <- 1 - plogis(r_[1]+r_[2]*y1+r_[3]*y2_s)         
          ft <- py1 * pr10 * fy2 * pr2_0
          hY2 <- dnorm(y2_s, mean=hmu2, sd=hsd2)
          wa[ss] <- ft / (qY1 * hY2 + 1e-10)
        }
        all_w <- c(all_w, wa)
        for (ss in 1:My1) all_s[[length(all_s)+1]] <- list(Y_1=y1, Y_2=y2s[ss])
      }
      all_w <- all_w/sum(all_w)
      for (idx in seq_along(all_s))
        wd[[length(wd)+1]] <- list(X=x, Y_1=all_s[[idx]]$Y_1, R_1=r1,
                                   Y_2=all_s[[idx]]$Y_2, R_2=r2, w=all_w[idx])
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
  
  alpha <- coef(glm(Y_1 ~ X, data=df, family=binomial(), weights=w))
  
  fit_phi <- lm(Y_2 ~ X + Y_1 + X:Y_1, data=df, weights=w)
  res_w <- residuals(fit_phi) * sqrt(df$w)
  sigma_y <- sqrt(sum(res_w^2) / sum(df$w))
  phi <- list(beta = coef(fit_phi), sigma = sigma_y)
  
  gamma <- coef(glm(R_1 ~ X + Y_1, data=df, family=binomial(), weights=w))
  
  rho <- coef(glm(R_2 ~ Y_1 + Y_2, data=df, family=binomial(), weights=w))   
  
  list(alpha=alpha, phi=phi, gamma=gamma, rho=rho, imp=imp)
}

# Calculate the observed log-likelihood
compute_observed_loglik_M_2 <- function(data, params) {
  a_ <- params$alpha; p_ <- params$phi
  g_ <- params$gamma; r_ <- params$rho
  ll <- 0; n_quad <- 20
  
  for (i in 1:nrow(data)) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    y1o <- data$Y_1_obs[i]; y2o <- data$Y_2_obs[i]
    
    if (!is.na(y1o) && !is.na(y2o)) {
      py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1o==1, py1, 1-py1)
      pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1o); pr1 <- ifelse(r1==1, pr1, 1-pr1)
      mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1o+p_$beta[4]*x*y1o
      lp <- log(py1+1e-10) + log(pr1+1e-10) +
        dnorm(y2o, mean=mu_y2, sd=p_$sigma, log=TRUE) +
        log(plogis(r_[1]+r_[2]*y1o+r_[3]*y2o)+1e-10)        
      
    } else if (!is.na(y1o) && is.na(y2o)) {
      py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1o==1, py1, 1-py1)
      pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1o); pr1 <- ifelse(r1==1, pr1, 1-pr1)
      mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1o+p_$beta[4]*x*y1o
      y2r <- seq(mu_y2-4*p_$sigma, mu_y2+4*p_$sigma, length.out=n_quad)
      dy2 <- y2r[2]-y2r[1]; ps <- 0
      for (y2v in y2r) {
        fy2 <- dnorm(y2v, mean=mu_y2, sd=p_$sigma)
        pr2_0 <- 1 - plogis(r_[1]+r_[2]*y1o+r_[3]*y2v)          
        ps <- ps + fy2 * pr2_0
      }
      ps <- ps * dy2
      lp <- log(py1+1e-10) + log(pr1+1e-10) + log(ps+1e-10)
      
    } else {
      ps <- 0
      for (y1 in 0:1) {
        py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1==1, py1, 1-py1)
        pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1); pr1 <- ifelse(r1==1, pr1, 1-pr1)
        mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1+p_$beta[4]*x*y1
        
        if (!is.na(y2o)) {
          fy2 <- dnorm(y2o, mean=mu_y2, sd=p_$sigma)
          pr2 <- plogis(r_[1]+r_[2]*y1+r_[3]*y2o)               
          ps <- ps + py1 * pr1 * fy2 * pr2
        } else {
          y2r <- seq(mu_y2-4*p_$sigma, mu_y2+4*p_$sigma, length.out=n_quad)
          dy2 <- y2r[2]-y2r[1]; y1_int <- 0
          for (y2v in y2r) {
            fy2 <- dnorm(y2v, mean=mu_y2, sd=p_$sigma)
            pr2_0 <- 1 - plogis(r_[1]+r_[2]*y1+r_[3]*y2v)       
            y1_int <- y1_int + fy2 * pr2_0
          }
          y1_int <- y1_int * dy2
          ps <- ps + py1 * pr1 * y1_int
        }
      }
      lp <- log(ps+1e-10)
    }
    ll <- ll + lp
  }
  ll
}


# Calculate the individual log-likelihood
compute_individual_loglik_M_2 <- function(data, params) {
  a_ <- params$alpha; p_ <- params$phi
  g_ <- params$gamma; r_ <- params$rho
  n <- nrow(data); ill <- numeric(n); n_quad <- 20
  
  for (i in 1:n) {
    x <- data$X[i]; r1 <- data$R_1[i]; r2 <- data$R_2[i]
    y1o <- data$Y_1_obs[i]; y2o <- data$Y_2_obs[i]
    
    if (!is.na(y1o) && !is.na(y2o)) {
      py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1o==1, py1, 1-py1)
      pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1o); pr1 <- ifelse(r1==1, pr1, 1-pr1)
      mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1o+p_$beta[4]*x*y1o
      ill[i] <- log(py1+1e-10) + log(pr1+1e-10) +
        dnorm(y2o, mean=mu_y2, sd=p_$sigma, log=TRUE) +
        log(plogis(r_[1]+r_[2]*y1o+r_[3]*y2o)+1e-10)
      
    } else if (!is.na(y1o) && is.na(y2o)) {
      py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1o==1, py1, 1-py1)
      pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1o); pr1 <- ifelse(r1==1, pr1, 1-pr1)
      mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1o+p_$beta[4]*x*y1o
      y2r <- seq(mu_y2-4*p_$sigma, mu_y2+4*p_$sigma, length.out=n_quad)
      dy2 <- y2r[2]-y2r[1]; ps <- 0
      for (y2v in y2r) ps <- ps + dnorm(y2v,mu_y2,p_$sigma)*(1-plogis(r_[1]+r_[2]*y1o+r_[3]*y2v))
      ill[i] <- log(py1+1e-10) + log(pr1+1e-10) + log(ps*dy2+1e-10)
      
    } else {
      ps <- 0
      for (y1 in 0:1) {
        py1 <- plogis(a_[1]+a_[2]*x); py1 <- ifelse(y1==1, py1, 1-py1)
        pr1 <- plogis(g_[1]+g_[2]*x+g_[3]*y1); pr1 <- ifelse(r1==1, pr1, 1-pr1)
        mu_y2 <- p_$beta[1]+p_$beta[2]*x+p_$beta[3]*y1+p_$beta[4]*x*y1
        if (!is.na(y2o)) {
          ps <- ps + py1*pr1*dnorm(y2o,mu_y2,p_$sigma)*plogis(r_[1]+r_[2]*y1+r_[3]*y2o)
        } else {
          y2r <- seq(mu_y2-4*p_$sigma, mu_y2+4*p_$sigma, length.out=n_quad)
          dy2 <- y2r[2]-y2r[1]; y1_int <- 0
          for (y2v in y2r) y1_int <- y1_int+dnorm(y2v,mu_y2,p_$sigma)*(1-plogis(r_[1]+r_[2]*y1+r_[3]*y2v))
          ps <- ps + py1*pr1*y1_int*dy2
        }
      }
      ill[i] <- log(ps+1e-10)
    }
  }
  ill
}

# EM algorithm
EM_algorithm_M_2 <- function(data, max_iter=500, tol=1e-4, M=100) {
  params <- initialize_parameters_M_2(data, M=M)
  lh <- numeric(max_iter)
  for (it in 1:max_iter) {
    wd <- E_step_M_2(data, params)
    params <- M_step_M_2(wd, params$imp)
    lh[it] <- compute_observed_loglik_M_2(data, params)
    if (it > 1 && abs(lh[it]-lh[it-1]) < tol) {
      lh <- lh[1:it]
      return(list(params=params, loglik_history=lh, converged=TRUE))
    }
  }
  list(params=params, loglik_history=lh, converged=FALSE)
}

# Calculate true parameters
compute_true_parameters_M_2 <- function(data) {
  fit_y2 <- lm(Y_2_true ~ X + Y_1_true + X:Y_1_true, data=data)
  list(
    alpha = coef(glm(Y_1_true ~ X, data=data, family=binomial())),
    phi   = list(beta = coef(fit_y2), sigma = sigma(fit_y2)),
    gamma = coef(glm(R_1 ~ X + Y_1_true, data=data, family=binomial())),
    rho   = coef(glm(R_2 ~ Y_1_true + Y_2_true, data=data, family=binomial()))
  )
}

#  params conversion 
params_to_theta_M_2 <- function(params) {
  c(params$alpha, params$phi$beta, params$phi$sigma, params$gamma, params$rho)
}

theta_to_params_M_2 <- function(theta) {
  list(
    alpha = theta[1:2],
    phi   = list(beta = theta[3:6], sigma = theta[7]),
    gamma = theta[8:10],
    rho   = theta[11:13]
  )
}

loglik_theta_M_2 <- function(theta, data) {
  compute_observed_loglik_M_2(data, theta_to_params_M_2(theta))
}

individual_loglik_theta_M_2 <- function(theta, data) {
  compute_individual_loglik_M_2(data, theta_to_params_M_2(theta))
}
