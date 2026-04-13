#================================================================
# Sat Mar 21 16:27:03 2026
#================================================================

options(future.globals.maxSize = +Inf, warn = -1, digits = 5)

#================================================================
# load R packages
#================================================================

if (!require("Sim.DiffProc")) { install.packages("Sim.DiffProc"); library("Sim.DiffProc") }
if (!require("yuima"))        { install.packages("yuima");        library("yuima")        }
if (!require("ttutils"))      { install.packages("ttutils");      library("ttutils")      }
if (!require("EQL"))          { install.packages("EQL");          library("EQL")          }
if (!require("gsl"))          { install.packages("gsl");          library("gsl")          }
if (!require("ggsci"))        { install.packages("ggsci");        library("ggsci")        }
if (!require("furrr"))        { install.packages("furrr");        library("furrr")        }
if (!require("Deriv"))        { install.packages("Deriv");        library("Deriv")        }
if (!require("zoo"))          { install.packages("zoo");          library("zoo")          }
if (!require("ggplot2"))      { install.packages("ggplot2");      library("ggplot2")      }
if (!require("tidyr"))        { install.packages("tidyr");        library("tidyr")        }
if (!require("dplyr"))        { install.packages("dplyr");        library("dplyr")        }
if (!require("scales"))       { install.packages("scales");       library("scales")       }
if (!require("MASS"))         { install.packages("MASS");         library("MASS")         }

#================================================================
# Simulation model Eq.(1)
#================================================================

Sim_mod <- function(x0, gt, k = NULL, Alpha = 1, Beta = 1, Gamma = 1,
                    H = 0.75, T = 1, N = 1000, M = 100, ncpu = 12) {
  
  if (H < 0.5 || H >= 1) stop("'H' must be in '(0, 1)'")
  if (Gamma <= 0)        stop("'Gamma' must be numeric > 0")
  if (!is.expression(gt)) stop("'g(t)' must be an expression in 't'")
  
  if (is.null(k)) {
    Gt <- gt
  } else {
    if (!ttutils::isInteger(k)) stop("k must be a non-negative integer.")
    Gt <- Deriv::Deriv(gt, "t", nderiv = k)
  }
  
  samp     <- yuima::setSampling(Terminal = T, n = N)
  temps    <- samp@grid[[1]]
  G_t_vec  <- sapply(temps, function(ti) eval(Gt, envir = list(t = ti)))
  G_t_0    <- eval(Gt, envir = list(t = 0))
  det_part <- G_t_vec^Alpha * (x0 * G_t_0^(-Alpha) + Beta * temps)
  cl <- parallel::makeCluster(ncpu)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, varlist = c("H", "T", "N"), envir = environment())
  parallel::clusterEvalQ(cl, {
    m     <- 2 * N
    k_seq <- 0:N
    acv   <- 0.5 * ((k_seq+1)^(2*H) - 2*k_seq^(2*H) + abs(k_seq-1)^(2*H))
    cr    <- c(acv[1:N], acv[N+1], acv[N:2])  
    .lam  <- pmax(Re(fft(cr)), 0)             
    .sl   <- sqrt(.lam)                        
    .sl2  <- sqrt(.lam / 2)                    
    .idx  <- 2:N                               
    rfbm <- function() {
      Z1 <- rnorm(m)
      Z2 <- rnorm(m)
      W  <- complex(length.out = m)
      W[1]              <- .sl[1]       * Z1[1]
      W[N + 1]          <- .sl[N + 1]   * Z1[N + 1]
      W[.idx]           <- .sl2[.idx]   * complex(real = Z1[.idx], imaginary = Z2[.idx])
      W[m + 2 - .idx]   <- Conj(W[.idx])
      fgn <- Re(fft(W, inverse = TRUE))[1:N] / sqrt(m)
      c(0, cumsum(fgn) * (T / N)^H)
    }
  })
  idx_list <- parallel::splitIndices(M, ncpu)
  W_h_list <- parallel::parLapply(cl, idx_list, function(idx) {
    vapply(seq_along(idx), function(i) rfbm(), numeric(N + 1))
  })
  W_h  <- do.call(cbind, W_h_list)
  X    <- det_part + Gamma * G_t_vec^Alpha * W_h
  name <- if (M > 1) paste0("X", 1:M) else "X"
  X    <- ts(X, start = 0, deltat = samp@delta, names = name)
  return(invisible(X))
}

#================================================================
# Exact moments: Lemma 1
#================================================================

HO_mom <- function(x0, gt, k = NULL, Alpha = 1, Beta = 1, Gamma = 1,
                   H = 0.75, p = 1, t) {
  
  if (!inherits(gt, "expression")) stop("'gt' must be an expression in 't'")
  if (H < 0.5 | H >= 1) stop("'H' must be in '(1/2,1)'")
  if (p <= 0)            stop("'p' must be numeric > 0'")
  if (Gamma <= 0)        stop("'Gamma' must be numeric > 0'")
  
  t <- as.numeric(t)
  is_single_t <- length(t) == 1
  if (all(t == 0)) return(x0^p)
  
  if (is.null(k)) {
    Gt <- gt
  } else {
    if (!ttutils::isInteger(k)) stop("k must be a non-negative integer.")
    Gt <- Deriv::Deriv(gt, "t", nderiv = k)
  }
  
  ## Pre-compute G_t for ALL t in one pass
  G_t_vec <- vapply(t, function(ti) eval(Gt, envir = list(t = ti)), numeric(1))
  G_t_0   <- eval(Gt, envir = list(t = 0))
  if (length(G_t_0) == 0) stop("Evaluation of g(t) at t=0 failed.")
  
  ## Separate t=0 and t>0
  idx_pos          <- t > 0
  result           <- numeric(length(t))
  result[!idx_pos] <- x0^p
  if (!any(idx_pos)) return(if (is_single_t) result else invisible(result))
  t_p  <- t[idx_pos]
  Gt_p <- G_t_vec[idx_pos]
  A_vec <- 2^(-0.5*p) * ((Gamma * Gt_p^Alpha) / 1i)^p * t_p^(H*p)
  B_vec <- 0.5 * sqrt(2) * 1i * ((G_t_0^(-Alpha) * x0 + Beta * t_p) / Gamma) * t_p^(-H)
  if (ttutils::isInteger(p)) {
    Herm_vec <- sapply(B_vec, function(b) EQL::hermite(x = b, n = p, prob = FALSE))
  } else {
    B2_vec <- Re(B_vec^2)
    ## Pre-compute p-dependent constants ONCE
    c1 <- 2^p * sqrt(pi) / gamma(0.5 - 0.5*p)
    c2 <- 2^(p+1) * sqrt(pi) / gamma(-0.5*p)
    
    Herm_vec <- vapply(seq_along(t_p), function(i) {
      M1 <- gsl::hyperg_1F1(-0.5*p,       0.5, B2_vec[i])
      M2 <- gsl::hyperg_1F1(0.5 - 0.5*p,  1.5, B2_vec[i])
      c1 * M1 - c2 * B_vec[i] * M2
    }, complex(1))
  }
  
  result[idx_pos] <- Re(A_vec * Herm_vec)
  return(if (is_single_t) result else invisible(result))
}

#================================================================
# Theoretical weighted fractional entropy H_w^(p)(X_t)
#================================================================

WFE_theor <- function(x0, gt, k = NULL, Alpha = 1, Beta = 1,
                      Gamma = 1, H = 0.75, p = 0, t,
                      methods   = c("cor1", "cor2", "cor3"),
                      mom_cache = NULL) {
  
  if (!inherits(gt, "expression")) stop("'gt' must be an expression in 't'")
  if (H < 0.5 | H >= 1) stop("'H' must be in '(1/2,1)'")
  if (p < 0)           stop("'p' must be numeric >= 0")
  if (Gamma <= 0)      stop("'Gamma' must be numeric > 0")
  if (any(c("cor2", "cor3") %in% methods) && !ttutils::isInteger(p)) {
    warning("Corollaries 2 and 3 require p in N. They will be skipped for p = ", p)
    methods <- setdiff(methods, c("cor2", "cor3"))
  }
  
  t <- as.numeric(t)
  if (any(t <= 0)) {
    warning("H_w^(p)(X_t) is defined for t > 0 only. Values t <= 0 removed.")
    t <- t[t > 0]
  }
  
  if (is.null(k)) {
    Gt <- gt
  } else {
    if (!ttutils::isInteger(k)) stop("k must be a non-negative integer.")
    Gt <- Deriv::Deriv(gt, "t", nderiv = k)
  }
  G_t   <- function(tv) vapply(tv, function(ti) eval(Gt, envir = list(t = ti)), numeric(1))
  G_t_0 <- G_t(0)
  G_t_vec <- G_t(t)
  mu_vec  <- G_t_vec^Alpha * (x0 * G_t_0^(-Alpha) + Beta * t)
  s2_vec  <- Gamma^2 * G_t_vec^(2 * Alpha) * t^(2 * H)
  HS_vec  <- log(Gamma * sqrt(2 * pi * exp(1))) + Alpha * log(G_t_vec) + H * log(t)
  get_mom_vec <- function(q) {
    if (q == 0) return(rep(1, length(t)))
    q_char <- as.character(q)
    if (!is.null(mom_cache) && q_char %in% names(mom_cache))
      return(as.numeric(mom_cache[[q_char]]))
    HO_mom(x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
           Gamma=Gamma, H=H, p=q, t=t)
  }
  
  result <- list()
  
  if ("cor1" %in% methods) {
    Ep_vec  <- get_mom_vec(p)
    Ep1_vec <- get_mom_vec(p + 1)
    Cov_p   <- Ep1_vec - mu_vec * Ep_vec
    result[["cor1"]] <- (HS_vec + p/2) * Ep_vec -
      (mu_vec / (2 * s2_vec)) * Cov_p
  }
  
  if ("cor2" %in% methods) {
    Ep_vec   <- get_mom_vec(p)
    Ep_1_vec <- if (p == 0) rep(0, length(t)) else get_mom_vec(p - 1)
    result[["cor2"]] <- (HS_vec + p/2) * Ep_vec -
      (p/2) * mu_vec * Ep_1_vec
  }
  
  if ("cor3" %in% methods) {
    Ep_vec   <- get_mom_vec(p)
    Ep_2_vec <- if (p <= 1) rep(0, length(t)) else get_mom_vec(p - 2)
    result[["cor3"]] <- HS_vec * Ep_vec +
      (p * (p - 1) / 2) * s2_vec * Ep_2_vec
  }
  return(invisible(result))
}

#================================================================
# Empirical weighted fractional entropy H_w^(p)(X_t)
#================================================================

WFE_emp <- function(X = NULL, x0, gt, k = NULL, Alpha = 1, Beta = 1,
                    Gamma = 1, H = 0.75, p = 0,
                    T = 1, N = 1000, M = 100, ncpu = 12,
                    method    = c("kernel", "gaussian"),
                    bw        = "SJ",
                    kde_cache = NULL) {
  
  if (!inherits(gt, "expression")) stop("'gt' must be an expression in 't'")
  if (H < 0.5 | H >= 1) stop("'H' must be in '(1/2,1)'")
  if (p < 0)           stop("'p' must be numeric >= 0")
  if (Gamma <= 0)      stop("'Gamma' must be numeric > 0")
  
  method <- match.arg(method)
  
  if (is.null(k)) {
    Gt <- gt
  } else {
    if (!ttutils::isInteger(k)) stop("k must be a non-negative integer.")
    Gt <- Deriv::Deriv(gt, "t", nderiv = k)
  }
  G_t   <- function(tv) vapply(tv, function(ti) eval(Gt, envir = list(t = ti)), numeric(1))
  G_t_0 <- G_t(0)
  mu_X  <- function(ti) G_t(ti)^Alpha * (x0 * G_t_0^(-Alpha) + Beta * ti)
  s2_X  <- function(ti) Gamma^2 * G_t(ti)^(2 * Alpha) * ti^(2 * H)
  
  samp      <- yuima::setSampling(Terminal = T, n = N)
  temps     <- samp@grid[[1]]
  temps_pos <- temps[temps > 0]
  
  if (method == "kernel") {
    if (is.null(kde_cache)) {
      if (is.null(X)) {
        X <- Sim_mod(x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
                     Gamma=Gamma, H=H, T=T, N=N, M=M, ncpu=ncpu)
      }
      idx_pos   <- which(as.numeric(time(X)) > 0)
      kde_cache <- lapply(idx_pos, function(i) {
        x_i <- as.numeric(X[i, ])
        if (var(x_i) < .Machine$double.eps || length(unique(x_i)) < 2) return(NULL)
        kde <- tryCatch(
          density(x_i, bw = bw, n = 512),
          error = function(e) tryCatch(
            density(x_i, bw = "nrd0", n = 512),
            error = function(e2) NULL
          )
        )
        if (is.null(kde)) return(NULL)
        p_hat <- approx(kde$x, kde$y, xout = x_i, rule = 2)$y
        list(x = x_i, log_p = -log(pmax(p_hat, .Machine$double.eps)))
      })
    }

    WFE_t <- vapply(seq_along(kde_cache), function(i) {
      if (is.null(kde_cache[[i]])) return(NA_real_)
      mean(kde_cache[[i]]$x^p * kde_cache[[i]]$log_p)
    }, numeric(1))
    
    if (any(is.na(WFE_t)))
      WFE_t <- zoo::na.approx(WFE_t, na.rm = FALSE)   
  } else { 
    WFE_t <- vapply(temps_pos, function(ti) {
      mu_i  <- mu_X(ti)
      sd_i  <- sqrt(s2_X(ti))
      x_i   <- rnorm(M, mean = mu_i, sd = sd_i)
      p_hat <- dnorm(x_i, mean = mu_i, sd = sd_i)
      p_hat <- pmax(p_hat, .Machine$double.eps)
      mean(x_i^p * (-log(p_hat)))
    }, numeric(1))
  } 
  return(invisible(WFE_t))
}
