##================================================================
## Fri Apr 10 20:51:00 2026
##================================================================

##================================================================
## Parameters of model Eq.(1)
##================================================================

x0    <- 1
gt    <- expression(4*t + 1)
k     <- 0
Alpha <- -1/2
Beta  <- 2
Gamma <- 1/4

##================================================================
## Model Eq.(1) displayed
##================================================================

g_k_1      <- Deriv(gt, "t", nderiv = k + 1)[[1]]
g_k        <- Deriv(gt, "t", nderiv = k)[[1]]
drift_raw  <- bquote(.(Alpha) * (.(g_k_1) / .(g_k)) * x + .(Beta) * (.(g_k))^.(Alpha))
diff_raw   <- bquote(.(Gamma) * (.(g_k))^.(Alpha))
drift_simp <- Simplify(drift_raw)
diff_simp  <- Simplify(diff_raw)
drift_str  <- paste(deparse(drift_simp), collapse = " ")
diff_str   <- paste(deparse(diff_simp), collapse = " ")
drift_final <- gsub("\\bx\\b", "X(t)", drift_str)
diff_final  <- gsub("\\bx\\b", "X(t)", diff_str)
cat(sprintf("--- Model FSDE ---\ndX(t) = [ %s ] dt + [ %s ] dW(t)^H\n\n",
            drift_final, diff_final))

##================================================================
## Parameters of simulation
##================================================================

T    <- 10
N    <- 1000
M    <- 100000
ncpu <- 12
bw   <- "nrd0"

##================================================================
## Change p_vals and H_vals as needed
##================================================================

p_vals <- seq(0, 4, by = 1)
H_vals <- c(0.5,0.6,0.7,0.8,0.9,0.99)

method_th  <- "cor1"
method_emp <- "kernel"

##================================================================
## Time grid
##================================================================

samp     <- yuima::setSampling(Terminal = T, n = N)
time_vec <- samp@grid[[1]]
t_pos    <- time_vec[-1]

##================================================================
## Scaling function
##================================================================

scaled <- TRUE
scale_mp <- function(x, p_val) {
  if (!scaled || p_val == 0) return(x)
  sign(x) * abs(x)^(1 / p_val)
}

##================================================================
## Compute WFE_theor() and WFE_emp() for each (H, p)
##================================================================

set.seed(123)
results_list <- list()
for (h in H_vals) {
  cat("-> H =", h, "\n")
  
  ## ---- Step 1: Simulate ONCE per H ----
  if (method_emp == "kernel") {
    cat("   [1/3] Simulating", M, "trajectories...\n")
    X_sim <- Sim_mod(x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
                     Gamma=Gamma, H=h, T=T, N=N, M=M, ncpu=ncpu)
    
    ## ---- Step 2: KDE cache built ONCE per H ----
    cat("   [2/3] Building KDE cache...\n")
    idx_pos     <- which(as.numeric(time(X_sim)) > 0)
    kde_cache_h <- lapply(idx_pos, function(i) {
      x_i <- as.numeric(X_sim[i, ])
      if (var(x_i) < .Machine$double.eps || length(unique(x_i)) < 2) return(NULL)
      kde <- tryCatch(
        density(x_i, bw = bw, n = 512),
        error = function(e) tryCatch(
          density(x_i, bw = bw, n = 512),
          error = function(e2) NULL
        )
      )
      if (is.null(kde)) return(NULL)
      p_hat <- approx(kde$x, kde$y, xout = x_i, rule = 2)$y
      list(x = x_i, log_p = -log(pmax(p_hat, .Machine$double.eps)))
    })
  } else {
    kde_cache_h <- NULL
  }

  ## ---- Step 3: Moments pre-computed ONCE per H ----
  cat("   [3/3] Pre-computing moments (q = 0 to", max(p_vals) + 1, ")...\n")
  all_q       <- 0:(max(p_vals) + 1)
  mom_cache_h <- setNames(
    lapply(all_q, function(q) {
      if (q == 0) return(rep(1, length(t_pos)))
      HO_mom(x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
             Gamma=Gamma, H=h, p=q, t=t_pos)
    }),
    as.character(all_q)
  )
  for (p in p_vals) {
    cat("   -> p =", p, "\n")
    WFE_th_list <- WFE_theor(
      x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
      Gamma=Gamma, H=h, p=p, t=t_pos,
      methods   = method_th,
      mom_cache = mom_cache_h
    )
    WFE_em <- WFE_emp(
      kde_cache = kde_cache_h,
      x0=x0, gt=gt, k=k, Alpha=Alpha, Beta=Beta,
      Gamma=Gamma, H=h, p=p,
      T=T, N=N, M=M, ncpu=ncpu,
      method = method_emp,
      bw     = bw
    )
    results_list[[paste0("H_", h, "_p_", p, "_th")]] <- data.frame(
      time   = t_pos,
      value  = scale_mp(WFE_th_list[[method_th]], p),
      p      = p,
      H      = h,
      type   = method_th,
      source = "Exact"
    )
    results_list[[paste0("H_", h, "_p_", p, "_emp")]] <- data.frame(
      time   = t_pos,
      value  = scale_mp(WFE_em, p),
      p      = p,
      H      = h,
      type   = method_emp,
      source = "Estimated"
    )
  }
}

##================================================================
## Combine all results
##================================================================

df_all <- bind_rows(results_list) %>%
  filter(is.finite(value)) %>%
  mutate(
    H_lab = paste0("bold(H == '", H, "')"),
    p_fac = factor(p)
  )

##================================================================
## ggplot
##================================================================

p_unique    <- sort(unique(df_all$p))
make_label  <- function(p_val) paste0("p = ", as.character(MASS::fractions(p_val)))
labels_expr <- sapply(p_unique, make_label)
Color <- ggsci::pal_lancet("lanonc", alpha = 0.5)(length(p_unique))
#Color        <- ggsci::pal_jama(alpha = 0.65)(length(p_unique))
names(Color) <- as.character(p_unique)

ltype_map <- c("Exact" = "dashed", "Estimated" = "solid")

plot_WFE <- ggplot() +
  geom_line(
    data = df_all %>% filter(source == "Estimated"),
    aes(x        = time,
        y        = value,
        group    = interaction(p_fac, type),
        color    = factor(p),
        linetype = source),
    linewidth = 1
  ) +
  geom_line(
    data = df_all %>% filter(source == "Exact"),
    aes(x        = time,
        y        = value,
        group    = interaction(p_fac, type),
        linetype = source),
    color     = "black",
    linewidth = 0.8
  ) +
  scale_linetype_manual(name = NULL, values = ltype_map) +
  scale_color_manual(name = NULL, values = Color, labels = labels_expr) +
  scale_x_continuous(
    #trans   = "log1p",
    breaks = pretty_breaks(n = 12),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 10),
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    x = "Time",
    y = bquote("Dynamic Weighted Fractional Entropy")
  ) +
  facet_wrap(
    ~ H_lab,
    labeller = label_parsed,
    nrow = 2,
    ncol = 3
  ) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.spacing    = unit(0, "lines"),
    legend.position      = "top",
    legend.justification = "centre", 
    legend.box           = "horizontal",
    legend.spacing.x     = unit(5, "cm"),
    legend.key.width     = unit(1.8, "cm"),
    legend.text          = element_text(size = 14, face = "bold.italic"),
    legend.margin        = margin(t = 0, r = 0, b = -10, l = 0),
    legend.box.spacing   = unit(0.5, "cm"),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text        = element_text(),
    plot.margin      = unit(c(0, 0.05, 0, 0), "cm"),
    strip.text       = element_text(size = 15, face = "bold.italic",
                                    margin = margin(1, 1, 1, 1)),
    strip.background = element_rect(fill = "#f0f0f0",
                                    color = "black", linewidth = 0.75)
  ) +
  guides(
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1)),
    color    = guide_legend(
      order    = 2,
      nrow     = 1,
      theme    = theme(legend.justification = "right"),
      keywidth = unit(0.2, "cm"),
      override.aes = list(linewidth = 4, linetype = "solid")
    )
  )

plot_WFE

