##================================================================
## Sat Apr 18 22:18:39 2026
##================================================================

##================================================================
## Parameters of model Eq.(1)
##================================================================

x0    <- 10
gt    <- expression(exp(theta*t))
k     <- 0
Alpha <- -1/2
Beta  <- 4
Gamma <- 2
theta <- 1/4

##================================================================
## Model Eq.(1) displayed
##================================================================

g_k_1      <- Deriv(gt, "t", nderiv = k + 1)[[1]]
g_k        <- Deriv(gt, "t", nderiv = k)[[1]]
drift_raw  <- bquote(.(Alpha) * (.(g_k_1) / .(g_k)) * x + .(Beta) * (.(g_k))^.(Alpha))
diff_raw   <- bquote(.(Gamma) * (.(g_k))^.(Alpha))
drift_simp <- Simplify(drift_raw)
diff_simp  <- Simplify(diff_raw)
drift_str  <- paste(deparse(drift_simp),collapse = " ")
diff_str   <- paste(deparse(diff_simp),collapse = " ")
drift_final <- gsub("\\bx\\b", "X(t)", drift_str)
diff_final  <- gsub("\\bx\\b", "X(t)", diff_str)
cat(sprintf("--- Model FSDE ---\ndX(t) = [ %s ] dt + [ %s ] dW(t)^H\n\n",
            drift_final, diff_final))

##================================================================
## Parameters of simulation
##================================================================

T    <- 50
N    <- 5000
M    <- 100000
ncpu <- 12
bw   <- "nrd0"

##================================================================
## Change p_vals and H_vals as needed
##================================================================

p_vals <- seq(0,3.5,by=1/2)
H_vals <- 0.75

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

scaled <- FALSE
scale_mp <- function(x, p_val) {
  if (!scaled || p_val == 0) return(x)
  sign(x) * abs(x)^(1 / p_val)
}

##================================================================
## Compute WFE_theor() and WFE_emp() for each (H, p)
##================================================================

set.seed(1234)
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
names(Color) <- as.character(p_unique)
ltype_map <- c("Exact" = "dashed", "Estimated" = "solid")
p_label_map <- setNames(
  paste0("p == ", sapply(p_unique, function(v) as.character(MASS::fractions(v)))),
  as.character(p_unique)
)

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
    ~ p_fac,
    scales   = "free_y",
    labeller = as_labeller(p_label_map, label_parsed),
    nrow = 4,
    ncol = 2
  ) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.spacing    = unit(0, "lines"),
    legend.position      = "top",
    legend.justification = "left", 
    legend.box           = "horizontal",
    legend.spacing.x     = unit(0.75, "cm"),
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
    color    = "none"
  )

plot_WFE

##================================================================
## computing critical time and maximum weighted fractional entropy
##================================================================

H_vals <- c(0.5,0.6,0.75,0.8,0.9,0.95)
p_vals <- seq(0,5,by=0.1)
grid <- expand.grid(H = H_vals, p = p_vals, KEEP.OUT.ATTRS = FALSE)

t_star_table <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  h <- grid$H[i]
  p <- grid$p[i]
  cat(sprintf("H = %.2f | p = %s\n", h, p))
  ## Objectif
  obj <- function(t) {
    val <- tryCatch(
      WFE_theor(
        x0 = x0, gt = gt, k = k,
        Alpha = Alpha, Beta = Beta, Gamma = Gamma,
        H = h, p = p, t = t,
        methods = "cor1"
      )[["cor1"]],
      error = function(e) -Inf
    )
    if (is.null(val) || !is.finite(val)) -Inf else val
  }
  ## Optimisation on (0, T]
  opt <- tryCatch(
    optimize(obj, interval = c(1e-3, T), maximum = TRUE),
    error = function(e) list(maximum = NA_real_, objective = NA_real_)
  )
  data.frame(H       = h,
             p       = p,
             t_star  = round(opt$maximum,  6),
             WFE_max = round(opt$objective, 6))
  
})) %>% dplyr::arrange(H, p)

print(t_star_table, digits = 5)

##================================================================
## Plot Critical time
##================================================================

H_unique  <- sort(unique(t_star_table$H))
Color_H   <- ggsci::pal_d3("category20", alpha = 0.5)(length(H_unique))
names(Color_H) <- as.character(H_unique)
H_labels  <- paste0("H = ", H_unique)
p_breaks  <- sort(unique(t_star_table$p))
p_labels  <- sapply(p_breaks, function(v) as.character(MASS::fractions(v)))

plot_t_stars <- ggplot(t_star_table, aes(x     = p,
                                         y     = t_star,
                                         color = factor(H),
                                         group = factor(H))) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(size = 1) +
  scale_color_manual(name   = NULL,
                     values = Color_H,
                     labels = H_labels) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 12),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 16),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(x = expression(bold("Weighting order"~italic(p))),
       y = expression(bold("Critical time"~italic(t)[italic(p)]^"*"))) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.spacing    = unit(0, "lines"),
    legend.position      = "top",
    legend.justification = "left", 
    legend.box           = "horizontal",
    legend.spacing.x     = unit(0.75, "cm"),
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
    color = guide_legend(
      nrow         = 1,
      override.aes = list(linewidth = 2, linetype = "solid")
    )
  )

##================================================================
## Plot maximum weighted fractional entropy
##================================================================

plot_WFE_max <- ggplot(t_star_table, aes(x     = p,
                         y     = WFE_max,
                         color = factor(H),
                         group = factor(H))) +
  geom_line(linewidth = 1, linetype = "solid") +
  geom_point(size = 1) +
  scale_color_manual(name   = NULL,
                     values = Color_H,
                     labels = H_labels) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 12),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    trans   = "log10",
    breaks  = trans_breaks("log10", function(x) 10^x, n = 6),
    labels  = trans_format("log10", math_format(10^.x))
  )+
  labs(
    x = expression(bold("Weighting order" ~ italic(p))),
    y = expression(
      max[t %in% group("[", list(0, T), "]")] ~
        bold(H)[italic(w)]^{"(p)"} * (X[t])
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid       = element_blank(),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
    legend.position      = "top",
    legend.justification = "centre",
    legend.box           = "horizontal",
    legend.spacing.x     = unit(0.75, "cm"),
    legend.key.width     = unit(1.8, "cm"),
    legend.text          = element_text(size = 14, face = "bold.italic"),
    legend.margin        = margin(t = 0, r = 0, b = -10, l = 0),
    legend.box.spacing   = unit(0.5, "cm"),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text        = element_text(size = 10),
    plot.margin      = unit(c(0, 0.05, 0, 0), "cm")
  ) +
  guides(
    color = guide_legend(
      nrow         = 1,
      override.aes = list(linewidth = 2, linetype = "solid")
    )
  )

##================================================================
## Plots Critical time and maximum weighted fractional entropy
## with cowplot
##================================================================

legend_H <- get_legend(
  plot_t_stars +
    theme(legend.position      = "top",
          legend.justification = "center")
)
p1 <- plot_t_stars + theme(legend.position = "none",axis.title.y= element_text(face = "bold", margin = margin(r = 0)))
p2 <- plot_WFE_max + theme(legend.position = "none",axis.title.y= element_text(face = "bold", margin = margin(r = 0)))
plots_row <- plot_grid(p1, p2, nrow = 1, align = "v", axis = "tb")

plot_grid(
  legend_H,
  plots_row,
  ncol        = 1,
  rel_heights = c(0.05, 1)
)
