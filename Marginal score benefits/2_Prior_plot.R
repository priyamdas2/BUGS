# ============================================================
# Conceptual 1x2 prior-mechanism figure for BUGS vs RHS
#
# Panel A: Guidance multiplier exp(eta * z_j^*) as a function of z_j^*
# Panel B: Effective prior variance \tilde{\kappa}_j^2 as a function of lambda_j^2
#
# No data are needed.
# ============================================================

rm(list = ls())
setwd("U:/BUGS/Marginal score benefits")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

# -----------------------------
# User-chosen hyperparameters
# -----------------------------
eta  <- 1
tau2 <- 0.05
c2   <- 1.0

# Representative clipped-guidance levels for panel B
z_low  <- -1.5
z_mid  <-  0.0
z_high <-  1.5

# Range shown in panel A
Cz_plot <- 2

# -----------------------------
# Legend labels
# -----------------------------
label_bugs_A <- "BUGS (η = 1)"
label_rhs_A  <- "RHS (η = 0)"

label_low  <- "BUGS (η = 1, z* = -1.5)"
label_mid  <- "BUGS (η = 1, z* = 0) / RHS (η = 0)"
label_high <- "BUGS (η = 1, z* = 1.5)"

# -----------------------------
# Colors
# -----------------------------
col_rhs   <- "#1F4E79"
col_low   <- "#B22222"
col_mid   <- "#1F4E79"
col_high  <- "#B22222"
col_bugsA <- col_high

# -----------------------------
# Functions
# -----------------------------
w_fun <- function(z, eta) {
  exp(eta * z)
}

kappa_bugs <- function(lambda2, z, eta, tau2, c2) {
  (c2 * tau2 * lambda2 * exp(eta * z)) /
    (c2 + tau2 * lambda2 * exp(eta * z))
}

# -----------------------------
# Data for Panel A
# -----------------------------
df_left <- tibble(
  z = seq(-Cz_plot, Cz_plot, length.out = 500)
) %>%
  mutate(
    BUGS = w_fun(z, eta),
    RHS  = 1
  ) %>%
  pivot_longer(
    cols = c(BUGS, RHS),
    names_to = "method",
    values_to = "weight"
  ) %>%
  mutate(
    method = factor(method, levels = c("BUGS", "RHS"))
  )

# -----------------------------
# Data for Panel B
# -----------------------------
lambda_grid <- exp(seq(log(1e-3), log(1e3), length.out = 600))

df_right <- bind_rows(
  tibble(
    lambda2 = lambda_grid,
    value   = kappa_bugs(lambda_grid, z = z_low, eta = eta, tau2 = tau2, c2 = c2),
    curve   = label_low
  ),
  tibble(
    lambda2 = lambda_grid,
    value   = kappa_bugs(lambda_grid, z = z_mid, eta = eta, tau2 = tau2, c2 = c2),
    curve   = label_mid
  ),
  tibble(
    lambda2 = lambda_grid,
    value   = kappa_bugs(lambda_grid, z = z_high, eta = eta, tau2 = tau2, c2 = c2),
    curve   = label_high
  )
) %>%
  mutate(
    curve = factor(curve, levels = c(label_low, label_mid, label_high))
  )

# -----------------------------
# Plot theme
# -----------------------------
theme_prof <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.border = element_rect(linewidth = 0.6, colour = "black"),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.minor = element_blank(),
      
      # Larger equation-related labels/titles
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 11, colour = "black"),
      plot.title = element_text(face = "bold", size = 19, hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5, colour = "black"),
      
      legend.position = c(0.03, 0.97),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.background = element_rect(
        fill = alpha("white", 0.90),
        colour = "grey70",
        linewidth = 0.3
      ),
      legend.key = element_blank(),
      legend.key.height = unit(0.48, "cm"),
      legend.text = element_text(size = 10.8),
      legend.title = element_blank(),
      
      plot.margin = margin(t = 8, r = 10, b = 6, l = 8)
    )
}

# -----------------------------
# Panel A: Guidance multiplier
# -----------------------------
p_left <- ggplot(
  df_left,
  aes(x = z, y = weight, color = method, linetype = method)
) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "grey65") +
  geom_line(linewidth = 1.10) +
  scale_color_manual(
    values = c("BUGS" = col_bugsA, "RHS" = col_rhs),
    breaks = c("BUGS", "RHS"),
    labels = c("BUGS" = label_bugs_A, "RHS" = label_rhs_A)
  ) +
  scale_linetype_manual(
    values = c("BUGS" = "solid", "RHS" = "solid"),
    breaks = c("BUGS", "RHS"),
    labels = c("BUGS" = label_bugs_A, "RHS" = label_rhs_A)
  ) +
  scale_x_continuous(
    limits = c(-Cz_plot, Cz_plot),
    breaks = seq(-Cz_plot, Cz_plot, by = 1)
  ) +
  scale_y_continuous(
    breaks = c(0.2, 0.5, 1, 2, 5, 10),
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Guidance multiplier",
    subtitle = expression(exp(eta %.% tilde(z)[j]^"*")),
    x = expression(tilde(z)[j]^"*"),
    y = "Multiplier"
  ) +
  coord_cartesian(
    ylim = c(exp(-eta * Cz_plot), exp(eta * Cz_plot))
  ) +
  theme_prof()

# -----------------------------
# Panel B: Effective prior variance
# -----------------------------
p_right <- ggplot(
  df_right,
  aes(x = lambda2, y = value, color = curve, linetype = curve)
) +
  geom_line(linewidth = 1.10) +
  scale_x_log10(
    breaks = 10^(-3:3),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_color_manual(
    values = setNames(
      c(col_low, col_mid, col_high),
      c(label_low, label_mid, label_high)
    ),
    breaks = c(label_low, label_mid, label_high),
    labels = c(label_low, label_mid, label_high)
  ) +
  scale_linetype_manual(
    values = setNames(
      c("dotted", "solid", "dotdash"),
      c(label_low, label_mid, label_high)
    ),
    breaks = c(label_low, label_mid, label_high),
    labels = c(label_low, label_mid, label_high)
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    title = "Effective prior variance",
    subtitle = expression(tilde(kappa)[j]^2),
    x = expression(lambda[j]^2),
    y = expression(tilde(kappa)[j]^2)
  ) +
  theme_prof()

# -----------------------------
# Combine
# -----------------------------
p_prior_panel <- p_left + p_right +
  plot_layout(ncol = 2, widths = c(1, 1))

p_prior_panel

# -----------------------------
# Save
# -----------------------------
ggsave(
  "Results/figures/Figure_prior_mechanism_1x2.pdf",
  p_prior_panel,
  width = 11.6,
  height = 5.4,
  device = cairo_pdf
)

ggsave(
  "Results/figures/Figure_prior_mechanism_1x2.png",
  p_prior_panel,
  width = 11.6,
  height = 5.4,
  dpi = 700
)