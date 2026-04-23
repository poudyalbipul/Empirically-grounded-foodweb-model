library(readxl)
library(tidyverse)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
alldata <- read_excel("parameter_sensitivity.xlsx")

# ── 2. Define column groups and row slices ────────────────────────────────────
alpha_cols <- c("beta_f","beta_sigma", "beta_theta",
                "beta_ai","beta_ak", "beta_x")
row_slices <- list(1:11, 12:22, 23:33, 34:44)

# Build parameter list: alpha sets only
parameter_list <- lapply(row_slices, \(r) alldata[r, c("Normalized_x", alpha_cols)])

# ── 3. Normalize: subtract row-6 mean, divide by per-column SD ───────────────
normalize_alpha <- function(df) {
  mu <- mean(unlist(df[6, alpha_cols]), na.rm = TRUE)
  sigma <- sd(unlist(df[6, alpha_cols]), na.rm = TRUE)
  df[alpha_cols] <- lapply(df[alpha_cols], \(x) (x- mu)/ sigma)
  df
}

parameter_list <- lapply(parameter_list, normalize_alpha)

# ── 4. Reshape to long format ─────────────────────────────────────────────────
group_labels <- c("Species Richness", "Biomass", "Mean Trophic Level", "Stability")

combined_alpha <- map2_dfr(parameter_list, group_labels, \(df, g) {
  df |>
    pivot_longer(-Normalized_x, names_to = "Variable", values_to = "Value") |>
    mutate(
      group = factor(g, levels = group_labels),
      Value = if (g == "Stability") -Value else Value
    )
})

# ── 5. Colors, linetypes, shared theme ───────────────────────────────────────
alpha_colors <- c(
  "beta_ai" = "#E31A1C",
  "beta_ak" = "#1F78B4",
  "beta_f"  = "#33A02C",
  "beta_sigma"  = "#FF7F00",  # orange
  "beta_theta"  = "#6A3D9A",  # purple
  "beta_x"  = "#B15928"   # brown
)

linetype_scale <- scale_linetype_manual(
  values = c(
    "Species Richness"   = "dashed",
    "Biomass"            = "solid",
    "Mean Trophic Level" = "dotted",
    "Stability"          = "dotdash"
  )
)

base_theme <- list(
  geom_line(linewidth = 0.5),
  geom_hline(yintercept = c(-2.576, 2.576), linewidth = 0.8),
  geom_vline(xintercept = 0,              linewidth = 0.8),
  coord_cartesian(ylim = c(-35, 35)),
  linetype_scale,
  theme_bw(),
  theme(
    legend.background = element_rect(fill = NA, color = NA),
    legend.box        = "vertical",
    legend.text       = element_text(size = 14),
    legend.title      = element_text(size = 14),
    axis.title        = element_text(size = 14)
  )
)

# ── 6. Plot ───────────────────────────────────────────────────────────────────
plot_alpha <- ggplot(combined_alpha, aes(Normalized_x, Value, color = Variable, linetype = group)) +
  base_theme +
  scale_color_manual(
    values = alpha_colors,
    labels = c("beta_f" = expression(beta[f]), "beta_ai" = expression(beta[ai]), 
               "beta_ak" = expression(beta[ak]),
               "beta_sigma" = expression(beta[sigma]), "beta_theta" = expression(beta[theta]), 
               "beta_x" = expression(beta[x])
    )
  ) +
  labs(
    x       = "Change from reference",
    y       = "",
    title   = "B. All parameters (slopes)",
    color   = "",
    linetype = "Community metrics"
  ) +
  guides(
    color    = guide_legend(order = 1, ncol=2, override.aes = list(linewidth = 1)),
    linetype = "none"
  ) +
  theme(legend.position = c(0.65, 0.3))

plot_alpha
saveRDS(plot_alpha, "plot_betaslopes.rds")
