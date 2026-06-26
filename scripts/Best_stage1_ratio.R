# =============================================================================
# STAGE 1 RATIO SENSITIVITY — find best S1 ratio for Stage 2 trend balance
# =============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("STAGE 1 RATIO SENSITIVITY CHECK\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

s1_ratios_to_test <- c(20, 30, 50, 100)

s1_sensitivity <- map_df(s1_ratios_to_test, function(r1) {
  cat(sprintf("\n--- Testing S1 ratio: 1:%d ---\n", r1))
  
  scheme_results <- map_df(english_schemes, function(s) {
    scheme_treated <- data_england %>%
      filter(treat_indicator == 1, scheme == s)
    scheme_data <- bind_rows(scheme_treated, control_pool)
    
    scheme_clean <- winsorise_and_log_s1(scheme_data, stage1_vars,
                                         log_transform_s1, log_nozero_s1)
    scheme_clean <- winsorise_and_log_s2(scheme_clean, stage2_levels,
                                         log_transform_s2_levels)
    
    s1v <- check_vars(scheme_clean, stage1_vars_log, paste("S1 sens", s))
    s2v <- check_vars(scheme_clean, stage2_vars_log, paste("S2 sens", s))
    
    # --- Stage 1 at ratio r1 ---
    formula_s1 <- reformulate(s1v, response = "treat_indicator")
    m_s1 <- tryCatch(
      matchit(formula_s1, data = scheme_clean, method = "nearest",
              distance = "mahalanobis", ratio = r1, replace = TRUE),
      error = function(e) NULL
    )
    if (is.null(m_s1)) {
      cat(sprintf("  S1 FAILED for scheme %s at ratio %d\n", s, r1))
      return(tibble(scheme = s, s1_ratio = r1,
                    n_s1_controls = NA_integer_,
                    best_s2_ratio = NA_integer_,
                    best_trend_smd = NA_real_,
                    best_mean_smd  = NA_real_))
    }
    
    mm_s1       <- m_s1$match.matrix
    treated_s1  <- scheme_clean[as.integer(rownames(mm_s1)), , drop = FALSE] %>%
      mutate(treat_indicator = 1L)
    ctrl_idx    <- unique(as.integer(mm_s1[!is.na(mm_s1)]))
    controls_s1 <- scheme_clean[ctrl_idx, , drop = FALSE] %>%
      mutate(treat_indicator = 0L)
    n_s1_controls <- nrow(controls_s1)
    
    # --- Prep Stage 2 data ---
    s2_raw      <- bind_rows(treated_s1, controls_s1) %>%
      select(-any_of(c("weights", "subclass", "distance")))
    treated_ref <- s2_raw %>% filter(treat_indicator == 1)
    
    for (v in intersect(stage2_trends, names(s2_raw))) {
      q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
      q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
      s2_raw[[v]] <- pmin(pmax(s2_raw[[v]], q_lo), q_hi)
    }
    for (v in intersect(stage2_levels, names(s2_raw))) {
      q_lo <- quantile(treated_ref[[v]], 0.01, na.rm = TRUE)
      q_hi <- quantile(treated_ref[[v]], 0.99, na.rm = TRUE)
      s2_raw[[paste0("log1p_", v)]] <-
        log1p(pmax(pmin(pmax(s2_raw[[v]], q_lo), q_hi), 0))
    }
    
    s2v_clean   <- check_vars(s2_raw, s2v, paste("S2 sens", s))
    formula_s2  <- reformulate(s2v_clean, response = "treat_indicator")
    
    # --- Stage 2: test ratios 1:10, pick best total_trend_smd ---
    n_t <- sum(s2_raw$treat_indicator == 1)
    n_c <- sum(s2_raw$treat_indicator == 0)
    s2_ratios <- 1:min(10, floor(n_c / n_t))
    
    trend_var <- intersect(s2v_clean, stage2_trends)[1]
    
    ratio_scan <- map_df(s2_ratios, function(r2) {
      m <- tryCatch(
        matchit(formula_s2, data = s2_raw, method = "nearest",
                distance = "mahalanobis", ratio = r2, replace = TRUE),
        error = function(e) NULL
      )
      if (is.null(m)) return(tibble(r2 = r2, trend_smd = NA_real_,
                                    mean_smd  = NA_real_))
      bt <- bal.tab(m, un = FALSE, stats = "mean.diffs")$Balance
      tibble(
        r2        = r2,
        trend_smd = if (trend_var %in% rownames(bt))
          abs(bt[trend_var, "Diff.Adj"]) else NA_real_,
        mean_smd  = mean(abs(bt$Diff.Adj), na.rm = TRUE)
      )
    })
    
    best <- ratio_scan %>%
      filter(!is.na(trend_smd)) %>%
      arrange(trend_smd, mean_smd) %>%
      slice(1)
    
    cat(sprintf("  Scheme %-20s | S1 pool: %4d controls | best S2 ratio: 1:%d | trend SMD: %.4f\n",
                s, n_s1_controls, best$r2, best$trend_smd))
    
    tibble(
      scheme         = s,
      s1_ratio       = r1,
      n_s1_controls  = n_s1_controls,
      best_s2_ratio  = best$r2,
      best_trend_smd = best$trend_smd,
      best_mean_smd  = best$mean_smd
    )
  })
  
  scheme_results
})

# =============================================================================
# SUMMARISE + PLOT SENSITIVITY RESULTS
# =============================================================================

cat("\n=== STAGE 1 RATIO SENSITIVITY SUMMARY ===\n\n")

# Per-scheme view
s1_sens_wide <- s1_sensitivity %>%
  select(scheme, s1_ratio, best_trend_smd) %>%
  pivot_wider(names_from = s1_ratio,
              values_from = best_trend_smd,
              names_prefix = "S1_")

print(s1_sens_wide)

# Overall: mean trend SMD across schemes per S1 ratio
s1_sens_summary <- s1_sensitivity %>%
  group_by(s1_ratio) %>%
  summarise(
    mean_trend_smd = round(mean(best_trend_smd, na.rm = TRUE), 4),
    max_trend_smd  = round(max(best_trend_smd,  na.rm = TRUE), 4),
    mean_pool_size = round(mean(n_s1_controls,   na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  arrange(mean_trend_smd)

cat("\nAggregated across schemes:\n")
print(s1_sens_summary)

best_s1_ratio <- s1_sens_summary %>%
  slice(1) %>%
  pull(s1_ratio)

cat(sprintf("\n>>> Recommended S1 ratio: 1:%d (lowest mean trend SMD = %.4f across schemes)\n",
            best_s1_ratio,
            s1_sens_summary %>% filter(s1_ratio == best_s1_ratio) %>% pull(mean_trend_smd)))

# Plot
p_s1_sens <- s1_sensitivity %>%
  filter(!is.na(best_trend_smd)) %>%
  ggplot(aes(x = s1_ratio, y = best_trend_smd,
             colour = scheme, group = scheme)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.10, linetype = "dashed", colour = "#888888") +
  geom_hline(yintercept = 0.05, linetype = "dotted", colour = "#555555") +
  scale_x_continuous(breaks = s1_ratios_to_test) +
  labs(
    title    = "Stage 1 ratio sensitivity — effect on Stage 2 trend balance",
    subtitle = "Best achievable total trend |SMD| (across S2 ratios 1–10) by S1 pool size",
    x        = "Stage 1 matching ratio (1:k)",
    y        = "Best Stage 2 total trend |SMD|",
    colour   = "Scheme",
    caption  = "Dashed = 0.10 threshold; dotted = 0.05 threshold"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "fig_s1_ratio_sensitivity.png"),
       p_s1_sens, width = 13, height = 7, dpi = 300, bg = "white")

cat("\nSaved: fig_s1_ratio_sensitivity.png\n")

# Save table
saveRDS(s1_sensitivity,
        here("data", "processed", "OA_s1_ratio_sensitivity.rds"))
write_csv(s1_sensitivity,
          here("output", "diagnostics", "pooled", "OA_s1_ratio_sensitivity.csv"))