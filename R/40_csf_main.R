#!/usr/bin/env Rscript
# =============================================================================
# 40_csf_main.R — Causal Survival Forest メイン解析
# grf::causal_survival_forest + RMST
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(grf)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(survival)
})

# setup ロード（here() でプロジェクトルートを参照）
project_root <- here::here()
setwd(project_root)

source("R/00_setup_survival.R")
source("R/utils/cf_helpers.R")   # get_gate, cochran_q_test, get_rate_df, get_vi 等を流用
source("R/utils/var_labels.R")  # 図の変数ラベルを読者向け表示名に置換

set.seed(CSF_SEED)

# =============================================================================
# 0. データ読み込み
# =============================================================================
cat("\n=== Loading data ===\n")
df <- read_csv(SURV_CSV, show_col_types = FALSE)
cat(sprintf("  n=%d, W=1(PTH)=%d, W=0(BP)=%d\n",
            nrow(df), sum(df$W == 1), sum(df$W == 0)))

# 施設クラスタリング
df$clusters_int <- as.integer(factor(df$SITE_CODE))

# =============================================================================
# 1. W.hat 推定（回帰フォレスト、cross-fit）
# =============================================================================
cat("\n=== Estimating W.hat via regression_forest ===\n")

# 共変量行列（W.hat 推定用: TRIAL_ID 含む）
covar_valid <- intersect(COVARIATES_INTEGRATED, names(df))
X_all <- as.matrix(df[, covar_valid])
W_all <- as.numeric(df$W)
cl_all <- df$clusters_int

rf_w <- regression_forest(
  X_all, W_all,
  clusters     = cl_all,
  num.trees    = 2000,
  seed         = CSF_SEED
)
W_hat_all <- predict(rf_w)$predictions
W_hat_all <- pmax(pmin(W_hat_all, 0.95), 0.05)  # 極端値をクリップ

cat(sprintf("  W.hat: min=%.3f, mean=%.3f, max=%.3f\n",
            min(W_hat_all), mean(W_hat_all), max(W_hat_all)))
cat(sprintf("  J04 W.hat: mean=%.3f; J05 W.hat: mean=%.3f\n",
            mean(W_hat_all[df$TRIAL_ID == 4]),
            mean(W_hat_all[df$TRIAL_ID == 5])))

# positivity 診断
n_low  <- sum(W_hat_all < 0.05)
n_high <- sum(W_hat_all > 0.95)
cat(sprintf("  Positivity: n(W.hat<0.05)=%d, n(W.hat>0.95)=%d\n", n_low, n_high))

df$W_hat <- W_hat_all

# W.hat 分布を保存
write_csv(
  data.frame(TRIAL_ID = df$TRIAL_ID, W = df$W, W_hat = df$W_hat),
  table_path_surv("ALL", "what_distribution")
)

# =============================================================================
# 2. アウトカムループ
# =============================================================================
results <- list()

for (out in SURV_OUTCOMES) {
  code      <- out$code
  time_col  <- out$time_col
  event_col <- out$event_col
  label     <- out$label

  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("Outcome: %s — %s\n", code, label))
  cat(strrep("=", 60), "\n")

  # ---- 完全観測のみ（time & event 両方 non-NA）----
  keep <- !is.na(df[[time_col]]) & !is.na(df[[event_col]]) & df[[time_col]] > 0
  df_sub   <- df[keep, ]
  n_sub    <- nrow(df_sub)
  n_event  <- sum(df_sub[[event_col]] == 1)
  cat(sprintf("  Valid n=%d, events=%d (%.1f%%)\n",
              n_sub, n_event, n_event / n_sub * 100))

  if (n_event < 20) {
    cat("  Skipped: too few events\n")
    next
  }

  X       <- as.matrix(df_sub[, covar_valid])
  Y_time  <- as.numeric(df_sub[[time_col]])
  Y_event <- as.numeric(df_sub[[event_col]])
  W_sub   <- as.numeric(df_sub$W)
  W_hat_sub <- as.numeric(df_sub$W_hat)
  cl_sub  <- as.integer(factor(df_sub$SITE_CODE))

  # ---- 2a. CSF フィット ----
  cat(sprintf("  Fitting CSF (target=%s, horizon=%.1f, trees=%d)...\n",
              CSF_TARGET, HORIZON, CSF_NUM_TREES))
  t_start <- proc.time()
  csf <- causal_survival_forest(
    X             = X,
    Y             = Y_time,
    W             = W_sub,
    D             = Y_event,
    W.hat         = W_hat_sub,
    target        = CSF_TARGET,
    horizon       = HORIZON,
    num.trees     = CSF_NUM_TREES,
    min.node.size = CSF_MIN_NODE,
    honesty       = TRUE,
    seed          = CSF_SEED,
    clusters      = cl_sub
  )
  elapsed <- (proc.time() - t_start)["elapsed"]
  cat(sprintf("  Fit time: %.1f sec\n", elapsed))

  # ---- 2b. CATE（OOB）----
  pred     <- predict(csf, estimate.variance = TRUE)
  tau_hat  <- pred$predictions          # 正値 = PTH benefit（RMST増加）
  tau_se   <- sqrt(pred$variance.estimates)

  # ---- 2c. ATE ----
  ate_obj  <- average_treatment_effect(csf)
  ate_est  <- ate_obj["estimate"]
  ate_se   <- ate_obj["std.err"]
  ate_lo   <- ate_est - 1.96 * ate_se
  ate_hi   <- ate_est + 1.96 * ate_se
  ate_p    <- 2 * pnorm(-abs(ate_est / ate_se))
  cat(sprintf("  ATE(RMST diff) = %.4f yr [%.4f, %.4f], p=%.3f\n",
              ate_est, ate_lo, ate_hi, ate_p))

  # ---- 2d. GATE（tau_hat で四分位）----
  # ★ CSF: benefit = 正のCATE → get_gate は tau_hat をそのまま渡す
  gate_df  <- get_gate(csf, tau_hat, probs = GATE_PROBS)
  q_test   <- cochran_q_test(gate_df)
  cat(sprintf("  Cochran Q=%.3f (df=%d, p=%.3f)\n",
              q_test$Q, q_test$df, q_test$p_value))

  # ---- 2e. RATE / AUTOC ----
  # ★ ranking_score = +tau_hat（binary CF とは逆符号）
  rate_df  <- get_rate_df(csf, tau_hat, ranking_score = tau_hat)
  cat(sprintf("  AUTOC=%.4f (p=%.3f)\n",
              rate_df$estimate[1], rate_df$p_value[1]))

  # ---- 2f. BLP ----
  blp_raw  <- best_linear_projection(csf, A = X)
  blp_df   <- data.frame(
    variable = rownames(blp_raw),
    estimate = blp_raw[, "Estimate"],
    std_err  = blp_raw[, "Std. Error"],
    t_stat   = blp_raw[, "t value"],
    p_value  = blp_raw[, "Pr(>|t|)"]
  )
  blp_df <- blp_df[order(blp_df$p_value, na.last = TRUE), ]

  # ---- 2g. Variable Importance ----
  vi_df    <- get_vi(csf, covar_valid)

  # ---- 2h. 保存 ----
  ate_out <- data.frame(
    outcome = code, label = label,
    estimate = ate_est, std_err = ate_se,
    ci_lo = ate_lo, ci_hi = ate_hi, p_value = ate_p,
    n = n_sub, n_event = n_event
  )
  cate_out <- data.frame(
    ID = df_sub$ID, W = W_sub, TRIAL_ID = df_sub$TRIAL_ID,
    time = Y_time, event = Y_event,
    tau_hat = tau_hat, tau_se = tau_se,
    tau_lo = tau_hat - 1.96 * tau_se,
    tau_hi = tau_hat + 1.96 * tau_se
  )

  write_csv(ate_out,              table_path_surv(code, "ate"))
  write_csv(cate_out,             table_path_surv(code, "cate"))
  write_csv(gate_df,              table_path_surv(code, "gate"))
  write_csv(rate_df,              table_path_surv(code, "rate"))
  write_csv(blp_df,               table_path_surv(code, "blp"))
  write_csv(vi_df,                table_path_surv(code, "vi"))
  write_csv(data.frame(Q = q_test$Q, df = q_test$df, p_value = q_test$p_value),
            table_path_surv(code, "cochranq"))

  # ---- 2i. 図: CATE 分布 ----
  p_cate <- ggplot(data.frame(tau = tau_hat, W = factor(W_sub, labels = c("BP","PTH"))),
                   aes(x = tau, fill = W)) +
    geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(title = paste0(code, ": CATE distribution (RMST diff, yr)"),
         x = "CATE (RMST difference, yr)", y = "Count") +
    theme_bw()
  ggsave(figure_path_surv(code, "cate_distribution"),
         p_cate, width = 7, height = 4, dpi = 150)

  # ---- 2j. 図: GATE ----
  gate_df$group_f <- factor(gate_df$group, levels = rev(gate_df$group))
  p_gate <- ggplot(gate_df, aes(x = group_f, y = estimate,
                                ymin = ci_lo, ymax = ci_hi)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(color = "steelblue", size = 0.8) +
    coord_flip() +
    labs(title = paste0(code, ": GATE (RMST diff by CATE quartile)"),
         x = "CATE Quartile", y = "GATE: RMST difference (yr)") +
    theme_bw()
  ggsave(figure_path_surv(code, "gate"),
         p_gate, width = 6, height = 4, dpi = 150)

  # ---- 2k. 図: BLP forest ----
  blp_plot <- head(blp_df, 15)
  # 生のデータ項目名を読者向け表示ラベルに置換（manuscript と表記統一）
  blp_plot$variable <- relabel_vars(blp_plot$variable)
  blp_plot$variable <- factor(blp_plot$variable,
                               levels = rev(blp_plot$variable))
  p_blp <- ggplot(blp_plot, aes(x = variable, y = estimate,
                                 ymin = estimate - 1.96 * std_err,
                                 ymax = estimate + 1.96 * std_err,
                                 color = p_value < 0.05)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(size = 0.7) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "tomato"),
                       name = "p<0.05") +
    labs(title = paste0(code, ": BLP — treatment effect modifiers"),
         x = NULL, y = "BLP coefficient (β)") +
    theme_bw()
  ggsave(figure_path_surv(code, "blp_forest"),
         p_blp, width = 8, height = 5, dpi = 150)

  # ---- 結果保存 ----
  results[[code]] <- list(
    code = code, csf = csf, df_sub = df_sub,
    tau_hat = tau_hat, tau_se = tau_se,
    ate = ate_out, gate = gate_df, q_test = q_test,
    rate = rate_df, blp = blp_df, vi = vi_df,
    n = n_sub, n_event = n_event
  )

  cat(sprintf("  Saved tables/figures for %s\n", code))
}

# =============================================================================
# 3. サマリーテーブル出力
# =============================================================================
cat("\n=== Summary ===\n")
summary_rows <- lapply(results, function(r) {
  data.frame(
    outcome   = r$code,
    label     = r$ate$label,
    n         = r$n,
    n_event   = r$n_event,
    ATE       = round(r$ate$estimate, 4),
    ATE_lo    = round(r$ate$ci_lo, 4),
    ATE_hi    = round(r$ate$ci_hi, 4),
    ATE_p     = round(r$ate$p_value, 3),
    AUTOC     = round(r$rate$estimate[1], 4),
    AUTOC_p   = round(r$rate$p_value[1], 3),
    CochranQ  = round(r$q_test$Q, 3),
    CochranQ_p = round(r$q_test$p_value, 3),
    HTE_sig   = ifelse(
      (r$rate$p_value[1] < 0.05) + (r$q_test$p_value < 0.05) >= 2,
      "✅ HTE", "partial/none"
    )
  )
})
summary_df <- do.call(rbind, summary_rows)
write_csv(summary_df, file.path(OUT_DIR, "CSF_summary.csv"))
print(summary_df)

cat("\n=== 40_csf_main.R 完了 ===\n")
cat(sprintf("出力先: %s\n", OUT_DIR))
