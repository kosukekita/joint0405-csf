#!/usr/bin/env Rscript
# =============================================================================
# 41_csf_sensitivity.R — CSF 感度分析
# 1. 骨折歴除外モデル（COVARIATES_INTEGRATED_NOHX）
# 2. horizon 変動（1.5 / 2.0 / 2.5 年）
# 3. J05 単独 CSF（真の RCT 内 HTE）
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

project_root <- here::here()
setwd(project_root)

source("R/00_setup_survival.R")
source("R/utils/cf_helpers.R")

set.seed(CSF_SEED)

cat("\n=== Loading data ===\n")
df <- read_csv(SURV_CSV, show_col_types = FALSE)
df$clusters_int <- as.integer(factor(df$SITE_CODE))
cat(sprintf("  n=%d, W=1(PTH)=%d, W=0(BP)=%d\n",
            nrow(df), sum(df$W == 1), sum(df$W == 0)))

# =============================================================================
# 共通: W.hat 推定（COVARIATES_INTEGRATED）
# =============================================================================
cat("\n=== Estimating W.hat ===\n")
covar_full <- intersect(COVARIATES_INTEGRATED, names(df))
X_all      <- as.matrix(df[, covar_full])
W_all      <- as.numeric(df$W)
cl_all     <- df$clusters_int

rf_w <- regression_forest(X_all, W_all, clusters = cl_all,
                          num.trees = 2000, seed = CSF_SEED)
W_hat_all <- pmax(pmin(predict(rf_w)$predictions, 0.95), 0.05)
df$W_hat  <- W_hat_all
cat(sprintf("  W.hat: mean=%.3f (J04=%.3f, J05=%.3f)\n",
            mean(W_hat_all),
            mean(W_hat_all[df$TRIAL_ID == 4]),
            mean(W_hat_all[df$TRIAL_ID == 5])))

# 結果格納
sens_results <- list()

# =============================================================================
# ヘルパー関数: 1つの CSF フィットと評価を実行
# =============================================================================
run_one_csf <- function(df_sub, covar_names, horizon_val, label, tag) {
  keep <- !is.na(df_sub$T_VF) & !is.na(df_sub$D_VF) & df_sub$T_VF > 0
  df_s <- df_sub[keep, ]
  n_evt <- sum(df_s$D_VF == 1)
  cat(sprintf("  [%s] n=%d, events=%d, horizon=%.1f\n",
              tag, nrow(df_s), n_evt, horizon_val))
  if (n_evt < 20) {
    cat("  Skipped: too few events\n")
    return(NULL)
  }

  covar_valid <- intersect(covar_names, names(df_s))
  X      <- as.matrix(df_s[, covar_valid])
  Y_time <- as.numeric(df_s$T_VF)
  Y_evt  <- as.numeric(df_s$D_VF)
  W_sub  <- as.numeric(df_s$W)
  W_hat  <- as.numeric(df_s$W_hat)
  cl_sub <- as.integer(factor(df_s$SITE_CODE))

  csf <- causal_survival_forest(
    X = X, Y = Y_time, W = W_sub, D = Y_evt,
    W.hat = W_hat, target = CSF_TARGET, horizon = horizon_val,
    num.trees = CSF_NUM_TREES, min.node.size = CSF_MIN_NODE,
    honesty = TRUE, seed = CSF_SEED, clusters = cl_sub
  )

  tau_hat <- predict(csf)$predictions
  ate_obj <- average_treatment_effect(csf)
  ate_est <- ate_obj["estimate"]; ate_se <- ate_obj["std.err"]
  ate_p   <- 2 * pnorm(-abs(ate_est / ate_se))

  gate_df <- get_gate(csf, tau_hat, probs = GATE_PROBS)
  q_test  <- cochran_q_test(gate_df)
  rate_df <- get_rate_df(csf, tau_hat, ranking_score = tau_hat)

  cat(sprintf("  ATE=%.4f (p=%.3f), AUTOC=%.4f (p=%.3f), CochranQ=%.3f (p=%.3f)\n",
              ate_est, ate_p,
              rate_df$estimate[1], rate_df$p_value[1],
              q_test$Q, q_test$p_value))

  list(
    tag = tag, label = label, horizon = horizon_val,
    n = nrow(df_s), n_event = n_evt,
    ATE = ate_est, ATE_se = ate_se, ATE_p = ate_p,
    ATE_lo = ate_est - 1.96 * ate_se, ATE_hi = ate_est + 1.96 * ate_se,
    AUTOC = rate_df$estimate[1], AUTOC_p = rate_df$p_value[1],
    CochranQ = q_test$Q, CochranQ_p = q_test$p_value,
    HTE_sig = ifelse(
      (rate_df$p_value[1] < 0.05) + (q_test$p_value < 0.05) >= 2,
      "HTE", "partial/none"
    ),
    gate = gate_df, rate = rate_df
  )
}

# =============================================================================
# 感度分析 1: 骨折歴除外（NOHX）— E1, E3 × horizon=2.0
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("Sensitivity 1: No fracture history covariates (NOHX)\n")
cat(strrep("=", 60), "\n")

for (outcome in list(
  list(time = "T_VF", event = "D_VF", lab = "E1-NOHX"),
  list(time = "T_E3", event = "D_E3", lab = "E3-NOHX")
)) {
  df_tmp <- df
  df_tmp$T_VF <- df[[outcome$time]]
  df_tmp$D_VF <- df[[outcome$event]]
  res <- run_one_csf(df_tmp, COVARIATES_INTEGRATED_NOHX,
                     HORIZON, outcome$lab, outcome$lab)
  if (!is.null(res)) sens_results[[outcome$lab]] <- res
}

# =============================================================================
# 感度分析 2: horizon 変動（E1 のみ）
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("Sensitivity 2: Horizon variation (E1)\n")
cat(strrep("=", 60), "\n")

for (h in c(1.5, 2.5)) {
  tag <- sprintf("E1-horizon%.1f", h)
  df_tmp <- df
  df_tmp$T_VF <- df$T_VF
  df_tmp$D_VF <- df$D_VF
  res <- run_one_csf(df_tmp, COVARIATES_INTEGRATED, h,
                     sprintf("E1 horizon=%.1f", h), tag)
  if (!is.null(res)) sens_results[[tag]] <- res
}

# =============================================================================
# 感度分析 3: J05 単独（真の RCT 内 HTE）
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("Sensitivity 3: J05 only (true RCT, W.hat=0.5 fixed)\n")
cat(strrep("=", 60), "\n")

df_j05 <- df[df$TRIAL_ID == 5, ]
cat(sprintf("  J05: n=%d, W=1=%d, W=0=%d\n",
            nrow(df_j05), sum(df_j05$W == 1), sum(df_j05$W == 0)))

# J05 単独では W.hat=0.5 固定（真のRCT）
df_j05$W_hat <- 0.5

for (outcome in list(
  list(time = "T_VF", event = "D_VF", lab = "E1-J05only"),
  list(time = "T_E3", event = "D_E3", lab = "E3-J05only")
)) {
  df_tmp <- df_j05
  df_tmp$T_VF <- df_j05[[outcome$time]]
  df_tmp$D_VF <- df_j05[[outcome$event]]
  res <- run_one_csf(df_tmp, COVARIATES_INTEGRATED,
                     HORIZON, outcome$lab, outcome$lab)
  if (!is.null(res)) sens_results[[outcome$lab]] <- res
}

# =============================================================================
# サマリーテーブル保存
# =============================================================================
cat("\n=== Sensitivity Analysis Summary ===\n")

summary_rows <- lapply(sens_results, function(r) {
  data.frame(
    tag       = r$tag,
    label     = r$label,
    horizon   = r$horizon,
    n         = r$n,
    n_event   = r$n_event,
    ATE       = round(r$ATE,  4),
    ATE_lo    = round(r$ATE_lo, 4),
    ATE_hi    = round(r$ATE_hi, 4),
    ATE_p     = round(r$ATE_p, 3),
    AUTOC     = round(r$AUTOC, 4),
    AUTOC_p   = round(r$AUTOC_p, 3),
    CochranQ  = round(r$CochranQ, 3),
    CochranQ_p = round(r$CochranQ_p, 3),
    HTE_sig   = r$HTE_sig,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
write_csv(summary_df, file.path(TABLE_DIR, "CSF_sensitivity_summary.csv"))
print(summary_df)

# メイン解析結果（参照用）との比較表も保存
main_ref <- data.frame(
  tag = c("E1-main", "E3-main"),
  label = c("E1 main (horizon=2.0, full covariates)",
            "E3 main (horizon=2.0, full covariates)"),
  horizon = 2.0,
  n = c(2586, 2608), n_event = c(182, 210),
  ATE = c(-0.0080, -0.0042),
  ATE_lo = c(-0.0222, -0.0202), ATE_hi = c(0.0062, 0.0118),
  ATE_p = c(0.271, 0.604),
  AUTOC = c(-0.0208, -0.0217),
  AUTOC_p = c(0.025, 0.012),
  CochranQ = c(5.185, 8.015),
  CochranQ_p = c(0.159, 0.046),
  HTE_sig = c("partial/none", "HTE"),
  stringsAsFactors = FALSE
)
write_csv(
  bind_rows(main_ref, summary_df),
  file.path(TABLE_DIR, "CSF_sensitivity_vs_main.csv")
)

cat("\n=== 41_csf_sensitivity.R 完了 ===\n")
cat(sprintf("出力先: %s\n", TABLE_DIR))
