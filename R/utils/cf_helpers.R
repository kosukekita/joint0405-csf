#!/usr/bin/env Rscript
# =============================================================================
# utils/cf_helpers.R — Causal Forest 共通ヘルパー関数
# =============================================================================

# -----------------------------------------------------------------------------
# データ準備: NA 除去、行列変換
# -----------------------------------------------------------------------------
prepare_data <- function(df, outcome_col, covariates) {
  # 列が存在するか確認
  missing_cols <- setdiff(covariates, names(df))
  if (length(missing_cols) > 0) {
    warning("Missing columns (will be excluded): ", paste(missing_cols, collapse = ", "))
    covariates <- intersect(covariates, names(df))
  }

  if (!outcome_col %in% names(df)) {
    stop("Outcome column not found: ", outcome_col)
  }

  # 解析対象行: アウトカム + 全共変量 + W が揃っている
  required_cols <- c(outcome_col, "W", covariates)
  complete_rows <- complete.cases(df[, required_cols])
  df_valid <- df[complete_rows, ]

  n_total  <- nrow(df)
  n_valid  <- nrow(df_valid)
  n_dropped <- n_total - n_valid
  cat(sprintf("  [%s] N total=%d, complete=%d (dropped=%d, %.1f%%)\n",
              outcome_col, n_total, n_valid, n_dropped, n_dropped / n_total * 100))

  X <- as.matrix(df_valid[, covariates])
  Y <- as.numeric(df_valid[[outcome_col]])
  W <- as.numeric(df_valid[["W"]])

  # 施設コード（クラスタリング用）
  clusters <- NULL
  if ("SITE_CODE" %in% names(df_valid)) {
    clusters <- as.integer(factor(df_valid[["SITE_CODE"]]))
  }

  list(
    X           = X,
    Y           = Y,
    W           = W,
    df_valid    = df_valid,
    covariates  = covariates,
    clusters    = clusters,
    n_valid     = n_valid
  )
}

# -----------------------------------------------------------------------------
# Causal Forest 実行（W.hat=0.5 固定、RCT 設計）
# -----------------------------------------------------------------------------
fit_cf <- function(X, Y, W, clusters = NULL,
                   num.trees  = CF_NUM_TREES,
                   min.node.size = CF_MIN_NODE,
                   seed       = CF_SEED,
                   use_clusters = TRUE) {
  cl_arg <- if (use_clusters && !is.null(clusters)) clusters else NULL

  cf <- causal_forest(
    X             = X,
    Y             = Y,
    W             = W,
    W.hat         = CF_W_HAT,        # 0.5 固定: RCT の既知割付確率
    tune.parameters = "all",
    num.trees     = num.trees,
    min.node.size = min.node.size,
    honesty       = TRUE,
    seed          = seed,
    clusters      = cl_arg
  )
  cf
}

# -----------------------------------------------------------------------------
# CATE 予測（point estimate + variance）
# -----------------------------------------------------------------------------
get_cate <- function(cf, X = NULL) {
  if (is.null(X)) {
    pred <- predict(cf, estimate.variance = TRUE)
  } else {
    pred <- predict(cf, newdata = X, estimate.variance = TRUE)
  }
  data.frame(
    cate    = pred$predictions,
    cate_var = pred$variance.estimates,
    cate_se  = sqrt(pred$variance.estimates),
    cate_lo  = pred$predictions - 1.96 * sqrt(pred$variance.estimates),
    cate_hi  = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
  )
}

# -----------------------------------------------------------------------------
# ATE 抽出
# -----------------------------------------------------------------------------
get_ate <- function(cf) {
  ate <- average_treatment_effect(cf, target.sample = "all")
  data.frame(
    estimate = ate["estimate"],
    std_err  = ate["std.err"],
    ci_lo    = ate["estimate"] - 1.96 * ate["std.err"],
    ci_hi    = ate["estimate"] + 1.96 * ate["std.err"],
    p_value  = 2 * pnorm(-abs(ate["estimate"] / ate["std.err"]))
  )
}

# -----------------------------------------------------------------------------
# GATE（四分位別 ATE）
# -----------------------------------------------------------------------------
get_gate <- function(cf, tau_hat, probs = GATE_PROBS) {
  breaks <- quantile(tau_hat, probs = probs, na.rm = TRUE)
  q_labels <- paste0("Q", seq_len(length(probs) - 1))
  groups   <- cut(tau_hat, breaks = breaks, labels = q_labels,
                  include.lowest = TRUE)

  gate_list <- lapply(q_labels, function(q) {
    idx <- which(groups == q)
    if (length(idx) < 5) return(NULL)
    ate_q <- average_treatment_effect(cf, subset = idx)
    data.frame(
      group    = q,
      n        = length(idx),
      estimate = ate_q["estimate"],
      std_err  = ate_q["std.err"],
      ci_lo    = ate_q["estimate"] - 1.96 * ate_q["std.err"],
      ci_hi    = ate_q["estimate"] + 1.96 * ate_q["std.err"]
    )
  })

  gate_df <- do.call(rbind, Filter(Negate(is.null), gate_list))
  rownames(gate_df) <- NULL
  gate_df
}

# -----------------------------------------------------------------------------
# Cochran's Q（GATE 間の異質性検定）
# -----------------------------------------------------------------------------
cochran_q_test <- function(gate_df) {
  est <- gate_df$estimate
  se  <- gate_df$std_err
  w   <- 1 / se^2
  w_mean <- sum(w * est) / sum(w)
  Q   <- sum(w * (est - w_mean)^2)
  df  <- nrow(gate_df) - 1
  p   <- pchisq(Q, df = df, lower.tail = FALSE)
  data.frame(Q = Q, df = df, p_value = p)
}

# -----------------------------------------------------------------------------
# RATE / AUTOC / QINI
# -----------------------------------------------------------------------------
get_rate <- function(cf, tau_hat, target = "AUTOC") {
  rate_obj <- rank_average_treatment_effect(cf, tau_hat, target = target)
  list(
    estimate = rate_obj$estimate,
    std_err  = rate_obj$std.err,
    ci_lo    = rate_obj$estimate - 1.96 * rate_obj$std.err,
    ci_hi    = rate_obj$estimate + 1.96 * rate_obj$std.err,
    p_value  = 2 * pnorm(-abs(rate_obj$estimate / rate_obj$std.err))
  )
}

get_rate_df <- function(cf, tau_hat, ranking_score = NULL) {
  # ranking_score: ランキングに使うスコア。
  # 骨折アウトカム(benefit_dir="neg_cate")では -tau_hat を渡すこと。
  if (is.null(ranking_score)) ranking_score <- tau_hat
  autoc <- get_rate(cf, ranking_score, "AUTOC")
  qini  <- get_rate(cf, ranking_score, "QINI")
  data.frame(
    metric   = c("AUTOC", "QINI"),
    estimate = c(autoc$estimate, qini$estimate),
    std_err  = c(autoc$std_err, qini$std_err),
    ci_lo    = c(autoc$ci_lo, qini$ci_lo),
    ci_hi    = c(autoc$ci_hi, qini$ci_hi),
    p_value  = c(autoc$p_value, qini$p_value)
  )
}

# -----------------------------------------------------------------------------
# RATE 並べ替え検定（n_perm 回 W をシャッフルして帰無分布を生成）
# negate=TRUE: 骨折アウトカム用（負の CATE = benefit → -tau_hat でランキング）
# -----------------------------------------------------------------------------
rate_permutation_test <- function(X, Y, W, tau_hat, negate = FALSE,
                                  clusters = NULL,
                                  n_perm = 500, seed = CF_SEED,
                                  num.trees = 1000) {
  set.seed(seed)
  obs_cf <- causal_forest(X, Y, W, W.hat = CF_W_HAT, num.trees = num.trees,
                          seed = seed, clusters = clusters)
  obs_rs <- if (negate) -tau_hat else tau_hat
  obs_autoc <- rank_average_treatment_effect(obs_cf, obs_rs)$estimate

  null_autoc <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    W_perm <- sample(W)
    cf_perm <- causal_forest(X, Y, W_perm, W.hat = CF_W_HAT,
                             num.trees = num.trees, seed = seed + i,
                             clusters = clusters)
    tau_perm <- predict(cf_perm)$predictions
    # 帰無分布も同じ方向でランキング
    rs_perm  <- if (negate) -tau_perm else tau_perm
    null_autoc[i] <- rank_average_treatment_effect(cf_perm, rs_perm)$estimate
  }

  p_val <- mean(null_autoc >= obs_autoc)
  list(
    observed_autoc = obs_autoc,
    null_autoc     = null_autoc,
    p_value        = p_val
  )
}

# -----------------------------------------------------------------------------
# Variable Importance
# -----------------------------------------------------------------------------
get_vi <- function(cf, covar_names) {
  vi <- variable_importance(cf)
  df <- data.frame(
    variable   = covar_names,
    importance = as.numeric(vi)
  )
  df[order(df$importance, decreasing = TRUE), ]
}

# -----------------------------------------------------------------------------
# Policy Tree（深さ 2）
# -----------------------------------------------------------------------------
get_policy_tree <- function(X, tau_hat, depth = 2) {
  pt <- policy_tree(X, cbind(-tau_hat, tau_hat), depth = depth)
  pt
}

# -----------------------------------------------------------------------------
# ストーリー判定（BLP 結果 + Story Mapping から）
# -----------------------------------------------------------------------------
determine_story <- function(blp_df, story_mapping = STORY_MAPPING) {
  sig_vars <- blp_df[!is.na(blp_df$p_value) & blp_df$p_value < 0.05, ]
  if (nrow(sig_vars) == 0) return("Undetermined")

  for (story in names(story_mapping)) {
    if (any(sig_vars$variable %in% story_mapping[[story]])) return(story)
  }
  "Other"
}

# -----------------------------------------------------------------------------
# フルパイプライン（1アウトカム）
# -----------------------------------------------------------------------------
run_cf_pipeline <- function(df, outcome_def, covariates,
                             use_clusters = TRUE, save_outputs = TRUE) {
  code  <- outcome_def$code
  col   <- outcome_def$col
  label <- outcome_def$label

  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("Outcome: %s — %s\n", code, label))
  cat(strrep("=", 60), "\n")

  # 1. データ準備
  prep <- prepare_data(df, col, covariates)
  if (prep$n_valid < 50) {
    cat("  Skipped: insufficient sample size\n")
    return(NULL)
  }

  # 2. CF フィット
  cat("  Fitting Causal Forest (W.hat=0.5, trees=", CF_NUM_TREES, ")...\n", sep = "")
  cf <- fit_cf(prep$X, prep$Y, prep$W, clusters = prep$clusters,
               use_clusters = use_clusters)

  # 3. CATE
  cate_df <- get_cate(cf)
  tau_hat  <- cate_df$cate

  # 4. ATE
  ate_df <- get_ate(cf)
  cat(sprintf("  ATE = %.4f (95%% CI: %.4f, %.4f), p=%.3f\n",
              ate_df$estimate, ate_df$ci_lo, ate_df$ci_hi, ate_df$p_value))

  # 5. Calibration test
  calib <- test_calibration(cf)

  # 6. BLP
  blp_raw <- best_linear_projection(cf, A = prep$X)
  blp_df  <- data.frame(
    variable  = rownames(blp_raw),
    estimate  = blp_raw[, "Estimate"],
    std_err   = blp_raw[, "Std. Error"],
    t_stat    = blp_raw[, "t value"],
    p_value   = blp_raw[, "Pr(>|t|)"]
  )
  blp_df <- blp_df[order(blp_df$p_value, na.last = TRUE), ]

  # 7. Variable Importance
  vi_df <- get_vi(cf, prep$covariates)

  # 8. GATE
  gate_df <- get_gate(cf, tau_hat)
  q_test  <- cochran_q_test(gate_df)
  cat(sprintf("  GATE Cochran Q=%.3f (df=%d, p=%.3f)\n",
              q_test$Q, q_test$df, q_test$p_value))

  # 9. RATE
  rate_df <- get_rate_df(cf, tau_hat)
  cat(sprintf("  AUTOC=%.4f (p=%.3f)\n",
              rate_df$estimate[1], rate_df$p_value[1]))

  # 10. Policy Tree
  pt <- tryCatch(
    get_policy_tree(prep$X, tau_hat, depth = 2),
    error = function(e) { cat("  PolicyTree error:", conditionMessage(e), "\n"); NULL }
  )

  # 11. ストーリー判定
  story <- determine_story(blp_df)
  cat(sprintf("  Story: %s\n", story))

  # ファイル出力
  if (save_outputs) {
    # CATE with ID
    cate_out <- cbind(
      data.frame(ID = prep$df_valid$ID, W = prep$W, Y = prep$Y),
      cate_df
    )
    write_csv(cate_out,    table_path(code, "cate"))
    write_csv(ate_df,      table_path(code, "ate"))
    write_csv(blp_df,      table_path(code, "blp"))
    write_csv(vi_df,       table_path(code, "vi"))
    write_csv(gate_df,     table_path(code, "gate"))
    write_csv(rate_df,     table_path(code, "rate"))
    write_csv(
      data.frame(
        beta_baseline  = calib["mean.forest.prediction",       "Estimate"],
        beta_het       = calib["differential.forest.prediction", "Estimate"],
        p_het          = calib["differential.forest.prediction", "Pr(>t)"]
      ),
      table_path(code, "calibration")
    )
    cat(sprintf("  Tables saved to %s/R_table_%s_*.csv\n", TABLE_DIR, code))
  }

  list(
    code      = code,
    label     = label,
    cf        = cf,
    prep      = prep,
    tau_hat   = tau_hat,
    ate_df    = ate_df,
    calib     = calib,
    blp_df    = blp_df,
    vi_df     = vi_df,
    gate_df   = gate_df,
    q_test    = q_test,
    rate_df   = rate_df,
    pt        = pt,
    story     = story
  )
}
